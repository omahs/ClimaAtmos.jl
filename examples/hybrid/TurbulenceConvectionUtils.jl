module TurbulenceConvectionUtils

using LinearAlgebra, StaticArrays
import ClimaAtmos
import ClimaAtmos.Parameters as CAP
import ClimaCore as CC
import ClimaCore.Geometry as CCG
import ClimaCore.Operators as CCO
import ClimaCore.Geometry: ⊗
import OrdinaryDiffEq as ODE
import ClimaAtmos.TurbulenceConvection
import ClimaAtmos.TurbulenceConvection as TC

import UnPack
import Logging
import TerminalLoggers

const ca_dir = pkgdir(ClimaAtmos)

include(joinpath(ca_dir, "tc_driver", "Cases.jl"))
import .Cases

include(joinpath(ca_dir, "tc_driver", "dycore.jl"))
include(joinpath(ca_dir, "tc_driver", "Surface.jl"))
include(joinpath(ca_dir, "tc_driver", "initial_conditions.jl"))
include(joinpath(ca_dir, "tc_driver", "generate_namelist.jl"))
import .NameList

function get_aux(edmf, Y, ::Type{FT}) where {FT}
    fspace = axes(Y.f)
    cspace = axes(Y.c)
    aux_cent_fields = TC.FieldFromNamedTuple(cspace, cent_aux_vars, FT, edmf)
    aux_face_fields = TC.FieldFromNamedTuple(fspace, face_aux_vars, FT, edmf)
    aux = CC.Fields.FieldVector(cent = aux_cent_fields, face = aux_face_fields)
    return aux
end

function get_edmf_cache(Y, namelist, param_set, parsed_args)
    tc_params = CAP.turbconv_params(param_set)
    Ri_bulk_crit = namelist["turbulence"]["EDMF_PrognosticTKE"]["Ri_crit"]
    case = Cases.get_case(namelist)
    FT = CC.Spaces.undertype(axes(Y.c))
    forcing =
        Cases.ForcingBase(case, FT; Cases.forcing_kwargs(case, namelist)...)
    radiation = Cases.RadiationBase(case, FT)
    surf_ref_state = Cases.surface_ref_state(case, tc_params, namelist)
    surf_params =
        Cases.surface_params(case, surf_ref_state, tc_params; Ri_bulk_crit)
    precip_name = TC.parse_namelist(
        namelist,
        "microphysics",
        "precipitation_model";
        default = "None",
        valid_options = ["None", "cutoff", "clima_1m"],
    )
    # TODO: move to grid mean model
    precip_model = if precip_name == "None"
        TC.NoPrecipitation()
    elseif precip_name == "cutoff"
        TC.CutoffPrecipitation()
    elseif precip_name == "clima_1m"
        TC.Clima1M()
    else
        error("Invalid precip_name $(precip_name)")
    end
    edmf = TC.EDMFModel(FT, namelist, precip_model, parsed_args)
    @show "EDMFModel: \n$(summary(edmf))"
    return (;
        edmf,
        case,
        forcing,
        radiation,
        surf_params,
        param_set,
        surf_ref_state,
        aux = get_aux(edmf, Y, FT),
        precip_model,
    )
end

function tc_column_state(prog, p, tendencies, colidx)
    prog_cent_column = CC.column(prog.c, colidx)
    prog_face_column = CC.column(prog.f, colidx)
    aux_cent_column = CC.column(p.edmf_cache.aux.cent, colidx)
    aux_face_column = CC.column(p.edmf_cache.aux.face, colidx)
    tends_cent_column = CC.column(tendencies.c, colidx)
    tends_face_column = CC.column(tendencies.f, colidx)
    prog_column =
        CC.Fields.FieldVector(cent = prog_cent_column, face = prog_face_column)
    aux_column =
        CC.Fields.FieldVector(cent = aux_cent_column, face = aux_face_column)
    tends_column = CC.Fields.FieldVector(
        cent = tends_cent_column,
        face = tends_face_column,
    )

    return TC.State(prog_column, aux_column, tends_column, p, colidx)
end

function tc_column_state(prog, p, tendencies::Nothing, colidx)
    prog_cent_column = CC.column(prog.c, colidx)
    prog_face_column = CC.column(prog.f, colidx)
    aux_cent_column = CC.column(p.edmf_cache.aux.cent, colidx)
    aux_face_column = CC.column(p.edmf_cache.aux.face, colidx)
    prog_column =
        CC.Fields.FieldVector(cent = prog_cent_column, face = prog_face_column)
    aux_column =
        CC.Fields.FieldVector(cent = aux_cent_column, face = aux_face_column)
    tends_column = nothing

    return TC.State(prog_column, aux_column, tends_column, p, colidx)
end


function init_tc!(Y, p, param_set, namelist)
    (; edmf_cache, Δt) = p
    (; edmf, param_set, surf_ref_state, surf_params, forcing, radiation, case) =
        edmf_cache
    tc_params = CAP.turbconv_params(param_set)

    FT = eltype(edmf)
    N_up = TC.n_updrafts(edmf)

    CC.Fields.bycolumn(axes(Y.c)) do colidx
        # `nothing` goes into State because OrdinaryDiffEq.jl owns tendencies.
        state = tc_column_state(Y, p, nothing, colidx)

        grid = TC.Grid(state)
        FT = eltype(grid)
        t = FT(0)
        compute_ref_state!(state, grid, tc_params; ts_g = surf_ref_state)

        Cases.initialize_profiles(case, grid, tc_params, state)
        set_thermo_state_pθq!(state, grid, edmf.moisture_model, tc_params)
        set_grid_mean_from_thermo_state!(tc_params, state, grid)
        assign_thermo_aux!(state, grid, edmf.moisture_model, tc_params)
        Cases.initialize_forcing(case, forcing, grid, state, tc_params)
        Cases.initialize_radiation(case, radiation, grid, state, tc_params)
        initialize_edmf(edmf, grid, state, surf_params, tc_params, t, case)
    end
end


function sgs_flux_tendency!(Yₜ, Y, p, t)
    (; edmf_cache, Δt) = p
    (; edmf, param_set, case, surf_params, radiation, forcing, precip_model) =
        edmf_cache
    tc_params = CAP.turbconv_params(param_set)

    # TODO: write iterator for this
    CC.Fields.bycolumn(axes(Y.c)) do colidx
        state = tc_column_state(Y, p, Yₜ, colidx)
        grid = TC.Grid(state)

        set_thermo_state_peq!(
            state,
            grid,
            edmf.moisture_model,
            edmf.compressibility_model,
            tc_params,
        )
        assign_thermo_aux!(state, grid, edmf.moisture_model, tc_params)

        aux_gm = TC.center_aux_grid_mean(state)

        surf = get_surface(surf_params, grid, state, t, tc_params)

        TC.affect_filter!(edmf, grid, state, tc_params, surf, t)

        # Update aux / pre-tendencies filters. TODO: combine these into a function that minimizes traversals
        # Some of these methods should probably live in `compute_tendencies`, when written, but we'll
        # treat them as auxiliary variables for now, until we disentangle the tendency computations.
        Cases.update_forcing(case, grid, state, t, tc_params)
        Cases.update_radiation(radiation, grid, state, t, tc_params)

        TC.update_aux!(edmf, grid, state, surf, tc_params, t, Δt)

        TC.compute_precipitation_sink_tendencies(
            precip_model,
            edmf,
            grid,
            state,
            tc_params,
            Δt,
        )
        TC.compute_precipitation_advection_tendencies(
            precip_model,
            edmf,
            grid,
            state,
            tc_params,
        )

        TC.compute_turbconv_tendencies!(edmf, grid, state, tc_params, surf, Δt)

        # TODO: incrementally disable this and enable proper grid mean terms
        compute_gm_tendencies!(
            edmf,
            grid,
            state,
            surf,
            radiation,
            forcing,
            tc_params,
        )
    end
end

function compute_up_tendencies_testing!(
    edmf::TC.EDMFModel,
    grid::TC.Grid,
    state::TC.State,
    param_set::TC.APS,
    surf::TC.SurfaceBase,
)
    N_up = TC.n_updrafts(edmf)
    kc_surf = TC.kc_surface(grid)
    kf_surf = TC.kf_surface(grid)
    FT = TC.float_type(state)

    aux_up = TC.center_aux_updrafts(state)
    aux_en = TC.center_aux_environment(state)
    aux_en_f = TC.face_aux_environment(state)
    aux_up_f = TC.face_aux_updrafts(state)
    prog_up = TC.center_prog_updrafts(state)
    prog_up_f = TC.face_prog_updrafts(state)
    tendencies_up = TC.center_tendencies_updrafts(state)
    tendencies_up_f = TC.face_tendencies_updrafts(state)
    prog_gm = TC.center_prog_grid_mean(state)
    aux_gm_f = TC.face_aux_grid_mean(state)
    ρ_c = prog_gm.ρ
    ρ_f = aux_gm_f.ρ
    au_lim = edmf.max_area

    @inbounds for i in 1:N_up
        aux_up_i = aux_up[i]
        @. aux_up_i.entr_turb_dyn = aux_up_i.entr_sc + aux_up_i.frac_turb_entr
        @. aux_up_i.detr_turb_dyn = aux_up_i.detr_sc + aux_up_i.frac_turb_entr
    end

    UB = TC.CCO.UpwindBiasedProductC2F(
        bottom = TC.CCO.SetValue(FT(0)),
        top = TC.CCO.SetValue(FT(0)),
    )
    Ic = TC.CCO.InterpolateF2C()

    wvec = TC.CC.Geometry.WVector
    ∇c = TC.CCO.DivergenceF2C()
    w_bcs =
        (; bottom = TC.CCO.SetValue(wvec(FT(0))), top = TC.CCO.SetValue(wvec(FT(0))))
    LBF = TC.CCO.LeftBiasedC2F(; bottom = TC.CCO.SetValue(FT(0)))

    # Solve for updraft area fraction
    @inbounds for i in 1:N_up
        aux_up_i = aux_up[i]
        w_up = aux_up_f[i].w
        a_up = aux_up_i.area
        q_tot_up = aux_up_i.q_tot
        q_tot_en = aux_en.q_tot
        θ_liq_ice_en = aux_en.θ_liq_ice
        θ_liq_ice_up = aux_up_i.θ_liq_ice
        entr_turb_dyn = aux_up_i.entr_turb_dyn
        detr_turb_dyn = aux_up_i.detr_turb_dyn
        θ_liq_ice_tendency_precip_formation =
            aux_up_i.θ_liq_ice_tendency_precip_formation
        qt_tendency_precip_formation = aux_up_i.qt_tendency_precip_formation

        ρarea = prog_up[i].ρarea
        ρaθ_liq_ice = prog_up[i].ρaθ_liq_ice
        ρaq_tot = prog_up[i].ρaq_tot

        tends_ρarea = tendencies_up[i].ρarea
        tends_ρaθ_liq_ice = tendencies_up[i].ρaθ_liq_ice
        tends_ρaq_tot = tendencies_up[i].ρaq_tot

        @. tends_ρarea =
            -∇c(wvec(LBF(Ic(w_up) * ρarea))) # +
            # (ρarea * Ic(w_up) * entr_turb_dyn) -
            # (ρarea * Ic(w_up) * detr_turb_dyn)

    #     @. tends_ρaθ_liq_ice =
    #         -∇c(wvec(LBF(Ic(w_up) * ρaθ_liq_ice))) +
    #         (ρarea * Ic(w_up) * entr_turb_dyn * θ_liq_ice_en) -
    #         (ρaθ_liq_ice * Ic(w_up) * detr_turb_dyn) +
    #         (ρ_c * θ_liq_ice_tendency_precip_formation)

    #     @. tends_ρaq_tot =
    #         -∇c(wvec(LBF(Ic(w_up) * ρaq_tot))) +
    #         (ρarea * Ic(w_up) * entr_turb_dyn * q_tot_en) -
    #         (ρaq_tot * Ic(w_up) * detr_turb_dyn) +
    #         (ρ_c * qt_tendency_precip_formation)

    #     if edmf.moisture_model isa NonEquilibriumMoisture

    #         q_liq_up = aux_up_i.q_liq
    #         q_ice_up = aux_up_i.q_ice
    #         q_liq_en = aux_en.q_liq
    #         q_ice_en = aux_en.q_ice

    #         ql_tendency_noneq = aux_up_i.ql_tendency_noneq
    #         qi_tendency_noneq = aux_up_i.qi_tendency_noneq
    #         ql_tendency_precip_formation = aux_up_i.ql_tendency_precip_formation
    #         qi_tendency_precip_formation = aux_up_i.qi_tendency_precip_formation

    #         ρaq_liq = prog_up[i].ρaq_liq
    #         ρaq_ice = prog_up[i].ρaq_ice

    #         tends_ρaq_liq = tendencies_up[i].ρaq_liq
    #         tends_ρaq_ice = tendencies_up[i].ρaq_ice

    #         @. tends_ρaq_liq =
    #             -∇c(wvec(LBF(Ic(w_up) * ρaq_liq))) +
    #             (ρarea * Ic(w_up) * entr_turb_dyn * q_liq_en) -
    #             (ρaq_liq * Ic(w_up) * detr_turb_dyn) +
    #             (ρ_c * (ql_tendency_precip_formation + ql_tendency_noneq))

    #         @. tends_ρaq_ice =
    #             -∇c(wvec(LBF(Ic(w_up) * ρaq_ice))) +
    #             (ρarea * Ic(w_up) * entr_turb_dyn * q_ice_en) -
    #             (ρaq_ice * Ic(w_up) * detr_turb_dyn) +
    #             (ρ_c * (qi_tendency_precip_formation + qi_tendency_noneq))

    #         tends_ρaq_liq[kc_surf] = 0
    #         tends_ρaq_ice[kc_surf] = 0
        end

    #     # prognostic entr/detr
    #     if edmf.entr_closure isa PrognosticNoisyRelaxationProcess
    #         c_gen_stoch = edmf.entr_closure.c_gen_stoch
    #         mean_entr = aux_up[i].ε_nondim
    #         mean_detr = aux_up[i].δ_nondim
    #         ε_λ = c_gen_stoch[3]
    #         δ_λ = c_gen_stoch[4]
    #         tends_ε_nondim = tendencies_up[i].ε_nondim
    #         tends_δ_nondim = tendencies_up[i].δ_nondim
    #         ε_nondim = prog_up[i].ε_nondim
    #         δ_nondim = prog_up[i].δ_nondim
    #         @. tends_ε_nondim = ε_λ * (mean_entr - ε_nondim)
    #         @. tends_δ_nondim = δ_λ * (mean_detr - δ_nondim)
    #     end
    #     tends_ρarea[kc_surf] = 0
    #     tends_ρaθ_liq_ice[kc_surf] = 0
    #     tends_ρaq_tot[kc_surf] = 0
    # end

    # # Solve for updraft velocity

    # # We know that, since W = 0 at z = 0, BCs for entr, detr,
    # # and buoyancy should not matter in the end
    # zero_bcs = (; bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0)))
    # I0f = CCO.InterpolateC2F(; zero_bcs...)
    # adv_bcs =
    #     (; bottom = CCO.SetValue(wvec(FT(0))), top = CCO.SetValue(wvec(FT(0))))
    # LBC = CCO.LeftBiasedF2C(; bottom = CCO.SetValue(FT(0)))
    # ∇f = CCO.DivergenceC2F(; adv_bcs...)

    # @inbounds for i in 1:N_up
    #     a_up_bcs = a_up_boundary_conditions(surf, edmf, i)
    #     Iaf = CCO.InterpolateC2F(; a_up_bcs...)
    #     ρaw = prog_up_f[i].ρaw
    #     tends_ρaw = tendencies_up_f[i].ρaw
    #     nh_pressure = aux_up_f[i].nh_pressure
    #     a_up = aux_up[i].area
    #     w_up = aux_up_f[i].w
    #     w_en = aux_en_f.w
    #     entr_w = aux_up[i].entr_turb_dyn
    #     detr_w = aux_up[i].detr_turb_dyn
    #     buoy = aux_up[i].buoy

    #     @. tends_ρaw = -(∇f(wvec(LBC(ρaw * w_up))))
    #     @. tends_ρaw +=
    #         (ρaw * (I0f(entr_w) * w_en - I0f(detr_w) * w_up)) +
    #         (ρ_f * Iaf(a_up) * I0f(buoy)) +
    #         nh_pressure
    #     tends_ρaw[kf_surf] = 0
    # end

    return nothing
end

function compute_turbconv_tendencies_testing!(
    edmf::TC.EDMFModel,
    grid::TC.Grid,
    state::TC.State,
    param_set::TC.APS,
    surf::TC.SurfaceBase,
    Δt::Real,
)
    compute_up_tendencies_testing!(edmf, grid, state, param_set, surf)
    # compute_en_tendencies!(edmf, grid, state, param_set, Val(:tke), Val(:ρatke))

    # if edmf.thermo_covariance_model isa PrognosticThermoCovariances
    #     compute_en_tendencies!(
    #         edmf,
    #         grid,
    #         state,
    #         param_set,
    #         Val(:Hvar),
    #         Val(:ρaHvar),
    #     )
    #     compute_en_tendencies!(
    #         edmf,
    #         grid,
    #         state,
    #         param_set,
    #         Val(:QTvar),
    #         Val(:ρaQTvar),
    #     )
    #     compute_en_tendencies!(
    #         edmf,
    #         grid,
    #         state,
    #         param_set,
    #         Val(:HQTcov),
    #         Val(:ρaHQTcov),
    #     )
    # end

    return nothing
end

function sgs_flux_tendency_testing!(Yₜ, Y, p, t)
    Yₜ .= 0 # sanity check
    (; edmf_cache, Δt) = p
    (;
        edmf,
        param_set,
        aux,
        case,
        surf_params,
        radiation,
        forcing,
        precip_model,
    ) = edmf_cache
    tc_params = CAP.turbconv_params(param_set)

    # TODO: write iterator for this
    for inds in TC.iterate_columns(Y.c)
        state = tc_column_state(Y, aux, Yₜ, inds...)
        grid = TC.Grid(state)

        set_thermo_state_peq!(
            state,
            grid,
            edmf.moisture_model,
            edmf.compressibility_model,
            tc_params,
        )
        assign_thermo_aux!(state, grid, edmf.moisture_model, tc_params)

        aux_gm = TC.center_aux_grid_mean(state)

        surf = get_surface(surf_params, grid, state, t, tc_params)

        TC.affect_filter!(edmf, grid, state, tc_params, surf, t)

        # Update aux / pre-tendencies filters. TODO: combine these into a function that minimizes traversals
        # Some of these methods should probably live in `compute_tendencies`, when written, but we'll
        # treat them as auxiliary variables for now, until we disentangle the tendency computations.
        # Cases.update_forcing(case, grid, state, t, tc_params)
        # Cases.update_radiation(radiation, grid, state, t, tc_params)

        TC.update_aux!(edmf, grid, state, surf, tc_params, t, Δt)

        # TC.compute_precipitation_sink_tendencies(
        #     precip_model,
        #     edmf,
        #     grid,
        #     state,
        #     tc_params,
        #     Δt,
        # )
        # TC.compute_precipitation_advection_tendencies(
        #     precip_model,
        #     edmf,
        #     grid,
        #     state,
        #     tc_params,
        # )

        compute_turbconv_tendencies_testing!(edmf, grid, state, tc_params, surf, Δt)

        # TODO: incrementally disable this and enable proper grid mean terms
        # compute_gm_tendencies!(
        #     edmf,
        #     grid,
        #     state,
        #     surf,
        #     radiation,
        #     forcing,
        #     tc_params,
        # )
    end
end

end # module
