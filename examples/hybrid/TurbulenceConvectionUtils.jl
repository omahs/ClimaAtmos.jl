module TurbulenceConvectionUtils

using LinearAlgebra
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

#####
##### No TurbulenceConvection scheme
#####

turbconv_cache(
    Y,
    turbconv_model::Nothing,
    precip_model,
    namelist,
    param_set,
    parsed_args,
) = (; turbconv_model)

sgs_flux_tendency!(Yₜ, Y, p, t, colidx, ::Nothing) = nothing

#####
##### EDMF
#####

function get_aux(edmf, Y, ::Type{FT}) where {FT}
    fspace = axes(Y.f)
    cspace = axes(Y.c)
    aux_cent_fields = TC.FieldFromNamedTuple(cspace, cent_aux_vars, FT, edmf)
    aux_face_fields = TC.FieldFromNamedTuple(fspace, face_aux_vars, FT, edmf)
    aux = CC.Fields.FieldVector(cent = aux_cent_fields, face = aux_face_fields)
    return aux
end

function turbconv_cache(
    Y,
    turbconv_model::TC.EDMFModel,
    precip_model,
    namelist,
    param_set,
    parsed_args,
)
    tc_params = CAP.turbconv_params(param_set)
    Ri_bulk_crit = namelist["turbulence"]["EDMF_PrognosticTKE"]["Ri_crit"]
    FT = CC.Spaces.undertype(axes(Y.c))
    test_consistency = parsed_args["test_edmf_consistency"]
    case = Cases.get_case(namelist["meta"]["casename"])
    thermo_params = CAP.thermodynamics_params(param_set)
    surf_thermo_state = Cases.surface_thermo_state(case, thermo_params)
    surf_params = Cases.surface_params(
        case,
        surf_thermo_state,
        thermo_params;
        Ri_bulk_crit,
    )
    edmf = turbconv_model
    ᶠspace_1 = axes(Y.f[CC.Fields.ColumnIndex((1, 1), 1)])
    ᶜspace_1 = axes(Y.c[CC.Fields.ColumnIndex((1, 1), 1)])
    logpressure_fun =
        CA.log_pressure_profile(ᶠspace_1, thermo_params, surf_thermo_state)
    ᶠz = CC.Fields.coordinate_field(ᶠspace_1).z
    ᶜz = CC.Fields.coordinate_field(ᶜspace_1).z
    ᶠp₀ = @. exp(logpressure_fun(ᶠz))
    ᶜp₀ = @. exp(logpressure_fun(ᶜz))
    @info "EDMFModel: \n$(summary(edmf))"
    cache = (;
        edmf,
        turbconv_model,
        ᶠp₀,
        ᶜp₀,
        case,
        test_consistency,
        surf_params,
        param_set,
        surf_thermo_state,
        aux = get_aux(edmf, Y, FT),
        precip_model,
    )
    return (; edmf_cache = cache, turbconv_model)
end

function init_tc!(Y, p, params)
    CC.Fields.bycolumn(axes(Y.c)) do colidx
        init_tc!(Y, p, params, colidx)
    end
end

function init_tc!(Y, p, params, colidx)

    (; edmf, surf_thermo_state, ᶠp₀, ᶜp₀, surf_params, case) = p.edmf_cache
    tc_params = CAP.turbconv_params(params)

    FT = eltype(edmf)
    # `nothing` goes into State because OrdinaryDiffEq.jl owns tendencies.
    state = TC.tc_column_state(Y, p, nothing, colidx)
    thermo_params = CAP.thermodynamics_params(params)

    grid = TC.Grid(state)
    FT = eltype(grid)
    t = FT(0)

    @. p.ᶜp[colidx] = ᶜp₀
    @. p.edmf_cache.aux.face.p[colidx] = ᶠp₀

    CA.compute_ref_density!(
        Y.c.ρ[colidx],
        p.ᶜp[colidx],
        thermo_params,
        surf_thermo_state,
    )
    CA.compute_ref_density!(
        p.edmf_cache.aux.face.ρ[colidx],
        p.edmf_cache.aux.face.p[colidx],
        thermo_params,
        surf_thermo_state,
    )

    # TODO: convert initialize_profiles to set prognostic state, not aux state
    Cases.initialize_profiles(case, grid, thermo_params, state)

    # Temporarily, we'll re-populate ρq_tot based on initial aux q_tot
    q_tot = p.edmf_cache.aux.cent.q_tot[colidx]
    @. Y.c.ρq_tot[colidx] = Y.c.ρ[colidx] * q_tot
    set_thermo_state_pθq!(Y, p, colidx)
    set_grid_mean_from_thermo_state!(thermo_params, state, grid)
    assign_thermo_aux!(state, grid, edmf.moisture_model, thermo_params)
    initialize_edmf(edmf, grid, state, surf_params, tc_params, t)
end


function sgs_flux_tendency!(Yₜ, Y, p, t, colidx, ::TC.EDMFModel)
    t < 100 && println("SGS: ‖Yₜ_aa($t)‖ = $(norm(Yₜ))")
    for prop_chain in CC.Fields.property_chains(Yₜ)
        field = CC.Fields.single_field(Yₜ, prop_chain)
        name = join(prop_chain, '.')
        t < 100 && println("‖Yₜ.$name($t)‖ = $(norm(parent(field)))")
    end

    (; edmf_cache, Δt, compressibility_model) = p
    (; edmf, param_set, surf_params, surf_thermo_state) = edmf_cache
    (; precip_model, test_consistency, ᶠp₀, ᶜp₀) = edmf_cache
    thermo_params = CAP.thermodynamics_params(param_set)
    tc_params = CAP.turbconv_params(param_set)
    state = TC.tc_column_state(Y, p, Yₜ, colidx)
    grid = TC.Grid(state)
    if test_consistency
        parent(state.aux.face) .= NaN
        parent(state.aux.cent) .= NaN
    end

    if compressibility_model isa CA.AnelasticFluid
        # TODO: how should this be computed if compressible?
        # pressure at cell centers have been populated, we need
        # pressure BCs
        @. p.edmf_cache.aux.face.p[colidx] = ᶠp₀
        CA.compute_ref_density!(
            p.edmf_cache.aux.face.ρ[colidx],
            p.edmf_cache.aux.face.p[colidx],
            thermo_params,
            surf_thermo_state,
        )
    end
    assign_thermo_aux!(state, grid, edmf.moisture_model, thermo_params)

    surf = get_surface(surf_params, grid, state, t, tc_params)

    TC.affect_filter!(edmf, grid, state, tc_params, surf, t)

    TC.update_aux!(edmf, grid, state, surf, tc_params, t, Δt)

    t < 100 && println("SGS: ‖Yₜ_a($t)‖ = $(norm(Yₜ))")
    for prop_chain in CC.Fields.property_chains(Yₜ)
        field = CC.Fields.single_field(Yₜ, prop_chain)
        name = join(prop_chain, '.')
        t < 100 && println("‖Yₜ.$name($t)‖ = $(norm(parent(field)))")
    end

    TC.compute_precipitation_sink_tendencies(
        precip_model,
        edmf.precip_fraction_model,
        grid,
        state,
        tc_params,
        Δt,
    )

    t < 100 && println("SGS: ‖Yₜ_b($t)‖ = $(norm(Yₜ))")
    for prop_chain in CC.Fields.property_chains(Yₜ)
        field = CC.Fields.single_field(Yₜ, prop_chain)
        name = join(prop_chain, '.')
        t < 100 && println("‖Yₜ.$name($t)‖ = $(norm(parent(field)))")
    end

    TC.compute_turbconv_tendencies!(edmf, grid, state, tc_params, surf, Δt)

    t < 100 && println("SGS: ‖Yₜ_c($t)‖ = $(norm(Yₜ))")
    for prop_chain in CC.Fields.property_chains(Yₜ)
        field = CC.Fields.single_field(Yₜ, prop_chain)
        name = join(prop_chain, '.')
        t < 100 && println("‖Yₜ.$name($t)‖ = $(norm(parent(field)))")
    end

    # TODO: incrementally disable this and enable proper grid mean terms
    compute_gm_tendencies!(edmf, grid, state, surf, tc_params)

    t < 100 && println("SGS: ‖Yₜ_d($t)‖ = $(norm(Yₜ))")
    for prop_chain in CC.Fields.property_chains(Yₜ)
        field = CC.Fields.single_field(Yₜ, prop_chain)
        name = join(prop_chain, '.')
        t < 100 && println("‖Yₜ.$name($t)‖ = $(norm(parent(field)))")
    end

    return nothing
end

end # module
