import ClimaCore.DataLayouts as DL
import ClimaCore.Fields
import ClimaComms
import ClimaCore as CC
import ClimaCore.Spaces
import OrdinaryDiffEq as ODE
import ClimaAtmos.Parameters as CAP
import DiffEqCallbacks as DEQ
import ClimaCore: InputOutput

import ClimaAtmos.TurbulenceConvection as TC

function get_callbacks(parsed_args, simulation, model_spec, params)
    FT = eltype(params)
    (; dt) = simulation

    callback_filters = call_every_n_steps(affect_filter!; skip_first = true)
    tc_callbacks =
        call_every_n_steps(turb_conv_affect_filter!; skip_first = true)

    additional_callbacks = if !isnothing(model_spec.radiation_model)
        # TODO: better if-else criteria?
        dt_rad = if parsed_args["config"] == "column"
            dt
        else
            FT(time_to_seconds(parsed_args["dt_rad"]))
        end
        (call_every_dt(rrtmgp_model_callback!, dt_rad),)
    else
        ()
    end

    if !isnothing(model_spec.turbconv_model) # && startswith(parsed_args["ode_algo"], "ODE.")
        additional_callbacks = (additional_callbacks..., tc_callbacks)
    end
    if model_spec.moisture_model isa EquilMoistModel &&
       parsed_args["apply_moisture_filter"]
        additional_callbacks = (additional_callbacks..., callback_filters)
    end

    dt_save_to_disk = time_to_seconds(parsed_args["dt_save_to_disk"])
    dt_save_restart = time_to_seconds(parsed_args["dt_save_restart"])

    dss_cb = if startswith(parsed_args["ode_algo"], "ODE.")
        call_every_n_steps(dss_callback)
    else
        nothing
    end
    save_to_disk_callback = if dt_save_to_disk == Inf
        nothing
    else
        call_every_dt(save_to_disk_func, dt_save_to_disk)
    end

    save_restart_callback = if dt_save_restart == Inf
        nothing
    else
        call_every_dt(save_restart_func, dt_save_restart)
    end
    return ODE.CallbackSet(
        dss_cb,
        save_to_disk_callback,
        save_restart_callback,
        additional_callbacks...,
        # DEQ.FunctionCallingCallback(
        #     (Y, t, integrator) -> begin
        #         @info "t = $t:"

        #         edmf_cache = integrator.p.edmf_cache
        #         a_min = edmf_cache.edmf.minimum_area
        #         ρarea = Spaces.column(Y.c.turbconv.up.:(1).ρarea, 1, 1, 1)
        #         ρaw = Spaces.column(Y.f.turbconv.up.:(1).ρaw, 1, 1, 1)
        #         ρ = Spaces.column(Y.c.ρ, 1, 1, 1)
        #         w = Spaces.column(edmf_cache.aux.face.turbconv.up.:(1).w, 1, 1, 1)
        #         ρ_f = Spaces.column(edmf_cache.aux.face.ρ, 1, 1, 1)

        #         # The actual values of SetValue BCs don't matter because they
        #         # don't appear in the Jacobian, so they can all be set to 0.
        #         ∇c = Operators.DivergenceF2C()
        #         ∇f = Operators.DivergenceC2F(;
        #             bottom = Operators.SetValue(Geometry.WVector(FT(0))),
        #             top = Operators.SetValue(Geometry.WVector(FT(0))),
        #         )
        #         LBC = Operators.LeftBiasedF2C(;
        #             bottom = Operators.SetValue(FT(0)),
        #         )
        #         LBF = Operators.LeftBiasedC2F(;
        #             bottom = Operators.SetValue(FT(0)),
        #         )
        #         Ic = Operators.InterpolateF2C()
        #         If = Operators.InterpolateC2F(;
        #             bottom = Operators.SetValue(FT(0)),
        #             top = Operators.SetValue(FT(0)),
        #         )

        #         ∇c_stencil = Operators.Operator2Stencil(∇c)
        #         ∇f_stencil = Operators.Operator2Stencil(∇f)
        #         LBC_stencil = Operators.Operator2Stencil(LBC)
        #         LBF_stencil = Operators.Operator2Stencil(LBF)
        #         Ic_stencil = Operators.Operator2Stencil(Ic)
        #         If_stencil = Operators.Operator2Stencil(If)
        #         compose = Operators.ComposeStencils()

        #         @. w = ifelse(
        #             If(ρarea / ρ) >= a_min,
        #             max(ρaw / (ρ_f * If(ρarea / ρ)), 0),
        #             0,
        #         ) # update w

        #         # w = ifelse(
        #         #     If(ρarea / ρ_c) >= a_min,
        #         #     max(ρaw / (ρ_f * If(ρarea / ρ_c)), 0),
        #         #     0
        #         # ) = ifelse(
        #         #     If(ρarea / ρ_c) >= a_min,
        #         #     ifelse(
        #         #         ρaw / (ρ_f * If(ρarea / ρ_c)) > 0,
        #         #         ρaw / (ρ_f * If(ρarea / ρ_c)),
        #         #         0,
        #         #     ),
        #         #     0,
        #         # )

        #         # ∂w/∂ρarea = ifelse(
        #         #     If(ρarea / ρ_c) >= a_min,
        #         #     ifelse(
        #         #         ρaw / (ρ_f * If(ρarea / ρ_c)) > 0,
        #         #         ∂(ρaw / (ρ_f * If(ρarea / ρ_c)))/∂ρarea,
        #         #         0,
        #         #     ),
        #         #     0,
        #         # )
        #         # ∂(ρaw / (ρ_f * If(ρarea / ρ_c)))/∂ρarea =
        #         #     -ρaw / (ρ_f * If(ρarea / ρ_c))^2 * ρ_f *
        #         #     ∂(If(ρarea / ρ_c))/∂ρarea
        #         # ∂(If(ρarea / ρ_c))/∂ρarea = If_stencil(1 / ρ_c)

        #         # ∂w/∂ρaw = ifelse(
        #         #     If(ρarea / ρ_c) >= a_min,
        #         #     ifelse(
        #         #         ρaw / (ρ_f * If(ρarea / ρ_c)) > 0,
        #         #         1 / (ρ_f * If(ρarea / ρ_c)),
        #         #         0,
        #         #     ),
        #         #     0,
        #         # )

        #         @info "∂ρareaₜ/∂ρarea"

        #         exact_block = exact_column_jacobian_block(
        #             TCU.sgs_flux_tendency_testing!,
        #             Y,
        #             integrator.p,
        #             t,
        #             1, # i
        #             1, # j
        #             1, # h
        #             (:c, :turbconv, :up, :1, :ρarea),
        #             (:c, :turbconv, :up, :1, :ρarea),
        #         )
        #         display(exact_block[1:10, 1:10])

        #         # ρareaₜ = -∇c(wvec(LBF(Ic(w) * ρarea)))
        #         # ∂ρareaₜ/∂ρarea =
        #         #     -∇c_stencil(wvec(1)) * LBF_stencil(1) * ∂(Ic(w) * ρarea)/∂ρarea
        #         # ∂(Ic(w) * ρarea)/∂ρarea =
        #         #     Ic_stencil(1) * ∂w/∂ρarea * ρarea + Ic(w)
        #         approx_block = @. compose(
        #             -∇c_stencil(Geometry.WVector(one(ρaw))),
        #             compose(
        #                 LBF_stencil(one(ρarea)),
        #                 compose(
        #                     Ic_stencil(one(w)),
        #                     ifelse(
        #                         If(ρarea / ρ) >= a_min,
        #                         ifelse(
        #                             ρaw / (ρ_f * If(ρarea / ρ)) > 0,
        #                             -ρaw / (ρ_f * If(ρarea / ρ))^2 * ρ_f,
        #                             0,
        #                         ),
        #                         0,
        #                     ) * If_stencil(1 / ρ),
        #                 ) * ρarea + Ic(w),
        #             ),
        #         )
        #         display(matrix_column(approx_block, axes(ρarea), 1, 1, 1)[1:10, 1:10])

        #         @info "∂ρaθ_liq_iceₜ/∂ρaθ_liq_ice and ∂ρaq_totₜ/∂ρaq_tot"

        #         exact_block = exact_column_jacobian_block(
        #             TCU.sgs_flux_tendency_testing!,
        #             Y,
        #             integrator.p,
        #             t,
        #             1, # i
        #             1, # j
        #             1, # h
        #             (:c, :turbconv, :up, :1, :ρaθ_liq_ice),
        #             (:c, :turbconv, :up, :1, :ρaθ_liq_ice),
        #         )
        #         display(exact_block[1:10, 1:10])
        #         exact_block = exact_column_jacobian_block(
        #             TCU.sgs_flux_tendency_testing!,
        #             Y,
        #             integrator.p,
        #             t,
        #             1, # i
        #             1, # j
        #             1, # h
        #             (:c, :turbconv, :up, :1, :ρaq_tot),
        #             (:c, :turbconv, :up, :1, :ρaq_tot),
        #         )
        #         display(exact_block[1:10, 1:10])

        #         # ρaq_totₜ = -∇c(wvec(LBF(Ic(w) * ρaq_tot)))
        #         # ∂ρaq_totₜ/∂ρaq_tot =
        #         #     -∇c_stencil(wvec(1)) * LBF_stencil(Ic(w))
        #         approx_block = @. compose(
        #             -∇c_stencil(Geometry.WVector(one(ρaw))),
        #             LBF_stencil(Ic(w)),
        #         )
        #         display(matrix_column(approx_block, axes(ρarea), 1, 1, 1)[1:10, 1:10])

        #         @info "∂ρawₜ/∂ρaw"

        #         exact_block = exact_column_jacobian_block(
        #             TCU.sgs_flux_tendency_testing!,
        #             Y,
        #             integrator.p,
        #             t,
        #             1, # i
        #             1, # j
        #             1, # h
        #             (:f, :turbconv, :up, :1, :ρaw),
        #             (:f, :turbconv, :up, :1, :ρaw),
        #         )
        #         display(exact_block[1:10, 1:10])

        #         # ρawₜ = -(∇f(wvec(LBC(ρaw * w))))
        #         # ∂ρawₜ/ρaw =
        #         #     -∇f_stencil(wvec(1)) * LBC_stencil(∂(ρaw * w)/∂ρaw)
        #         # ∂(ρaw * w)/∂ρaw = w + ρaw * ∂w/∂ρaw
        #         approx_block = @. compose(
        #             -∇f_stencil(Geometry.WVector(one(ρarea))),
        #             LBC_stencil(
        #                 w + ρaw * ifelse(
        #                     If(ρarea / ρ) >= a_min,
        #                     ifelse(
        #                         ρaw / (ρ_f * If(ρarea / ρ)) > 0,
        #                         1 / (ρ_f * If(ρarea / ρ)),
        #                         0,
        #                     ),
        #                     0,
        #                 ),
        #             ),
        #         )
        #         display(matrix_column(approx_block, axes(ρaw), 1, 1, 1)[1:10, 1:10])

        #         if t >= 500
        #             error("STOPPING")
        #         end
        #     end,
        #     func_start = true,
        # ),
    )
end

function call_every_n_steps(f!, n = 1; skip_first = false, call_at_end = false)
    previous_step = Ref(0)
    return ODE.DiscreteCallback(
        (u, t, integrator) ->
            (previous_step[] += 1) % n == 0 ||
                (call_at_end && t == integrator.sol.prob.tspan[2]),
        f!;
        initialize = (cb, u, t, integrator) -> skip_first || f!(integrator),
        save_positions = (false, false),
    )
end

function call_every_dt(f!, dt; skip_first = false, call_at_end = false)
    next_t = Ref{typeof(dt)}()
    affect! = function (integrator)
        t = integrator.t
        t_end = integrator.sol.prob.tspan[2]
        while t >= next_t[]
            f!(integrator)
            next_t[] =
                (call_at_end && t < t_end) ? min(t_end, next_t[] + dt) :
                next_t[] + dt
        end
    end
    return ODE.DiscreteCallback(
        (u, t, integrator) -> t >= next_t[],
        affect!;
        initialize = (cb, u, t, integrator) -> begin
            skip_first || f!(integrator)
            t_end = integrator.sol.prob.tspan[2]
            next_t[] =
                (call_at_end && t < t_end) ? min(t_end, t + dt) : t + dt
        end,
        save_positions = (false, false),
    )
end

function affect_filter!(Y::Fields.FieldVector)
    @. Y.c.ρq_tot = max(Y.c.ρq_tot, 0)
    return nothing
end

function dss_callback(integrator)
    Y = integrator.u
    ghost_buffer = integrator.p.ghost_buffer
    @nvtx "dss callback" color = colorant"yellow" begin
        Spaces.weighted_dss_start!(Y.c, ghost_buffer.c)
        Spaces.weighted_dss_start!(Y.f, ghost_buffer.f)
        Spaces.weighted_dss_internal!(Y.c, ghost_buffer.c)
        Spaces.weighted_dss_internal!(Y.f, ghost_buffer.f)
        Spaces.weighted_dss_ghost!(Y.c, ghost_buffer.c)
        Spaces.weighted_dss_ghost!(Y.f, ghost_buffer.f)
    end
    # ODE.u_modified!(integrator, false) # TODO: try this
end


function affect_filter!(integrator)
    (; apply_moisture_filter) = integrator.p
    affect_filter!(integrator.u)
    # We're lying to OrdinaryDiffEq.jl, in order to avoid
    # paying for an additional tendency call, which is required
    # to support supplying a continuous representation of the solution.
    ODE.u_modified!(integrator, false)
end

function turb_conv_affect_filter!(integrator)
    (; edmf_cache, Δt) = p
    (; edmf, param_set, aux, case, surf_params) = edmf_cache
    t = integrator.t
    Y = integrator.u
    tc_params = CAP.turbconv_params(param_set)

    Fields.bycolumn(axes(Y.c)) do colidx
        state = TCU.tc_column_state(Y, p, nothing, colidx)
        grid = TC.Grid(state)
        surf = TCU.get_surface(surf_params, grid, state, t, tc_params)
        TC.affect_filter!(edmf, grid, state, tc_params, surf, t)
    end

    # We're lying to OrdinaryDiffEq.jl, in order to avoid
    # paying for an additional `∑tendencies!` call, which is required
    # to support supplying a continuous representation of the
    # solution.
    ODE.u_modified!(integrator, false)
end

function save_to_disk_func(integrator)

    (; t, u, p) = integrator
    (; output_dir) = p.simulation
    Y = u

    if :ᶜS_ρq_tot in propertynames(p)
        (; ᶜts, ᶜp, ᶜS_ρq_tot, ᶜ3d_rain, ᶜ3d_snow, params, ᶜK) = p
    else
        (; ᶜts, ᶜp, params, ᶜK) = p
    end

    thermo_params = CAP.thermodynamics_params(params)
    cm_params = CAP.microphysics_params(params)

    ᶜuₕ = Y.c.uₕ
    ᶠw = Y.f.w
    # kinetic
    @. ᶜK = norm_sqr(C123(ᶜuₕ) + C123(ᶜinterp(ᶠw))) / 2

    # thermo state
    thermo_state!(ᶜts, Y, params, ᶜinterp, ᶜK)
    @. ᶜp = TD.air_pressure(thermo_params, ᶜts)
    ᶜT = @. TD.air_temperature(thermo_params, ᶜts)
    ᶜθ = @. TD.dry_pottemp(thermo_params, ᶜts)

    # vorticity
    curl_uh = @. curlₕ(Y.c.uₕ)
    ᶜvort = Geometry.WVector.(curl_uh)
    Spaces.weighted_dss!(ᶜvort)

    dry_diagnostic = (;
        pressure = ᶜp,
        temperature = ᶜT,
        potential_temperature = ᶜθ,
        kinetic_energy = ᶜK,
        vorticity = ᶜvort,
    )

    # cloudwater (liquid and ice), watervapor and RH for moist simulation
    if :ρq_tot in propertynames(Y.c)

        ᶜq = @. TD.PhasePartition(thermo_params, ᶜts)
        ᶜcloud_liquid = @. ᶜq.liq
        ᶜcloud_ice = @. ᶜq.ice
        ᶜwatervapor = @. TD.vapor_specific_humidity(ᶜq)
        ᶜRH = @. TD.relative_humidity(thermo_params, ᶜts)

        moist_diagnostic = (;
            cloud_liquid = ᶜcloud_liquid,
            cloud_ice = ᶜcloud_ice,
            water_vapor = ᶜwatervapor,
            relative_humidity = ᶜRH,
        )
        # precipitation
        if :ᶜS_ρq_tot in propertynames(p)

            @. ᶜS_ρq_tot =
                Y.c.ρ * CM.Microphysics0M.remove_precipitation(
                    cm_params,
                    TD.PhasePartition(thermo_params, ᶜts),
                )

            # rain vs snow
            @. ᶜ3d_rain = ifelse(ᶜT >= FT(273.15), ᶜS_ρq_tot, FT(0))
            @. ᶜ3d_snow = ifelse(ᶜT < FT(273.15), ᶜS_ρq_tot, FT(0))

            col_integrated_rain =
                vertical∫_col(ᶜ3d_rain) ./ FT(CAP.ρ_cloud_liq(params))
            col_integrated_snow =
                vertical∫_col(ᶜ3d_snow) ./ FT(CAP.ρ_cloud_liq(params))

            moist_diagnostic = (
                moist_diagnostic...,
                precipitation_removal = ᶜS_ρq_tot,
                column_integrated_rain = col_integrated_rain,
                column_integrated_snow = col_integrated_snow,
            )
        end
    else
        moist_diagnostic = NamedTuple()
    end

    if :edmf_cache in propertynames(p) && p.simulation.is_debugging_tc

        tc_cent(p) = p.edmf_cache.aux.cent.turbconv
        tc_face(p) = p.edmf_cache.aux.face.turbconv
        turbulence_convection_diagnostic = (;
            bulk_up_area = tc_cent(p).bulk.area,
            bulk_up_h_tot = tc_cent(p).bulk.h_tot,
            bulk_up_buoyancy = tc_cent(p).bulk.buoy,
            bulk_up_q_tot = tc_cent(p).bulk.q_tot,
            bulk_up_q_liq = tc_cent(p).bulk.q_liq,
            bulk_up_q_ice = tc_cent(p).bulk.q_ice,
            bulk_up_temperature = tc_cent(p).bulk.T,
            bulk_up_cloud_fraction = tc_cent(p).bulk.cloud_fraction,
            bulk_up_e_tot_tendency_precip_formation = tc_cent(
                p,
            ).bulk.e_tot_tendency_precip_formation,
            bulk_up_qt_tendency_precip_formation = tc_cent(
                p,
            ).bulk.qt_tendency_precip_formation,
            env_w = tc_cent(p).en.w,
            env_area = tc_cent(p).en.area,
            env_q_tot = tc_cent(p).en.q_tot,
            env_q_liq = tc_cent(p).en.q_liq,
            env_q_ice = tc_cent(p).en.q_ice,
            env_theta_liq_ice = tc_cent(p).en.θ_liq_ice,
            env_theta_virt = tc_cent(p).en.θ_virt,
            env_theta_dry = tc_cent(p).en.θ_dry,
            env_e_tot = tc_cent(p).en.e_tot,
            env_e_kin = tc_cent(p).en.e_kin,
            env_h_tot = tc_cent(p).en.h_tot,
            env_RH = tc_cent(p).en.RH,
            env_s = tc_cent(p).en.s,
            env_temperature = tc_cent(p).en.T,
            env_buoyancy = tc_cent(p).en.buoy,
            env_cloud_fraction = tc_cent(p).en.cloud_fraction,
            env_TKE = tc_cent(p).en.tke,
            env_Hvar = tc_cent(p).en.Hvar,
            env_QTvar = tc_cent(p).en.QTvar,
            env_HQTcov = tc_cent(p).en.HQTcov,
            env_e_tot_tendency_precip_formation = tc_cent(
                p,
            ).en.e_tot_tendency_precip_formation,
            env_qt_tendency_precip_formation = tc_cent(
                p,
            ).en.qt_tendency_precip_formation,
            env_Hvar_rain_dt = tc_cent(p).en.Hvar_rain_dt,
            env_QTvar_rain_dt = tc_cent(p).en.QTvar_rain_dt,
            env_HQTcov_rain_dt = tc_cent(p).en.HQTcov_rain_dt,
            face_bulk_w = tc_face(p).bulk.w,
            face_env_w = tc_face(p).en.w,
        )
    else
        turbulence_convection_diagnostic = NamedTuple()
    end

    if vert_diff
        (; dif_flux_uₕ, dif_flux_energy, dif_flux_ρq_tot) = p
        vert_diff_diagnostic = (;
            sfc_flux_momentum = dif_flux_uₕ,
            sfc_flux_energy = dif_flux_energy,
            sfc_evaporation = dif_flux_ρq_tot,
        )
    else
        vert_diff_diagnostic = NamedTuple()
    end

    if !isnothing(model_spec.radiation_model)
        (; face_lw_flux_dn, face_lw_flux_up, face_sw_flux_dn, face_sw_flux_up) =
            p.rrtmgp_model
        rad_diagnostic = (;
            lw_flux_down = RRTMGPI.array2field(FT.(face_lw_flux_dn), axes(Y.f)),
            lw_flux_up = RRTMGPI.array2field(FT.(face_lw_flux_up), axes(Y.f)),
            sw_flux_down = RRTMGPI.array2field(FT.(face_sw_flux_dn), axes(Y.f)),
            sw_flux_up = RRTMGPI.array2field(FT.(face_sw_flux_up), axes(Y.f)),
        )
        if model_spec.radiation_model isa
           RRTMGPI.AllSkyRadiationWithClearSkyDiagnostics
            (;
                face_clear_lw_flux_dn,
                face_clear_lw_flux_up,
                face_clear_sw_flux_dn,
                face_clear_sw_flux_up,
            ) = p.rrtmgp_model
            rad_clear_diagnostic = (;
                clear_lw_flux_down = RRTMGPI.array2field(
                    FT.(face_clear_lw_flux_dn),
                    axes(Y.f),
                ),
                clear_lw_flux_up = RRTMGPI.array2field(
                    FT.(face_clear_lw_flux_up),
                    axes(Y.f),
                ),
                clear_sw_flux_down = RRTMGPI.array2field(
                    FT.(face_clear_sw_flux_dn),
                    axes(Y.f),
                ),
                clear_sw_flux_up = RRTMGPI.array2field(
                    FT.(face_clear_sw_flux_up),
                    axes(Y.f),
                ),
            )
        else
            rad_clear_diagnostic = NamedTuple()
        end
    else
        rad_diagnostic = NamedTuple()
        rad_clear_diagnostic = NamedTuple()
    end

    diagnostic = merge(
        dry_diagnostic,
        moist_diagnostic,
        vert_diff_diagnostic,
        rad_diagnostic,
        rad_clear_diagnostic,
        turbulence_convection_diagnostic,
    )

    day = floor(Int, t / (60 * 60 * 24))
    sec = floor(Int, t % (60 * 60 * 24))
    @info "Saving diagnostics to HDF5 file on day $day second $sec"
    output_file = joinpath(output_dir, "day$day.$sec.hdf5")
    hdfwriter = InputOutput.HDF5Writer(output_file, comms_ctx)
    InputOutput.HDF5.write_attribute(hdfwriter.file, "time", t) # TODO: a better way to write metadata
    InputOutput.write!(hdfwriter, Y, "Y")
    InputOutput.write!(
        hdfwriter,
        Fields.FieldVector(; pairs(diagnostic)...),
        "diagnostics",
    )
    Base.close(hdfwriter)
    return nothing
end

function save_restart_func(integrator)
    (; t, u, p) = integrator
    (; output_dir) = p.simulation
    Y = u
    day = floor(Int, t / (60 * 60 * 24))
    sec = floor(Int, t % (60 * 60 * 24))
    @info "Saving restart file to HDF5 file on day $day second $sec"
    mkpath(joinpath(output_dir, "restart"))
    output_file = joinpath(output_dir, "restart", "day$day.$sec.hdf5")
    hdfwriter = InputOutput.HDF5Writer(output_file, comms_ctx)
    InputOutput.HDF5.write_attribute(hdfwriter.file, "time", t) # TODO: a better way to write metadata
    InputOutput.write!(hdfwriter, Y, "Y")
    Base.close(hdfwriter)
    return nothing
end
