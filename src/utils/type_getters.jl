using Dates: DateTime, @dateformat_str
using NCDatasets
using Dierckx
using DiffEqBase
using ImageFiltering
using Interpolations
import ClimaCore: InputOutput, Meshes, Spaces
import ClimaAtmos.RRTMGPInterface as RRTMGPI
import ClimaAtmos as CA
import LinearAlgebra
import ClimaCore.Fields
import OrdinaryDiffEq as ODE
import ClimaTimeSteppers as CTS
import DiffEqCallbacks as DEQ

function get_atmos(::Type{FT}, params, turbconv_params) where {FT}

    # should this live in the radiation model?

    moisture_model = get_moisture_model(params)
    precip_model = get_precipitation_model(params)
    radiation_mode = get_radiation_mode(params, FT)
    forcing_type = get_forcing_type(params)
    surface_scheme = get_surface_scheme(FT, params)

    diffuse_momentum =
        !(forcing_type isa HeldSuarezForcing) && !isnothing(surface_scheme)

    edmfx_adv_test = params["edmfx_adv_test"]
    @assert edmfx_adv_test in (false, true)

    edmfx_entr_detr = params["edmfx_entr_detr"]
    @assert edmfx_entr_detr in (false, true)

    model_config = get_model_config(params)
    vert_diff = get_vertical_diffusion_model(diffuse_momentum, params, FT)
    atmos = AtmosModel(;
        moisture_model,
        model_config,
        coupling = get_coupling_type(params),
        perf_mode = get_perf_mode(params),
        energy_form = get_energy_form(params, vert_diff),
        radiation_mode,
        subsidence = get_subsidence_model(params, radiation_mode, FT),
        ls_adv = get_large_scale_advection_model(params, FT),
        edmf_coriolis = get_edmf_coriolis(params, FT),
        edmfx_adv_test,
        edmfx_entr_detr,
        precip_model,
        forcing_type,
        turbconv_model = get_turbconv_model(
            FT,
            moisture_model,
            precip_model,
            params,
            turbconv_params,
        ),
        surface_scheme,
        non_orographic_gravity_wave = get_non_orographic_gravity_wave_model(
            params,
            model_config,
            FT,
        ),
        orographic_gravity_wave = get_orographic_gravity_wave_model(params, FT),
        hyperdiff = get_hyperdiffusion_model(params, FT),
        vert_diff,
        viscous_sponge = get_viscous_sponge_model(params, FT),
        rayleigh_sponge = get_rayleigh_sponge_model(params, FT),
    )

    @info "AtmosModel: \n$(summary(atmos))"
    return atmos
end

function get_numerics(params)
    # wrap each upwinding mode in a Val for dispatch
    numerics = (;
        energy_upwinding = Val(Symbol(params["energy_upwinding"])),
        tracer_upwinding = Val(Symbol(params["tracer_upwinding"])),
        density_upwinding = Val(Symbol(params["density_upwinding"])),
        edmfx_upwinding = Val(Symbol(params["edmfx_upwinding"])),
        apply_limiter = params["apply_limiter"],
        bubble = params["bubble"],
    )
    @info "numerics" numerics...

    return numerics
end

function get_spaces(params, ca_phys_params, comms_ctx)

    FT = eltype(ca_phys_params)
    z_elem = Int(params["z_elem"])
    z_max = FT(params["z_max"])
    dz_bottom = FT(params["dz_bottom"])
    dz_top = FT(params["dz_top"])
    topography = params["topography"]
    bubble = params["bubble"]

    @assert topography in ("NoWarp", "DCMIP200", "Earth", "Agnesi", "Schar")
    if topography == "DCMIP200"
        warp_function = topography_dcmip200
    elseif topography == "Agnesi"
        warp_function = topography_agnesi
    elseif topography == "Schar"
        warp_function = topography_schar
    elseif topography == "NoWarp"
        warp_function = nothing
    elseif topography == "Earth"
        data_path = joinpath(topo_elev_dataset_path(), "ETOPO1_coarse.nc")
        earth_spline = NCDataset(data_path) do data
            zlevels = data["elevation"][:]
            lon = data["longitude"][:]
            lat = data["latitude"][:]
            # Apply Smoothing
            smooth_degree = 15
            esmth = imfilter(zlevels, Kernel.gaussian(smooth_degree))
            linear_interpolation(
                (lon, lat),
                esmth,
                extrapolation_bc = (Periodic(), Flat()),
            )
        end
        @info "Generated interpolation stencil"
        warp_function = generate_topography_warp(earth_spline)
    end
    @info "Topography" topography


    h_elem = params["h_elem"]
    radius = CAP.planet_radius(ca_phys_params)
    center_space, face_space = if params["config"] == "sphere"
        nh_poly = params["nh_poly"]
        quad = Spaces.Quadratures.GLL{nh_poly + 1}()
        horizontal_mesh = cubed_sphere_mesh(; radius, h_elem)
        h_space =
            make_horizontal_space(horizontal_mesh, quad, comms_ctx, bubble)
        z_stretch = if params["z_stretch"]
            Meshes.GeneralizedExponentialStretching(dz_bottom, dz_top)
        else
            Meshes.Uniform()
        end
        if warp_function == nothing
            make_hybrid_spaces(h_space, z_max, z_elem, z_stretch)
        else
            make_hybrid_spaces(
                h_space,
                z_max,
                z_elem,
                z_stretch;
                surface_warp = warp_function,
            )
        end
    elseif params["config"] == "column" # single column
        @warn "perturb_initstate flag is ignored for single column configuration"
        FT = eltype(ca_phys_params)
        Δx = FT(1) # Note: This value shouldn't matter, since we only have 1 column.
        quad = Spaces.Quadratures.GL{1}()
        horizontal_mesh = periodic_rectangle_mesh(;
            x_max = Δx,
            y_max = Δx,
            x_elem = 1,
            y_elem = 1,
        )
        if bubble
            @warn "Bubble correction not compatible with single column configuration. It will be switched off."
            bubble = false
        end
        h_space =
            make_horizontal_space(horizontal_mesh, quad, comms_ctx, bubble)
        z_stretch = if params["z_stretch"]
            Meshes.GeneralizedExponentialStretching(dz_bottom, dz_top)
        else
            Meshes.Uniform()
        end
        make_hybrid_spaces(h_space, z_max, z_elem, z_stretch)
    elseif params["config"] == "box"
        FT = eltype(ca_phys_params)
        nh_poly = params["nh_poly"]
        quad = Spaces.Quadratures.GLL{nh_poly + 1}()
        x_elem = Int(params["x_elem"])
        x_max = FT(params["x_max"])
        y_elem = Int(params["y_elem"])
        y_max = FT(params["y_max"])
        horizontal_mesh = periodic_rectangle_mesh(;
            x_max = x_max,
            y_max = y_max,
            x_elem = x_elem,
            y_elem = y_elem,
        )
        h_space =
            make_horizontal_space(horizontal_mesh, quad, comms_ctx, bubble)
        z_stretch = if params["z_stretch"]
            Meshes.GeneralizedExponentialStretching(dz_bottom, dz_top)
        else
            Meshes.Uniform()
        end
        make_hybrid_spaces(
            h_space,
            z_max,
            z_elem,
            z_stretch;
            surface_warp = warp_function,
        )
    elseif params["config"] == "plane"
        FT = eltype(ca_phys_params)
        nh_poly = params["nh_poly"]
        quad = Spaces.Quadratures.GLL{nh_poly + 1}()
        x_elem = Int(params["x_elem"])
        x_max = FT(params["x_max"])
        horizontal_mesh =
            periodic_line_mesh(; x_max = x_max, x_elem = x_elem)
        h_space =
            make_horizontal_space(horizontal_mesh, quad, comms_ctx, bubble)
        z_stretch = if params["z_stretch"]
            Meshes.GeneralizedExponentialStretching(dz_bottom, dz_top)
        else
            Meshes.Uniform()
        end
        make_hybrid_spaces(
            h_space,
            z_max,
            z_elem,
            z_stretch;
            surface_warp = warp_function,
        )
    end
    return (;
        center_space,
        face_space,
        horizontal_mesh,
        quad,
        z_max,
        z_elem,
        z_stretch,
    )
end

# get_state(simulation, params, spaces, ca_phys_params, atmos)
function get_state(simulation, args...)
    if simulation.restart
        return get_state_restart(comms_ctx)
    else
        return get_state_fresh_start(args...)
    end
end

function get_spaces_restart(Y)
    center_space = axes(Y.c)
    face_space = axes(Y.f)
    hspace = Spaces.horizontal_space(center_space)
    horizontal_mesh = hspace.topology.mesh
    quad = horizontal_mesh.ne + 1
    vertical_mesh = Spaces.vertical_topology(center_space).mesh
    z_max = vertical_mesh.domain.coord_max.z
    z_elem = length(vertical_mesh.faces) - 1
    return (; center_space, face_space, horizontal_mesh, quad, z_max, z_elem)
end

function get_state_restart(comms_ctx)
    @assert haskey(ENV, "RESTART_FILE")
    reader = InputOutput.HDF5Reader(ENV["RESTART_FILE"], comms_ctx)
    Y = InputOutput.read_field(reader, "Y")
    t_start = InputOutput.HDF5.read_attribute(reader.file, "time")
    return (Y, t_start)
end

function get_initial_condition(params)
    if isnothing(params["turbconv_case"])
        if params["initial_condition"] in [
            "DryBaroclinicWave",
            "MoistBaroclinicWave",
            "DecayingProfile",
            "MoistBaroclinicWaveWithEDMF",
            "DryAdiabaticProfileEDMFX",
        ]
            return getproperty(ICs, Symbol(params["initial_condition"]))(
                params["perturb_initstate"],
            )
        elseif params["initial_condition"] in [
            "IsothermalProfile",
            "Bomex",
            "AgnesiHProfile",
            "DryDensityCurrentProfile",
            "ScharProfile",
        ]
            return getproperty(ICs, Symbol(params["initial_condition"]))()
        else
            error("Unknown `initial_condition`: $(params["initial_condition"])")
        end
    else
        # turbconv_case is also used for surface fluxes for TRMM and ARM cases.
        # I don't want to change that right now, so I'm leaving the
        # EDMF logic as is. This should be obsolete soon.
        return getproperty(ICs, Symbol(params["turbconv_case"]))()
    end
end

is_explicit_CTS_algo_type(alg_or_tableau) =
    alg_or_tableau <: CTS.ERKAlgorithmName

is_imex_CTS_algo_type(alg_or_tableau) =
    alg_or_tableau <: CTS.IMEXARKAlgorithmName

is_implicit_type(::typeof(ODE.IMEXEuler)) = true
is_implicit_type(alg_or_tableau) =
    alg_or_tableau <: Union{
        ODE.OrdinaryDiffEqImplicitAlgorithm,
        ODE.OrdinaryDiffEqAdaptiveImplicitAlgorithm,
    } || is_imex_CTS_algo_type(alg_or_tableau)

is_ordinary_diffeq_newton(::typeof(ODE.IMEXEuler)) = true
is_ordinary_diffeq_newton(alg_or_tableau) =
    alg_or_tableau <: Union{
        ODE.OrdinaryDiffEqNewtonAlgorithm,
        ODE.OrdinaryDiffEqNewtonAdaptiveAlgorithm,
    }

is_imex_CTS_algo(::CTS.IMEXAlgorithm) = true
is_imex_CTS_algo(::DiffEqBase.AbstractODEAlgorithm) = false

is_implicit(::ODE.OrdinaryDiffEqImplicitAlgorithm) = true
is_implicit(::ODE.OrdinaryDiffEqAdaptiveImplicitAlgorithm) = true
is_implicit(ode_algo) = is_imex_CTS_algo(ode_algo)

is_rosenbrock(::ODE.Rosenbrock23) = true
is_rosenbrock(::ODE.Rosenbrock32) = true
is_rosenbrock(::DiffEqBase.AbstractODEAlgorithm) = false
use_transform(ode_algo) =
    !(is_imex_CTS_algo(ode_algo) || is_rosenbrock(ode_algo))

additional_integrator_kwargs(::DiffEqBase.AbstractODEAlgorithm) = (;
    adaptive = false,
    progress = isinteractive(),
    progress_steps = isinteractive() ? 1 : 1000,
)
additional_integrator_kwargs(::CTS.DistributedODEAlgorithm) = (;
    kwargshandle = ODE.KeywordArgSilent, # allow custom kwargs
    adjustfinal = true,
    # TODO: enable progress bars in ClimaTimeSteppers
)

is_cts_algo(::DiffEqBase.AbstractODEAlgorithm) = false
is_cts_algo(::CTS.DistributedODEAlgorithm) = true

jacobi_flags(::TotalEnergy) =
    (; ∂ᶜ𝔼ₜ∂ᶠ𝕄_mode = :no_∂ᶜp∂ᶜK, ∂ᶠ𝕄ₜ∂ᶜρ_mode = :exact)
jacobi_flags(::PotentialTemperature) =
    (; ∂ᶜ𝔼ₜ∂ᶠ𝕄_mode = :exact, ∂ᶠ𝕄ₜ∂ᶜρ_mode = :exact)

function jac_kwargs(ode_algo, Y, energy_form)
    if is_implicit(ode_algo)
        W = SchurComplementW(
            Y,
            use_transform(ode_algo),
            jacobi_flags(energy_form),
        )
        if use_transform(ode_algo)
            return (; jac_prototype = W, Wfact_t = Wfact!)
        else
            return (; jac_prototype = W, Wfact = Wfact!)
        end
    else
        return NamedTuple()
    end
end

#=
    ode_configuration(Y, params, atmos)

Returns the ode algorithm
=#
function ode_configuration(Y, params, atmos)
    FT = Spaces.undertype(axes(Y.c))
    ode_name = params["ode_algo"]
    alg_or_tableau = if startswith(ode_name, "ODE.")
        @warn "apply_limiter flag is ignored for OrdinaryDiffEq algorithms"
        getproperty(ODE, Symbol(split(ode_name, ".")[2]))
    else
        getproperty(CTS, Symbol(ode_name))
    end
    @info "Using ODE config: `$alg_or_tableau`"

    if is_explicit_CTS_algo_type(alg_or_tableau)
        return CTS.ExplicitAlgorithm(alg_or_tableau())
    elseif !is_implicit_type(alg_or_tableau)
        return alg_or_tableau()
    elseif is_ordinary_diffeq_newton(alg_or_tableau)
        if params["max_newton_iters"] == 1
            error("OridinaryDiffEq requires at least 2 Newton iterations")
        end
        # κ like a relative tolerance; its default value in ODE is 0.01
        nlsolve = ODE.NLNewton(;
            κ = params["max_newton_iters"] == 2 ? Inf : 0.01,
            max_iter = params["max_newton_iters"],
        )
        return alg_or_tableau(; linsolve = linsolve!, nlsolve)
    elseif is_imex_CTS_algo_type(alg_or_tableau)
        newtons_method = CTS.NewtonsMethod(;
            max_iters = params["max_newton_iters"],
            krylov_method = if params["use_krylov_method"]
                CTS.KrylovMethod(;
                    jacobian_free_jvp = CTS.ForwardDiffJVP(;
                        step_adjustment = FT(params["jvp_step_adjustment"]),
                    ),
                    forcing_term = if params["use_dynamic_krylov_rtol"]
                        α = FT(params["eisenstat_walker_forcing_alpha"])
                        CTS.EisenstatWalkerForcing(; α)
                    else
                        CTS.ConstantForcing(FT(params["krylov_rtol"]))
                    end,
                )
            else
                nothing
            end,
            convergence_checker = if params["use_newton_rtol"]
                norm_condition =
                    CTS.MaximumRelativeError(FT(params["newton_rtol"]))
                CTS.ConvergenceChecker(; norm_condition)
            else
                nothing
            end,
        )
        return CTS.IMEXAlgorithm(alg_or_tableau(), newtons_method)
    else
        return alg_or_tableau(; linsolve = linsolve!)
    end
end

function get_integrator(args, kwargs)
    @time "Define integrator" integrator = ODE.init(args...; kwargs...)
    return integrator
end

thermo_state_type(::DryModel, ::Type{FT}) where {FT} = TD.PhaseDry{FT}
thermo_state_type(::EquilMoistModel, ::Type{FT}) where {FT} = TD.PhaseEquil{FT}
thermo_state_type(::NonEquilMoistModel, ::Type{FT}) where {FT} =
    TD.PhaseNonEquil{FT}


function get_callbacks(params, simulation, atmos, ca_phys_params)
    FT = eltype(ca_phys_params)
    (; dt) = simulation

    tc_callbacks =
        call_every_n_steps(turb_conv_affect_filter!; skip_first = true)
    flux_accumulation_callback = call_every_n_steps(
        flux_accumulation!;
        skip_first = true,
        call_at_end = true,
    )

    additional_callbacks =
        if atmos.radiation_mode isa RRTMGPI.AbstractRRTMGPMode
            # TODO: better if-else criteria?
            dt_rad = if params["config"] == "column"
                dt
            else
                FT(time_to_seconds(params["dt_rad"]))
            end
            (call_every_dt(rrtmgp_model_callback!, dt_rad),)
        else
            ()
        end

    if atmos.turbconv_model isa TC.EDMFModel
        additional_callbacks = (additional_callbacks..., tc_callbacks)
    end

    if params["check_conservation"]
        additional_callbacks =
            (flux_accumulation_callback, additional_callbacks...)
    end

    dt_save_to_disk = time_to_seconds(params["dt_save_to_disk"])
    dt_save_restart = time_to_seconds(params["dt_save_restart"])

    dss_cb = if startswith(params["ode_algo"], "ODE.")
        call_every_n_steps(dss_callback)
    else
        nothing
    end
    save_to_disk_callback = if dt_save_to_disk == Inf
        nothing
    elseif simulation.restart
        call_every_dt(save_to_disk_func, dt_save_to_disk; skip_first = true)
    else
        call_every_dt(save_to_disk_func, dt_save_to_disk)
    end

    save_restart_callback = if dt_save_restart == Inf
        nothing
    else
        call_every_dt(save_restart_func, dt_save_restart)
    end

    gc_callback = if CA.is_distributed(simulation.comms_ctx)
        call_every_n_steps(
            gc_func,
            parse(Int, get(ENV, "CLIMAATMOS_GC_NSTEPS", "1000")),
            skip_first = true,
        )
    else
        nothing
    end

    return ODE.CallbackSet(
        dss_cb,
        save_to_disk_callback,
        save_restart_callback,
        gc_callback,
        additional_callbacks...,
    )
end


function get_cache(
    Y,
    params,
    ca_phys_params,
    spaces,
    atmos,
    numerics,
    simulation,
    initial_condition,
)
    _default_cache = default_cache(
        Y,
        params,
        ca_phys_params,
        atmos,
        spaces,
        numerics,
        simulation,
    )
    merge(
        _default_cache,
        additional_cache(
            Y,
            _default_cache,
            params,
            ca_phys_params,
            atmos,
            simulation.dt,
            initial_condition,
        ),
    )
end

function get_simulation(::Type{FT}, params, comms_ctx) where {FT}

    job_id = if isnothing(params["job_id"])
        s = argparse_settings()
        job_id_from_parsed_args(s, params)
    else
        params["job_id"]
    end
    default_output = haskey(ENV, "CI") ? job_id : joinpath("output", job_id)
    out_dir = params["output_dir"]
    output_dir = isnothing(out_dir) ? default_output : out_dir
    mkpath(output_dir)

    sim = (;
        comms_ctx,
        is_debugging_tc = params["debugging_tc"],
        output_dir,
        restart = haskey(ENV, "RESTART_FILE"),
        job_id,
        dt = FT(time_to_seconds(params["dt"])),
        start_date = DateTime(params["start_date"], dateformat"yyyymmdd"),
        t_end = FT(time_to_seconds(params["t_end"])),
    )
    n_steps = floor(Int, sim.t_end / sim.dt)
    @info(
        "Time info:",
        dt = params["dt"],
        t_end = params["t_end"],
        floor_n_steps = n_steps,
    )

    return sim
end

function args_integrator(params, Y, p, tspan, ode_algo, callback)
    (; atmos, simulation) = p
    (; dt) = simulation
    dt_save_to_sol = time_to_seconds(params["dt_save_to_sol"])

    @time "Define ode function" func = if params["split_ode"]
        implicit_func = ODE.ODEFunction(
            implicit_tendency!;
            jac_kwargs(ode_algo, Y, atmos.energy_form)...,
            tgrad = (∂Y∂t, Y, p, t) -> (∂Y∂t .= 0),
        )
        if is_cts_algo(ode_algo)
            CTS.ClimaODEFunction(;
                T_lim! = limited_tendency!,
                T_exp! = remaining_tendency!,
                T_imp! = implicit_func,
                # Can we just pass implicit_tendency! and jac_prototype etc.?
                lim! = limiters_func!,
                dss!,
            )
        else
            ODE.SplitFunction(implicit_func, remaining_tendency!)
        end
    else
        remaining_tendency! # should be total_tendency!
    end
    problem = ODE.ODEProblem(func, Y, tspan, p)
    saveat = if dt_save_to_sol == Inf
        tspan[2]
    elseif tspan[2] % dt_save_to_sol == 0
        dt_save_to_sol
    else
        [tspan[1]:dt_save_to_sol:tspan[2]..., tspan[2]]
    end # ensure that tspan[2] is always saved
    args = (problem, ode_algo)
    kwargs = (; saveat, callback, dt, additional_integrator_kwargs(ode_algo)...)
    return (args, kwargs)
end

import ClimaComms, Logging, NVTX
function get_comms_context(device = ClimaComms.CPUDevice())
    comms_ctx = ClimaComms.context(device)
    ClimaComms.init(comms_ctx)
    if ClimaComms.iamroot(comms_ctx)
        Logging.global_logger(Logging.ConsoleLogger(stderr, Logging.Info))
    else
        Logging.global_logger(Logging.NullLogger())
    end
    if comms_ctx isa ClimaComms.SingletonCommsContext
        @info "Setting up single-process ClimaAtmos run"
    else
        @info "Setting up distributed ClimaAtmos run" nprocs =
            ClimaComms.nprocs(comms_ctx)
    end
    if NVTX.isactive()
        # makes output on buildkite a bit nicer
        if ClimaComms.iamroot(comms_ctx)
            atexit() do
                println("--- Saving profiler information")
            end
        end
    end

    return comms_ctx
end
