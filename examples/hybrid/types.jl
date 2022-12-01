using Dates: DateTime, @dateformat_str
using NCDatasets
using Dierckx
import ClimaCore: InputOutput
import ClimaAtmos.RRTMGPInterface as RRTMGPI
import ClimaAtmos as CA
import ClimaAtmos:
    DryModel,
    EquilMoistModel,
    NonEquilMoistModel,
    CompressibleFluid,
    AnelasticFluid,
    PotentialTemperature,
    TotalEnergy,
    InternalEnergy,
    Microphysics0Moment,
    HeldSuarezForcing,
    BulkSurfaceScheme,
    MoninObukhovSurface,
    SingleColumnModel,
    SphericalModel,
    BoxModel

import ClimaCore: InputOutput

function get_atmos(::Type{FT}, parsed_args, namelist) where {FT}
    # should this live in the radiation model?
    non_orographic_gravity_wave = parsed_args["non_orographic_gravity_wave"]
    @assert non_orographic_gravity_wave in (true, false)

    moisture_model = CA.moisture_model(parsed_args)
    precip_model = CA.precipitation_model(parsed_args)
    radiation_mode = CA.radiation_mode(parsed_args, FT)

    atmos = CA.AtmosModel(;
        moisture_model,
        model_config = CA.model_config(parsed_args),
        coupling = CA.coupling_type(parsed_args),
        perf_mode = CA.perf_mode(parsed_args),
        energy_form = CA.energy_form(parsed_args),
        radiation_mode,
        subsidence = CA.subsidence_model(parsed_args, radiation_mode, FT),
        ls_adv = CA.large_scale_advection_model(parsed_args, FT),
        edmf_coriolis = CA.edmf_coriolis(parsed_args, FT),
        precip_model,
        forcing_type = CA.forcing_type(parsed_args),
        turbconv_model = CA.turbconv_model(
            FT,
            moisture_model,
            precip_model,
            parsed_args,
            namelist,
        ),
        compressibility_model = CA.compressibility_model(parsed_args),
        surface_scheme = CA.surface_scheme(FT, parsed_args),
        non_orographic_gravity_wave,
    )

    return atmos
end

function get_numerics(parsed_args)
    # wrap each upwinding mode in a Val for dispatch
    numerics = (;
        energy_upwinding = Val(Symbol(parsed_args["energy_upwinding"])),
        tracer_upwinding = Val(Symbol(parsed_args["tracer_upwinding"])),
        apply_limiter = parsed_args["apply_limiter"],
    )
    @info "numerics" numerics...

    return numerics
end

function get_simulation(::Type{FT}, parsed_args) where {FT}

    job_id = if isnothing(parsed_args["job_id"])
        (s, default_parsed_args) = parse_commandline()
        job_id_from_parsed_args(s, parsed_args)
    else
        parsed_args["job_id"]
    end
    default_output = haskey(ENV, "CI") ? job_id : joinpath("output", job_id)
    output_dir = parse_arg(parsed_args, "output_dir", default_output)
    mkpath(output_dir)

    sim = (;
        is_distributed = haskey(ENV, "CLIMACORE_DISTRIBUTED"),
        is_debugging_tc = parsed_args["debugging_tc"],
        output_dir,
        restart = haskey(ENV, "RESTART_FILE"),
        job_id,
        dt = FT(time_to_seconds(parsed_args["dt"])),
        start_date = DateTime(parsed_args["start_date"], dateformat"yyyymmdd"),
        t_end = FT(time_to_seconds(parsed_args["t_end"])),
    )
    n_steps = floor(Int, sim.t_end / sim.dt)
    @info(
        "Time info:",
        dt = parsed_args["dt"],
        t_end = parsed_args["t_end"],
        floor_n_steps = n_steps,
    )

    return sim
end

function get_spaces(parsed_args, params, comms_ctx)

    FT = eltype(params)
    z_elem = Int(parsed_args["z_elem"])
    z_max = FT(parsed_args["z_max"])
    dz_bottom = FT(parsed_args["dz_bottom"])
    dz_top = FT(parsed_args["dz_top"])
    topography = parsed_args["topography"]

    if topography == "DCMIP200"
        warp_function = CA.topography_dcmip200
    elseif topography == "NoWarp"
        warp_function = nothing
    end
    @assert topography in ("NoWarp", "DCMIP200")
    @info "Topography" topography

    h_elem = parsed_args["h_elem"]
    radius = CAP.planet_radius(params)
    center_space, face_space = if parsed_args["config"] == "sphere"
        nh_poly = parsed_args["nh_poly"]
        quad = Spaces.Quadratures.GLL{nh_poly + 1}()
        horizontal_mesh = cubed_sphere_mesh(; radius, h_elem)
        h_space = make_horizontal_space(horizontal_mesh, quad, comms_ctx)
        z_stretch = if parsed_args["z_stretch"]
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
    elseif parsed_args["config"] == "column" # single column
        @warn "perturb_initstate flag is ignored for single column configuration"
        FT = eltype(params)
        Δx = FT(1) # Note: This value shouldn't matter, since we only have 1 column.
        quad = Spaces.Quadratures.GL{1}()
        horizontal_mesh = periodic_rectangle_mesh(;
            x_max = Δx,
            y_max = Δx,
            x_elem = 1,
            y_elem = 1,
        )
        h_space = make_horizontal_space(horizontal_mesh, quad, comms_ctx)
        z_stretch = if parsed_args["z_stretch"]
            Meshes.GeneralizedExponentialStretching(dz_bottom, dz_top)
        else
            Meshes.Uniform()
        end
        make_hybrid_spaces(h_space, z_max, z_elem, z_stretch)
    elseif parsed_args["config"] == "box"
        FT = eltype(params)
        nh_poly = parsed_args["nh_poly"]
        quad = Spaces.Quadratures.GLL{nh_poly + 1}()
        x_elem = Int(parsed_args["x_elem"])
        x_max = FT(parsed_args["x_max"])
        y_elem = Int(parsed_args["y_elem"])
        y_max = FT(parsed_args["y_max"])
        horizontal_mesh = periodic_rectangle_mesh(;
            x_max = x_max,
            y_max = y_max,
            x_elem = x_elem,
            y_elem = y_elem,
        )
        h_space = make_horizontal_space(horizontal_mesh, quad, comms_ctx)
        z_stretch = if parsed_args["z_stretch"]
            Meshes.GeneralizedExponentialStretching(dz_bottom, dz_top)
        else
            Meshes.Uniform()
        end
        make_hybrid_spaces(h_space, z_max, z_elem, z_stretch)
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

# get_state(simulation, parsed_args, spaces, params, atmos)
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

import ClimaAtmos.InitialConditions as ICs

function get_state_fresh_start(parsed_args, spaces, params, atmos)
    (; center_space, face_space) = spaces
    FT = eltype(params)
    t_start = FT(0)

    center_initial_condition = if is_baro_wave(parsed_args)
        ICs.center_initial_condition_baroclinic_wave
    elseif parsed_args["config"] == "sphere"
        ICs.center_initial_condition_3d
    elseif parsed_args["config"] == "column"
        ICs.center_initial_condition_column
    elseif parsed_args["config"] == "box"
        ICs.center_initial_condition_box
    end
    perturb_initstate = parsed_args["perturb_initstate"]

    Y = ICs.init_state(
        center_initial_condition,
        ICs.face_initial_condition,
        center_space,
        face_space,
        params,
        atmos,
        perturb_initstate,
    )
    return (Y, t_start)
end

is_implicit(::Type{T<:ODE.OrdinaryDiffEqImplicitAlgorithm}) where {T} = true
is_implicit(::Type{T<:ODE.OrdinaryDiffEqAdaptiveImplicitAlgorithm}) where {T} = true
is_implicit(::Type{T<:Any}) where {T} = is_imex_CTS_algo(T)
is_implicit(x::Function) = is_implicit(typeof(x()))
is_implicit(x) = is_implicit(typeof(x))

is_imex_CTS_algo(::Type{<:Any}) = false
is_imex_CTS_algo(::Type{<:CTS.IMEXARKAlgorithm}) = true
is_imex_CTS_algo(x) = is_imex_CTS_algo(typeof(x))
is_imex_CTS_algo(x::Function) = is_imex_CTS_algo(typeof(x()))
is_rosenbrock(::Type{<:ODE.Rosenbrock23}) = true
is_rosenbrock(::Type{<:ODE.Rosenbrock32}) = false
is_rosenbrock(::Type{<:Any}) = false
is_rosenbrock(x::Function) = is_rosenbrock(typeof(x()))
is_rosenbrock(x) = is_rosenbrock(typeof(x))

is_ordinary_diffeq_newton(x::Function) = is_ordinary_diffeq_newton(typeof(x()))
is_ordinary_diffeq_newton(x) = is_ordinary_diffeq_newton(typeof(x))
is_ordinary_diffeq_newton(::Type{T}) where {T} =
    T <: Union{
        ODE.OrdinaryDiffEqNewtonAlgorithm,
        ODE.OrdinaryDiffEqNewtonAdaptiveAlgorithm,
    }

is_CTS_algo(::Type{<:CTS.DistributedODEAlgorithm}) = true
is_CTS_algo(::Type{<:Any}) = false
is_CTS_algo(x::Function) = is_CTS_algo(typeof(x()))
is_CTS_algo(x) = is_CTS_algo(typeof(x))
is_ordinary_diffeq_algo(x) = !is_CTS_algo(x)

use_transform(ode_algo_type) =
    !(is_imex_CTS_algo(ode_algo_type) || is_rosenbrock(ode_algo_type))

function jac_kwargs(ode_algorithm, Y, energy_form)
    if is_implicit(ode_algorithm)
        W = CA.SchurComplementW(
            Y,
            use_transform(ode_algorithm),
            jacobi_flags(energy_form),
        )
        if use_transform(ode_algorithm)
            return (; jac_prototype = W, Wfact_t = CA.Wfact!)
        else
            return (; jac_prototype = W, Wfact = CA.Wfact!)
        end
    else
        return NamedTuple()
    end
end

import OrdinaryDiffEq as ODE
import ClimaTimeSteppers as CTS
#=
(; ode_algo_type, alg_kwargs) =
    ode_config(Y, parsed_args, atmos)
=#
function ode_configuration(Y, parsed_args, atmos)
    _is_ordinary_diffeq_algo = startswith(parsed_args["ode_algo"], "ODE.")
    ode_algo_type = if _is_ordinary_diffeq_algo
        @warn "apply_limiter flag is ignored for OrdinaryDiffEq algorithms"
        getproperty(ODE, Symbol(split(parsed_args["ode_algo"], ".")[2]))
    else
        getproperty(CTS, Symbol(parsed_args["ode_algo"]))
    end
    @assert is_CTS_algo(ode_algo_type) ≠ _is_ordinary_diffeq_algo

    if !is_implicit(ode_algo_type)
        return ode_algo_type()
    elseif is_ordinary_diffeq_newton(ode_algo_type)
        if parsed_args["max_newton_iters"] == 1
            error("OridinaryDiffEq requires at least 2 Newton iterations")
        end
        # κ like a relative tolerance; its default value in ODE is 0.01
        nlsolve = ODE.NLNewton(;
            κ = parsed_args["max_newton_iters"] == 2 ? Inf : 0.01,
            max_iter = parsed_args["max_newton_iters"],
        )
        return ode_algo_type(;linsolve = CA.linsolve!, nlsolve)
    elseif is_imex_CTS_algo(ode_algo_type)
        newtons_method = NewtonsMethod(;
            max_iters = parsed_args["max_newton_iters"],
            krylov_method = if parsed_args["use_krylov_method"]
                KrylovMethod(;
                    jacobian_free_jvp = ForwardDiffJVP(;
                        step_adjustment = FT(
                            parsed_args["jvp_step_adjustment"],
                        ),
                    ),
                    forcing_term = if parsed_args["use_dynamic_krylov_rtol"]
                        α = FT(parsed_args["eisenstat_walker_forcing_alpha"])
                        EisenstatWalkerForcing(; α)
                    else
                        ConstantForcing(FT(parsed_args["krylov_rtol"]))
                    end,
                )
            else
                nothing
            end,
            convergence_checker = if parsed_args["use_newton_rtol"]
                norm_condition =
                    MaximumRelativeError(FT(parsed_args["newton_rtol"]))
                ConvergenceChecker(; norm_condition)
            else
                nothing
            end,
        )
        return ode_algo_type(;newtons_method)
    else
        return ode_algo_type(;linsolve = CA.linsolve!)
    end
end

function args_integrator(parsed_args, Y, p, tspan, ode_config, callback)
    (; alg_kwargs, ode_algo) = ode_config
    (; dt) = p.simulation
    FT = eltype(tspan)
    dt_save_to_sol = time_to_seconds(parsed_args["dt_save_to_sol"])
    show_progress_bar = isinteractive()

    @time "Define problem" problem = if parsed_args["split_ode"]
        remaining_func =
            is_ordinary_diffeq_algo(ode_algo) ? ForwardEulerODEFunction(remaining_tendency_increment!) :
            remaining_tendency!
        ODE.SplitODEProblem(
            ODE.ODEFunction(
                implicit_tendency!;
                jac_kwargs(ode_algo, Y, atmos.energy_form)...,
                tgrad = (∂Y∂t, Y, p, t) -> (∂Y∂t .= FT(0)),
            ),
            remaining_func,
            Y,
            tspan,
            p,
        )
    else
        ODE.ODEProblem(remaining_tendency!, Y, tspan, p)
    end
    if startswith(parsed_args["ode_algo"], "ODE.")
        ode_algo = ode_algo_type(; alg_kwargs...)
        integrator_kwargs = (;
            adaptive = false,
            progress = show_progress_bar,
            progress_steps = isinteractive() ? 1 : 1000,
        )
    else
        ode_algo = ode_algo_type(alg_kwargs...)
        integrator_kwargs = (;
            kwargshandle = KeywordArgSilent, # allow custom kwargs
            adjustfinal = true,
            # TODO: enable progress bars in ClimaTimeSteppers
        )
    end
    saveat = if dt_save_to_sol == Inf
        tspan[2]
    elseif tspan[2] % dt_save_to_sol == 0
        dt_save_to_sol
    else
        [tspan[1]:dt_save_to_sol:tspan[2]..., tspan[2]]
    end # ensure that tspan[2] is always saved
    args = (problem, ode_algo)
    kwargs = (; saveat, callback, dt, integrator_kwargs...)
    return (args, kwargs)
end

function get_integrator(args, kwargs)
    @time "Define integrator" integrator = ODE.init(args...; kwargs...)
    return integrator
end
