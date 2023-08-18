import DiffEqCallbacks
#####
##### Callback helpers
#####

function call_every_n_steps(f!, n = 1; skip_first = false, call_at_end = false)
    previous_step = Ref(0)
    @assert n ≠ Inf "Adding callback that never gets called!"
    cb! = AtmosCallback(f!, EveryNSteps(n))
    return ODE.DiscreteCallback(
        (u, t, integrator) ->
            (previous_step[] += 1) % n == 0 ||
                (call_at_end && t == integrator.sol.prob.tspan[2]),
        cb!;
        initialize = (cb, u, t, integrator) -> skip_first || cb!(integrator),
        save_positions = (false, false),
    )
end

function call_every_dt(
    f!,
    dt;
    skip_first = false,
    call_at_end = false,
    initialize = nothing,
)
    cb! = AtmosCallback(f!, EveryΔt(dt))
    @assert dt ≠ Inf "Adding callback that never gets called!"
    next_t = Ref{typeof(dt)}()
    affect! = function (integrator)
        cb!(integrator)

        t = integrator.t
        t_end = integrator.sol.prob.tspan[2]
        next_t[] = max(t, next_t[] + dt)
        if call_at_end
            next_t[] = min(next_t[], t_end)
        end
    end
    return ODE.DiscreteCallback(
        (u, t, integrator) -> t >= next_t[],
        affect!;
        initialize = if isnothing(initialize)
            (cb, u, t, integrator) -> begin
                skip_first || cb!(integrator)
                t_end = integrator.sol.prob.tspan[2]
                next_t[] =
                    (call_at_end && t < t_end) ? min(t_end, t + dt) : t + dt
            end
        else
            initialize
        end,
        save_positions = (false, false),
    )
end

callback_from_affect(x::AtmosCallback) = x
function callback_from_affect(affect!)
    for p in propertynames(affect!)
        x = getproperty(affect!, p)
        if x isa AtmosCallback
            return x
        elseif x isa DiffEqCallbacks.SavedValues
            return x
        end
    end
    error("Callback not found in $(affect!)")
end
function atmos_callbacks(cbs::ODE.CallbackSet)
    all_cbs = [cbs.continuous_callbacks..., cbs.discrete_callbacks...]
    callback_objs = map(cb -> callback_from_affect(cb.affect!), all_cbs)
    filter!(x -> !(x isa DiffEqCallbacks.SavedValues), callback_objs)
    return callback_objs
end

n_measured_calls(integrator) = n_measured_calls(integrator.callback)
n_measured_calls(cbs::ODE.CallbackSet) =
    map(x -> x.n_measured_calls, atmos_callbacks(cbs))

n_expected_calls(integrator) = n_expected_calls(
    integrator.callback,
    integrator.dt,
    integrator.sol.prob.tspan,
)
n_expected_calls(cbs::ODE.CallbackSet, dt, tspan) =
    map(x -> n_expected_calls(x, dt, tspan), atmos_callbacks(cbs))

n_steps_per_cycle(integrator) =
    n_steps_per_cycle(integrator.callback, integrator.dt)
function n_steps_per_cycle(cbs::ODE.CallbackSet, dt)
    nspc = n_steps_per_cycle_per_cb(cbs, dt)
    return isempty(nspc) ? 1 : lcm(nspc)
end

n_steps_per_cycle_per_cb(integrator) =
    n_steps_per_cycle_per_cb(integrator.callback, integrator.dt)

function n_steps_per_cycle_per_cb(cbs::ODE.CallbackSet, dt)
    return map(atmos_callbacks(cbs)) do cb
        cbf = callback_frequency(cb)
        if cbf isa EveryΔt
            Int(ceil(cbf.Δt / dt))
        elseif cbf isa EveryNSteps
            cbf.n
        else
            error("Uncaught case")
        end
    end
end
