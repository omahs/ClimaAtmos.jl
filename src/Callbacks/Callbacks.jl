module Callbacks

using DiffEqCallbacks
using UnPack
using JLD2

export cfl_cb,
       generate_callback,
       AbstractCallback,
       CFLCallback,
       JLD2Callback

abstract type AbstractCallback end
function generate_callback(::AbstractCallback)
    return nothing
end

###
### CFLCallback 
###
struct CFLCallback <: AbstractCallback 
    interval::Number
end
function f_test!(integrator)
    var = integrator.u.swm
    u = getproperty(integrator.u.swm, :u)
    # Default callback simulations out of scope
    # # Default callback simulations out of scope
    u = var.u 
    c = var.c
    h = var.h
    t = integrator.t
    return nothing
end
function generate_callback(c::CFLCallback; kwargs...)
    return PeriodicCallback(f_test!,c.interval; initial_affect=false, kwargs...)
end

###
### JLD2 Callback
### 
struct JLD2Callback <: AbstractCallback 
    filename :: String
    interval :: Number
end
function g_test!(integrator)
    mkpath("./TestOutput/")
    jldsave(joinpath("./TestOutput/", "Test.jl"); 
            h = integrator.u.swm.h,
            c = integrator.u.swm.c,
            u = integrator.u.swm.u,
            t = integrator.t)
    return nothing
end
function write_output!(::JLD2Callback)
    return g_test!
end
function generate_callback(c::JLD2Callback; kwargs...)
    return PeriodicCallback(write_output!(JLD2Callback(c.filename,1)), 
                            c.interval; 
                            initial_affect=false, 
                            kwargs...)
end

end
