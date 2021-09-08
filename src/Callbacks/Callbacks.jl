module Callbacks
using DiffEqCallbacks
using UnPack
using JLD2

export cfl_cb,
       AbstractCallback,
       CFLCallback,
       generate_callback,
       JLD2Callback

abstract type AbstractCallback end


###
### CFLCallback 
###
struct CFLCallback <: AbstractCallback 
    interval::Number
end
function generate_callback(::AbstractCallback)
    return nothing
end
function f_test!(integrator)
    var = integrator.sol.u
    #TODO Unpack from known `Model` variables
    #Requires generalisation 
    u = var.swm.u 
    c = var.swm.c
    h = var.swm.h
    t = sol.t
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
    @show "YESSSSSSSSSSSSS"
    jldsave(joinpath("./TestOutput/", "Test.jl"); 
            u = integrator.sol.u.swm.u, 
            c = integrator.sol.u.swm.c, 
            h = integrator.sol.u.swm.h, 
            t = integrator.sol.t)
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
