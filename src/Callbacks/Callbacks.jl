module Callbacks

using DiffEqCallbacks
using UnPack
using JLD2

export cfl_cb,
       generate_callback,
       AbstractCallback,
       CFLInfo,
       JLD2Output

abstract type AbstractCallback end
function generate_callback(::AbstractCallback)
    return nothing
end

# CFL Information Callback
struct CFLInfo <: AbstractCallback 
    interval::Number
end
function (F::CFLInfo)(integrator)
    # Simulation out of scope , model need not be
    var = integrator.u.swm
    u = getproperty(integrator.u.swm, :u)
    u = var.u 
    c = var.c
    h = var.h
    t = integrator.t
    return nothing
end
function generate_callback(F::CFLInfo; kwargs...)
    return PeriodicCallback(F,F.interval; initial_affect=false, kwargs...)
end

# JLD2 Callback
struct JLD2Output <: AbstractCallback 
    filedir  :: String
    filename :: String
    interval :: Number
end
function (F::JLD2Output)(integrator)
    u = integrator.u.swm.u
    h = integrator.u.swm.h 
    c = integrator.u.swm.c
    t = integrator.t
    mkpath(F.filedir)
    jldsave(joinpath("./", F.filedir, F.filename*".jld2"), 
           u=u,
           h=h,
           c=c,
           t=t)
    return nothing
end
function generate_callback(F::JLD2Output; kwargs...)
    return PeriodicCallback(F,F.interval; initial_affect=false, kwargs...)
end
end
