module Callbacks
# module begin 
#
using DiffEqCallbacks
using UnPack
using JLD2

export cfl_cb,
       generate_callback,
       AbstractCallback,
       CFLInfo,
       JLD2Output

"""
    AbstractCallback
    Abstract type for callback definitions. 
"""
abstract type AbstractCallback end
function generate_callback(::AbstractCallback)
    return nothing
end

# CFL Information Callback
"""
    CFLInfo <: AbstractCallback
    Container for CFL information callback
"""
struct CFLInfo <: AbstractCallback 
    interval::Number
end
function (F::CFLInfo)(integrator)
    # Solution variables
    var = integrator.u.swm
    u = getproperty(integrator.u.swm, :u)
    u = var.u 
    c = var.c
    h = var.h
    # CFL computation 
    t = integrator.t
    return nothing
end

"""
    generate_callback(F::CFLInfo; kwargs...)

    Creates a PeriodicCallback object that calculates and prints CFL 
    information
"""
function generate_callback(F::CFLInfo; kwargs...)
    return PeriodicCallback(F,F.interval; initial_affect=false, kwargs...)
end

# JLD2 Callback
"""
    JLD2Output <: AbstractCallback
    Container for JLD2 output callback (writes to disk)
"""
struct JLD2Output <: AbstractCallback 
    filedir  :: String
    filename :: String
    interval :: Number
end
function (F::JLD2Output)(integrator)
    # Solution variables
    u = integrator.u.swm.u
    h = integrator.u.swm.h 
    c = integrator.u.swm.c
    t = integrator.t
    # Create directory
    mkpath(F.filedir)
    # Save data
    jldsave(joinpath("./", F.filedir, F.filename*".jld2"), 
           u=u,
           h=h,
           c=c,
           t=t)
    return nothing
end
"""
    generate_callback(F::JLD2Output; kwargs...)

    Creates a PeriodicCallback object that extracts solution
    variables from the integrator and stores them in a jld2 file. 
"""
function generate_callback(F::JLD2Output; kwargs...)
    return PeriodicCallback(F,F.interval; initial_affect=false, kwargs...)
end

# module end
end
