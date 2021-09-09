
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
    # CFL computation 
    dt = integrator.dt
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
    model :: AbstractModel
    filedir  :: String
    filename :: String
    interval :: Number
end
function (F::JLD2Output)(integrator)
    # Create directory
    mkpath(joinpath(@__DIR__,F.filedir))
    # Save data
    jldsave(joinpath(@__DIR__, F.filedir, F.filename*".jld2"), 
            integrator=integrator)
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
