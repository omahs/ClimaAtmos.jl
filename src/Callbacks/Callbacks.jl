module Callbacks

export cfl_cb,
       AbstractCallback,
       CFLCallback,
       generate_callback

abstract type AbstractCallback end

struct CFLCallback <: AbstractCallback 
    interval::Number
end

function generate_callback(::AbstractCallback)
    return nothing
end

function generate_callback(c::CFLCallback; kwargs...)
    return PeriodicCallback(f,c.interval; initial_affect=true, kwargs...)
end

end
