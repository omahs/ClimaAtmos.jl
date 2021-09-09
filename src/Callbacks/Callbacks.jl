module Callbacks
# module begin 
#
using DiffEqCallbacks
using UnPack
using JLD2

using ClimaAtmos.Models: AbstractModel

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

include("callback.jl")

end
