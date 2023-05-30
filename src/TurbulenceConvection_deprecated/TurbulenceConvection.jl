module TurbulenceConvection

# Import ClimaAtmos types
import ..DryModel
import ..EquilMoistModel
import ..NonEquilMoistModel
import ..RadiationDYCOMS_RF01
import ..RadiationTRMM_LBA
import ..AbstractPrecipitationModel
import ..Microphysics0Moment
import ..Microphysics1Moment

import ..compute_kinetic!

import ClimaCore as CC
import ClimaCore.Geometry as CCG
import ClimaCore.Operators as Operators
import ClimaCore.Fields as Fields
import ClimaCore.Spaces as Spaces
import ClimaCore.Geometry: ⊗, _norm_sqr
import ClimaCore.Operators as CCO
import LinearAlgebra as LA
import DocStringExtensions
import StaticArrays as SA
import StatsBase
import Dierckx
import LambertW
import Thermodynamics as TD
import Distributions
import FastGaussQuadrature
import CloudMicrophysics as CM
import CloudMicrophysics.Microphysics0M as CM0
import CloudMicrophysics.Microphysics1M as CM1
import Random

const liq_type = CM.CommonTypes.LiquidType()
const ice_type = CM.CommonTypes.IceType()
const rain_type = CM.CommonTypes.RainType()
const snow_type = CM.CommonTypes.SnowType()

include("Parameters.jl")
import .Parameters as TCP
const APS = TCP.AbstractTurbulenceConvectionParameters

Base.broadcastable(param_set::APS) = tuple(param_set)

#=
    debug_state(state, code_location::String)

A simple function for debugging the entire state,
specifically for when quantities that should remain
positive-definite become negative.

=#
function debug_state(state, code_location::String)
    prog_gm = center_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    prog_gm_f = face_prog_grid_mean(state)
    aux_gm_f = face_aux_grid_mean(state)

    prog_up = center_prog_updrafts(state)
    aux_up = center_aux_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    aux_up_f = face_aux_updrafts(state)

    prog_en = center_prog_environment(state)
    aux_en = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)

    ######
    ###### Positive-definite variables
    ######

    vars_positive = [
        vec(prog_gm_f.u₃),
        vec(prog_up[1].ρarea),
        vec(prog_up_f[1].w),
        vec(aux_en.area),
        vec(aux_en.θ_liq_ice),
    ]
    vars = vars_positive
    vars_conds = map(v -> any(v .< 0), vars)

    if any(vars_conds)
        @show code_location
        for (i, vc, v) in zip(1:length(vars), vars_conds, vars)
            vc || continue
            @show i, v
        end
        @show vars_conds
        error("Negative state for positive-definite field(s)")
    end

    ######
    ###### All listed variables
    ######
    vars = vars_positive
    vars_conds = map(v -> any(isnan.(v)) || any(isinf.(v)), vars)

    if any(vars_conds)
        @show code_location
        for (i, vc, v) in zip(1:length(vars), vars_conds, vars)
            vc || continue
            @show i, v
        end
        @show vars_conds
        error("Nan/Inf state for field(s)")
    end
end

include("Grid.jl")
include("dycore_api.jl")
include("Fields.jl")
include("types.jl")
include("Operators.jl")

include("turbulence_functions.jl")
include("utility_functions.jl")
include("variables.jl")
include("update_aux.jl")
include("EDMF_functions.jl")
include("thermodynamics.jl")
include("closures/microphysics.jl")
include("closures/perturbation_pressure.jl")
include("closures/entr_detr.jl")
include("closures/mixing_length.jl")
include("closures/buoyancy_gradients.jl")

end
