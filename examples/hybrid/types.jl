abstract type AbstractMoistureModel end
struct DryModel <: AbstractMoistureModel end
struct EquilMoistModel <: AbstractMoistureModel end
struct NonEquilMoistModel <: AbstractMoistureModel end

abstract type AbstractCompressibilityModel end
struct CompressibleFluid <: AbstractCompressibilityModel end
struct AnelasticFluid <: AbstractCompressibilityModel end

function moisture_model(parsed_args)
    moisture_name = parsed_args["moist"]
    @assert moisture_name in ("dry", "equil", "nonequil")
    return if moisture_name == "dry"
        DryModel()
    elseif moisture_name == "equil"
        EquilMoistModel()
    elseif moisture_name == "nonequil"
        NonEquilMoistModel()
    end
end

abstract type AbstractEnergyFormulation end
struct PotentialTemperature <: AbstractEnergyFormulation end
struct TotalEnergy <: AbstractEnergyFormulation end
struct InternalEnergy <: AbstractEnergyFormulation end

function energy_form(parsed_args)
    energy_name = parse_arg(parsed_args, "energy_name", "rhoe")
    @assert energy_name in ("rhoe", "rhoe_int", "rhotheta")
    return if energy_name == "rhoe"
        TotalEnergy()
    elseif energy_name == "rhoe_int"
        InternalEnergy()
    elseif energy_name == "rhotheta"
        PotentialTemperature()
    end
end


function radiation_model(parsed_args)
    radiation_name = parsed_args["rad"]
    @assert radiation_name in
            (nothing, "clearsky", "gray", "allsky", "allskywithclear")
    return if radiation_name == "clearsky"
        RRTMGPI.ClearSkyRadiation()
    elseif radiation_name == "gray"
        RRTMGPI.GrayRadiation()
    elseif radiation_name == "allsky"
        RRTMGPI.AllSkyRadiation()
    elseif radiation_name == "allskywithclear"
        RRTMGPI.AllSkyRadiationWithClearSkyDiagnostics()
    else
        nothing
    end
end

abstract type AbstractMicrophysicsModel end
struct Microphysics0Moment <: AbstractMicrophysicsModel end

function microphysics_model(parsed_args)
    microphysics_name = parsed_args["microphy"]
    @assert microphysics_name in (nothing, "0M")
    return if microphysics_name == nothing
        nothing
    elseif microphysics_name == "0M"
        Microphysics0Moment()
    end
end

abstract type AbstractForcing end
struct HeldSuarezForcing <: AbstractForcing end

function forcing_type(parsed_args)
    forcing = parsed_args["forcing"]
    @assert forcing in (nothing, "held_suarez")
    return if forcing == nothing
        nothing
    elseif forcing == "held_suarez"
        HeldSuarezForcing()
    end
end

function turbconv_model(FT, parsed_args, namelist)
    turbconv = parsed_args["turbconv"]
    precip_model = nothing
    @assert turbconv in (nothing, "edmf")
    return if turbconv == "edmf"
        TC.EDMFModel(FT, namelist, precip_model)
    else
        nothing
    end
end

Base.broadcastable(x::AbstractMoistureModel) = Ref(x)
Base.broadcastable(x::AbstractEnergyFormulation) = Ref(x)
Base.broadcastable(x::AbstractMicrophysicsModel) = Ref(x)
Base.broadcastable(x::AbstractForcing) = Ref(x)


function get_model_spec(::Type{FT}, parsed_args, namelist) where {FT}

    model_spec = (;
        moisture_model = moisture_model(parsed_args),
        energy_form = energy_form(parsed_args),
        radiation_model = radiation_model(parsed_args),
        microphysics_model = microphysics_model(parsed_args),
        forcing_type = forcing_type(parsed_args),
        turbconv_model = turbconv_model(FT, parsed_args, namelist),
    )

    return model_spec
end

function get_numerics(parsed_args)

    numerics = (;
        upwinding_mode = Symbol(
            parse_arg(parsed_args, "upwinding", "third_order"),
        )
    )
    @assert numerics.upwinding_mode in (:none, :first_order, :third_order)

    return numerics
end

function get_simulation(parsed_args)

    sim = (; is_distributed = haskey(ENV, "CLIMACORE_DISTRIBUTED"))

    return sim
end
