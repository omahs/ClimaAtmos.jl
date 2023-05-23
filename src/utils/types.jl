abstract type AbstractMoistureModel end
struct DryModel <: AbstractMoistureModel end
struct EquilMoistModel <: AbstractMoistureModel end
struct NonEquilMoistModel <: AbstractMoistureModel end

abstract type AbstractEnergyFormulation end
struct PotentialTemperature <: AbstractEnergyFormulation end
struct TotalEnergy <: AbstractEnergyFormulation end

abstract type AbstractPrecipitationModel end
struct NoPrecipitation <: AbstractPrecipitationModel end
struct Microphysics0Moment <: AbstractPrecipitationModel end
struct Microphysics1Moment <: AbstractPrecipitationModel end

abstract type AbstractModelConfig end
struct SingleColumnModel <: AbstractModelConfig end
struct SphericalModel <: AbstractModelConfig end
struct BoxModel <: AbstractModelConfig end
struct PlaneModel <: AbstractModelConfig end

abstract type AbstractHyperdiffusion end
Base.@kwdef struct ClimaHyperdiffusion{FT} <: AbstractHyperdiffusion
    κ₄::FT
    divergence_damping_factor::FT
end

abstract type AbstractVerticalDiffusion end
Base.@kwdef struct VerticalDiffusion{DM, FT} <: AbstractVerticalDiffusion
    C_E::FT
end
diffuse_momentum(::VerticalDiffusion{DM}) where {DM} = DM
diffuse_momentum(::Nothing) = false

abstract type AbstractSponge end
Base.@kwdef struct ViscousSponge{FT} <: AbstractSponge
    zd::FT
    κ₂::FT
end

Base.@kwdef struct RayleighSponge{FT} <: AbstractSponge
    zd::FT
    α_uₕ::FT
    α_w::FT
end

abstract type AbstractGravityWave end
Base.@kwdef struct NonOrographyGravityWave{FT} <: AbstractGravityWave
    source_pressure::FT = 31500
    damp_pressure::FT = 85
    source_height::FT = 15000
    Bw::FT = 1.0
    Bn::FT = 1.0
    dc::FT = 0.6
    cmax::FT = 99.6
    c0::FT = 0
    nk::FT = 1
    cw::FT = 40.0
    cw_tropics::FT = 40.0
    cn::FT = 40.0
    Bt_0::FT = 0.0003
    Bt_n::FT = 0.0003
    Bt_s::FT = 0.0003
    Bt_eq::FT = 0.0003
    ϕ0_n::FT = 30
    ϕ0_s::FT = -30
    dϕ_n::FT = 5
    dϕ_s::FT = -5
end

Base.@kwdef struct OrographicGravityWave{FT, S} <: AbstractGravityWave
    γ::FT = 0.4
    ϵ::FT = 0.0
    β::FT = 0.5
    h_frac::FT = 0.1
    ρscale::FT = 1.2
    L0::FT = 80e3
    a0::FT = 0.9
    a1::FT = 3.0
    Fr_crit::FT = 0.7
    topo_info::S = "gfdl_restart"
end

abstract type AbstractForcing end
struct HeldSuarezForcing <: AbstractForcing end
struct Subsidence{T} <: AbstractForcing
    prof::T
end
# TODO: is this a forcing?
struct LargeScaleAdvection{PT, PQ}
    prof_dTdt::PT # Set large-scale cooling
    prof_dqtdt::PQ # Set large-scale drying
end

struct EDMFCoriolis{U, V, FT}
    prof_ug::U
    prof_vg::V
    coriolis_param::FT
end

struct EDMFX{N, FT}
    a_half::FT # WARNING: this should never be used outside of divide_by_ρa
end
EDMFX{N}(a_half::FT) where {N, FT} = EDMFX{N, FT}(a_half)
n_mass_flux_subdomains(::EDMFX{N}) where {N} = N
n_mass_flux_subdomains(::Any) = 0

abstract type AbstractSurfaceThermoState end
struct GCMSurfaceThermoState <: AbstractSurfaceThermoState end

# Define broadcasting for types
Base.broadcastable(x::AbstractSurfaceThermoState) = tuple(x)
Base.broadcastable(x::AbstractMoistureModel) = tuple(x)
Base.broadcastable(x::AbstractEnergyFormulation) = tuple(x)
Base.broadcastable(x::AbstractPrecipitationModel) = tuple(x)
Base.broadcastable(x::AbstractForcing) = tuple(x)
Base.broadcastable(x::EDMFX) = tuple(x)

Base.@kwdef struct RadiationDYCOMS_RF01{FT}
    "Large-scale divergence"
    divergence::FT = 3.75e-6
    alpha_z::FT = 1.0
    kappa::FT = 85.0
    F0::FT = 70.0
    F1::FT = 22.0
end
import AtmosphericProfilesLibrary as APL

struct RadiationTRMM_LBA{R}
    rad_profile::R
    function RadiationTRMM_LBA(::Type{FT}) where {FT}
        rad_profile = APL.TRMM_LBA_radiation(FT)
        return new{typeof(rad_profile)}(rad_profile)
    end
end

# TODO: remove AbstractPerformanceMode and all subtypes
# This is temporarily needed to investigate performance of
# our handling of tracers.
abstract type AbstractPerformanceMode end
struct PerfExperimental <: AbstractPerformanceMode end
struct PerfStandard <: AbstractPerformanceMode end
Base.broadcastable(x::AbstractPerformanceMode) = tuple(x)

Base.@kwdef struct AtmosModel{
    MC,
    PEM,
    MM,
    EF,
    PM,
    F,
    S,
    RM,
    LA,
    EC,
    EAT,
    EED,
    TCM,
    SS,
    NOGW,
    OGW,
    HD,
    VD,
    VS,
    RS,
}
    model_config::MC = nothing
    perf_mode::PEM = nothing
    moisture_model::MM = nothing
    energy_form::EF = nothing
    precip_model::PM = nothing
    forcing_type::F = nothing
    subsidence::S = nothing
    radiation_mode::RM = nothing
    ls_adv::LA = nothing
    edmf_coriolis::EC = nothing
    edmfx_adv_test::EAT = nothing
    edmfx_entr_detr::EED = nothing
    turbconv_model::TCM = nothing
    surface_scheme::SS = nothing
    non_orographic_gravity_wave::NOGW = nothing
    orographic_gravity_wave::OGW = nothing
    hyperdiff::HD = nothing
    vert_diff::VD = nothing
    viscous_sponge::VS = nothing
    rayleigh_sponge::RS = nothing
end

Base.broadcastable(x::AtmosModel) = tuple(x)

function Base.summary(io::IO, atmos::AtmosModel)
    pns = string.(propertynames(atmos))
    buf = maximum(length.(pns))
    keys = propertynames(atmos)
    vals = repeat.(" ", map(s -> buf - length(s) + 2, pns))
    bufs = (; zip(keys, vals)...)
    print(io, '\n')
    for pn in propertynames(atmos)
        prop = getproperty(atmos, pn)
        # Skip some data:
        prop isa Bool && continue
        prop isa NTuple && continue
        prop isa Int && continue
        prop isa Float64 && continue
        prop isa Float32 && continue
        s = string(
            "  ", # needed for some reason
            getproperty(bufs, pn),
            '`',
            string(pn),
            '`',
            "::",
            '`',
            typeof(prop),
            '`',
            '\n',
        )
        print(io, s)
    end
end
