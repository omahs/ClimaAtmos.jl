#####
##### Precipitation models
#####

import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics as CM
import Thermodynamics as TD
import ClimaCore.Spaces as Spaces
import ClimaCore.Operators as Operators
import ClimaCore.Fields as Fields

#####
##### No Precipitation
#####

precipitation_cache(Y, precip_model::NoPrecipitation) = (; precip_model)
precipitation_tendency!(Yₜ, Y, p, t, colidx, ::NoPrecipitation) = nothing

#####
##### 0-Moment without sgs scheme
#####

function precipitation_cache(Y, precip_model::Microphysics0Moment)
    FT = Spaces.undertype(axes(Y.c))
    return (;
        precip_model,
        ᶜS_ρq_tot = similar(Y.c, FT),
        ᶜλ = similar(Y.c, FT),
        ᶜ3d_rain = similar(Y.c, FT),
        ᶜ3d_snow = similar(Y.c, FT),
        col_integrated_rain = similar(Fields.level(Y.c.ρ, 1), FT),
        col_integrated_snow = similar(Fields.level(Y.c.ρ, 1), FT),
    )
end

function compute_precipitation_cache!(
    Y,
    p,
    colidx,
    ::Microphysics0Moment,
    ::Nothing,
)
    (; ᶜts, ᶜS_ρq_tot, params) = p
    cm_params = CAP.microphysics_params(params)
    thermo_params = CAP.thermodynamics_params(params)
    @. ᶜS_ρq_tot[colidx] =
        Y.c.ρ[colidx] * CM.Microphysics0M.remove_precipitation(
            cm_params,
            TD.PhasePartition(thermo_params, ᶜts[colidx]),
        )
end

function compute_precipitation_cache!(
    Y,
    p,
    colidx,
    ::Microphysics0Moment,
    ::TC.EDMFModel,
)
    (; ᶜS_ρq_tot) = p
    qt_tendency_precip_formation_en =
        p.edmf_cache.aux.cent.turbconv.en.qt_tendency_precip_formation[colidx]
    qt_tendency_precip_formation_bulk =
        p.edmf_cache.aux.cent.turbconv.bulk.qt_tendency_precip_formation[colidx]

    @. ᶜS_ρq_tot[colidx] =
        Y.c.ρ[colidx] *
        (qt_tendency_precip_formation_bulk + qt_tendency_precip_formation_en)
end

function precipitation_tendency!(
    Yₜ,
    Y,
    p,
    t,
    colidx,
    precip_model::Microphysics0Moment,
)
    FT = Spaces.undertype(axes(Y.c))
    (;
        ᶜts,
        ᶜΦ,
        ᶜT,
        ᶜ3d_rain,
        ᶜ3d_snow,
        ᶜS_ρq_tot,
        ᶜλ,
        col_integrated_rain,
        col_integrated_snow,
        params,
        turbconv_model,
    ) = p # assume ᶜts has been updated
    thermo_params = CAP.thermodynamics_params(params)
    cm_params = CAP.microphysics_params(params)
    compute_precipitation_cache!(Y, p, colidx, precip_model, turbconv_model)
    @. Yₜ.c.ρq_tot[colidx] += ᶜS_ρq_tot[colidx]
    @. Yₜ.c.ρ[colidx] += ᶜS_ρq_tot[colidx]

    # update precip in cache for coupler's use
    # 3d rain and snow
    @. ᶜT[colidx] = TD.air_temperature(thermo_params, ᶜts[colidx])
    @. ᶜ3d_rain[colidx] =
        ifelse(ᶜT[colidx] >= FT(273.15), ᶜS_ρq_tot[colidx], FT(0))
    @. ᶜ3d_snow[colidx] =
        ifelse(ᶜT[colidx] < FT(273.15), ᶜS_ρq_tot[colidx], FT(0))
    Operators.column_integral_definite!(
        col_integrated_rain[colidx],
        ᶜ3d_rain[colidx],
    )
    Operators.column_integral_definite!(
        col_integrated_snow[colidx],
        ᶜ3d_snow[colidx],
    )

    @. col_integrated_rain[colidx] =
        col_integrated_rain[colidx] / CAP.ρ_cloud_liq(params)
    @. col_integrated_snow[colidx] =
        col_integrated_snow[colidx] / CAP.ρ_cloud_liq(params)

    # liquid fraction
    @. ᶜλ[colidx] = TD.liquid_fraction(thermo_params, ᶜts[colidx])

    if :ρe_tot in propertynames(Y.c)
        @. Yₜ.c.ρe_tot[colidx] +=
            ᶜS_ρq_tot[colidx] * (
                ᶜλ[colidx] *
                TD.internal_energy_liquid(thermo_params, ᶜts[colidx]) +
                (1 - ᶜλ[colidx]) *
                TD.internal_energy_ice(thermo_params, ᶜts[colidx]) +
                ᶜΦ[colidx]
            )
    end
    return nothing
end


#####
##### 1-Moment coupled to sgs
#####
function precipitation_cache(Y, precip_model::Microphysics1Moment)
    FT = Spaces.undertype(axes(Y.c))

    return (;
        precip_model,
        ᶜS_ρe_tot = similar(Y.c, FT),
        ᶜS_ρq_tot = similar(Y.c, FT),
        ᶜS_ρq_rai = similar(Y.c, FT),
        ᶜS_ρq_sno = similar(Y.c, FT),
        ᶜS_qr_evap = similar(Y.c, FT),
        #ᶜS_qs_sub_dep = similar(Y.c, FT),
        #ᶜS_qs_melt = similar(Y.c, FT),
    )
end

function compute_precipitation_cache!(
    Y,
    p,
    colidx,
    ::Microphysics1Moment,
    ::TC.EDMFModel,
)
    (;
        params,
        ᶜts,
        ᶜΦ,
        ᶜT,
        Δt,
        ᶜS_ρe_tot,
        ᶜS_ρq_tot,
        ᶜS_ρq_rai,
        ᶜS_ρq_sno,
        ᶜS_qr_evap,
        #ᶜS_qs_sub_dep,
        #ᶜS_qs_melt
    ) = p

    FT = Spaces.undertype(axes(Y.c))

    th_p = CAP.thermodynamics_params(params)
    mp_p = CAP.microphysics_params(params)
    rain_type = CM.CommonTypes.RainType()
    snow_type = CM.CommonTypes.SnowType()

    # SGS precipitation sources from updrafts and environment
    e_tot_tendency_precip_formation_en =
        p.edmf_cache.aux.cent.turbconv.en.e_tot_tendency_precip_formation[colidx]
    e_tot_tendency_precip_formation_bulk =
        p.edmf_cache.aux.cent.turbconv.bulk.e_tot_tendency_precip_formation[colidx]
    qt_tendency_precip_formation_en =
        p.edmf_cache.aux.cent.turbconv.en.qt_tendency_precip_formation[colidx]
    qt_tendency_precip_formation_bulk =
        p.edmf_cache.aux.cent.turbconv.bulk.qt_tendency_precip_formation[colidx]
    qr_tendency_precip_formation_en =
        p.edmf_cache.aux.cent.turbconv.en.qr_tendency_precip_formation[colidx]
    qr_tendency_precip_formation_bulk =
        p.edmf_cache.aux.cent.turbconv.bulk.qr_tendency_precip_formation[colidx]
    qs_tendency_precip_formation_en =
        p.edmf_cache.aux.cent.turbconv.en.qs_tendency_precip_formation[colidx]
    qs_tendency_precip_formation_bulk =
        p.edmf_cache.aux.cent.turbconv.bulk.qs_tendency_precip_formation[colidx]

    # Grid scale precipitation sinks
    # (evaporation, melting, deposition/sublimation)
    ρ_c = Y.c.ρ[colidx]
    ρq_rai = Y.c.ρq_rai[colidx]
    ρq_sno = Y.c.ρq_sno[colidx]
    @. ᶜT[colidx] = TD.air_temperature(th_p, ᶜts[colidx])

    @. ᶜS_qr_evap[colidx] = -min(
        max(ρq_rai / ρ_c / Δt, FT(0)),
        CM1.evaporation_sublimation(
            mp_p,
            rain_type,
            TD.PhasePartition(th_p, ᶜts[colidx]),
            ρq_rai / ρ_c,
            ρ_c,
            ᶜT[colidx],
        ),
    )
    #@info "RHS_evap" ᶜS_qr_evap
    #@info "ρq_rai" ρq_rai
    #@info "ρ" ρ_c
    #@info "T" ᶜT

    #@. ᶜS_qs_melt[colidx] = -min(
    #    ρq_sno / ρ_c / Δt,
    #    CM1.snow_melt(
    #        mp_p,
    #        ρq_sno / ρ_c,
    #        ρ_c,
    #        TD.air_temperature(th_p, ᶜts[colidx]),
    #    )
    #)
    #@. ᶜS_qs_sub_dep[colidx] = CM1.evaporation_sublimation(
    #    mp_p,
    #    snow_type,
    #    TD.PhasePartition(th_p, ᶜts[colidx]),
    #    ρq_sno / ρ_c,
    #    ρ_c,
    #    TD.air_temperature(th_p, ᶜts[colidx]),
    #)
    #@. ᶜS_qs_sub_dep[colidx] = ifelse(
    #    ᶜS_qs_sub_dep[colidx] > FT(0),
    #    min(
    #        TD.vapor_specific_humidity(th_p, ᶜts[colidx]) / Δt,
    #        ᶜS_qs_sub_dep[colidx]
    #    ),
    #    -min(ρq_sno / ρ_c / Δt, -ᶜS_qs_sub_dep[colidx])
    #)

    # Total microphysics source terms
    @. ᶜS_ρe_tot[colidx] =
        Y.c.ρ[colidx] * (
            e_tot_tendency_precip_formation_bulk +
            e_tot_tendency_precip_formation_en
            #-
            #ᶜS_qr_evap[colidx] * (
            #    TD.internal_energy_liquid(th_p, ᶜts[colidx]) +
            #    ᶜΦ[colidx]
            #)
            #-
            #ᶜS_qs_sub_dep[colidx] * (
            #    TD.internal_energy_ice(th_p, ᶜts[colidx]) +
            #    ᶜΦ[colidx]
            #) +
            #ᶜS_qs_melt[colidx] * TD.latent_heat_fusion(th_p, ᶜts[colidx])
        )
    @. ᶜS_ρq_tot[colidx] =
        Y.c.ρ[colidx] * (
            qt_tendency_precip_formation_bulk +
            qt_tendency_precip_formation_en
            #-
            #ᶜS_qr_evap[colidx]
            #-
            #ᶜS_qs_sub_dep[colidx]
        )
    @. ᶜS_ρq_rai[colidx] =
        Y.c.ρ[colidx] * (
            qr_tendency_precip_formation_bulk +
            qr_tendency_precip_formation_en
            #+
            #ᶜS_qr_evap[colidx]
            #-
            #ᶜS_qs_melt[colidx]
        )
    @. ᶜS_ρq_sno[colidx] =
        Y.c.ρ[colidx] * (
            qs_tendency_precip_formation_bulk +
            qs_tendency_precip_formation_en
            #+
            #ᶜS_qs_sub_dep[colidx] +
            #ᶜS_qs_melt[colidx]
        )
end

"""
Computes the rain and snow advection (down) tendency
"""
function precipitation_advection_tendency!(
    Yₜ,
    Y,
    p,
    colidx,
    ::Microphysics1Moment,
)
    FT = Spaces.undertype(axes(Y.c))
    (; params) = p

    ρ_c = Y.c.ρ[colidx]

    # helper to calculate the rain velocity
    # TODO: assuming w_gm = 0
    # TODO: verify translation

    ρq_rai = Y.c.ρq_rai[colidx]
    ρq_sno = Y.c.ρq_sno[colidx]

    RB = Operators.RightBiasedC2F(; top = Operators.SetValue(FT(0)))
    ᶜdivᵥ = Operators.DivergenceF2C(; bottom = Operators.Extrapolate())
    wvec = Geometry.WVector
    microphys_params = CAP.microphysics_params(params)
    rain_type = CM.CommonTypes.RainType()
    snow_type = CM.CommonTypes.SnowType()

    # TODO - some positivity limiters are needed

    # TODO: need to add horizontal advection + vertical velocity of air

    # TODO: use correct advection operators
    @. Yₜ.c.ρq_rai[colidx] += ᶜdivᵥ(
        wvec(
            RB(
                ρq_rai * CM1.terminal_velocity(
                    microphys_params,
                    rain_type,
                    ρ_c,
                    ρq_rai / ρ_c,
                ),
            ),
        ),
    )

    @. Yₜ.c.ρq_sno[colidx] += ᶜdivᵥ(
        wvec(
            RB(
                ρq_sno * CM1.terminal_velocity(
                    microphys_params,
                    snow_type,
                    ρ_c,
                    ρq_sno / ρ_c,
                ),
            ),
        ),
    )
    return nothing
end

function precipitation_tendency!(
    Yₜ,
    Y,
    p,
    t,
    colidx,
    precip_model::Microphysics1Moment,
)
    (; ᶜS_ρe_tot, ᶜS_ρq_tot, ᶜS_ρq_rai, ᶜS_ρq_sno) = p
    compute_precipitation_cache!(
        Y,
        p,
        colidx,
        precip_model,
        p.atmos.turbconv_model,
    )

    @. Yₜ.c.ρ[colidx] += ᶜS_ρq_tot[colidx]
    @. Yₜ.c.ρq_tot[colidx] += ᶜS_ρq_tot[colidx]
    @. Yₜ.c.ρq_rai[colidx] += ᶜS_ρq_rai[colidx]
    @. Yₜ.c.ρq_sno[colidx] += ᶜS_ρq_sno[colidx]

    if :ρe_tot in propertynames(Y.c)
        @. Yₜ.c.ρe_tot[colidx] += ᶜS_ρe_tot[colidx]
    else
        error(
            "1-moment microphysics can only be coupled to ρe_tot energy variable",
        )
    end

    precipitation_advection_tendency!(Yₜ, Y, p, colidx, precip_model)
    return nothing
end
