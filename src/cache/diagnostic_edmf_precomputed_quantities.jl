#####
##### Precomputed quantities
#####
import Thermodynamics as TD
import ClimaCore: Spaces, Fields, RecursiveApply

function set_diagnostic_edmfx_draft_quantities_level!(
    thermo_params,
    K,
    ts,
    ρ,
    uₕ,
    local_geometry_level,
    u³,
    local_geometry_halflevel,
    h_tot,
    q_tot,
    p,
    Φ,
)
    FT = eltype(thermo_params)
    @. K =
        (
            dot(
                C123(uₕ, local_geometry_level),
                CT123(uₕ, local_geometry_level),
            ) +
            dot(
                C123(u³, local_geometry_halflevel),
                CT123(u³, local_geometry_halflevel),
            ) +
            2 * dot(
                CT123(uₕ, local_geometry_level),
                C123(u³, local_geometry_halflevel),
            )
        ) / 2
    @. ts =
        TD.PhaseEquil_phq(thermo_params, p, h_tot - K - Φ, q_tot, 8, FT(0.0003))
    @. ρ = TD.air_density(thermo_params, ts)
    return nothing
end

function set_diagnostic_edmfx_env_quantities_level!(
    ρ_level,
    ρaʲs_level,
    u³_halflevel,
    u³ʲs_halflevel,
    u³⁰_halflevel,
    turbconv_model,
)
    @. u³⁰_halflevel = divide_by_ρa(
        ρ_level * u³_halflevel - mapreduce(*, +, ρaʲs_level, u³ʲs_halflevel),
        ρ_level,
        ρ_level * u³_halflevel,
        ρ_level,
        turbconv_model,
    )
    return nothing
end

function surface_flux_tke(
    turbconv_params,
    ρ_int,
    u_int,
    ustar,
    interior_local_geometry,
    surface_local_geometry,
)
    c_d = TCP.tke_diss_coeff(turbconv_params)
    c_m = TCP.tke_ed_coeff(turbconv_params)
    k_star² = TCP.tke_surf_scale(turbconv_params)
    speed = Geometry._norm(
        CA.CT12(u_int, interior_local_geometry),
        interior_local_geometry,
    )
    c3_unit = C3(unit_basis_vector_data(C3, surface_local_geometry))
    return ρ_int * (1 - c_d * c_m * k_star²^2) * ustar^2 * speed * c3_unit
end

"""
    set_diagnostic_edmf_precomputed_quantities!(Y, p, t)

Updates the precomputed quantities stored in `p` for diagnostic edmfx.
"""
function set_diagnostic_edmf_precomputed_quantities!(Y, p, t)
    (; moisture_model, turbconv_model) = p.atmos
    FT = eltype(Y)
    n = n_mass_flux_subdomains(turbconv_model)
    ᶜz = Fields.coordinate_field(Y.c).z
    ᶜdz = Fields.Δz_field(axes(Y.c))
    (; params) = p
    (; dt) = p.simulation
    (; ᶜp, ᶜΦ, ᶜρ_ref, ᶠu³, ᶜu, ᶜts, ᶜh_tot) = p
    (; q_tot) = p.ᶜspecific
    (; ustar, obukhov_length, buoyancy_flux, ρ_flux_h_tot, ρ_flux_q_tot) =
        p.sfc_conditions
    (;
        ᶜρaʲs,
        ᶠu³ʲs,
        ᶜuʲs,
        ᶜKʲs,
        ᶜh_totʲs,
        ᶜq_totʲs,
        ᶜtsʲs,
        ᶜρʲs,
        ᶜentr_detrʲs,
        ᶜnh_pressureʲs,
    ) = p
    (; ᶠu³⁰, ᶜu⁰, ᶜtke⁰, ᶜlinear_buoygrad, ᶜshear², ᶜmixing_length) = p
    (; ᶜK_h, ᶜK_u, ρatke_flux) = p
    thermo_params = CAP.thermodynamics_params(params)
    ᶜlg = Fields.local_geometry_field(Y.c)

    @. ᶜtke⁰ = Y.c.sgs⁰.ρatke / Y.c.ρ

    ᶜ∇Φ³ = p.ᶜtemp_CT3
    @. ᶜ∇Φ³ = CT3(ᶜgradᵥ(ᶠinterp(ᶜΦ)))
    @. ᶜ∇Φ³ += CT3(gradₕ(ᶜΦ))

    ρaʲu³ʲ_data = p.temp_data_level
    u³ʲ_datau³ʲ_data = p.temp_data_level_2
    ρaʲu³ʲ_datah_tot = ρaʲu³ʲ_dataq_tot = p.temp_data_level_3

    ρ_int_level = Fields.field_values(Fields.level(Y.c.ρ, 1))
    uₕ_int_level = Fields.field_values(Fields.level(Y.c.uₕ, 1))
    u³_int_halflevel = Fields.field_values(Fields.level(ᶠu³, half))
    h_tot_int_level = Fields.field_values(Fields.level(ᶜh_tot, 1))
    q_tot_int_level = Fields.field_values(Fields.level(q_tot, 1))

    p_int_level = Fields.field_values(Fields.level(ᶜp, 1))
    Φ_int_level = Fields.field_values(Fields.level(ᶜΦ, 1))

    local_geometry_int_level =
        Fields.field_values(Fields.level(Fields.local_geometry_field(Y.c), 1))
    local_geometry_int_halflevel = Fields.field_values(
        Fields.level(Fields.local_geometry_field(Y.f), half),
    )
    z_int_level =
        Fields.field_values(Fields.level(Fields.coordinate_field(Y.c).z, 1))
    z_sfc_halflevel =
        Fields.field_values(Fields.level(Fields.coordinate_field(Y.f).z, half))
    buoyancy_flux_sfc_halflevel = Fields.field_values(buoyancy_flux)
    ρ_flux_h_tot_sfc_halflevel = Fields.field_values(ρ_flux_h_tot)
    ρ_flux_q_tot_sfc_halflevel = Fields.field_values(ρ_flux_q_tot)
    ustar_sfc_halflevel = Fields.field_values(ustar)
    obukhov_length_sfc_halflevel = Fields.field_values(obukhov_length)

    # boundary condition
    for j in 1:n
        ᶜρaʲ = ᶜρaʲs.:($j)
        ᶠu³ʲ = ᶠu³ʲs.:($j)
        ᶜKʲ = ᶜKʲs.:($j)
        ᶜh_totʲ = ᶜh_totʲs.:($j)
        ᶜtsʲ = ᶜtsʲs.:($j)
        ᶜρʲ = ᶜρʲs.:($j)
        ᶜq_totʲ = ᶜq_totʲs.:($j)

        ρaʲ_int_level = Fields.field_values(Fields.level(ᶜρaʲ, 1))
        u³ʲ_int_halflevel = Fields.field_values(Fields.level(ᶠu³ʲ, half))
        Kʲ_int_level = Fields.field_values(Fields.level(ᶜKʲ, 1))
        h_totʲ_int_level = Fields.field_values(Fields.level(ᶜh_totʲ, 1))
        q_totʲ_int_level = Fields.field_values(Fields.level(ᶜq_totʲ, 1))
        tsʲ_int_level = Fields.field_values(Fields.level(ᶜtsʲ, 1))
        ρʲ_int_level = Fields.field_values(Fields.level(ᶜρʲ, 1))

        @. u³ʲ_int_halflevel = CT3(
            Geometry.WVector($(FT(0.01)), local_geometry_int_halflevel),
            local_geometry_int_halflevel,
        )
        @. h_totʲ_int_level = sgs_scalar_first_interior_bc(
            z_int_level - z_sfc_halflevel,
            ρ_int_level,
            h_tot_int_level,
            buoyancy_flux_sfc_halflevel,
            ρ_flux_h_tot_sfc_halflevel,
            ustar_sfc_halflevel,
            obukhov_length_sfc_halflevel,
            local_geometry_int_halflevel,
        )
        @. q_totʲ_int_level = sgs_scalar_first_interior_bc(
            z_int_level - z_sfc_halflevel,
            ρ_int_level,
            q_tot_int_level,
            buoyancy_flux_sfc_halflevel,
            ρ_flux_q_tot_sfc_halflevel,
            ustar_sfc_halflevel,
            obukhov_length_sfc_halflevel,
            local_geometry_int_halflevel,
        )

        set_diagnostic_edmfx_draft_quantities_level!(
            thermo_params,
            Kʲ_int_level,
            tsʲ_int_level,
            ρʲ_int_level,
            uₕ_int_level,
            local_geometry_int_level,
            u³ʲ_int_halflevel,
            local_geometry_int_halflevel,
            h_totʲ_int_level,
            q_totʲ_int_level,
            p_int_level,
            Φ_int_level,
        )

        @. ρaʲ_int_level = ρʲ_int_level * turbconv_model.a_int
    end

    ρaʲs_int_level = Fields.field_values(Fields.level(ᶜρaʲs, 1))
    u³ʲs_int_halflevel = Fields.field_values(Fields.level(ᶠu³ʲs, half))
    u³⁰_int_halflevel = Fields.field_values(Fields.level(ᶠu³⁰, half))
    set_diagnostic_edmfx_env_quantities_level!(
        ρ_int_level,
        ρaʲs_int_level,
        u³_int_halflevel,
        u³ʲs_int_halflevel,
        u³⁰_int_halflevel,
        turbconv_model,
    )

    # integral
    for i in 2:Spaces.nlevels(axes(Y.c))
        ρ_level = Fields.field_values(Fields.level(Y.c.ρ, i))
        uₕ_level = Fields.field_values(Fields.level(Y.c.uₕ, i))
        u³_halflevel = Fields.field_values(Fields.level(ᶠu³, i - half))
        h_tot_level = Fields.field_values(Fields.level(ᶜh_tot, i))
        q_tot_level = Fields.field_values(Fields.level(q_tot, i))
        p_level = Fields.field_values(Fields.level(ᶜp, i))
        Φ_level = Fields.field_values(Fields.level(ᶜΦ, i))
        local_geometry_level = Fields.field_values(
            Fields.level(Fields.local_geometry_field(Y.c), i),
        )
        local_geometry_halflevel = Fields.field_values(
            Fields.level(Fields.local_geometry_field(Y.f), i - half),
        )
        ∂x∂ξ_level = local_geometry_level.∂x∂ξ.components.data
        end_index = fieldcount(eltype(∂x∂ξ_level)) # This will be 4 in 2D and 9 in 3D.
        ∂x³∂ξ³_level = ∂x∂ξ_level.:($end_index)

        ρ_ref_prev_level = Fields.field_values(Fields.level(ᶜρ_ref, i - 1))
        ∇Φ³_prev_level = Fields.field_values(Fields.level(ᶜ∇Φ³, i - 1))
        ∇Φ³_data_prev_level = ∇Φ³_prev_level.components.data.:1
        ρ_prev_level = Fields.field_values(Fields.level(Y.c.ρ, i - 1))
        u³⁰_prev_halflevel =
            Fields.field_values(Fields.level(ᶠu³⁰, i - 1 - half))
        u³⁰_data_prev_halflevel = u³⁰_prev_halflevel.components.data.:1
        h_tot_prev_level = Fields.field_values(Fields.level(ᶜh_tot, i - 1))
        q_tot_prev_level = Fields.field_values(Fields.level(q_tot, i - 1))
        ts_prev_level = Fields.field_values(Fields.level(ᶜts, i - 1))
        p_prev_level = Fields.field_values(Fields.level(ᶜp, i - 1))
        z_prev_level = Fields.field_values(Fields.level(ᶜz, i - 1))

        local_geometry_prev_level = Fields.field_values(
            Fields.level(Fields.local_geometry_field(Y.c), i - 1),
        )
        local_geometry_prev_halflevel = Fields.field_values(
            Fields.level(Fields.local_geometry_field(Y.f), i - 1 - half),
        )

        for j in 1:n
            ᶜρaʲ = ᶜρaʲs.:($j)
            ᶠu³ʲ = ᶠu³ʲs.:($j)
            ᶜKʲ = ᶜKʲs.:($j)
            ᶜh_totʲ = ᶜh_totʲs.:($j)
            ᶜtsʲ = ᶜtsʲs.:($j)
            ᶜρʲ = ᶜρʲs.:($j)
            ᶜq_totʲ = ᶜq_totʲs.:($j)
            ᶜentr_detrʲ = ᶜentr_detrʲs.:($j)
            ᶜnh_pressureʲ = ᶜnh_pressureʲs.:($j)

            ρaʲ_level = Fields.field_values(Fields.level(ᶜρaʲ, i))
            u³ʲ_halflevel = Fields.field_values(Fields.level(ᶠu³ʲ, i - half))
            Kʲ_level = Fields.field_values(Fields.level(ᶜKʲ, i))
            h_totʲ_level = Fields.field_values(Fields.level(ᶜh_totʲ, i))
            q_totʲ_level = Fields.field_values(Fields.level(ᶜq_totʲ, i))
            tsʲ_level = Fields.field_values(Fields.level(ᶜtsʲ, i))
            ρʲ_level = Fields.field_values(Fields.level(ᶜρʲ, i))

            ρaʲ_prev_level = Fields.field_values(Fields.level(ᶜρaʲ, i - 1))
            u³ʲ_prev_halflevel =
                Fields.field_values(Fields.level(ᶠu³ʲ, i - 1 - half))
            u³ʲ_data_prev_halflevel = u³ʲ_prev_halflevel.components.data.:1
            h_totʲ_prev_level =
                Fields.field_values(Fields.level(ᶜh_totʲ, i - 1))
            q_totʲ_prev_level =
                Fields.field_values(Fields.level(ᶜq_totʲ, i - 1))
            ρʲ_prev_level = Fields.field_values(Fields.level(ᶜρʲ, i - 1))
            tsʲ_prev_level = Fields.field_values(Fields.level(ᶜtsʲ, i - 1))
            entr_detrʲ_prev_level =
                Fields.field_values(Fields.level(ᶜentr_detrʲ, i - 1))
            nh_pressureʲ_prev_level =
                Fields.field_values(Fields.level(ᶜnh_pressureʲ, i - 1))
            scale_height =
                CAP.R_d(params) * CAP.T_surf_ref(params) / CAP.grav(params)

            @. entr_detrʲ_prev_level = pi_groups_entr_detr(
                params,
                p.atmos.edmfx_entr_detr,
                z_prev_level,
                p_prev_level,
                ρ_prev_level,
                buoyancy_flux_sfc_halflevel,
                ρaʲ_prev_level / ρʲ_prev_level,
                get_physical_w(
                    u³ʲ_prev_halflevel,
                    local_geometry_prev_halflevel,
                ),
                TD.relative_humidity(thermo_params, tsʲ_prev_level),
                ᶜbuoyancy(params, ρ_ref_prev_level, ρʲ_prev_level),
                get_physical_w(
                    u³⁰_prev_halflevel,
                    local_geometry_prev_halflevel,
                ),
                TD.relative_humidity(thermo_params, ts_prev_level),
                ᶜbuoyancy(params, ρ_ref_prev_level, ρ_prev_level),
                dt,
            )

            @. nh_pressureʲ_prev_level = ᶠupdraft_nh_pressure(
                params,
                p.atmos.edmfx_nh_pressure,
                local_geometry_prev_halflevel,
                -∇Φ³_prev_level * (ρʲ_prev_level - ρ_prev_level) /
                ρʲ_prev_level,
                u³ʲ_prev_halflevel,
                u³⁰_prev_halflevel,
                scale_height,
            )

            nh_pressureʲ_data_prev_level =
                nh_pressureʲ_prev_level.components.data.:1

            @. ρaʲu³ʲ_data =
                (1 / local_geometry_halflevel.J) * (
                    local_geometry_prev_halflevel.J *
                    ρaʲ_prev_level *
                    u³ʲ_data_prev_halflevel
                )

            @. ρaʲu³ʲ_data +=
                (1 / local_geometry_halflevel.J) * (
                    local_geometry_prev_level.J *
                    ρaʲ_prev_level *
                    (entr_detrʲ_prev_level.entr - entr_detrʲ_prev_level.detr)
                )

            # Using constant exponents in broadcasts allocate, so we use
            # local_geometry_halflevel.J * local_geometry_halflevel.J instead.
            # See ClimaCore.jl issue #1126.
            @. u³ʲ_datau³ʲ_data =
                (
                    1 /
                    (local_geometry_halflevel.J * local_geometry_halflevel.J)
                ) * (
                    local_geometry_prev_halflevel.J *
                    local_geometry_prev_halflevel.J *
                    u³ʲ_data_prev_halflevel *
                    u³ʲ_data_prev_halflevel
                )

            @. u³ʲ_datau³ʲ_data -=
                (
                    1 /
                    (local_geometry_halflevel.J * local_geometry_halflevel.J)
                ) * (
                    local_geometry_prev_level.J *
                    local_geometry_prev_level.J *
                    2 *
                    (
                        ∇Φ³_data_prev_level * (ρʲ_prev_level - ρ_prev_level) /
                        ρʲ_prev_level
                    )
                )

            @. u³ʲ_datau³ʲ_data +=
                (
                    1 /
                    (local_geometry_halflevel.J * local_geometry_halflevel.J)
                ) * (
                    local_geometry_prev_level.J *
                    local_geometry_prev_level.J *
                    2 *
                    (
                        entr_detrʲ_prev_level.entr * u³⁰_data_prev_halflevel -
                        entr_detrʲ_prev_level.entr * u³ʲ_data_prev_halflevel
                    )
                )

            @. u³ʲ_datau³ʲ_data -=
                (
                    1 /
                    (local_geometry_halflevel.J * local_geometry_halflevel.J)
                ) * (
                    local_geometry_prev_level.J *
                    local_geometry_prev_level.J *
                    2 *
                    nh_pressureʲ_data_prev_level
                )

            minimum_value = FT(1e-6)
            @. u³ʲ_halflevel = ifelse(
                (
                    (
                        u³ʲ_datau³ʲ_data <
                        (minimum_value / (∂x³∂ξ³_level * ∂x³∂ξ³_level))
                    ) | (ρaʲu³ʲ_data < (minimum_value / ∂x³∂ξ³_level))
                ),
                CT3(0),
                CT3(sqrt(max(0, u³ʲ_datau³ʲ_data))),
            )
            @. ρaʲ_level = ifelse(
                (
                    (
                        u³ʲ_datau³ʲ_data <
                        (minimum_value / (∂x³∂ξ³_level * ∂x³∂ξ³_level))
                    ) | (ρaʲu³ʲ_data < (minimum_value / ∂x³∂ξ³_level))
                ),
                0,
                ρaʲu³ʲ_data / sqrt(max(0, u³ʲ_datau³ʲ_data)),
            )

            @. ρaʲu³ʲ_datah_tot =
                (1 / local_geometry_halflevel.J) * (
                    local_geometry_prev_halflevel.J *
                    ρaʲ_prev_level *
                    u³ʲ_data_prev_halflevel *
                    h_totʲ_prev_level
                )
            @. ρaʲu³ʲ_datah_tot +=
                (1 / local_geometry_halflevel.J) * (
                    local_geometry_prev_level.J *
                    ρaʲ_prev_level *
                    (
                        entr_detrʲ_prev_level.entr * h_tot_prev_level -
                        entr_detrʲ_prev_level.detr * h_totʲ_prev_level
                    )
                )
            @. h_totʲ_level = ifelse(
                (
                    (
                        u³ʲ_datau³ʲ_data <
                        (minimum_value / (∂x³∂ξ³_level * ∂x³∂ξ³_level))
                    ) | (ρaʲu³ʲ_data < (minimum_value / ∂x³∂ξ³_level))
                ),
                h_tot_level,
                ρaʲu³ʲ_datah_tot / ρaʲu³ʲ_data,
            )

            @. ρaʲu³ʲ_dataq_tot =
                (1 / local_geometry_halflevel.J) * (
                    local_geometry_prev_halflevel.J *
                    ρaʲ_prev_level *
                    u³ʲ_data_prev_halflevel *
                    q_totʲ_prev_level
                )
            @. ρaʲu³ʲ_dataq_tot +=
                (1 / local_geometry_halflevel.J) * (
                    local_geometry_prev_level.J *
                    ρaʲ_prev_level *
                    (
                        entr_detrʲ_prev_level.entr * q_tot_prev_level -
                        entr_detrʲ_prev_level.detr * q_totʲ_prev_level
                    )
                )
            @. q_totʲ_level = ifelse(
                (
                    (
                        u³ʲ_datau³ʲ_data <
                        (minimum_value / (∂x³∂ξ³_level * ∂x³∂ξ³_level))
                    ) | (ρaʲu³ʲ_data < (minimum_value / ∂x³∂ξ³_level))
                ),
                q_tot_level,
                ρaʲu³ʲ_dataq_tot / ρaʲu³ʲ_data,
            )

            set_diagnostic_edmfx_draft_quantities_level!(
                thermo_params,
                Kʲ_level,
                tsʲ_level,
                ρʲ_level,
                uₕ_level,
                local_geometry_level,
                u³ʲ_halflevel,
                local_geometry_halflevel,
                h_totʲ_level,
                q_totʲ_level,
                p_level,
                Φ_level,
            )
        end

        ρaʲs_level = Fields.field_values(Fields.level(ᶜρaʲs, i))
        u³ʲs_halflevel = Fields.field_values(Fields.level(ᶠu³ʲs, i - half))
        u³⁰_halflevel = Fields.field_values(Fields.level(ᶠu³⁰, i - half))
        set_diagnostic_edmfx_env_quantities_level!(
            ρ_level,
            ρaʲs_level,
            u³_halflevel,
            u³ʲs_halflevel,
            u³⁰_halflevel,
            turbconv_model,
        )
    end

    # set values for the top level
    i_top = Spaces.nlevels(axes(Y.c))
    u³⁰_halflevel = Fields.field_values(Fields.level(ᶠu³⁰, i_top + half))
    @. u³⁰_halflevel = CT3(0)
    for j in 1:n
        ᶠu³ʲ = ᶠu³ʲs.:($j)
        ᶜuʲ = ᶜuʲs.:($j)
        ᶠu³ʲ = ᶠu³ʲs.:($j)
        ᶜentr_detrʲ = ᶜentr_detrʲs.:($j)
        u³ʲ_halflevel = Fields.field_values(Fields.level(ᶠu³ʲ, i_top + half))
        @. u³ʲ_halflevel = CT3(0)
        entr_detrʲ_level = Fields.field_values(Fields.level(ᶜentr_detrʲ, i_top))
        fill!(entr_detrʲ_level, RecursiveApply.rzero(eltype(entr_detrʲ_level)))
        @. ᶜuʲ = C123(Y.c.uₕ) + ᶜinterp(C123(ᶠu³ʲ))
    end

    @. ᶜu⁰ = C123(Y.c.uₕ) + ᶜinterp(C123(ᶠu³⁰))

    @. ᶜlinear_buoygrad = buoyancy_gradients(
        params,
        moisture_model,
        EnvBuoyGrad(
            BuoyGradMean(),
            TD.air_temperature(thermo_params, ᶜts),                           # t_sat
            TD.vapor_specific_humidity(thermo_params, ᶜts),                   # qv_sat
            q_tot,                                                            # qt_sat
            TD.dry_pottemp(thermo_params, ᶜts),                               # θ_sat
            TD.liquid_ice_pottemp(thermo_params, ᶜts),                        # θ_liq_ice_sat
            projected_vector_data(
                C3,
                ᶜgradᵥ(ᶠinterp(TD.virtual_pottemp(thermo_params, ᶜts))),
                ᶜlg,
            ),                                                                 # ∂θv∂z_unsat
            projected_vector_data(C3, ᶜgradᵥ(ᶠinterp(q_tot)), ᶜlg),            # ∂qt∂z_sat
            projected_vector_data(
                C3,
                ᶜgradᵥ(ᶠinterp(TD.liquid_ice_pottemp(thermo_params, ᶜts))),
                ᶜlg,
            ),                                                                 # ∂θl∂z_sat
            ᶜp,                                                                # p
            ifelse(TD.has_condensate(thermo_params, ᶜts), 1, 0),               # en_cld_frac
            Y.c.ρ,                                                             # ρ
        ),
    )

    # TODO: This is not correct with topography or with grid stretching
    ᶠu⁰ = p.ᶠtemp_C123
    @. ᶠu⁰ = C123(ᶠinterp(Y.c.uₕ)) + C123(ᶠu³⁰)
    ct3_unit = p.ᶜtemp_CT3
    @. ct3_unit = CT3(Geometry.WVector(FT(1)), ᶜlg)
    @. ᶜshear² = norm_sqr(adjoint(CA.ᶜgradᵥ(ᶠu⁰)) * ct3_unit)
    ᶜprandtl_nvec = p.ᶜtemp_scalar
    @. ᶜprandtl_nvec = turbulent_prandtl_number(
        params,
        obukhov_length,
        ᶜlinear_buoygrad,
        ᶜshear²,
    )
    ᶜtke_exch = p.ᶜtemp_scalar_2
    @. ᶜtke_exch = 0
    # using ᶜu⁰ would be more correct, but this is more consistent with the
    # TKE equation, where using ᶜu⁰ results in allocation
    for j in 1:n
        @. ᶜtke_exch +=
            ᶜρaʲs.:($$j) * ᶜentr_detrʲs.:($$j).detr / Y.c.ρ *
            (1 / 2 * norm_sqr(ᶜinterp(ᶠu³⁰) - ᶜinterp(ᶠu³ʲs.:($$j))) - ᶜtke⁰)
    end

    sfc_tke = Fields.level(ᶜtke⁰, 1)
    @. ᶜmixing_length = mixing_length(
        params,
        ustar,
        ᶜz,
        ᶜdz,
        max(sfc_tke, 0),
        ᶜlinear_buoygrad,
        max(ᶜtke⁰, 0),
        obukhov_length,
        ᶜshear²,
        ᶜprandtl_nvec,
        ᶜtke_exch,
    )

    turbconv_params = CAP.turbconv_params(params)
    c_m = TCP.tke_ed_coeff(turbconv_params)
    @. ᶜK_u = c_m * ᶜmixing_length * sqrt(max(ᶜtke⁰, 0))
    # TODO: add Prantdl number
    @. ᶜK_h = ᶜK_u / ᶜprandtl_nvec

    ρatke_flux_values = Fields.field_values(ρatke_flux)
    ρ_int_values = Fields.field_values(Fields.level(Y.c.ρ, 1))
    u_int_values = Fields.field_values(Fields.level(ᶜu, 1))
    ustar_values = Fields.field_values(ustar)
    int_local_geometry_values =
        Fields.field_values(Fields.level(Fields.local_geometry_field(Y.c), 1))
    sfc_local_geometry_values = Fields.field_values(
        Fields.level(Fields.local_geometry_field(Y.f), half),
    )
    @. ρatke_flux_values = surface_flux_tke(
        turbconv_params,
        ρ_int_values,
        u_int_values,
        ustar_values,
        int_local_geometry_values,
        sfc_local_geometry_values,
    )

    return nothing
end
