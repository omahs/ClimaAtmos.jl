#####
##### State debugging tools
#####

import ClimaCore.Fields as Fields

function debug_state_generic!(
    state,
    var::Union{Fields.FieldVector, NamedTuple};
    name = "",
)
    for pn in propertynames(var)
        pfx = isempty(name) ? "" : "$name."
        debug_state_generic!(state, getproperty(var, pn); name = "$pfx$pn")
    end
end

function debug_state_generic!(
    state,
    var::Union{Fields.FieldVector, NamedTuple},
    colidx,
    name = "",
)
    for pn in propertynames(var)
        pfx = isempty(name) ? "" : "$name."
        debug_state_generic!(
            state,
            getproperty(var, pn),
            colidx;
            name = "$pfx$pn",
        )
    end
end

debug_state_generic!(state, var::Number; name = "") = nothing
debug_state_generic!(state, var::AbstractString; name = "") = nothing
debug_state_generic!(state, var::Bool; name = "") = nothing
debug_state_generic!(state, var::Nothing; name = "") = nothing
debug_state_generic!(state, var::Any; name = "") = nothing # TODO: should we try to catch more types?

debug_state_generic!(state, var::Number, colidx; name = "") = nothing
debug_state_generic!(state, var::AbstractString, colidx; name = "") = nothing
debug_state_generic!(state, var::Bool, colidx; name = "") = nothing
debug_state_generic!(state, var::Nothing, colidx; name = "") = nothing
debug_state_generic!(state, var::Any, colidx; name = "") = nothing # TODO: should we try to catch more types?

debug_state_generic!(state, var::Fields.Field, colidx; name = "") =
    debug_state_column_field!(state, var[colidx]; name)
debug_state_generic!(state, var::Fields.Field; name = "") =
    debug_state_field!(state, var; name)

debug_state_field!(state, ::Nothing; name = "") = nothing
debug_state_field!(state, ::Nothing, colidx; name = "") = nothing

debug_state_field!(state, prog::Fields.Field, colidx; name = "") =
    debug_state_column_field!(state, prog[colidx]; name)
debug_state_field!(state, prog::Fields.Field; name = "") =
    debug_state_full_field!(state, prog; name)

function debug_state_full_field!(state, prog::Fields.Field; name = "")
    isbad(x) = isnan(x) || isinf(x)
    (; msg) = state
    for prop_chain in Fields.property_chains(prog)
        var = Fields.single_field(prog, prop_chain)
        nan = any(isnan.(parent(var)))
        inf = any(isinf.(parent(var)))
        any(isbad.(parent(var))) || continue
        pfx = isempty(name) ? "" : "$name."
        push!(
            msg,
            "-------------------- Bad data (nan=$nan, inf=$inf) in $(state.name).$pfx$prop_chain",
        )
        push!(msg, sprint(show, var))
    end
end

function debug_state_column_field!(state, prog::Fields.Field; name = "") # can we change this to prof::Fields.ColumnField ?
    isbad(x) = isnan(x) || isinf(x)
    (; msg) = state
    for prop_chain in Fields.property_chains(prog)
        var = Fields.single_field(prog, prop_chain)
        nan = any(isnan.(parent(var)))
        inf = any(isinf.(parent(var)))
        pfx = isempty(name) ? "" : "$name."
        any(isbad.(parent(var))) || continue
        push!(
            msg,
            "-------------------- Bad data (nan=$nan, inf=$inf) in $(state.name).$pfx$prop_chain",
        )
        push!(msg, sprint(show, var))
    end
end

"""
    debug_state(t
        [, colidx];                                            # colidx is optional
        Yₜ::Union{Fields.Field, Fields.FieldVector} = nothing, # Yₜ is optional
        Y::Union{Fields.Field, Fields.FieldVector} = nothing,  # Y is optional
        p::Any = nothing,                                      # p is optional
    )

Helper function for debugging `NaN`s and `Inf`s.

To avoid jumbled printed messages, it's recommended to use this
feature with `ClimaCore.enable_threading() = false`.

## Example
```julia
function precomputed_quantities!(Y, p, t, colidx)
    ᶜuₕ = Y.c.uₕ
    ᶠu₃ = Y.f.u₃
    (; ᶜu_bar, ᶜK, ᶜts, ᶜp, ca_phys_params, thermo_dispatcher) = p

    @. ᶜu_bar[colidx] = C123(ᶜuₕ[colidx]) + C123(ᶜinterp(ᶠu₃[colidx]))
    @. ᶜK[colidx] = norm_sqr(ᶜu_bar[colidx]) / 2
    thermo_params = CAP.thermodynamics_params(ca_phys_params)

    CA.debug_state(t, colidx; Y, p) # Debug Y and p state here!

    CA.thermo_state!(Y, p, ᶜinterp, colidx)
    @. ᶜp[colidx] = TD.air_pressure(thermo_params, ᶜts[colidx])
    return nothing
end
```
"""
debug_state(t; kwargs...) = debug_state(t, (); kwargs...)
debug_state(t, colidx::Fields.ColumnIndex; kwargs...) =
    debug_state(t, (colidx,); kwargs...)

function debug_state(t, colidx; Yₜ = nothing, Y = nothing, p = nothing)
    states = Dict()
    states["Yₜ"] = (; msg = String[], name = "Yₜ")
    states["Y"] = (; msg = String[], name = "Y")
    states["p"] = (; msg = String[], name = "p")
    debug_state_generic!(states["Yₜ"], Yₜ, colidx...)
    debug_state_generic!(states["Y"], Y, colidx...)
    debug_state_generic!(states["p"], p, colidx...)
    if !all(map(x -> isempty(states[x].msg), collect(keys(states))))
        for key in keys(states)
            for msg in states[key].msg
                println(msg)
            end
        end
        if colidx == ()
            error("Bad state at time $t")
        else
            error("Bad state at time $t in column $colidx")
        end
    end
    return nothing
end

#####
##### Recursive function for filling auxiliary state with NaNs
#####

function fill_with_nans_generic!(var::Union{Fields.FieldVector, NamedTuple})
    for pn in propertynames(var)
        fill_with_nans_generic!(getproperty(var, pn))
    end
end

function fill_with_nans_generic!(
    state,
    var::Union{Fields.FieldVector, NamedTuple},
    colidx,
)
    for pn in propertynames(var)
        fill_with_nans_generic!(getproperty(var, pn), colidx)
    end
end

fill_with_nans_generic!(var::Number) = nothing
fill_with_nans_generic!(var::AbstractString) = nothing
fill_with_nans_generic!(var::Bool) = nothing
fill_with_nans_generic!(var::Nothing) = nothing
fill_with_nans_generic!(var::Any) = nothing # TODO: should we try to catch more types?

fill_with_nans_generic!(var::Number, colidx) = nothing
fill_with_nans_generic!(var::AbstractString, colidx) = nothing
fill_with_nans_generic!(var::Bool, colidx) = nothing
fill_with_nans_generic!(var::Nothing, colidx) = nothing
fill_with_nans_generic!(var::Any, colidx) = nothing # TODO: should we try to catch more types?

fill_with_nans_generic!(var::Fields.Field) = fill_with_nans_field!(var)

fill_with_nans_field!(::Nothing) = nothing
fill_with_nans_field!(::Nothing, colidx) = nothing
function fill_with_nans_field!(prog::Fields.Field)
    parent(prog) .= NaN
end

"""
    fill_with_nans!(p)

Fill a data structure's `Field`s / `FieldVector`s with NaNs.
"""
fill_with_nans!(p) = fill_with_nans_generic!(p)
