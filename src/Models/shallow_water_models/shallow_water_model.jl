"""
    ShallowWaterModel <: AbstractModel
"""
Base.@kwdef struct ShallowWaterModel{DT<:AbstractHorizontalDomain,PT} <: AbstractModel
    domain::DT
    parameters::PT
end
prognostic_state_names(::ShallowWaterModel) = (:h, :u, :c)
diagnostic_state_names(::ShallowWaterModel) = nothing 

"""
    default_initial_conditions(model::ShallowWaterModel)
"""
function default_initial_conditions(model::ShallowWaterModel)
    function_space = make_function_space(model.domain)
    @unpack x1, x2 = Fields.coordinate_field(function_space)

    # function that initilizates the model state locally
    # to zero fields everywhere
    init_func(_...) = (
        prognostic = (
            h = zero(Float64),
            u = Geometry.Cartesian12Vector(zero(Float64), zero(Float64)),
            c = zero(Float64),    
        ),
    )

    return init_func.(x1, x2)
end

"""
    make_ode_function(model::ShallowWaterModel)
"""
function make_ode_function(model::ShallowWaterModel)
    function rhs!(dY, Y, _, t)
        @unpack D₄, g = model.parameters
        Yp = Y.prognostic
        dYp = dY.prognostic

        # function space
        function_space = axes(Y)

        # instantiate operators
        sdiv  = Operators.Divergence()
        wdiv  = Operators.WeakDivergence()
        sgrad = Operators.Gradient()
        wgrad = Operators.WeakGradient()
        scurl = Operators.Curl()
        wcurl = Operators.WeakCurl()

        # compute hyperviscosity first because it requires a direct stiffness summation
        @. dYp.u =
            wgrad(sdiv(Yp.u)) -
            Geometry.Cartesian12Vector(wcurl(Geometry.Covariant3Vector(scurl(Yp.u))))
        Spaces.weighted_dss!(dYp)
        @. dYp.u =
            -D₄ * (
                wgrad(sdiv(dYp.u)) -
                Geometry.Cartesian12Vector(wcurl(Geometry.Covariant3Vector(scurl(dYp.u))))
            )

        # add in advection terms
        J = Fields.Field(function_space.local_geometry.J, function_space)
        @. begin
            dYp.h = -wdiv(Yp.h * Yp.u)
            dYp.u +=
                -sgrad(g * Yp.h + norm(Yp.u)^2 / 2) +
                Geometry.Cartesian12Vector(J * (Yp.u × scurl(Yp.u)))
            dYp.c += -wdiv(Yp.c * Yp.u)
        end
        Spaces.weighted_dss!(dYp)

        return dY
    end

    return rhs!
end