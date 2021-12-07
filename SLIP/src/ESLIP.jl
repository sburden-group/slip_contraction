module ESLIP
import ..SimTools
import ..Model

struct Params
    slip_params::Model.Params
    Kxdot
    KEP
end

function stance_control(x::Vector{T},Ebar::T,p::Params) where T<:Real
    E = Model.stance_energy(x,p.slip_params)
    -p.KEP*(x[3]*x[7]+x[2]*x[8])/sqrt(x[3]^2+x[4]^2)*(E-Ebar)
end

"""
Retuns the stance dynamics accordint to xdot = f(x,u). 
"""
function stance_dynamics(x::Vector{T},Ebar::T,p::Params) where T<:Real
    u = stance_control(x,Ebar,p) 
    return Model.stance_dynamics(x,u,p.slip_params)
end

"""
Retuns the flight dynamics accordint to xdot = f(x). 
"""
function flight_dynamics(x::Vector{T}, p::Params) where T<:Real
    return Model.flight_dynamics(x,p.slip_params)
end

function touchdown_angle(x::Vector{T},θbar::T,xdot_bar::T,p::Params) where T<:Real
    θbar + p.Kxdot*(xdot_bar-x[5])
end
"""
Returns the value of the guard for a touchdown event
"""
function stance_guard(x::Vector{T},θbar::T,xdot_bar::T,p::Params) where T<:Real
    θ = touchdown_angle(x,θbar,xdot_bar,p)
    return Model.stance_guard(x,θ,T(0.),p.slip_params)
end

"""
Computes the saltation matrix for a touchdown event
"""
function stance_saltation(x::Vector{T},Ebar::T,θbar::T,xdot_bar::T,p::Params) where T<:Real
    θ = touchdown_angle(x,θbar,xdot_bar,p)
    u = stance_control(x,Ebar,p)
    return Model.stance_saltation(x,θ,T(0.),u,p.slip_params)
end

function flight_guard(x::Vector{T},p::Params) where T<:Real
    Model.flight_guard(x,p.slip_params)
end

function flight_reset(x::Vector{T}, p::Params) where T<:Real
    Model.flight_reset(x, p.slip_params)
end

function flight_saltation(x::Vector{T}, Ebar::T, p::Params) where T<:Real
    u = stance_control(x,Ebar,p)
    model.flight_saltation(x,u,p.slip_params)
end

end