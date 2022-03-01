module Model
using LinearAlgebra
using ForwardDiff
const g = 9.81  # gravity in m/s^2
const n_stance = 4     # number dimension of SLIP configuration manifold
const n_flight = 2     # number dimension of SLIP configuration manifold

struct Params
    m::Float64 # mass
    k::Float64 # spring stiffness
    l0::Float64 # resting leg length
end
const default_params = Params(74.5,4000,1.0)

"""
Stance coordinates
x[1] = body world frame horizontal
x[2] = body world frame vertical
x[3] = body foot frame horizontal
x[4] = body foot frame vertical
x[5:8] = d/dt x[1:4]
"""

"""
Flight coordinates
x[1] = body world frame horizontal
x[2] = body world frame vertical
x[3:4] = d/dt x[1:2]
"""
function stance_energy(x::Vector{T},p::Params) where T<:Real
    .5*p.m*(x[5]^2+x[6]^2)+p.m*g*x[2]+.5*p.k*(sqrt(x[3]^2+x[4]^2)-p.l0)^2
end

"""
Retuns the stance dynamics accordint to xdot = f(x,u). 
"""
function stance_dynamics(x::Vector{T},u::T,p::Params) where T<:Real
    xdot = zeros(T,2n_stance)
    xdot[1:n_stance] = x[n_stance+1:2n_stance]
    l = sqrt(x[3]^2+x[4]^2)
    xdot[5] = xdot[7] = (p.k*(p.l0-l)+u)*x[3]/l/p.m
    xdot[6] = xdot[8] = (p.k*(p.l0-l)+u)*x[4]/l/p.m - g
    return xdot
end

function flight_energy(x::Vector{T},p::Params) where T<:Real
    .5*p.m*(x[3]^2+x[4]^2)+p.m*g*x[2]
end

"""
Retuns the flight dynamics accordint to xdot = f(x). 
"""
function flight_dynamics(x::Vector{T}, p::Params) where T<:Real
    xdot = zeros(T,2n_flight)
    xdot[1:2] = x[3:4]
    xdot[end] = -g
    return xdot
end

"""
Returns the value of the guard for a touchdown event
"""
function stance_guard(x::Vector{T},θ::T,p::Params) where T<:Real
    return x[2]-p.l0*cos(θ)
end

function stance_reset(x::Vector{T},θ::T,p::Params) where T<:Real
    xnew = zeros(T,2n_stance)
    xnew[1:2] = x[1:2]
    xnew[3] = -p.l0*sin(θ)
    xnew[4] = x[2]
    xnew[5:6] = xnew[7:8] = x[3:4]
    return xnew
end

"""
Computes the saltation matrix for a touchdown event
"""
function stance_saltation(x::Vector{T},θ::T,u::T,p::Params) where T<:Real
    f1 = flight_dynamics(x,p)
    R(x) = stance_reset(x,θ,p)
    f2 = stance_dynamics(R(x),u,p)
    DR = ForwardDiff.jacobian(R,x)
    g(x) = stance_guard(x,eltype(x)(θ),p)
    Dg = ForwardDiff.gradient(g,x)
    return  DR - ((f2 - DR*f1)*Dg')/(Dg'*f1)
end


function flight_guard(x::Vector{T},p::Params) where T<:Real
    p.l0-x[3]^2-x[4]^2
end

function flight_reset(x::Vector{T},p::Params) where T<:Real
    xnew = zeros(T,4)
    xnew[1:2] = x[1:2]
    xnew[3:4] = x[5:6]
    return xnew
end

function flight_saltation(x::Vector{T},u::T,p::Params) where T<:Real
    f1 = stance_dynamics(x,u,p)
    R(x) = flight_reset(x,p)
    f2 = flight_dynamics(R(x),p)
    DR = ForwardDiff.jacobian(R,x)
    g(x) = flight_guard(x,p)
    Dg = ForwardDiff.gradient(g,x)
    return  DR - ((f2 - DR*f1)*Dg')/(Dg'*f1)
end


end