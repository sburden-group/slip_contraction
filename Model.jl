module Model
using LinearAlgebra
using ForwardDiff
const g = 9.81  # gravity in m/s^2
const n = 2     # number dimension of SLIP configuration manifold

struct Params
    m::Float64 # mass
    k::Float64 # spring stiffness
    l0::Float64 # resting leg length
end
const default_params = Params(74.5,4000,1.0)

"""
coordinates
x[1] = body foot frame horizontal
x[2] = body foot frame vertical
x[3:4] = d/dt x[1:4]
"""

function energy(x::Vector{T},p::Params) where T<:Real
    .5*p.m*(x[3]^2+x[4]^2)+p.m*g*x[2]+.5*p.k*(sqrt(x[1]^2+x[2]^2)-p.l0)^2
end

"""
Retuns the stance dynamics accordint to xdot = f(x,u). 
"""
function dynamics(x::Vector{T},u::T,p::Params) where T<:Real
    xdot = zeros(T,2n)
    xdot[1:n] = x[n+1:2n]
    l = sqrt(x[1]^2+x[2]^2)
    xdot[3] = (p.k*(p.l0-l)+u)*x[1]/l/p.m
    xdot[4] = (p.k*(p.l0-l)+u)*x[2]/l/p.m - g
    return xdot
end


end