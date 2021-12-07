using LinearAlgebra
using NLopt
using Printf
import SLIP: Model, PeriodicSolutions

"""
This file demonstrates how the routines in SLIP.Model can be used to
optimize the leg touchdown angle to maximize contraction in the 
direction of energy error.
"""

# nominal trajectory information: state at touchdown
p = Model.default_params
x0,t = PeriodicSolutions.forward_running(1.5,17*pi/180,p)

# Energy state
Ex = p.m*x0[3]^2/2
Ey = (p.m*Model.g*x0[2]+.5*p.m*x0[4]^2)

# The saltation matrix at touchdown
Ξ(x,θ,u) = Model.stance_saltation(x,θ,eltype(x)(0.),u,p)

# Example problem:
# suppose that energy is perturbed
ΔEx = +Ex/2
ΔEy = 0.

# from touchdown angle, determine state at touchdown
function touchdown_velocity(ΔEx,ΔEy,qtd)
    [sqrt(2*(Ex+ΔEx)/p.m),-sqrt(2/p.m*(Ey+ΔEy-p.m*Model.g*qtd[2]))]
end

θtd = pi/8
qtd = [-p.l0*sin(θtd),p.l0*cos(θtd)]
vtd = touchdown_velocity(ΔEx,ΔEy,qtd)
xtd = vcat(qtd,vtd)
error = x0-xtd

# shows about 15% error expansion for this problem default_params
norm(Ξ(xtd,θtd,300.)*error)/norm(error) 

"""
Now we can define an optimization problem which, given a perturbation in energy,
determines the location and touchdown force that maximally contracts the error.
"""

using ForwardDiff
using NLopt

# energy perterbation
ΔEx = Ex/sqrt(2) # equivalent to 50% increase in horizontal speed
ΔEy = Ey/20

# optimization objective
function objective(x::Vector{T}) where T<:Real
    θ = x[1]
    u = x[2]
    qtd = [-p.l0*sin(θ),p.l0*cos(θ)]
    vtd = touchdown_velocity(T(ΔEx),T(ΔEy),qtd)
    state = Array{T}(vcat(qtd,vtd))
    error = Array{T}(x0)-state
    return norm(Ξ(state,θ,u)*error)/norm(error)
end

function f(x::Vector,grad::Vector)
    if length(grad) > 0
        grad[:] = ForwardDiff.gradient(objective,x)
    end
    return objective(x)
end

opt = NLopt.Opt(:LD_SLSQP,2)
opt.min_objective = f
opt.xtol_rel = 1e-8
opt.lower_bounds = [0.,0.]
opt.upper_bounds = [pi/4,p.m*Model.g/2]
initial_guess = [pi/8,300.]
@time (minf,minx,ret) = optimize(opt, initial_guess)
@printf("Result of optimization:\n")
@printf("|Ξδx|/|δ|x: %f, θ: %f\n", minf, minx[1])
