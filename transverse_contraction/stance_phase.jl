include("../SLIP.jl")
Model = SLIP.Model
ESLIP = SLIP.ESLIP
PeriodicSolutions = SLIP.PeriodicSolutions

using LinearAlgebra
using ForwardDiff
using Convex
using JuMP
using SCS

# model parameters
p = ESLIP.Params(Model.default_params,0.25,1.0)

# get samples of a stance period
x0,t = PeriodicSolutions.forward_running(.5,10*pi/180,p.slip_params)
Ebar = Model.energy(x0,p.slip_params)
θbar = atan(-x0[1]/x0[2])
xdot_bar = x0[3]

# data is in the named tuple "sim"
# pick a point in mid-stance
sim = ESLIP.sim_to_reset(x0,1e-3,Ebar,θbar,xdot_bar,p)
# x = sim.x[:,Int(floor(length(sim.t)/2))]
x = sim.x[:,50]

## now set up search for a contraction metric, which should be a convex problem...
# the data we need is f(x), Df(x)
f = ESLIP.dynamics(x,Ebar,p)
Df = ForwardDiff.jacobian(x->ESLIP.dynamics(x,eltype(x)(Ebar),p),x)
Q = f*f'
λ = 6

# definition of variables
W = Convex.Semidefinite(4)
ρ = Convex.Variable(1,Positive())
constraints = [W*Df'+Df*W+λ*W-ρ*Q<=-1e-3]
problem = Convex.maximize(0,constraints)
Convex.solve!(problem, Mosek.Optimizer)
problem.status


## something is wrong, let's sanity check these solutions by perturbing theorem
# we will perturb Ebar and θbar
using Random
using Distributions
μ = zeros(2); σ=[0.05,.3]
d = Distributions.Normal.(μ,σ)
y = x0[2]+.5*x0[4]^2/9.81 .+ rand(d[1],3)
xdot = x0[3] .+ rand(d[1],3)

init_conds = map(
    i->begin
        E = p.slip_params.m*(.5*xdot[i]^2+9.81*y[i])
        θ = ESLIP.touchdown_angle([0.,y[i],xdot[i],0.],θbar,xdot_bar,p)
        _x = -p.slip_params.l0*sin(θ)
        _y = p.slip_params.l0*cos(θ)
        [_x,_y,xdot[i],-sqrt(2*(9.81*(y[i]-_y)))] 
    end,
    1:3
)

sims = map(i->ESLIP.sim(init_conds[i],10.0,1e-3,Ebar,θbar,xdot_bar,p),1:3)
plt = plot()
for sim in sims
    plot!(plt,sim.x[2,end-2000:end])
end
plt