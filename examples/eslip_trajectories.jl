using Revise
using Distributions
using Plots
using ForwardDiff
import SLIP: Model, PeriodicSolutions, ESLIP, SimTools
##
p = ESLIP.Params(Model.default_params,0.25,20.)
##
x0,t = PeriodicSolutions.forward_running(.5,10*pi/180,p.slip_params)
Ebar = Model.flight_energy(x0,p.slip_params)
θbar = atan(-x0[1],x0[2])
xdot_bar = x0[3] # horizontal velocity of nominal gait during flight

# change initial condition to apex state centered at x[1] = 0
x0[1] = 0.
x0[2] += (x0[4]^2)/(2*Model.g)
x0[4] = 0

# create distribution for disturbances at apex height
μ = zeros(4)
σ = [0.,.05,.1,0.]
d = Normal.(μ,σ)

# sample 20 disturbances
N = 5
x = map(i->rand.(d)+x0,1:N)

# simulate disturbances for 10 seconds
sims = map(x->ESLIP.sim(x,2,10.,1e-4,Ebar,θbar,xdot_bar,p),x)

plt = plot(layout=(2,1));
for i = 1:N
    plot!(plt[1],[x[1] for x in sims[i][1]])
    plot!(plt[2],[x[2] for x in sims[i][1]])
end

## now for some contraction stuff
# nominal trajectory
xstar,jstar,t,events = ESLIP.hop(x0[2],x0[3],Ebar,θbar,xdot_bar,p)

stance_trajectory = xstar[events[1][end]:events[2][end]-1]

# too many states, need to downsample
len = length(stance_trajectory)
idx = 1:Int(floor(len/100)):Int(floor(len/100))*100
stance_trajectory = stance_trajectory[idx]

## now that we have the stance trajectory, we must calculate
# the dynamics and jacobian for each point
varidx = [3,4,7,8]
f = map(x->ESLIP.stance_dynamics(x,Ebar,p),stance_trajectory[varidx])
Df = map(x->ForwardDiff.jacobian(x->ESLIP.stance_dynamics(x,eltype(x)(Ebar),p),x),stance_trajectory)
