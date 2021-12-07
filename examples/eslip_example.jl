# using SLIP
import SLIP: Model, PeriodicSolutions, ESLIP

p = ESLIP.Params(Model.default_params,1.,1.)

x0,t = PeriodicSolutions.forward_running(1.0,20*pi/180,p.slip_params)
xdot_bar = x0[5]
Ebar = Model.flight_energy(x0,p.slip_params)
θbar = atan(-x0[3],x0[4])

# create an energy perterbation
δxapex = [0,.25,-.5,0] # a perturbation to the apex state

# Determine touchdown angle
θ = θbar+p.Kxdot*δxapex[3]
# Determine touchdown perturbation 
Ξ=ESLIP.stance_(x0+δx,Ebar,θ,xdot_bar,p)

# # forward simulation of ESLIP during stance
# f2(x,t) = ESLIP.stance_dynamics(x,Ebar,p)
# x = SimTools.flow(f2,x0+δx,t/2,1e-4)

# f2(x,t) = SLIP.stance_dynamics(x,0.,p.slip_params)
# x = SimTools.flow(f2,x0,t/2,1e-4)
