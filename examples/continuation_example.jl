using LinearAlgebra
using ForwardDiff
using NLsolve
import SLIP: Model, PeriodicSolutions 
# import .SLIP: Model, PeriodicSolutions

"""
This file demonstrates how to use two functions in SLIP.PeriodicSolutions
to find periodic orbits, and to generate new periodic orbits from given ones.

Periodic gaits can be determined by imposing certain symmetry conditions on the stance period,
namely, that the initial and final horizontal velocities are the same, the initial and final
vertical velocities are reversed, and the initial and final leg angle is reversed.

These symmetries are expressed in a boundary value problem.
New periodic gaits can be found from old ones by exploiting the fact that the 
jacobian of these boundary values is degenerate about solutions.
"""

# use PeriodicSolutions to get a starting running gait
p = Model.default_params
x0,t = PeriodicSolutions.forward_running(1.0,20*pi/180,p)

# evaluate the jacobian of the boundary value problem about this solution
f(y) = PeriodicSolutions.shooting_method(y[1:4],y[5],p)

"""
The next line uses automatic differentiation to compute the jacobian of the boundary
value problem. This could be done with finite differencing, but automatic differentiation
is correct to machine precision.

This is a big reason I use Julia. Differentiating this function is not trivial and 
would likely be impossible using existing Python libraries.
"""
Df(y) = ForwardDiff.jacobian(f,y)

# example, we compute the jacobian of the boundary value problem
y = vcat(x0,t)
A = Df(y)

# the singular vectors of 'A' give directions that we can perturb the 
# initial conditions and get a new periodic solution

U,S,V = svd(A)

# the perterbation dy will result in lower momentum during stance
dy = V[:,end-1] 

# we make sure to correct errors by resolving the boundary value problem
result = nlsolve(f,Df,y+5e-1*dy)
y = result.zero 
# this process of numerical continuation can be repeated
# to generate a family of periodic solutions