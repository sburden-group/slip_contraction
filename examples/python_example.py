
"""
An example of how to use PyJulia to interact with the SLIP module from python.
Note that Julia is compiled, and needs to be started in each Python interpreter instance.

This means that calling a script which executes Julia code will incur significant latency
due to Julia compiling local packages at startup. 

Running this code in an interactive Python session (REPL or Jupyter) will have much improved
latency.
"""

import julia
# run julia.install() if this is the first time

from julia import Pkg # need the package manager to activate this project

# place the relative path to slip_contraction root folder in the variable path
# for instance, if python interpreter is started in slip_contraction, then simply
# have path = "."
path = "."

# this sets the julia environment to this project, and loads all
# of the dependencies in the top-level Project.TOML, including the SLIP module
Pkg.activate(path)

# bring SLIP into scope
from julia import SLIP

# now we can call functions from SLIP

# for instance, this code generates some initial conditions that result in a periodic orbit
p = SLIP.Model.default_params
x0,t = SLIP.PeriodicSolutions.forward_running(1.0,20*3.14/180,p)
print("x0: %s, t: %f\n" % (x0,t))

# the above state corresponds to the flight mode, to evaluate the dynamics in stance,
# a reset needs to occur
x = SLIP.Model.stance_reset(x0,p)

# now we can calculate for instance the dynamics during stance

import numpy as np
td_angle = np.arctan2(-x[3],x[4]) # touchdown angle
u = 1.0 # hypothetical input force
xdot = SLIP.Model.stance_dynamics(x,u,p)

# or look for instance at the saltation matrix at this touchdown event
Xi = SLIP.Model.stance_saltation(x0,td_angle,0.,u,p)
print(Xi)