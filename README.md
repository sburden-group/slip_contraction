# About
This project exists to explore applications of infinitesimal contraction theory to a hybrid Spring Loaded Inverted Pendulum (SLIP) system.
It is written in Julia, with the bulk of the code in a local Julia package called SLIP. The root directory of this repsository is itself a Julia Project (environment).

# Requirements
This repository requires that Julia v1.7.0 is installed on the system. 

# Julia Usage
If Julia v1.7.0 is installed, then using this repository is quite simple. For instance, suppose that this repository is located on your filesystem at the relative path "./slip_contraction/".
To begin, start Julia, and execute the following commands in the REPL:

```julia
using Pkg
# the following line makes this project the active julia environment,
# and automatically handles fetching and installing all dependencies
Pkg.develop(path="slip_contraction")

using SLIP
```

Note that the SLIP package should be loaded in *develop* mode, so that local changes to the SLIP package will be reflected when Julia reloads. Alternatively, if you use the [Revise](https://github.com/timholy/Revise.jl) package, then changes to SLIP should be recompiled and loaded into the active Julia session automatically.

# Python compatibility and requirements
Included in this repository is a simple example of how the SLIP package be used within Python. This requires python3 and [PyJulia](https://pyjulia.readthedocs.io/en/latest/) to be installed. See the script examples/python_examples.py for more information.
