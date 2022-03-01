module PeriodicSolutions
import ..Model
import ..SimTools
using NLsolve
using LinearAlgebra
"""
This file demonstrates how to find periodic passive solutions to SLIP motion
by searching for trajectories with symmetric stance periods.

Periodic solutions are found by specifying the average stance horizontal speed
and leg touchdown angle. These parameters are used edo construct a boundary value
problem, which is approximately solved using a direct method. The solutions
yielded by the direct method are then refined using shooting to the desired precision.

"""

"""
Takes positional positional coordinates of SLIP (IN FOOT FRAME) during stance, and produces
a boundary value problem that enforces symmetry conditions at the boundaries.
"""
function direct_method( f::Function,    # dynamics in form xdot = f(x,t)
                        q::Vector{T},   # positional coordinates, a (2N) vector
                        θtd::T,         # touchdown angle
                        t0::T,          # initial time (for when f(x,t) is time varying)
                        h::T,           # time increment between position measurements
                        p::Model.Params) where T<:Real
    n = Int(length(q)/2)
    x = reshape(q,(n,2))
    d1 = zeros(T,(n-2,n))   # central difference approximation to d/dt
    d2 = zeros(T,(n-2,n))   # central difference approximation to d^2/dt^2
    t = t0.+Array(0:h:(n-1)*h)
    for i=1:n-2
        d1[i,i:i+2] = 1/(2h)*[-1,0,1] # fun fact, this is the mean value theorem
        d2[i,i:i+2] = 1/h^2*[1,-2,1]
    end
    vint = d1*x # approximate velocity at interior points
    aint = d2*x # approximate dynamics at interior points
    fint = reduce(vcat,map(i->f(vcat(x[i+1,:],vint[i,:]),t[i+1])',1:n-2))

    # calculate residuals
    vcat(
        vec(aint)-vec(fint),      # dynamics satisfied
        x[1,1]+p.l0*sin(θtd),     # initial horizontal position
        x[1,2]-p.l0*cos(θtd),     # initial vertical position
        x[end,1]-p.l0*sin(θtd),   # final horizontal position
        x[end,2]-p.l0*cos(θtd)    # final vertical position
    )
end


function shooting_method(θ::T,v::Vector{T},t::T,p::Model.Params) where T<:Real
    x0 = [-p.l0*sin(θ),p.l0*cos(θ),v[1],v[2]]
    xi = Model.stance_reset(x0,θ,p)

    f(x,t) = begin
        Model.stance_dynamics(x,T(0.),p)
    end

    # flows the system to time τ from x0 and returns the leg-length guard
    g(τ) = begin
        xf = SimTools.flow(f,xi,T(0),τ,T(1e-4))
        return sqrt(xf[3]^2+xf[4]^2)-p.l0
    end

    # solve for the time when stance ends
    result = nlsolve(t->g(t[1]),[t])

    # flow to end of stance
    xf = SimTools.flow(f,xi,T(0.),result.zero[1],T(1e-4))

    # # evauate boundary conditions
    vcat(
        θ + atan(-xf[3],xf[4]),
        v[1]-xf[7],
        v[2]+xf[8],
        t-result.zero[1]
    )
end

"""
Searches for a forward running gait given a specified average forward speed during stance,
and a touchdown angle for the leg. Returns the full initial state and the duration of the 
stance period.

Output needs to be sanity checked, since there are some solutions to this BVP that are physically
implausible.
"""
function forward_running(speed, θtd, p::Model.Params)
    f(x,t) = begin
        θ = atan(-x[1],x[2])
        y = Model.stance_reset(x,θ,p)
        Model.stance_dynamics(y,0.,p)[5:6]
    end

    # specify touchdown angle and average stance horizontal speed
    N = 10
    t = 2*p.l0*sin(θtd)/speed

    # initial guess for a solution that is moving forward, and downward at significant speed
    q1 = Array(range(-p.l0*sin(θtd),length=10,stop=p.l0*sin(θtd)))
    q2 = Array(range(p.l0*cos(θtd),length=5,step=-p.l0/10))
    q = hcat(q1,vcat(q2,reverse(q2)))

    # solve nonlinear equation
    result = nlsolve(q->direct_method(f,q,θtd,0.,t/N,p),vec(q))

    # extract initial conditions from this solution
    q = reshape(result.zero,(N,2))
    qi = q[1,:]
    vi = q[1:3,:]'*[-3/2,2,-1/2]*N/t

    # # refine solution using shooting method
    result = nlsolve(x->shooting_method(x[1],x[2:3],x[4],p),vcat(atan(-qi[1],qi[2]),vi,t);ftol=1e-16)
    x = result.zero[1:4]
    θ = result.zero[1]
    v = result.zero[2:3]
    t = result.zero[4]
    x = vcat(-p.l0*sin(θ),p.l0*cos(θ),v)
    return x, t
end


end