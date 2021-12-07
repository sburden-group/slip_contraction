module PeriodicSolutions
import ..Model
import ..SimTools
using NLsolve
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

"""
Integrates the stance dynamics given initial state *IN FOOT FRAME* until stance ends,
then evaluates the stance symmetry boundary conditions.
"""
function shooting_method(x0::Vector{T},t::T,p::Model.Params) where T<:Real
    f(x,t) = begin
        y = Model.stance_reset(x,p)
        Model.stance_dynamics(y,T(0.),p)[[1,2,5,6]]
    end

    # flows the system to time τ from x0 and returns the leg-length guard
    g(τ) = begin
        xf = SimTools.flow(f,x0,τ,T(1e-4))
        return xf[1]^2+xf[2]^2-p.l0
    end

    # solve for the time when stance ends
    result = nlsolve(x->g(x[1]),[t])

    # flow to end of stance
    xf = SimTools.flow(f,x0,result.zero[1],T(1e-4))

    # evauate boundary conditions
    vcat(
        x0[1]+xf[1],
        x0[2]-xf[2],
        x0[3]-xf[3],
        x0[4]+xf[4],
        t-result.zero[1]
    )
end

"""
Searches for a forward running gait given a specified average forward speed during stance,
and a touchdown angle for the leg. Returns the full initial state and the duration of the 
stance period.
"""
function forward_running(speed, θtd, p::Model.Params)
    f(x,t) = begin
        y = Model.stance_reset(x,p)
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
    result = nlsolve(x->shooting_method(x[1:4],x[5],p),vcat(qi,vi,t);ftol=1e-7,xtol=1e-6)
    x = result.zero[1:4]
    t = result.zero[5]
    return x, t
end

end