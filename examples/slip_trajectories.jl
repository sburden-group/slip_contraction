using Revise
import SLIP: Model, PeriodicSolutions, SimTools

p = Model.default_params

x0,t = PeriodicSolutions.forward_running(.5,10*pi/180,p)
θ = atan(-x0[1],x0[2])

# change initial condition to apex state centered at x[1] = 0
x0[1] = 0.
x0[2] += (x0[4]^2)/(2*Model.g)
x0[4] = 0.

# need code to simulate this system through multiple cycles

# function sim(x0::Vector{T},t::T,p::ESLIP.Params) where T<:Real

#     θtd = ESLIP.touchdown_angle(x,θbar,x0[3],p)
#     ytd = p.slip_params.l0 * cos(θtd)
#     t_reset = sqrt(2*x0[2]/Model.g)
#     tvec = Vector(0:1e-3:t_reset)
#     x_flight = SimTools.flow((x,t)->ESLIP.flight_dynamics(x,p),x0,tvec)
#     x_stance_reset = ESLIP.stance_reset(x_flight[:,end],p)
#     f(t) = SimTools
# end

function time_bisect(g,ti,tf;tol=1e-6)
    a = ti
    b = tf
    ga = g(a)
    gb = g(b)
    @assert(sign(ga) == -sign(gb))
    sign_a = sign(ga)
    c = (a+b)/2
    gc = g(c)
    count = 0
    while abs(gc)>tol
        count += 1
        if count > 100
            break
        end
        if sign_a*gc > 0
            a = c
        else
            b = c
        end
        c = (a+b)/2
        gc = g(c)
    end
    return c
end

# Do I just want to write a hybrid simulation method?
# Pass a time vector and initial state
# return simulation state at those times, including pre/post reset states
#TODO capture hybrid events in a separate data structure
function sim(x0::Vector{T},ti::T,tf::T,h::T) where T<:Real
    x = Array{Vector{T}}([])    # empty array to store state
    t = Vector{T}([])
    push!(x,x0)
    push!(t,ti)
    j = 1 # discrete state, 1=flight, 2=stance
    if length(x0) == 8
        j = 2
    end
    f = (x,t)->nothing
    g = (x)->nothing
    if j == 1
        f = (x,t)->Model.flight_dynamics(x,p)
        g = (x)->Model.stance_guard(x,eltype(x)(θ),p)
    else
        f = (x,t)->Model.stance_dynamics(x,p)
        g = (x)->Model.flight_guard(x,eltype(x)(θ),p)
    end

    while t[end] < tf
        if t[end]+h > tf
            h = tf-t[end]
        end
        _x = SimTools.heun_step(f,x[end],t[end],h)
        _t = t[end]+h
        if g(_x) < 0
            print("FOO")
            c = time_bisect(τ->g(SimTools.heun_step(f,x[end],t[end],τ)),0.,h)
            _x = SimTools.heun_step(f,x[end],t[end],c)
            _t = t[end]+c
            if j == 1
                j = 2
                _x = Model.stance_reset(_x,θ,p)
                f = (x,t)->Model.stance_dynamics(x,eltype(x)(0.),p)
                g = (x)->Model.flight_guard(x,p)
            else
                j = 1 
                _x = Model.flight_reset(_x,p)
                f = (x,t)->Model.flight_dynamics(x,p)
                g = (x)->Model.stance_guard(x,eltype(x)(θ),p)
            end
        end
        push!(x,_x)
        push!(t,_t)
    end
    return x,t
end