module ESLIP
using ForwardDiff
import ..SimTools
import ..Model

struct Params
    slip_params::Model.Params
    Kxdot
    KEP
end

function control(x::Vector{T},Ebar::T,p::Params) where T<:Real
    E = Model.energy(x,p.slip_params)
    -p.KEP*(x[1]*x[3]+x[2]*x[4])/sqrt(x[1]^2+x[2]^2)*(E-Ebar)
end

"""
Retuns the stance dynamics accordint to xdot = f(x,u). 
"""
function dynamics(x::Vector{T},Ebar::T,p::Params) where T<:Real
    u = control(x,Ebar,p) 
    return Model.dynamics(x,u,p.slip_params)
end

function touchdown_angle(x::Vector{T},θbar::T,xdot_bar::T,p::Params) where T<:Real
    θbar + p.Kxdot*(x[3]-xdot_bar)
end
"""
Returns the value of the guard for a liftoff event
"""
function guard(x::Vector{T},p::Params) where T<:Real
    return p.slip_params.l0^2-x[1]^2-x[2]^2
end
function reset(x::Vector{T},θbar::T,xdot_bar::T,p::Params) where T<:Real
    θ = touchdown_angle(x,θbar,xdot_bar,p)
    Ey = .5*p.slip_params.m*x[4]^2+p.slip_params.m*Model.g*x[2]
    ptd = [-p.slip_params.l0*sin(θ),p.slip_params.l0*cos(θ)]
    vtd = [x[3],-sqrt(2(Ey-p.slip_params.m*Model.g*ptd[2])/p.slip_params.m)]
    return [ptd...,vtd...]
end

"""
Computes the saltation matrix for a liftoff event
"""
function saltation(x::Vector{T},Ebar::T,θbar::T,xdot_bar::T,p::Params) where T<:Real
    u = stance_control(x,Ebar,p)
    xp = reset(x,θbar,xdot_bar,p)
    up = stance_control(xp,Ebar,p)
    F = Model.dynamics(x,u,p)
    Fp = Model.dynamics(xp,up,p)
    DR = ForwardDiff.jacobian(x->reset(x,eltype(x)(θbar),eltype(x)(xdot_bar),p),x)
    Dg = ForwardDiff.gradient(x->guard(x,p),x)
    return DR + dot(Fp-DR*F,Dg)/dot(Dg,F)
end

function time_bisect(g,ti,tf;tol=1e-6)
    a = ti
    b = tf
    ga = g(a)
    gb = g(b)
    if abs(ga) < tol
        return a
    elseif abs(gb) < tol
        return b
    end
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

function sim(x0::Vector{T},tf::T,h::T,Ebar::T,θbar::T,xdot_bar::T,p::Params) where T<:Real
    f = (x,t)->dynamics(x,Ebar,p)
    g = (x)->guard(x,p)
    r = (x)->reset(x,θbar,xdot_bar,p)

    x = Array{Vector{T}}([])
    t = Vector{T}([])
    events = Array([])

    push!(x,x0)
    push!(t,T(0.))

    while t[end] < tf
        _x = SimTools.heun_step(f,x[end],t[end],h)
        _t = t[end]+h
        if g(_x) < 0
            _t = time_bisect(t->g(SimTools.heun_step(f,x[end],t[end],t-t[end])),t[end],_t;tol=1e-16)
            xm = SimTools.heun_step(f,x[end],t[end],_t-t[end])
            _x = r(xm)
            push!(events,(_x,length(t)+1))
        end
        push!(x,_x)
        push!(t,_t)
    end
    return (x=reduce(hcat,x),t=t,events=events)
end

function sim_to_reset(x0::Vector{T},h::T,Ebar::T,θbar::T,xdot_bar::T,p::Params) where T<:Real
    f = (x,t)->dynamics(x,Ebar,p)
    g = (x)->guard(x,p)
    r = (x)->reset(x,θbar,xdot_bar,p)

    x = Array{Vector{T}}([])
    t = Vector{T}([])
    events = Array([])

    push!(x,x0)
    push!(t,T(0.))

    done = false
    while !done
        _x = SimTools.heun_step(f,x[end],t[end],h)
        _t = t[end]+h
        if g(_x) < 0
            _t = time_bisect(t->g(SimTools.heun_step(f,x[end],t[end],t-t[end])),t[end],_t;tol=1e-16)
            xm = SimTools.heun_step(f,x[end],t[end],_t-t[end])
            _x = r(xm)
            push!(events,(_x,length(t)+1))
            done = true
        end
        push!(x,_x)
        push!(t,_t)
    end
    return (x=reduce(hcat,x),t=t,events=events)
end

function hop(height::T,velocity::T,Ebar::T,θbar::T,xdot_bar::T,p::Params) where T<:Real

    x0 = [T(0.),height,velocity,T(0.)]
    x = Array{Vector{T}}([x0])
    t = Vector{T}([T(0.)])
    events = Array([])

    f = (x,t)->dynamics(x,Ebar,p)
    g = (x)->guard(x,p)
    r = (x)->reset(x,θbar,xdot_bar,p)
    done = false
    while !done
        _x = SimTools.heun_step(f,x[end],t[end],1e-4)
        _t = t[end]+1e-4
        if g(_x) < 0
            _t = time_bisect(t->g(SimTools.heun_step(f,x[end],t[end],t-t[end])),t[end],_t;tol=1e-16)
            xm = SimTools.heun_step(f,x[end],t[end],_t-t[end])
            _x = r(xm)
            push!(events,(xm,_x,length(t)+1))
            done = true
        end
        push!(x,_x)
        push!(t,_t)
    end
    return x=reduce(hcat,x),t=t,events=events
end

end