module ESLIP
import ..SimTools
import ..Model

"""
TODO: need to add some control leg so that the touchdown
guard is only evaluated after hitting apex height.
"""

struct Params
    slip_params::Model.Params
    Kxdot
    KEP
end

function stance_control(x::Vector{T},Ebar::T,p::Params) where T<:Real
    E = Model.stance_energy(x,p.slip_params)
    -p.KEP*(x[3]*x[7]+x[4]*x[8])/sqrt(x[3]^2+x[4]^2)*(E-Ebar)
end

"""
Retuns the stance dynamics accordint to xdot = f(x,u). 
"""
function stance_dynamics(x::Vector{T},Ebar::T,p::Params) where T<:Real
    u = stance_control(x,Ebar,p) 
    return Model.stance_dynamics(x,u,p.slip_params)
end

"""
Retuns the flight dynamics accordint to xdot = f(x). 
"""
function flight_dynamics(x::Vector{T}, p::Params) where T<:Real
    return Model.flight_dynamics(x,p.slip_params)
end

function touchdown_angle(x::Vector{T},θbar::T,xdot_bar::T,p::Params) where T<:Real
    θbar + p.Kxdot*(x[3]-xdot_bar)
end
"""
Returns the value of the guard for a touchdown event
"""
function stance_guard(x::Vector{T},θbar::T,xdot_bar::T,p::Params) where T<:Real
    θ = touchdown_angle(x,θbar,xdot_bar,p)
    return Model.stance_guard(x,θ,p.slip_params)
end

function stance_reset(x::Vector{T},θbar::T,xdot_bar::T,p::Params) where T<:Real
    θ = touchdown_angle(x,θbar,xdot_bar,p)
    Model.stance_reset(x,θ,p.slip_params)
end

"""
Computes the saltation matrix for a touchdown event
"""
function stance_saltation(x::Vector{T},Ebar::T,θbar::T,xdot_bar::T,p::Params) where T<:Real
    θ = touchdown_angle(x,θbar,xdot_bar,p)
    u = stance_control(x,Ebar,p)
    return Model.stance_saltation(x,θ,u,p.slip_params)
end

function flight_guard(x::Vector{T},p::Params) where T<:Real
    Model.flight_guard(x,p.slip_params)
end

function flight_reset(x::Vector{T}, p::Params) where T<:Real
    Model.flight_reset(x, p.slip_params)
end

function flight_saltation(x::Vector{T}, Ebar::T, p::Params) where T<:Real
    u = stance_control(x,Ebar,p)
    model.flight_saltation(x,u,p.slip_params)
end

function descent_guard(x::Vector{T}) where T<:Real
    return x[4]
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

function sim(x0::Vector{T},j0::Int,tf::T,h::T,Ebar::T,θbar::T,xdot_bar::T,p::Params) where T<:Real
    f = nothing
    g = nothing
    if j0 == 1
        f = (x,t)->flight_dynamics(x,p)
        g = (x)->descent_guard(x)
        r = (x)->(x,2)
    elseif j0 == 2
        f = (x,t)->flight_dynamics(x,p)
        g = (x)->stance_guard(x,θbar,xdot_bar,p)
        r = (x)->(stance_reset(x,θbar,xdot_bar,p),3)
    elseif j0 == 3
        f = (x,t)->stance_dynamics(x,Ebar,p)
        g = (x)->flight_guard(x,p)
        r = (x)->(flight_reset(x,p),1)
    else
        throw(ArgumentError("Invalid initial discrete state."))
    end

    x = Array{Vector{T}}([])
    j = Array{Int}([])
    t = Vector{T}([])
    events = Array([])

    push!(x,x0)
    push!(j,j0)
    push!(t,T(0.))

    _j = j0
    while t[end] < tf
        _x = SimTools.heun_step(f,x[end],t[end],h)
        _t = t[end]+h
        if g(_x) < 0
            _t = time_bisect(t->g(SimTools.heun_step(f,x[end],t[end],t-t[end])),t[end],_t;tol=1e-16)
            jm = _j
            xm = SimTools.heun_step(f,x[end],t[end],_t-t[end])
            _x,_j = r(xm)
            push!(events,(jm,_j,xm,_x,length(t)+1))
            if _j == 1
                f = (x,t)->flight_dynamics(x,p)
                g = (x)->descent_guard(x)
                r = (x)->(x,2)
            elseif _j == 2
                f = (x,t)->flight_dynamics(x,p)
                g = (x)->stance_guard(x,θbar,xdot_bar,p)
                r = (x)->(stance_reset(x,θbar,xdot_bar,p),3)
            elseif _j == 3
                f = (x,t)->stance_dynamics(x,Ebar,p)
                g = (x)->flight_guard(x,p)
                r = (x)->(flight_reset(x,p),1)
            end
        end
        push!(x,_x)
        push!(j,_j)
        push!(t,_t)
    end
    return x,j,t,events
end

function hop(height::T,velocity::T,Ebar::T,θbar::T,xdot_bar::T,p::Params) where T<:Real

    x0 = [T(0.),height,velocity,T(0.)]
    j0 = 2
    x = Array{Vector{T}}([x0])
    j = Array{Int}([j0])
    t = Vector{T}([T(0.)])
    events = Array([])
    _j = j0

    # sim to touchdown
    f = (x,t)->flight_dynamics(x,p)
    g = (x)->stance_guard(x,θbar,xdot_bar,p)
    r = (x)->(stance_reset(x,θbar,xdot_bar,p),3)
    done = false
    while !done
        _x = SimTools.heun_step(f,x[end],t[end],1e-4)
        _t = t[end]+1e-4
        if g(_x) < 0
            _t = time_bisect(t->g(SimTools.heun_step(f,x[end],t[end],t-t[end])),t[end],_t;tol=1e-16)
            jm = _j
            xm = SimTools.heun_step(f,x[end],t[end],_t-t[end])
            _x,_j = r(xm)
            push!(events,(jm,_j,xm,_x,length(t)+1))
            done = true
        end
        push!(x,_x)
        push!(j,_j)
        push!(t,_t)
    end

    # sim to liftoff
    f = (x,t)->stance_dynamics(x,Ebar,p)
    g = (x)->flight_guard(x,p)
    r = (x)->(flight_reset(x,p),1)
    done = false
    while !done
        _x = SimTools.heun_step(f,x[end],t[end],1e-4)
        _t = t[end]+1e-4
        if g(_x) < 0
            _t = time_bisect(t->g(SimTools.heun_step(f,x[end],t[end],t-t[end])),t[end],_t;tol=1e-16)
            jm = _j
            xm = SimTools.heun_step(f,x[end],t[end],_t-t[end])
            _x,_j = r(xm)
            push!(events,(jm,_j,xm,_x,length(t)+1))
            done = true
        end
        push!(x,_x)
        push!(j,_j)
        push!(t,_t)
    end

    # sim to descent
    f = (x,t)->flight_dynamics(x,p)
    g = (x)->descent_guard(x)
    r = (x)->(x,2)
    done = false
    while !done
        _x = SimTools.heun_step(f,x[end],t[end],1e-4)
        _t = t[end]+1e-4
        if g(_x) < 0
            _t = time_bisect(t->g(SimTools.heun_step(f,x[end],t[end],t-t[end])),t[end],_t;tol=1e-16)
            jm = _j
            xm = SimTools.heun_step(f,x[end],t[end],_t-t[end])
            _x,_j = r(xm)
            push!(events,(jm,_j,xm,_x,length(t)+1))
            done = true
        end
        push!(x,_x)
        push!(j,_j)
        push!(t,_t)
    end

    return x,j,t,events
end
end