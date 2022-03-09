
module SimTools
    """
    Computes an explicit Heun step to integrate the equation xdot = f(x,t)
    """
    function heun_step(f::Function,x::Vector{T},t::T,h::T) where T<:Real
        x_euler = x + h*f(x,t)
        x_heun = x + h/2*(f(x_euler,t+h)+f(x,t))
        return x_heun
    end
    """
    Flows x(ti)=x0 to x(tf) using Heun's method to integrate xdot = f(x,t)
    """
    function flow(f::Function,x0::Vector{T},ti::T,tf::T,h::T) where T<:Real
        t = ti
        x = x0
        while t < tf
            if t+h > tf 
                h = tf-t
            end
            x = heun_step(f,x,t,h)
            t = t + h
        end
        return x
    end
    """
    Flows x(t[1])=x0 to x(t[end]) using Heun's method to integrate xdot = f(x,t)
    """
    function flow(f::Function,x0::Vector{T},t::Vector{T}) where T<:Real
        X = zeros(T,(length(x0),length(t)))
        X[:,1] = x0
        for i=1:length(t)-1
            X[:,i+1] = flow(f,X[:,i],t[i],t[i+1],t[i+1]-t[i])
        end
        return X
    end
end