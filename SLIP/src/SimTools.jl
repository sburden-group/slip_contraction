
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
    Flows x(0)=x0 to x(t) using Heun's method to integrate xdot = f(x,t)
    """
    function flow(f::Function,x0::Vector{T},t::T,h::T) where T<:Real
        τ = T(0.)
        x = x0
        while τ < t
            if τ+h > t 
                h = t-τ
            end
            x = heun_step(f,x,τ,h)
            τ = τ + h
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
            X[:,i+1] = flow(f,X[:,i],t[i+1],t[i+1]-t[i])
        end
        return X
    end
end