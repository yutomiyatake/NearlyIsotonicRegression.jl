export nearly_isotonic_regression, nearly_isotonic_path
export neariso, neariso_path, neariso_Normal, neariso_path_Normal, neariso_AIC_Normal, neariso_AIC_value_Normal, neariso_Binomial, neariso_path_Binomial, neariso_AIC_Binomial, neariso_AIC_value_Binomial, neariso_Poisson, neariso_path_Poisson, neariso_AIC_Poisson, neariso_AIC_value_Poisson, neariso_Chisq, neariso_path_Chisq, neariso_AIC_Chisq, neariso_AIC_value_Chisq

mutable struct _NearIsoState
    a::Vector{Float64}
    β::Vector{Float64}
    m::Vector{Float64}
    c::Vector{Int}
    λ_tmp::Float64
    t_min::Float64
    idx::Vector{Int}
    K::Int
end

function init(v::_NearIsoState, x::Vector{Float64}, w::Vector{Float64}, n::Int)
    (v.β[1], v.a[1], v.c[1]) = (x[1], w[1], 1)
    v.K = 1
    @inbounds for l = 2:n
        if x[l] == x[l - 1]
            v.a[v.K] = v.a[v.K] + w[l]
            v.c[v.K] = v.c[v.K] + 1
        else
            v.K = v.K + 1
            (v.β[v.K], v.a[v.K], v.c[v.K]) = (x[l], w[l], 1)
        end
    end
    (v.β, v.a, v.c) = (v.β[1:v.K], v.a[1:v.K], v.c[1:v.K])
end

function compute_t_min(v::_NearIsoState)
    s = sign.(max.(v.β[1:v.K - 1] - v.β[2:v.K], zeros(v.K - 1)))
    v.m = ([0;s] - [s;0]) ./ v.a

    if v.K == 1
        v.t_min = Inf
    else
        # compute t_min
        t = (v.β[2:v.K] - v.β[1:v.K - 1]) ./
            (v.m[1:v.K - 1] - v.m[2:v.K]) .+ v.λ_tmp
        candidates = filter(value -> isfinite(value) && value > v.λ_tmp, t)

        if isempty(candidates)
            v.t_min = Inf
            v.idx = Int[]
        else
            v.t_min = minimum(candidates)
            tolerance = sqrt(eps(Float64)) * max(1.0, abs(v.t_min))
            v.idx = findall(value -> isapprox(value, v.t_min; atol=tolerance, rtol=0), t)
        end
    end

end

function update_a(v::_NearIsoState)
    for i = 1:length(v.idx)
        v.a[v.idx[end-i+1]] = v.a[v.idx[end-i+1]] + v.a[v.idx[end-i+1]+1]
    end
    v.a = v.a[setdiff(1:end, v.idx + ones(Int64, length(v.idx)))]
end

function update_c(v::_NearIsoState)
    for i = 1:length(v.idx)
        v.c[v.idx[end-i+1]] = v.c[v.idx[end-i+1]] + v.c[v.idx[end-i+1]+1]
    end
    v.c = v.c[setdiff(1:end, v.idx + ones(Int64, length(v.idx)))]
end

function update_β(v::_NearIsoState)
    v.β = v.β + (v.t_min - v.λ_tmp) * v.m
    v.β = v.β[setdiff(1:end, v.idx + ones(Int64, length(v.idx)))]
end

function update_λ_tmp(v::_NearIsoState)
    v.λ_tmp = v.t_min
end

function update_K(v::_NearIsoState)
    v.K = v.K - length(v.idx)
end

function _expand_fit(v::_NearIsoState)
    output = Vector{Float64}(undef, sum(v.c))
    offset = 0
    for i in 1:v.K
        count = v.c[i]
        fill!(view(output, (offset + 1):(offset + count)), v.β[i])
        offset += count
    end
    return output
end


"""
    neariso(x::Vector, λ, w::Vector) -> y, K

Perform weighted nearly isotonic regression 

# Arguments
- `x`: input vector
- `λ`: parameter
- `w`: weights, default to ones if not provided

# Keywords
- `var`: variances of the Gaussian distributions

# Outputs
- `y`: output vector, which is piecewise monotone
- `K`: the number of clusters

# Algorithm
- modified PAVA 
"""
function neariso(
    x::AbstractVector{<:Real},
    λ::Real,
    w::AbstractVector{<:Real}=ones(length(x));
    var=nothing,
)

    n = length(x)
    _require_same_length(x, w, "weights")
    _require_finite(x, "x")
    _require_positive(w, "w")
    isfinite(λ) && λ >= 0 || throw(DomainError(λ, "λ must be nonnegative and finite"))

    x_work = collect(Float64, x)
    w_work = collect(Float64, w)
    
    if var !== nothing && !isempty(var)
        _require_same_length(x, var, "var")
        _require_positive(var, "var")
        x_work ./= var
        w_work .*= var
    end

    n == 0 && return Float64[], 0


    # if λ == 0
    #     return x, param.K
    # end
        
    param = _NearIsoState(zeros(n), zeros(n), zeros(n), zeros(Int, n), 0.0, 0.0, Int[], 0)
    init(param,x_work,w_work,n)

    if λ == 0
        return x_work, param.K
    end

    while param.λ_tmp < λ
        compute_t_min(param)

        if λ < param.t_min
            param.β = param.β + (λ - param.λ_tmp) * param.m
            return _expand_fit(param), param.K
        end

        update_β(param)
        update_a(param)
        update_c(param)
        update_λ_tmp(param)
        update_K(param)

        if param.K == 1 || λ <= param.λ_tmp
            return _expand_fit(param), param.K
        end
    end

    return _expand_fit(param), param.K
end

"""
    neariso_Normal(x::Vector,λ,variance::Vector) -> y, K

Perform weighted nearly isotonic regression  (Normal)

# Arguments
- `x`: input vector
- `λ`: parameter
- `variance`: weights, default to ones if not provided

# Outputs
- `y`: output vector, which is piecewise monotone
- `K`: the number of clusters

"""
neariso_Normal(x::AbstractVector{<:Real}, λ::Real, variance::AbstractVector{<:Real}) =
    neariso(x, λ, 1 ./ variance)
neariso_Normal(x::AbstractVector{<:Real}, λ::Real) = neariso(x, λ)


"""
    neariso_Binomial(success::Vector, λ, trial::Vector) -> y, K

Perform weighted nearly isotonic regression  (Binomial)

# Arguments
- `success`: input vector consisting of the number of success
- `λ`: parameter
- `trials`: a vector consisting of the number of trials

# Outputs
- `y`: output vector, which is piecewise monotone
- `K`: the number of clusters

"""
function neariso_Binomial(
    success::AbstractVector{<:Real},
    λ::Real,
    trial::AbstractVector{<:Real},
)
    _validate_binomial(success, trial)
    return neariso(float.(success) ./ trial, λ, trial)
end

"""
    neariso_Poisson(x::Vector, λ) -> y, K

Perform weighted nearly isotonic regression  (Poisson)

# Arguments
- `x`: input vector consisting of the number of events
- `λ`: parameter

# Outputs
- `y`: output vector, which is piecewise monotone
- `K`: the number of clusters

"""
function neariso_Poisson(x::AbstractVector{<:Real}, λ::Real)
    _require_nonnegative(x, "x")
    return neariso(x, λ)
end

"""
    neariso_Chisq(x::Vector, λ, d::Vector) -> x, K

Perform weighted nearly isotonic regression  (Chisq)

# Arguments
- `x`: input vector 
- `λ`: parameter
- `d`: a vector consisting of the degrees of freedom

# Outputs
- `y`: output vector, which is piecewise monotone
- `K`: the number of clusters

"""
function neariso_Chisq(
    x::AbstractVector{<:Real},
    λ::Real,
    d::AbstractVector{<:Real},
)
    _validate_chisq(x, d)
    return neariso(2 .* float.(x) ./ d, λ, d ./ 2)
end


"""
    neariso_path(x::Vector, w::Vector) -> knot, K

λ-path of weighted nearly isotonic regression

# Arguments
- `x`: input vector
- `w`: weights, default to ones if not provided

# Outputs
- `knot`: checkpoints of the relaxation parameter
- `K` : the number of clusters at each checkpoint

# Algorithm
modified PAVA 
"""
function neariso_path(
    x::AbstractVector{<:Real},
    w::AbstractVector{<:Real}=ones(length(x)),
)

    n = length(x)
    _require_same_length(x, w, "weights")
    _require_finite(x, "x")
    _require_positive(w, "w")

    n == 0 && return Float64[0.0], Int[0]

    x_work = collect(Float64, x)
    w_work = collect(Float64, w)

    param = _NearIsoState(zeros(n), zeros(n), zeros(n), zeros(Int, n), 0.0, 0.0, Int[], 0)
    knot = Float64[]
    K = Int[]
    init(param,x_work,w_work,n)

    while true
        append!(knot, param.λ_tmp)
        append!(K, param.K)
        compute_t_min(param)

        if param.t_min == Inf
            return knot, K
        end

        update_β(param)
        update_a(param)
        update_c(param)
        update_λ_tmp(param)
        update_K(param)

    end
end

"""
    neariso_path_Normal(x::Vector, variance::Vector) -> knot, K

λ-path of weighted nearly isotonic regression (Normal)

# Arguments
- `x`: input vector
- `variance` a vector consisting of variances, default to ones if not provided

# Outputs
- `knot`: checkpoints of the relaxation parameter
- `K` : the number of clusters at each checkpoint
"""
neariso_path_Normal(x::AbstractVector{<:Real}, variance::AbstractVector{<:Real}) =
    neariso_path(x, 1 ./ variance)
neariso_path_Normal(x::AbstractVector{<:Real}) = neariso_path(x)

"""
    neariso_path_Binomial(success::Vector, trial::Vector) -> knot, K

λ-path of weighted nearly isotonic regression (Binomial)

# Arguments
- `success`: input vector consisting of the number of success
- `trials`: a vector consisting of the number of trials

# Outputs
- `knot`: checkpoints of the relaxation parameter
- `K` : the number of clusters at each checkpoint
"""
function neariso_path_Binomial(
    success::AbstractVector{<:Real},
    trial::AbstractVector{<:Real},
)
    _validate_binomial(success, trial)
    return neariso_path(float.(success) ./ trial, trial)
end


"""
    neariso_path_Poisson(x::Vector) -> knot, K

λ-path of weighted nearly isotonic regression (Poisson)

# Arguments
- `x`: input vector consisting of the number of events

# Outputs
- `knot`: checkpoints of the relaxation parameter
- `K` : the number of clusters at each checkpoint
"""
function neariso_path_Poisson(x::AbstractVector{<:Real})
    _require_nonnegative(x, "x")
    return neariso_path(x)
end


"""
    neariso_path_Chisq(x::Vector, d::Vector) -> knot, K

λ-path of weighted nearly isotonic regression (Chisq)

# Arguments
- `x`: input vector 
- `d`: a vector consisting of the degrees of freedom

# Outputs
- `knot`: checkpoints of the relaxation parameter
- `K` : the number of clusters at each checkpoint
"""
function neariso_path_Chisq(x::AbstractVector{<:Real}, d::AbstractVector{<:Real})
    _validate_chisq(x, d)
    return neariso_path(2 .* float.(x) ./ d, d ./ 2)
end



"""
    neariso_AIC_Normal(x::Vector, variance::Vector) -> knot, AIC

λ-path and AIC of weighted nearly isotonic regression (Normal)

# Arguments
- `x`: input vector
- `variance`: variances, default to ones if not provided

# Outputs
- `knot`: checkpoints of the relaxation parameter
- `AIC` : AIC at each checkpoint
"""
function neariso_AIC_Normal(x::AbstractVector{<:Real}, variance::AbstractVector{<:Real})
    knot, K = neariso_path_Normal(x,variance)
    AIC = 2.0 .* K
    n = length(x)

    for i=1:length(knot)
        η, _ = neariso_Normal(x,knot[i],variance)
        for j=1:n
            AIC[i] += -2 * (-(x[j] - η[j])^2 / (2 * variance[j]) - 1/2 * log(2 * π * variance[j]))
        end
    end

    return knot, AIC
end
neariso_AIC_Normal(x::AbstractVector{<:Real}) = neariso_AIC_Normal(x, ones(length(x)))

"""
    neariso_AIC_value_Normal(x::Vector, λ, variance::Vector) -> AIC

AIC value of weighted nearly isotonic regression (Normal)

# Arguments
- `x`: input vector
- `variance`: variances, default to ones if not provided

# Outputs
- `AIC` : AIC at λ
"""
function neariso_AIC_value_Normal(
    x::AbstractVector{<:Real},
    λ::Real,
    variance::AbstractVector{<:Real},
)
    η, K = neariso_Normal(x,λ,variance)
    AIC = 2*K
    n = length(x)

    for j=1:n
        AIC += -2 * (-(x[j] - η[j])^2 / (2 * variance[j]) - 1/2 * log(2 * π * variance[j]))
    end

    return AIC
end
neariso_AIC_value_Normal(x::AbstractVector{<:Real}, λ::Real) =
    neariso_AIC_value_Normal(x, λ, ones(length(x)))


"""
    neariso_AIC_Binomial(success::Vector, trial::Vector) -> knot, AIC

λ-path and AIC of weighted nearly isotonic regression (Binomial)

# Arguments
- `success`: input vector consisting of the number of success
- `trials`: a vector consisting of the number of trials

# Outputs
- `knot`: checkpoints of the relaxation parameter
- `AIC` : AIC at each checkpoint
"""
function neariso_AIC_Binomial(
    success::AbstractVector{<:Real},
    trial::AbstractVector{<:Real},
)
    knot, K = neariso_path_Binomial(success,trial)
    AIC = 2.0 .* K
    n = length(success)

    for i=1:length(knot)
        η, _ = neariso_Binomial(success,knot[i],trial)
        for j=1:n
            loglikelihood = _xlogy(success[j], η[j]) +
                _xlogy(trial[j] - success[j], 1 - η[j]) +
                log_binomial(trial[j], success[j])
            AIC[i] -= 2 * loglikelihood
        end
    end

    return knot, AIC
end

function log_binomial(n, k)
    return loggamma(n + 1) - loggamma(k + 1) - loggamma(n - k + 1)
end

_xlogy(x, y) = iszero(x) ? 0.0 : x * log(y)


"""
    neariso_AIC_value_Binomial(success::Vector, λ, trial::Vector) -> AIC

AIC value of weighted nearly isotonic regression (Binomial)

# Arguments
- `success`: input vector consisting of the number of success
- `trials`: a vector consisting of the number of trials

# Outputs
- `AIC` : AIC at λ
"""
function neariso_AIC_value_Binomial(
    success::AbstractVector{<:Real},
    λ::Real,
    trial::AbstractVector{<:Real},
)
    η, K = neariso_Binomial(success,λ,trial)
    AIC = 2*K
    n = length(success)

    for j=1:n
        loglikelihood = _xlogy(success[j], η[j]) +
            _xlogy(trial[j] - success[j], 1 - η[j]) +
            log_binomial(trial[j], success[j])
        AIC -= 2 * loglikelihood
    end

    return AIC
end

"""
    neariso_AIC_Poisson(x::Vector) -> knot, AIC

λ-path and AIC of weighted nearly isotonic regression (Poisson)

# Arguments
- `x`: input vector consisting of the number of events

# Outputs
- `knot`: checkpoints of the relaxation parameter
- `AIC` : AIC at each checkpoint
"""
function neariso_AIC_Poisson(x::AbstractVector{<:Real})
    knot, K = neariso_path_Poisson(x)
    AIC = 2.0 .* K
    n = length(x)

    for i=1:length(knot)
        η, _ = neariso_Poisson(x,knot[i])
        for j=1:n
            AIC[i] -= 2 * (_xlogy(x[j], η[j]) - η[j] - loggamma(x[j] + 1))
        end
    end

    return knot, AIC
end

"""
    neariso_AIC_value_Poisson(x::Vector, λ) -> AIC

AIC value of weighted nearly isotonic regression (poisson)

# Arguments
- `x`: input vector consisting of the number of events

# Outputs
- `AIC` : AIC at λ
"""
function neariso_AIC_value_Poisson(x::AbstractVector{<:Real}, λ::Real)
    η, K = neariso_Poisson(x,λ)
    AIC = 2*K
    n = length(x)

    for j=1:n
        AIC -= 2 * (_xlogy(x[j], η[j]) - η[j] - loggamma(x[j] + 1))
    end

    return AIC
end

"""
    neariso_AIC_Chisq(x::Vector, d::Vector) -> knot, AIC

λ-path and AIC of weighted nearly isotonic regression (Chisq)

# Arguments
- `x`: input vector 
- `d`: a vector consisting of the degrees of freedom

# Outputs
- `knot`: checkpoints of the relaxation parameter
- `AIC` : AIC at each checkpoint
"""
function neariso_AIC_Chisq(x::AbstractVector{<:Real}, d::AbstractVector{<:Real})
    _require_positive(x, "x")
    knot, K = neariso_path_Chisq(x,d)
    AIC = 2.0 .* K
    n = length(x)

    for i=1:length(knot)
        η, _ = neariso_Chisq(x,knot[i],d)
        for j=1:n
            loglikelihood = _xlogy(d[j] / 2 - 1, x[j]) - x[j] / η[j] -
                loggamma(d[j] / 2) - d[j] / 2 * log(η[j])
            AIC[i] -= 2 * loglikelihood
        end
    end

    return knot, AIC
end

"""
    neariso_AIC_value_Chisq(x::Vector, λ, d::Vector) -> AIC

AIC value of weighted nearly isotonic regression (Chisq)

# Arguments
- `x`: input vector 
- `d`: a vector consisting of the degrees of freedom

# Outputs
- `AIC` : AIC at λ
"""
function neariso_AIC_value_Chisq(
    x::AbstractVector{<:Real},
    λ::Real,
    d::AbstractVector{<:Real},
)
    _require_positive(x, "x")
    η, K = neariso_Chisq(x,λ,d)
    AIC = 2*K
    n = length(x)

    for j=1:n
        loglikelihood = _xlogy(d[j] / 2 - 1, x[j]) - x[j] / η[j] -
            loggamma(d[j] / 2) - d[j] / 2 * log(η[j])
        AIC -= 2 * loglikelihood
    end

    return AIC
end

"""
    nearly_isotonic_regression(x, λ[, w]; var=nothing) -> y, K

Descriptive alias for [`neariso`](@ref).
"""
function nearly_isotonic_regression(
    x::AbstractVector{<:Real},
    λ::Real,
    w::AbstractVector{<:Real}=ones(length(x));
    var=nothing,
)
    return neariso(x, λ, w; var=var)
end

"""
    nearly_isotonic_path(x[, w]) -> knot, K

Descriptive alias for [`neariso_path`](@ref).
"""
nearly_isotonic_path(x::AbstractVector{<:Real}) = neariso_path(x)
nearly_isotonic_path(x::AbstractVector{<:Real}, w::AbstractVector{<:Real}) =
    neariso_path(x, w)
