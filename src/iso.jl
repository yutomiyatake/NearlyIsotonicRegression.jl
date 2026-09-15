export isotonic_regression!, isotonic_regression
export iso!, iso, iso_Normal!, iso_Normal, iso_Binomial!, iso_Binomial, iso_Poisson!, iso_Poisson, iso_Chisq!, iso_Chisq

function _require_same_length(x::AbstractVector, y::AbstractVector, y_name::AbstractString)
    length(x) == length(y) ||
        throw(DimensionMismatch("Lengths of input vector and $(y_name) mismatch"))
    return nothing
end

function _require_finite(x::AbstractVector, name::AbstractString)
    all(isfinite, x) || throw(DomainError(x, "Every element of $(name) must be finite"))
    return nothing
end

function _require_positive(x::AbstractVector, name::AbstractString)
    all(value -> isfinite(value) && value > 0, x) ||
        throw(DomainError(x, "Every element of $(name) must be positive and finite"))
    return nothing
end

function _require_nonnegative(x::AbstractVector, name::AbstractString)
    all(value -> isfinite(value) && value >= 0, x) ||
        throw(DomainError(x, "Every element of $(name) must be nonnegative and finite"))
    return nothing
end

"""
    iso!(x::Vector, w::Union{Vector, Nothing}=nothing) -> x

Perform isotonic regression 

# Arguments
- `x`: input vector 
- `w`: weights, default to ones if not provided

# Outputs
- `x`: output vector, which is monotone

# Algorithm
- PAVA (Pool-Adjacent-Violators algorithm)
"""
function iso!(x::AbstractVector{<:AbstractFloat}, w::Union{Nothing,AbstractVector{<:Real}}=nothing)
    
    n = length(x)

    # Basic checks
    _require_finite(x, "x")
    if w !== nothing
        _require_same_length(x, w, "weights")
        _require_positive(w, "w")
    end

    if n <= 1
        return x
    end

  
    @inbounds begin
        n -= 1
        while true
            i = 1
            pooled = 0
            while i <= n
                k = i
                while k <= n && x[k] >= x[k+1]
                    k += 1
                end

                if x[i] != x[k]
                    numerator = 0.0
                    denominator = 0.0
                    for j in i : k
                        wj = w === nothing ? 1.0 : w[j]
                        numerator += x[j] * wj
                        denominator += wj
                    end

                    for j in i : k
                        x[j] = numerator / denominator
                    end
                    pooled = 1
                end
                i = k + 1
            end
            if pooled == 0
                break
            end
        end
    end
    return x
end

"""
    iso(x::Vector) -> y

Same as ```iso!``` with ```w=Nothing```, but allocates an output vector ```y```.
"""
iso(x::AbstractVector{<:Real}) = iso!(collect(float.(x)))

"""
    iso(x::Vector, w::Vector) -> y

Same as ```iso!```, but allocates an output vector ```y```.
"""
iso(x::AbstractVector{<:Real}, w::AbstractVector{<:Real}) = iso!(collect(float.(x)), w)


"""
    iso_Normal!(x::Vector, variance::Vector) -> x

Perform isotonic regression (Normal)

# Arguments
- `x`: input vector 
- `variance` a vector consisting of variances, default to ones if not provided

# Outputs
- `x`: output vector, which is monotone

"""
iso_Normal!(x::AbstractVector{<:AbstractFloat}, variance::AbstractVector{<:Real}) = iso!(x, 1 ./ variance)
iso_Normal!(x::AbstractVector{<:AbstractFloat}) = iso!(x)


"""
    iso_Normal(x::Vector, variance::Vector) -> y

Same as ```iso_Normal!```, but allocates an output vector ```y```.
"""
iso_Normal(x::AbstractVector{<:Real}, variance::AbstractVector{<:Real}) = iso(x, 1 ./ variance)
iso_Normal(x::AbstractVector{<:Real}) = iso(x)

function _validate_binomial(success::AbstractVector{<:Real}, trial::AbstractVector{<:Real})
    _require_same_length(success, trial, "trials")
    _require_positive(trial, "trial")
    _require_nonnegative(success, "success")
    all(isinteger, trial) || throw(DomainError(trial, "Every element of trial must be an integer"))
    all(isinteger, success) ||
        throw(DomainError(success, "Every element of success must be an integer"))
    all(success .<= trial) ||
        throw(DomainError(success, "Every element of success must not exceed trial"))
    return nothing
end

"""
    iso_Binomial!(success::Vector, trial::Vector) -> success

Perform isotonic regression (Binomial)

# Arguments
- `success`: input vector consisting of the number of success
- `trials`: a vector consisting of the number of trials

# Outputs
- `success`: output vector, which is monotone

"""
function iso_Binomial!(success::AbstractVector{<:AbstractFloat}, trial::AbstractVector{<:Real})
    _validate_binomial(success, trial)
    success ./= trial
    return iso!(success, trial)
end

"""
    iso_Binomial(success::Vector, trial::Vector) -> x

Same as ```iso_Binomial!```, but allocates an output vector ```x```.
"""
function iso_Binomial(success::AbstractVector{<:Real}, trial::AbstractVector{<:Real})
    _validate_binomial(success, trial)
    return iso(float.(success) ./ trial, trial)
end

"""
    iso_Poisson!(x::Vector) -> x

Perform isotonic regression (Poisson)

# Arguments
- `x`: input vector consisting of the number of events

# Outputs
- `x`: output vector, which is monotone

"""
function iso_Poisson!(x::AbstractVector{<:AbstractFloat})
    _require_nonnegative(x, "x")
    return iso!(x)
end

"""
    iso_Poisson(x::Vector) -> y

Same as ```iso_Poisson!```, but allocates an output vector ```y```.
"""
function iso_Poisson(x::AbstractVector{<:Real})
    _require_nonnegative(x, "x")
    return iso(x)
end

function _validate_chisq(x::AbstractVector{<:Real}, d::AbstractVector{<:Real})
    _require_same_length(x, d, "degrees of freedom")
    _require_nonnegative(x, "x")
    _require_positive(d, "d")
    return nothing
end

"""
    iso_Chisq!(x::Vector, d::Vector) -> x

Perform isotonic regression (Chisq)

# Arguments
- `x`: input vector 
- `d`: a vector consisting of the degrees of freedom

# Outputs
- `x`: output vector, which is monotone

"""
function iso_Chisq!(x::AbstractVector{<:AbstractFloat}, d::AbstractVector{<:Real})
    _validate_chisq(x, d)
    x .= 2 .* x ./ d
    return iso!(x, d ./ 2)
end

"""
    iso_Chisq(x::Vector, d::Vector) -> y

Same as ```iso_Chisq!```, but allocates an output vector ```y```.
"""
function iso_Chisq(x::AbstractVector{<:Real}, d::AbstractVector{<:Real})
    _validate_chisq(x, d)
    return iso(2 .* float.(x) ./ d, d ./ 2)
end

"""
    isotonic_regression!(x[, w]) -> x

Descriptive alias for [`iso!`](@ref). The input `x` is modified in place.
"""
isotonic_regression!(x::AbstractVector{<:AbstractFloat}) = iso!(x)
isotonic_regression!(x::AbstractVector{<:AbstractFloat}, w::AbstractVector{<:Real}) =
    iso!(x, w)

"""
    isotonic_regression(x[, w]) -> y

Descriptive alias for [`iso`](@ref). A new vector is returned.
"""
isotonic_regression(x::AbstractVector{<:Real}) = iso(x)
isotonic_regression(x::AbstractVector{<:Real}, w::AbstractVector{<:Real}) = iso(x, w)
