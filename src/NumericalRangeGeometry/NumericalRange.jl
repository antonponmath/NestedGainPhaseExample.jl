include("AbstractNumericalRange.jl")


"""
    NumericalRange = Union{AbstractNumericalRange, Number}

The numerical range is either a proper set or a complex number.
"""
const NumericalRange = Union{AbstractNumericalRange, Number}


"""
    NumericalRange(A::AbstractMatrix{<:Number})

Find the numerical range of a matrix.
"""
NumericalRange(A::AbstractMatrix{<:Number}) = begin
    n, m = size(A)
    n ≠ m && throw(ArgumentError("matrix should be square, got a $n×$m matrix"))
    n ≡ 1 && return A
    n ≡ 2 && return _numerical_range_2x2(A)
    throw(ArgumentError("only 2×2 matrices are supported for now"))
end


"""
    _numerical_range_2x2(A::AbstractMatrix{<:Number})

Compute the numerical range of a 2×2 matrix. It is a point, a segment, or an ellipse.
"""
function _numerical_range_2x2(A::AbstractMatrix{<:Number})
    size(A) ≡ (2, 2) || throw(ArgumentError("matrix must be 2×2"))
    λ1, λ2 = eigvals(A)

    center = (λ1 + λ2) / 2
    Δλ = λ1 - λ2
    majoraxis = (Δλ ≡ 0) ? 1.0 : Δλ / abs(Δλ)        # major axis direction

    trace = real(tr(A' * A))
    a = sqrt(max(0.0, trace - 2 * real(λ1' * λ2))) / 2    # semi-major axis length

    a ≡ 0 && return center

    semimajor = a * majoraxis
    b = sqrt(max(0.0, trace - abs2(λ1) - abs2(λ2))) / 2  # semi-minor axis length

    b ≡ 0 && return NumRangeSegment(center - semimajor, center + semimajor)

    ratio = min(1.0, b / a)

    return NumRangeEllipse(; center, semimajor, ratio)
end
