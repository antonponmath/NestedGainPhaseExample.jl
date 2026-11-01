"""
    AbstractNumericalRange

Numerical range when it is a proper set (not a single point).
"""
abstract type AbstractNumericalRange end


"""
    NumRangeSegment <: AbstractNumericalRange

Numerical range when it is a segment.

    NumRangeSegment(z1, z2)

Create a segment with endpoints z1 and z2.
"""
struct NumRangeSegment{T <: Real} <: AbstractNumericalRange
    z1::Complex{T}
    z2::Complex{T}
    function NumRangeSegment{T}(z1::Number, z2::Number) where {T <: Real}
        z1, z2 = (z1, z2) .|> Complex{T}
        (isfinite(z1) && isfinite(z2) && z1 ≠ z2) ||
            throw(ArgumentError("endpoints must be finite and not coincident, got $z1 and $z2"))
        return new(z1, z2)
    end
end

NumRangeSegment(z1::Number, z2::Number) = begin
    T = promote_type(((z1, z2) .|> typeof .|> real)...)
    return NumRangeSegment{T}(z1, z2)
end


"""
    NumRangeEllipse <: AbstractNumericalRange

Numerical range when it is an ellipse.

    NumRangeEllipse(; center, semimajor, ratio)

Create an ellipse:
- `center ∈ ℂ` is the center;
- `semimajor ∈ ℂ` is the semi-major axis;
- `ratio ∈ ℝ` is the ratio of the minor axis to the major one.
"""
struct NumRangeEllipse{T <: Real} <: AbstractNumericalRange
    center::Complex{T}
    semimajor::Complex{T}
    ratio::T
    function NumRangeEllipse{T}(;
            center::Number, semimajor::Number, ratio::Real
        ) where {T <: Real}
        center, semimajor = (center, semimajor) .|> Complex{T}
        ratio = T(ratio)
        isfinite(center) ||
            throw(ArgumentError("center point must be finite"))
        (isfinite(semimajor) && semimajor ≠ 0) ||
            throw(ArgumentError("semi-major axis must be finite and not zero, got $semimajor"))
        (isfinite(ratio) && 0 < ratio ≤ 1) ||
            throw(ArgumentError("axes ratio must be in (0, 1], got $ratio"))
        return new(center, semimajor, ratio)
    end
end

NumRangeEllipse(; center::Number, semimajor::Number, ratio::Real) = begin
    T = promote_type(((center, semimajor, ratio) .|> typeof .|> real)...)
    return NumRangeEllipse{T}(; center, semimajor, ratio)
end


function Base.in(z::Number, W::NumRangeSegment)
    return isreal((z - W.z1) / (W.z2 - W.z1))
end

function Base.in(z::Number, W::NumRangeEllipse)
    z = (z - W.center) / W.semimajor
    z = complex(real(z), imag(z) / W.ratio)
    return abs(z) ≤ 1
end
