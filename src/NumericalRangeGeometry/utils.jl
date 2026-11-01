"""
    boundary(W, n)

Compute `n` points of the boundary of the numerical range `W`.

    boundary(Φ)

Return the endpoints of the phase sector `Φ` as the tuple `(Φ.α, Φ.β)`.
"""
boundary(W::Number, ::Any = nothing) = W

boundary(W::NumRangeSegment, n::Integer = 100) =
    W.z1 + range(0, 1; length = n) * (W.z2 - W.z1)

boundary(W::NumRangeEllipse, n::Integer = 360) = begin
    t = range(0, 2π; length = n + 1)[1:(end - 1)]
    return @. W.center + W.semimajor * complex(cos(t), W.ratio * sin(t))
end

boundary(Φ::PhaseSector) = (Φ.α, Φ.β)


"""
    crawford_number(W)
    crawford_number(A)

Compute the Crawford number of the matrix `A` or numerical range `W`.
"""
crawford_number(W::NumericalRange) =
    0 ∈ W ? 0 : minimum(abs, boundary(W))

crawford_number(A::AbstractMatrix{<:Number}) =
    A |> NumericalRange |> crawford_number


"""
    PhaseSector(W)
    PhaseSector(A)

Compute the phase sector of the matrix `A` or numerical range `W`.
"""
PhaseSectors.PhaseSector(W::NumericalRange) =
    W |> boundary |> PhaseSector

PhaseSectors.PhaseSector(A::AbstractMatrix{<:Number}) =
    A |> NumericalRange |> PhaseSector
