_ensure_either(
    bound_const::AbstractVector{<:Real},
    bound_fn::Function,
    params::AbstractVector{<:Any},
    given::Real
) = minimum(
    [
        bound_fn(given, p...)
            for (c, p) in zip(bound_const, params)
            if given > c
    ]; init = Inf
)


"""
    ensure_either(bound1, bound2, params; given1)

For each `x1 ∈ given1` find the largest `x2` such that for each `i` at least one of the inequalities holds: `x1 ≤ bound1[i]` or `x2 ≤ bound2(x1, params[i]...)[i]`.

Return the vector of tuples `(x1, x2)`.
"""
ensure_either(
    bound1::AbstractVector{<:Real},
    bound2::Function,
    params::AbstractVector{<:Any}
    ;
    given1::AbstractVector{<:Real},
    given2::Nothing = nothing
) = [
    (x1, _ensure_either(bound1, bound2, params, x1)) for x1 in given1
]


"""
    ensure_either(bound1, bound2, params; given2)

For each `x2 ∈ given2` find the largest `x1` such that for each `i` at least one of the inequalities holds: `x1 ≤ bound1(x2, params[i]...)[i]` or `x2 ≤ bound2[i]`.

Return the vector of tuples `(x1, x2)`.
"""
ensure_either(
    bound1::Function,
    bound2::AbstractVector{<:Real},
    params::AbstractVector{<:Any}
    ;
    given1::Nothing = nothing,
    given2::AbstractVector{<:Real}
) = [
    (_ensure_either(bound2, bound1, params, x2), x2) for x2 in given2
]


"""
    xray(θ::ℝ, start::ℂ, direction::ℂ)

Find the intersection of two rays in ℂ:
- the first ray is {μ exp(jθ) : μ ≥ 0};
- the second ray starts at `start` and points in the direction of `direction`.

Returns the smallest μ₂ ≥ 0 for which there exists μ₁ ≥ 0 such that μ₁eʲᶿ = z + μ₂w where z = `start` and w = `direction`. If there is no intersection, returns `Inf`.
"""
function xray(θ::Real, start::Number, direction::Number)
    start ≡ 0 && return 0
    direction ≡ 0 && return Inf

    z, w = (start, direction) .* cis(-θ)
    if imag(w) ≡ 0
        imag(z) ≠ 0 && return Inf
        real(z) ≥ 0 && return 0
        real(w) ≤ 0 && return Inf
        return -real(z) / real(w)
    else
        μ = -imag(z) / imag(w)
        μ < 0 && return Inf
        return real(z) + μ * real(w) ≥ 0 ? μ : Inf
    end
end
