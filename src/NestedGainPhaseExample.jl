module NestedGainPhaseExample

export timed, precompute, compute_simple, compute_parallel, compute_inner_unity,
    compute_outer_unity, compute_sampling, plot_results

using ControlSystemsBase: feedback, isstable
using LinearAlgebra: opnorm

using PhaseSectors

include("NumericalRangeGeometry/NumericalRangeGeometry.jl")
using .NumericalRangeGeometry

include("Plotting.jl")
using .Plotting

include("Timing.jl")
using .Timing

include("helpers.jl")


function precompute(; G1, G2, ω)
    _G1 = G1.(im * ω)
    γ1 = _G1 .|> opnorm
    Φ1 = _G1 .|> PhaseSector |> dilate

    _G2 = G2.(im * ω)
    γ2 = _G2 .|> opnorm
    W2 = _G2 .|> NumericalRange
    Φ2 = W2 .|> PhaseSector |> dilate

    G1inv(s) = inv(G1(s))
    W1inv = G1inv.(im * ω) .|> NumericalRange
    c1inv = W1inv .|> crawford_number

    return (; γ1, Φ1, W1inv, c1inv, γ2, Φ2, W2)
end


"""
```
───◒─── μ₁G₁ ──┬──▶
   ̂│           │
   └── I + Δ ──┘
```
with ‖Δ‖ ≤ μ₂‖G₂‖
"""
function compute_simple(; γ1, Φ1, γ2, μ2_range)
    μ2_phase_compliant = map(Φ1, γ2) do Φ1, γ2
        sin(min(π / 2, dist(π, Φ1))) / γ2
    end

    μ1_gain_compliant = (μ2, γ1, γ2) -> 1 / (γ1 * (1 + μ2 * γ2))

    return ensure_either(
        μ1_gain_compliant,
        μ2_phase_compliant,
        tuple.(γ1, γ2);
        given2 = μ2_range
    )
end


"""
```
───◒──── μ₁G₁ ────┬──▶
   ̂│              │
   └── I + μ₂G₂ ──┘
```
"""
function compute_parallel(; γ1, Φ1, γ2, W2, μ2_range)
    μ2_phase_compliant = map(Φ1, W2) do Φ1, W2
        π ∈ Φ1 ? 0 :
            minimum(
                xray.(
                    boundary(π - Φ1),
                    1,
                    transpose(boundary(W2))
                )
            )
    end

    μ1_gain_compliant = (μ2, γ1, γ2) -> 1 / (γ1 * (1 + μ2 * γ2))

    return ensure_either(
        μ1_gain_compliant,
        μ2_phase_compliant,
        tuple.(γ1, γ2);
        given2 = μ2_range
    )
end


"""
```
──◒──◒── μ₁G₁ ──┬──┬─▶
  ̂│  ̂└──────────┘  │
  └───── μ₂G₂ ─────┘
```
"""
function compute_inner_unity(; γ1, Φ1, W1inv, c1inv, γ2, Φ2, μ1_range)
    μ1_inner_compliant = map(Φ1, γ1) do Φ1, γ1
        dist(π, Φ1) > 0 ? Inf : 1 / γ1
    end

    μ1_outer_phase_compliant = map(Φ1, Φ2, W1inv) do Φ1, Φ2, W1inv
        dist(π - Φ2, Φ1) ≡ 0 ? 0.0 :
            minimum(
                xray.(
                    boundary(Φ2 - π),
                    transpose(boundary(W1inv)),
                    1
                )
            )
    end

    γ_inner = (μ1, γ1, c1inv, Φ1) -> begin
        θ = dist(π, Φ1)
        return min(
            μ1 * γ1 ≥ 1 ? Inf :
                1 / (1 / (μ1 * γ1) - 1),
            θ ≡ 0 ? Inf :
                1 / (max(1, c1inv / μ1) * sin(min(π / 2, θ)))
        )
    end

    μ2_outer_gain_compliant = (μ1, γ1, c1inv, Φ1, γ2) ->
    1 / (γ2 * γ_inner(μ1, γ1, c1inv, Φ1))

    return ensure_either(
        μ1_outer_phase_compliant,
        μ2_outer_gain_compliant,
        tuple.(γ1, c1inv, Φ1, γ2);
        given1 = μ1_range[μ1_range .< minimum(μ1_inner_compliant)]
    )
end


"""
```
──◒──◒── μ₁G₁ ──┬──┬─▶
  ̂│  ̂└── μ₂G₂ ──┘  │
  └────────────────┘
```
"""
function compute_outer_unity(; γ1, Φ1, W1inv, c1inv, γ2, Φ2, W2, μ1_bounds, μ1μ2_range)
    μ1μ2_inner_compliant = map(Φ1, Φ2, γ1, γ2) do Φ1, Φ2, γ1, γ2
        dist(π - Φ1, Φ2) > 0 ? Inf : 1 / (γ1 * γ2)
    end

    μ1μ2_outer_phase_compliant = map(Φ1, W1inv, W2) do Φ1, W1inv, W2
        π ∈ Φ1 ? 0.0 : minimum(
                xray.(
                    π,
                    boundary(W1inv),
                    transpose(boundary(W2))
                )
            )
    end

    γ_inner = (μ1μ2, γ1, γ2, c1inv, Φ1, Φ2) -> begin
        θ = dist(π - Φ1, Φ2)
        return min(
            γ1 * μ1μ2 * γ2 ≥ 1 ? Inf :
                1 / (1 / γ1 - μ1μ2 * γ2),
            θ ≡ 0 ? Inf :
                1 / (c1inv * sin(min(π / 2, θ)))
        )
    end

    μ1_outer_gain_compliant = (μ1μ2, γ1, γ2, c1inv, Φ1, Φ2) ->
    1 / γ_inner(μ1μ2, γ1, γ2, c1inv, Φ1, Φ2)

    return [
        (min(μ1, μ1_bounds[2]), μ1μ2 / min(μ1, μ1_bounds[2]))
            for (μ1, μ1μ2) in ensure_either(
                μ1_outer_gain_compliant,
                μ1μ2_outer_phase_compliant,
                tuple.(γ1, γ2, c1inv, Φ1, Φ2);
                given2 = μ1μ2_range[μ1μ2_range .< minimum(μ1μ2_inner_compliant)]
            )
            if μ1 > 0
    ]
end


function compute_sampling(; G1ss, G2ss, μ1_range, μ2_range)
    return (
        μ1_range, μ2_range,
        [
            isstable(feedback(feedback(μ1 * G1ss), μ2 * G2ss))
                for μ1 in μ1_range, μ2 in μ2_range
        ],
    )
end

end # module NestedGainPhaseExample
