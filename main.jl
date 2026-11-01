using ControlSystemsBase: tf, pade, ss
using NestedGainPhaseExample


## SYSTEM DEFINITION

s = tf("s")

g(ω) = ω^2 / (s^2 + 0.2ω * s + ω^2)
h(ω) = ω / (s + ω)
d(τ) = pade(τ, 2)

function construct_system(p1, p2, f1, f2)
    P = [p1 0; 0 p2]
    F = [1 f1; f2 1]
    G = P * F
    Gss = ss(P) * ss(F)
    return G, Gss
end

G1, G1ss = construct_system(
    0.1 * g(10) * d(0.9),
    0.07 * g(20) * d(1),
    0.9 * h(1),
    0.8 * h(2)
)

G2, G2ss = construct_system(
    1.0 * g(2) * h(0.4),
    0.7 * g(1) * h(0.2),
    0.2,
    0.3
)


## SAMPLING SETUP

ω_bounds = (1.0e-2, 1.0e3)
ω = logrange(ω_bounds...; length = 1000)

Ns = 500

μ1_bounds = (10^(-0.1), 10^(1.1))
μ1_range = logrange(μ1_bounds...; length = Ns)

μ2_bounds = (10^(-0.2), 10^(1.1))
μ2_range = logrange(μ2_bounds...; length = Ns)

μ1μ2_range = logrange((μ1_bounds .* μ2_bounds)...; length = Ns)


## COMPUTE

(; γ1, Φ1, W1inv, c1inv, γ2, Φ2, W2) = timed(:PREPARE) do
    precompute(; G1, G2, ω)
end

dom_simple = timed(:SIMPLE) do
    compute_simple(; γ1, Φ1, γ2, μ2_range)
end

dom_parallel = timed(:PARALLEL) do
    compute_parallel(; γ1, Φ1, γ2, W2, μ2_range)
end

dom_inner_unity = timed(:INNER_UNITY) do
    compute_inner_unity(; γ1, Φ1, W1inv, c1inv, γ2, Φ2, μ1_range)
end

dom_outer_unity = timed(:OUTER_UNITY) do
    compute_outer_unity(; γ1, Φ1, W1inv, c1inv, γ2, Φ2, W2, μ1_bounds, μ1μ2_range)
end

dom_sampling = timed(:SAMPLING) do
    compute_sampling(; G1ss, G2ss, μ1_range, μ2_range)
end


## DRAWING

plot_results(;
    ω, γ1, γ2, Φ1, Φ2,
    μ1_bounds, μ2_bounds,
    dom_simple, dom_parallel, dom_inner_unity, dom_outer_unity, dom_sampling
)
