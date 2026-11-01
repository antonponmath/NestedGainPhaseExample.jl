module Plotting

export plot_results

using CairoMakie, TeXLayout

using PhaseSectors: PhaseSector, sectors!


const _COL = Dict(
    :G1_LINE => :firebrick3,
    :G2_LINE => :steelblue,
    :G1_FILL => (:firebrick3, 0.1),
    :G2_FILL => (:steelblue, 0.1),
    :SIMPLE => :black,
    :PARALLEL => :gray40,
    :INNER_UNITY => :deepskyblue3,
    :OUTER_UNITY => :firebrick1,
    :STABLE_FILL => :peachpuff,
    :UNSTABLE_FILL => :cadetblue4,
    :STABLE_TEXT => :chocolate4,
    :UNSTABLE_TEXT => :lightcyan2,
)

function plot_results(;
        ω::AbstractVector{<:Real},
        γ1::AbstractVector{<:Real},
        γ2::AbstractVector{<:Real},
        Φ1::AbstractVector{<:PhaseSector},
        Φ2::AbstractVector{<:PhaseSector},
        μ1_bounds::Tuple{Real, Real},
        μ2_bounds::Tuple{Real, Real},
        dom_simple::AbstractVector{<:Tuple{Real, Real}},
        dom_parallel::AbstractVector{<:Tuple{Real, Real}},
        dom_inner_unity::AbstractVector{<:Tuple{Real, Real}},
        dom_outer_unity::AbstractVector{<:Tuple{Real, Real}},
        dom_sampling::Tuple{AbstractVector{<:Real}, AbstractVector{<:Real}, AbstractMatrix{Bool}}
    )
    ω_bounds = (ω[1], ω[end])

    fonts_latex, fonts_makie = _tac_fonts()
    set_default_font_family!(fonts_latex)
    set_theme!(fonts = fonts_makie)

    fig = Figure()

    ax_gain = Axis(
        fig[1, 1], xscale = log10, yscale = log10,
        limits = (ω_bounds, (10^(-2), 10^(0.1))),
        title = "Gain",
        xlabel = L"\omega \text{ (rad/s)}",
        yticks = LogTicks(-2:0),
        yminorticks = IntervalsBetween(10), yminorgridvisible = true,
        xminorticks = IntervalsBetween(10), xminorgridvisible = true
    )
    ax_phase = Axis(
        fig[2, 1], xscale = log10,
        limits = (ω_bounds, (-1.1π, 1.1π)),
        title = "Phase Sector",
        xlabel = L"\omega \text{ (rad/s)}",
        yticks = (π * (-1:1), [L"[π]", L"[0]", L"[π]"]),
        xminorticks = IntervalsBetween(10), xminorgridvisible = true,
        yminorticks = IntervalsBetween(4), yminorgridvisible = true
    )
    ax_dom = Axis(
        fig[3, 1], xscale = log10, yscale = log10,
        limits = (μ1_bounds, μ2_bounds),
        title = "Stability Domain",
        xlabel = L"\mu_1 \text{ (log scale)}",
        ylabel = L"\mu_2 \text{ (log scale)}",
        xlabelpadding = -8,
        ylabelpadding = -8,
        xticks = logrange(1.0e0, 1.0e1; length = 2),
        yticks = logrange(1.0e0, 1.0e1; length = 2),
        xminorticks = IntervalsBetween(10), xminorgridvisible = true,
        yminorticks = IntervalsBetween(10), yminorgridvisible = true
    )

    ## GAIN

    lines!(
        ax_gain, ω, γ1;
        color = _COL[:G1_LINE],
        label = L"\Vert \mathrm\mathbf{G}_1 \Vert"
    )
    lines!(
        ax_gain, ω, γ2;
        color = _COL[:G2_LINE],
        label = L"\Vert \mathrm\mathbf{G}_2 \Vert"
    )

    axislegend(ax_gain; position = :lb)

    ## PHASE

    sectors!(
        ax_phase, ω, Φ1;
        linecolor = _COL[:G1_LINE],
        fillcolor = _COL[:G1_FILL]
    )
    sectors!(
        ax_phase, ω, Φ2;
        linecolor = _COL[:G2_LINE],
        fillcolor = _COL[:G2_FILL]
    )
    hlines!(
        ax_phase, [-π, π];
        color = :black,
        linewidth = 0.5
    )

    axislegend(
        ax_phase, [
            PolyElement(color = _COL[:G1_FILL], strokecolor = _COL[:G1_LINE], strokewidth = 2),
            PolyElement(color = _COL[:G2_FILL], strokecolor = _COL[:G2_LINE], strokewidth = 2),
        ], [
            L"\Phi(\mathrm\mathbf{G}_1)",
            L"\Phi(\mathrm\mathbf{G}_2)",
        ];
        position = :lt
    )

    ## DOMAIN

    bg_dom = heatmap!(
        ax_dom, dom_sampling...;
        colormap = [_COL[:UNSTABLE_FILL], _COL[:STABLE_FILL]],
        rasterize = 5
    )
    translate!(bg_dom, 0, 0, -100)
    text!(
        [tuple.(μ1_bounds, μ2_bounds)...];
        offset = [(25, 20), (-25, -20)],
        text = ["stable", "unstable"],
        font = :sansitalic,
        fontsize = 14,
        color = [_COL[:STABLE_TEXT], _COL[:UNSTABLE_TEXT]],
        align = [(:left, :bottom), (:right, :top)]
    )

    clip = (x, bounds) -> max(0.9 * bounds[1], min(x, 1.1 * bounds[2]))
    clip_μ = μ -> (clip(μ[1], μ1_bounds), clip(μ[2], μ2_bounds))

    lines_dom_simple = lines!(
        ax_dom, clip_μ.(dom_simple);
        color = _COL[:SIMPLE],
        linestyle = :dot
    )
    lines_dom_parallel = lines!(
        ax_dom, clip_μ.(dom_parallel);
        color = _COL[:PARALLEL]
    )
    lines_dom_inner_unity = lines!(
        ax_dom, clip_μ.(dom_inner_unity);
        color = _COL[:INNER_UNITY]
    )
    μ1, μ2 = dom_outer_unity[end]
    lines_dom_outer_unity = lines!(
        ax_dom, [dom_outer_unity; (μ1_bounds[1], μ1 * μ2 / μ1_bounds[1])];
        color = _COL[:OUTER_UNITY]
    )

    Legend(
        fig[4, 1], [
            lines_dom_simple,
            lines_dom_parallel,
            lines_dom_inner_unity,
            lines_dom_outer_unity,
        ], [
            L"\mu_1 \mathrm\mathbf{G}_1 \:⇌\: (\mathrm{I} + \mathrm\mathbf{Δ}) \text{ with } \Vert\mathrm\mathbf{Δ}\Vert \leq \mu_2 \Vert\mathrm\mathbf{G}_2\Vert",
            L"\mu_1 \mathrm\mathbf{G}_1 \:⇌\: (\mathrm{I} + \mu_2 \mathrm\mathbf{G}_2)",
            L"(\mu_1 \mathrm\mathbf{G}_1 \:⇌\: \mathrm{I}) \:⇌\: \mu_2 \mathrm\mathbf{G}_2",
            L"(\mu_1 \mathrm\mathbf{G}_1 \:⇌\: \mu_2 \mathrm\mathbf{G}_2) \:⇌\: \mathrm{I}",
        ],
        framevisible = false
    )

    rowsize!(fig.layout, 1, 130)
    rowsize!(fig.layout, 2, 160)
    rowsize!(fig.layout, 3, 180)
    rowsize!(fig.layout, 4, 70)
    colsize!(fig.layout, 1, 250)

    resize_to_layout!(fig)
    resize!(fig, 303, size(fig.scene)[2])

    return fig
end

function _tac_fonts()
    set_default_font_family!(:new_cm)
    cm = default_font_family().math

    set_default_font_family!(:termes)
    ff = default_font_family()

    fonts_latex = font_family(
        cm,
        regular = ff.regular,
        bold = ff.bold,
        italic = ff.italic,
        sans = ff.sans,
        monospace = ff.monospace
    )

    fonts_makie = (
        regular = ff.regular,
        bold = ff.bold,
        bolditalic = ff.bolditalic,
        italic = ff.italic,
        sans = ff.sans.regular,
        sansbold = ff.sans.bold,
        sansitalic = ff.sans.italic,
        sansbolditalic = ff.sans.bolditalic,
    )

    return fonts_latex, fonts_makie
end

end # module Plotting
