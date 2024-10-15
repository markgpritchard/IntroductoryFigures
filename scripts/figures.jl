
using DrWatson
@quickactivate :IntroductoryFigures

using CairoMakie, DifferentialEquations, MakieTeX 
using Makie.Colors

sirintroplot = let 
    sol1 = let 
        u0 = [ 1.0, 0.0 ]
        tspan = ( 0.0, 3.0 )
        p = [ 1 ]
        prob = ODEProblem(irmodel!, u0, tspan, p)
        solve(prob; saveat=0.01)
    end
    
    sol2 = let 
        u0 = [ 1.0, 0.0, 0.0 ]
        tspan = ( 0.0, 3.0 )
        p = [ 1 ]
        prob = ODEProblem(iirmodel!, u0, tspan, p)
        solve(prob; saveat=0.01)
    end
        
    td1 = TeXDocument(read(scriptsdir("tikz_sir.tex"), String))
    td2 = TeXDocument(read(scriptsdir("tikz_siir.tex"), String))

    fig = with_theme(theme_latexfonts()) do

        fig = Figure(; size=( 500, 350 ))

        ga = GridLayout(fig[1, 0:2])
        lt1 = LTeX(ga[1, 1], td1; tellwidth=false)
        lt2 = LTeX(ga[1, 2], td2; tellwidth=false)

        ax1a = Axis(fig[2, 1])
        ax1b = Axis(fig[3, 1])
        ax2a = Axis(fig[2, 2])
        ax2b = Axis(fig[3, 2])

        lines!(
            ax1a, sol1.t, [ sol1.u[t][1] for t ∈ eachindex(sol1.t) ]; 
            color=COLOURVECTOR[1]
        )
        lines!(
            ax1b, sol1.t, [ sol1.u[t][1] for t ∈ eachindex(sol1.t) ]; 
            color=COLOURVECTOR[1]
        )
        lines!(
            ax2a, sol2.t, [ sol2.u[t][1] + sol2.u[t][2] for t ∈ eachindex(sol2.t) ]; 
            color=COLOURVECTOR[1]
        )
        lines!(
            ax2b, sol2.t, 2 .* [ sol2.u[t][2] for t ∈ eachindex(sol2.t) ]; 
            color=COLOURVECTOR[1]
        )
        for ax ∈ [ ax1a, ax1b, ax2a, ax2b ]
            vlines!(ax, 1; color=:black, linestyle=:dot)
        end

        linkaxes!(ax1a, ax1b, ax2a, ax2b)
        formataxis!(ax1a; hidex=true)
        formataxis!(ax1b)
        formataxis!(ax2a; hidex=true, hidey=true)
        formataxis!(ax2b; hidey=true)

        Label(
            fig[2, 0], L"Proportion\\infectious$$"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(
            fig[3, 0], L"Generation \\ interval$$"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(fig[4, 1:2], L"Time, multiple of $\gamma$"; fontsize=11.84, tellwidth=false)

        colgap!(fig.layout, 1, 5)
        rowgap!(fig.layout, 3, 5)

        labelplots!(
            [ "A", "B", "C", "D" ], 
            [ ga, ga, fig.layout, fig.layout ]; 
            cols=[ 1, 2, 0, 2 ], rows=[ 1, 1, 2, 2 ]
        )
        fig
    end
    fig
end

safesave(plotsdir("sirintroplot.pdf"), sirintroplot)

