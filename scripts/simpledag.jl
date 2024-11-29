
using DrWatson

@quickactivate :IntroductoryFigures

using CairoMakie, MakieTeX 

simpledag = let 
    td1 = TeXDocument(read(scriptsdir("simpledag.tex"), String))
    fig = Figure(; size=( 500, 150 ))
    ga = GridLayout(fig[1, 1:2])
    lt1 = LTeX(ga[1, 1], td1; tellwidth=false)

    fig
end

safesave(plotsdir("simpledag.pdf"), simpledag)

