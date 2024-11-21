
using DrWatson
@quickactivate :IntroductoryFigures

using CairoMakie, DataFrames, DifferentialEquations, Distributions, PlotFormatting, Random, StochasticTransitionModels, Turing 

function sisrates(u, t)
    s, i = u
    n = s + i 
    return [
        0.48 * s * i / n,  # infection rate
        0.21 * i  # recovery rate
    ]
end

u0 = [ 145, 5 ]

sistransitionmatrix = [
    # s   i
     -1   1  # infection
      1  -1  # recovery
]

Random.seed!(1729)
mockdata = stochasticmodel(sisrates,  u0, 1:30, sistransitionmatrix)

function sis!(du, u, p, t)
    s, i = u 
    beta, gamma = p 
    du[1] = -beta * s * i / (s + i) + gamma * i 
    du[2] = beta * s * i / (s + i) - gamma * i
end 

p = [ 0.1, 0.1 ]

prob = ODEProblem(sis!, u0, ( 1, 30 ), p)
sol = solve(prob; saveat=1)

errordistances = zeros(30, 20)

for (i, b) ∈ enumerate(0.05:0.05:1.5), (j, g) ∈ enumerate(0.05:0.05:1.0)
    p = [ b, g ]
    sol = solve(prob; p, saveat=1)
    ev = 0 
    for t ∈ 1:30 
        ev += (mockdata[t, 2] - sol[t][2])^2
    end
    errordistances[i, j] = ev
end

a, b = findmin(errordistances) 
p = [ (0.05:0.05:1.5)[b[1]], (0.05:0.05:1.0)[b[2]] ]
println(p)
sol = solve(prob; p, saveat=1)
soldf = DataFrame(sol) 


@model function sisfit(
    data, prob;
    betaprior=Uniform(0.05, 1.5), gammaprior=Uniform(0.05, 1),
)
    beta ~ betaprior 
    gamma ~ gammaprior 
    sigma2 ~ Exponential(1)

    p = [ beta, gamma ]
    
    sol = solve(prob; p, saveat=1)

    for i ∈ 1:30
        data[i, 2] ~ Normal(sol[i][2], sigma2)
    end 
end

model = sisfit(mockdata, prob)

chain = sample(model, NUTS(), MCMCThreads(), 1000, 3)
chaindf = DataFrame(chain)

outputmat = zeros(30, 1000)
for j ∈ 1:1000
    p = [ chaindf.beta[j], chaindf.gamma[j] ]
    sol = solve(prob; p, saveat=1)
    for t ∈ 1:30
        outputmat[t, j] = sol[t][2]
    end
end

medians = zeros(30)
lci = zeros(30)
uci = zeros(30)
for t ∈ 1:30 
    lv, mv, uv = quantile(outputmat[t, :], [ 0.05, 0.5, 0.95 ])
    medians[t] = mv
    lci[t] = lv 
    uci[t] = uv 
end



fittingfig = with_theme(theme_latexfonts()) do
    bandcolour = ( COLOURVECTOR[1], 0.5 )
    fig = Figure(; size=( 500, 350 ))
    ga = GridLayout(fig[1, 1])
    ga2 = GridLayout(fig[1, 2])
    gb = GridLayout(fig[2, 1])
    gc = GridLayout(fig[2, 2])
    
    axa = Axis(ga[1, 1])
    cp = contourf!(axa, 0.05:0.05:1.5, 0.05:0.05:1.0, log.(errordistances))
    hlines!(axa, 0.21; color=:red, linestyle=:dot)
    vlines!(axa, 0.48; color=:red, linestyle=:dot)
    cb = Colorbar(ga[1, 2], cp)
    Label(ga[1, 0], L"$\gamma$"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(ga[2, 1], L"$\beta$"; fontsize=11.84, tellwidth=false)
    Label(ga[1, 3], "log square error"; fontsize=11.84, rotation=3π/2, tellheight=false)

    axa2 = Axis(ga2[1, 1])
    lines!(axa2, 1:30, soldf.value2; color=COLOURVECTOR[1])
    scatter!(axa2, 1:30, mockdata[:, 2]; color=:black, markersize=3)
    Label(ga2[1, 0], "prevalence"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(ga2[2, 1], "time"; fontsize=11.84, tellwidth=false)

    axsb1 = [ Axis(gb[i, 1]; xticks=WilkinsonTicks(3), yticks=WilkinsonTicks(3)) for i ∈ 1:2 ]
    axsb2 = [ Axis(gb[i, 3]; xticks=WilkinsonTicks(3), yticks=WilkinsonTicks(3)) for i ∈ 1:2 ]
    for ch ∈ 3:-1:1
        d = filter(:chain => x -> x == ch, chaindf)
        lines!(axsb1[1], d.iteration, d.beta; color=COLOURVECTOR[ch])
        lines!(axsb1[2], d.iteration, d.gamma; color=COLOURVECTOR[ch])
        density!(axsb2[1], d.beta; color=( :white, 0 ), strokecolor=COLOURVECTOR[ch], strokewidth = 1)
        density!(axsb2[2], d.gamma; color=( :white, 0 ), strokecolor=COLOURVECTOR[ch], strokewidth = 1)
    end
    vlines!(axsb2[1], 0.48; color=:red, linestyle=:dot)
    vlines!(axsb2[2], 0.21; color=:red, linestyle=:dot)
    Label(gb[1, 0], L"$\beta$"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(gb[2, 0], L"$\gamma$"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(gb[3, 1], "iteration"; fontsize=11.84, tellwidth=false)
    Label(gb[1:2, 2], "density"; fontsize=11.84, rotation=π/2, tellheight=false)
    #Label(gb[2, 3], L"$\beta$"; fontsize=11.84, tellwidth=false)
    #Label(gb[4, 3], L"$\gamma$"; fontsize=11.84, tellwidth=false)

    axc = Axis(gc[1, 1])
    lines!(axc, 1:30, medians; color=COLOURVECTOR[1])
    band!(axc, 1:30, lci, uci; color=bandcolour)
    scatter!(axc, 1:30, mockdata[:, 2]; color=:black, markersize=3)
    Label(gc[1, 0], "prevalence"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(gc[2, 1], "time"; fontsize=11.84, tellwidth=false)

    formataxis!(axa)
    formataxis!(axa2)
    formataxis!(cb)
    formataxis!(axsb1[1]; hidex=true)
    formataxis!(axsb1[2])
    formataxis!(axsb2)
    formataxis!(axc)

    for c ∈ [ 1, 2, 3 ] colgap!(ga, c, 5) end
    rowgap!(ga, 1, 5)
    colgap!(ga2, 1, 5)
    rowgap!(ga2, 1, 5)
    for c ∈ [ 1, 3 ] colgap!(gb, c, 5) end
    rowgap!(gb, 2, 5)
    colgap!(gc, 1, 5)
    rowgap!(gc, 1, 5)

    labelplots!([ "A", "B", "C", "D" ], [ ga, ga2, gb, gc ]; rows=1)

    fig 
end

safesave(plotsdir("fittingfig.pdf"), fittingfig)
