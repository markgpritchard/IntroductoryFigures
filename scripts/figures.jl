
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

        ga = GridLayout(fig[1, 1:2])
        gb = GridLayout(fig[2, 1])
        gc = GridLayout(fig[2, 2])

        lt1 = LTeX(ga[1, 1], td1; tellwidth=false)
        lt2 = LTeX(ga[1, 2], td2; tellwidth=false)

        ax1a = Axis(gb[1, 1])
        ax1b = Axis(gb[2, 1])
        ax2a = Axis(gc[1, 1])
        ax2b = Axis(gc[2, 1])

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
        formataxis!(ax2a; hidex=true)
        formataxis!(ax2b)

        for gl ∈ [ gb, gc ]
            Label(
                gl[1, 0], L"Proportion\\infectious$$"; 
                fontsize=11.84, rotation=π/2, tellheight=false
            )
            Label(
                gl[2, 0], L"Generation \\ interval$$"; 
                fontsize=11.84, rotation=π/2, tellheight=false
            )
            Label(gl[3, 1], L"Time, multiple of $\gamma$"; fontsize=11.84, tellwidth=false)
    
            colgap!(gl, 1, 5)
            rowgap!(gl, 2, 5)
        end


        labelplots!(
            [ "A", "B", "C", "D" ], [ ga, ga, gb, gc ]; 
            cols=[ 1, 2, 0, 0 ], rows=1
        )
        fig
    end
    fig
end

safesave(plotsdir("sirintroplot.pdf"), sirintroplot)
















using CSV, DataFrames, Dates

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# UK-wide data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df = CSV.read(
    datadir("exp_raw", "OxCGRT_fullwithnotes_GBR_v1.csv"), DataFrame
)

regioncodedf = DataFrame(
    :RegionCode => [ "UK_ENG", "UK_NIR", "UK_SCO", "UK_WAL" ],
    :RegionId => 1:4
)

# Government Response Index includes C1, C2, C3, C4, C5, C6, C7, C8, E1, E2, H1, H2, H3,
# H6, H7, H8 

# format date as `Date`
rename!(
    df, 
    Dict(
        :Date => "_date",  # rename old version of date 
        Symbol("C1E_School closing") => "C1E_Schoolclosing",
        Symbol("C2E_Workplace closing") => "C2E_Workplaceclosing",
        Symbol("C3E_Cancel public events") => "C3E_Cancelpublicevents",
        Symbol("C4E_Restrictions on gatherings") => "C4E_Restrictionsongatherings",
        Symbol("C5E_Close public transport") => "C5E_Closepublictransport",
        Symbol("C6E_Stay at home requirements") => "C6E_Stayathome",
        Symbol("C7E_Restrictions on internal movement") => 
            "C7E_Restrictionsoninternalmovement",
        Symbol("C8E_International travel controls") => "C8E_Internationaltravelcontrols",
        Symbol("E1_Income support") => "E1E_Incomesupport",
        Symbol("E2_Debt/contract relief") => "E2E_Debtcontractrelief",
        Symbol("H1_Public information campaigns") => "H1E_Publicinformationcampaigns",   
        Symbol("H2_Testing policy") => "H2E_Testingpolicy",
        Symbol("H3_Contact tracing") => "H3E_Contacttracing",
        Symbol("H6E_Facial Coverings") => "H6E_FacialCoverings",
        Symbol("H7_Vaccination policy") => "H7E_Vaccinationpolicy",
        Symbol("H8E_Protection of elderly people") => "H8E_Protectionofelderlypeople",
    )
) 
insertcols!(df, :Date => [ Date("$d", "yyyymmdd") for d ∈ df._date ])

# choose necessary columns 
select!(
    df,
    :RegionName,
    :RegionCode,
    :Jurisdiction,
    :Date,
    :C1E_Schoolclosing,
    :C2E_Workplaceclosing,
    :C3E_Cancelpublicevents,
    :C4E_Restrictionsongatherings,
    :C5E_Closepublictransport,
    :C6E_Stayathome,
    :C7E_Restrictionsoninternalmovement,
    :C8E_Internationaltravelcontrols,
    :E1E_Incomesupport,
    :E2E_Debtcontractrelief,
    :H1E_Publicinformationcampaigns,
    :H2E_Testingpolicy,
    :H3E_Contacttracing,
    :H6E_FacialCoverings,
    :H7E_Vaccinationpolicy,
    :H8E_Protectionofelderlypeople,
    :ConfirmedCases,
    :ConfirmedDeaths,
    :StringencyIndex_WeightedAverage,
    :GovernmentResponseIndex_WeightedAverage,
)

# H6E_FacialCoverings is missing for each nation on 2020-09-24, but == 3 for UK as a whole 
for i ∈ axes(df, 1)
    if df.Date[i] == Date("2020-09-24") df.H6E_FacialCoverings[i] = 3 end
end

# convert cumulative cases to daily new cases 
_calcdiff(a::Real, b::Real) = b - a 
_calcdiff(::Missing, b::Real) = b 
_calcdiff(::Any, ::Missing) = 0  # includes first case for each country

insertcols!(
    df,
    :NewConfirmedCases => [ 
        i <=7 ? 0 : max(_calcdiff(df.ConfirmedCases[i-7], df.ConfirmedCases[i]), 0) 
        for i ∈ axes(df, 1) 
    ],
)

# choose only "state-level" data 
filter!(:Jurisdiction => x -> x == "STATE_TOTAL", df)

leftjoin!(df, regioncodedf; on=:RegionCode)

df


## Plot cases 

casesplot = let 
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 500, 350 ))
        axs = [ 
            Axis(
                fig[i, 1]; 
                xticks=(
                    Dates.value.(
                        Date.(
                            [ "2020-01-01", "2021-01-01", "2022-01-01", "2023-01-01" ]
                        ) .- Date("2020-01-01")
                    ),
                    [ "2020", "2021", "2022", "2023" ]
                ),
                yticks=WilkinsonTicks(4)
            ) 
            for i ∈ 1:4 
        ]

        for region ∈ 1:4 
            inds1 = findall(x -> x == region, df.RegionId)
            inds2 = findall(x -> x == 7, dayofweek.(df.Date))
            inds = inds1[findall(x -> x ∈ inds2, inds1)]
            lines!(
                axs[region], 
                Dates.value.(df.Date[inds] .- Date("2020-01-01")), 
                df.NewConfirmedCases[inds] ./ 1000; 
                color=COLOURVECTOR[region],
            )            
            text!(
                axs[region], 0, maximum(df.NewConfirmedCases[inds] ./ 1000); 
                text=[ "England", "Northern Ireland", "Scotland", "Wales" ][region],
                align=( :left, :top ),
                fontsize=11.84
            )
            formataxis!(axs[region]; hidex=(region != 4))
        end

        Label(
            fig[1:4, 0], "Weekly recorded incidence, thousands"; 
            fontsize=11.84, rotation=π/2, tellheight=false,
        )
        Label(fig[5, 1], "Date"; fontsize=11.84, tellwidth=false)

        colgap!(fig.layout, 1, 5)
        rowgap!(fig.layout, 4, 5)
            
        fig
    end
    fig
end

safesave(plotsdir("casesplot.pdf"), casesplot)

function calcv(vec, mv)
    return [ 
        ismissing(v) || v == 0 ? 
            0 : 
            v == mv ? 
                0.05 : 
                0.025
        for v ∈ vec 
    ]
end

function plotinds(vec)
    plotinds = Int[ ]
    for (i, v) ∈ enumerate(vec)
        if v > 0 
            push!(plotinds, i)
        elseif i >= 2 && vec[i-1] > 0 
            push!(plotinds, i) 
        elseif i <= length(vec) - 1 && vec[i+1] > 0 
            push!(plotinds, i)  
        end
    end
    return plotinds 
end
 

interventionsplot = let 
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 500, 650 ))
        ga = GridLayout(fig[1, 1])
        ax = Axis(
            ga[1, 1];
            xticks=(
                Dates.value.(
                    Date.(
                        [ "2020-01-01", "2021-01-01", "2022-01-01", "2023-01-01" ]
                    ) .- Date("2020-01-01")
                ),
                [ "2020", "2021", "2022", "2023" ]
            ),
            yticks=(
                -0.5:-1:-15.5,
                [
                    "C1, School closure",
                    "C2, Workplace closure",
                    "C3, Cancel events",
                    "C4, Restrict gatherings",
                    "C5, Close public transport",
                    "C6, Stay at home",
                    "C7, Restrict internal movement",
                    "C8, International travel controls",
                    "E1, Income support",
                    "E2, Debt/contract relief",
                    "H1, Public information",
                    "H2, Testing policy",
                    "H3, Contact tracing",
                    "H6, Facial coverings",
                    "H7, Vaccination policy",
                    "H8, Protection of elderly people",
                ]
            )
        )

        for region ∈ 1:4 
            inds = findall(x -> x == region, df.RegionId)

            for (var, varmax, d) ∈ zip(
                [
                    :C1E_Schoolclosing,
                    :C2E_Workplaceclosing,
                    :C3E_Cancelpublicevents,
                    :C4E_Restrictionsongatherings,
                    :C5E_Closepublictransport,
                    :C6E_Stayathome,
                    :C7E_Restrictionsoninternalmovement,
                    :C8E_Internationaltravelcontrols,
                    :E1E_Incomesupport,
                    :E2E_Debtcontractrelief,
                    :H1E_Publicinformationcampaigns,
                    :H2E_Testingpolicy,
                    :H3E_Contacttracing,
                    :H6E_FacialCoverings,
                    :H7E_Vaccinationpolicy,
                    :H8E_Protectionofelderlypeople,
                ],
                [
                    3,  # :C1E_Schoolclosing,
                    3,  # :C2E_Workplaceclosing,
                    2,  # :C3E_Cancelpublicevents,
                    4,  # :C4E_Restrictionsongatherings,
                    2,  # :C5E_Closepublictransport,
                    3,  # :C6E_Stayathome,
                    2,  # :C7E_Restrictionsoninternalmovement,
                    4,  # :C8E_Internationaltravelcontrols,
                    2,  # :E1E_Incomesupport,
                    2,  # :E2E_Debtcontractrelief,
                    2,  # :H1E_Publicinformationcampaigns,
                    3,  # :H2E_Testingpolicy,
                    2,  # :H3E_Contacttracing,
                    4,  # :H6E_FacialCoverings,
                    5,  # :H7E_Vaccinationpolicy,
                    3,  # :H8E_Protectionofelderlypeople,
                ],
                -1:-1:-16
            )

                v = calcv(getproperty(df, var)[inds], varmax)
                pinds = plotinds(v)

                band!(
                    ax,
                    Dates.value.(df.Date[inds] .- Date("2020-01-01")),
                    (d + 0.875 - 0.15 * region) .+ v,
                    (d + 0.875 - 0.15 * region) .- v;
                    color=( COLOURVECTOR[region], 0.5 )
                )
                for sv = [ -v, v ]
                    lines!(
                        ax,
                        Dates.value.(df.Date[inds] .- Date("2020-01-01")),
                        [ 
                            i ∉ pinds ? 
                                missing : 
                                (d + 0.875 - 0.15 * region) + sv[i] 
                            for i ∈ eachindex(sv) 
                        ],
                        color=COLOURVECTOR[region]
                    )
                end
            end


        end

        formataxis!(ax)
        ylims!(ax, -15.9, 0)

        Label(fig[2, 1], "Date"; fontsize=11.84, tellwidth=false)

        leg = Legend(
            fig[0, 1],
            [ LineElement(color=COLOURVECTOR[i]) for i ∈ 1:4 ],
            [ "England", "Northern Ireland", "Scotland", "Wales" ]
        )

        formataxis!(leg)

        for r ∈ 1:2 rowgap!(fig.layout, r, 5) end
            
        fig
    end
    fig
end

safesave(plotsdir("interventionsplot.pdf"), interventionsplot)
