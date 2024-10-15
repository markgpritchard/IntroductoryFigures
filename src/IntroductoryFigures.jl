
module IntroductoryFigures

using CairoMakie 

include("odemodels.jl")
include("formatplots.jl")

## odemodels.jl
export iirmodel!, irmodel!
## formatplots.jl
export COLOURVECTOR
export formataxis!, labelplots!, setorigin!, setvalue!

end  # module IntroductoryFigures
