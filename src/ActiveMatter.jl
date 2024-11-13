module ActiveMatter

using StaticArrays

include("definitions.jl")
include("runsimulation.jl")
include("boundaryconditions.jl")
include("integrators.jl")
include("neighbours.jl")
include("interactions.jl")

end # module ActiveMatter
