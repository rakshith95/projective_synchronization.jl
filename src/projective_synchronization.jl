module projective_synchronization

import StructArrays, Distributions

using LinearAlgebra, StaticArrays, Statistics, SparseArrays, Graphs

#Common
include("common/datatypes.jl")

#Distances
include("distance_measures/distances.jl")

#Averaging
include("Averaging/spherical_mean.jl")
include("Averaging/weiszfeld.jl")
include("Averaging/least_squares.jl")

#Synchroniation
include("synchronization/iterative_synchronization.jl")

#Synthetic
include("synthetic/simulation.jl")

end # module
