module projective_synchronization

import StructArrays, Distributions, MAT, PlotlyJS, StatsBase, Random

using LinearAlgebra, StaticArrays, Statistics, SparseArrays, Graphs, Arpack, ProgressBars, TiledIteration

#Common
include("common/datatypes.jl")
include("common/losses.jl")

#Distances
include("distance_measures/distances.jl")

#Averaging
include("Averaging/spherical_mean.jl")
include("Averaging/weiszfeld.jl")
include("Averaging/least_squares.jl")

#Synchroniation
include("synchronization/iterative_synchronization.jl")
include("synchronization/global_synchronization.jl")
include("synchronization/weighted_synchronization.jl")
# include("synchronization/synchronization_matlab.jl")

#Synthetic
include("synthetic/simulation.jl")
include("synthetic/experiments.jl")
include("synthetic/averaging_experiments.jl")

#Plots
include("plots/plot_simulation.jl")

end # module
