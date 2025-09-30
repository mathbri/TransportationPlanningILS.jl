module TransportationPlanningILS

# Packages needed across the project

using Graphs
# using Graphs.LinAlg
using MetaGraphsNext
using CSV
using IterTools
using Random
using JuMP
using Gurobi
using HiGHS
# using CPLEX
using SparseArrays
# using Flux
# using InferOpt
using OhMyThreads
using Dates
# using ProfileView
using JLD2
using Statistics
using Logging

# Project files

include("constants.jl")
include("utils.jl")

# Structures 
include("Structures/Network.jl")
include("Structures/Commodity.jl")
include("Structures/Order.jl")
include("Structures/Bundle.jl")
include("Structures/Bin.jl")
include("Structures/TravelTime.jl")
include("Structures/TimeSpace.jl")
include("Structures/projectors.jl")
include("Structures/struct_utils.jl")
include("Structures/Instance.jl")
include("Structures/Solution.jl")
include("Structures/Perturbation.jl")

# Import and Export of data
include("Reading/read_instance.jl")
include("Reading/read_solution.jl")
include("Writing/write_solution.jl")

# Algorithms
include("Algorithms/bin_packing.jl")
include("Algorithms/bin_updating.jl")
include("Algorithms/solution_updating.jl")
# Greedy heuristic
include("Algorithms/Utils/greedy_utils.jl")
include("Algorithms/greedy.jl")
# Lower Bound computation and heuristic
include("Algorithms/Utils/lb_utils.jl")
include("Algorithms/lower_bound.jl")
# Local search heuristic
include("Algorithms/Utils/ls_utils.jl")
include("Algorithms/local_search.jl")
# Large Neighborhood Search
include("Algorithms/LNS/lns_utils.jl")
include("Algorithms/LNS/lns_milp.jl")
include("Algorithms/LNS/lns.jl")
# Miscellaneous heuristics (benchmarks and comparisons)
include("Algorithms/miscellaneous.jl")

include("run.jl")
include("julia_main.jl")

# Functions to be made public

export julia_main_test, julia_main_logged

end
