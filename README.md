# TransportationPlanningILS

Iterated Local Search for solving large-scale Shipper Transportation Planning Problems.

## Context

This package tackles large-scale shipper transportation planning problems, as defined in our paper

[Optimizing a Worldwide-Scale Shipper Transportation Planning in a Carmaker Inbound Supply Chain](https://arxiv.org/abs/2509.07576)

In our setting, routes can last several weeks, and the instances have 40 industrial sites, 100 platforms, 3000 customers, 700 000 commodities and a 26-week horizon.

## Solution pipeline

1. To create an STPP instance from a folder architecture, we call `read_instance`. This function is designed to browse multiple CSV files and combine them into an Instance object.

2. A lower bound can be computed on the given instance using a mixed-giant container relaxation, with `lower_bound!`.

3. An instance filtering procedure is used to fix some of the decision variables, in order to reduce the size of the instance to be actually optimized. It can be called with `lower_bound_filtering!`.  

4. The constructive heuristic is implemented in `greedy!` and the (large) local search procedure in `local_search!`.

5. The iterated local search can then be applied to improve the initial solution with `ILS!`.

## Reproducing our results

First, you need to clone the repository and open a Julia REPL at its root. Then, run the following commands:

```julia
using Pkg
Pkg.activate(".")
using TransportationPlanningILS
julia_main_test(instanceName="[instance name]", timeLimit=[time limit])
```

The industrial problem instances used in our numerical experiements are available at https://zenodo.org/records/17234091. 
For smooth integration, put the instances downloaded in the data directory of the repository.