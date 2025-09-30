# Compute path and cost for the greedy insertion of a bundle, not handling path admissibility
function greedy_path(
    solution::Solution,
    TTGraph::TravelTimeGraph,
    TSGraph::TimeSpaceGraph,
    bundle::Bundle,
    src::Int,
    dst::Int,
    CAPACITIES::Vector{Int};
    sorted::Bool=false,
    use_bins::Bool=true,
    opening_factor::Float64=1.0,
    current_cost::Bool=false,
    findSources::Bool=true,
)
    update_cost_matrix!(
        solution,
        TTGraph,
        TSGraph,
        bundle,
        CAPACITIES;
        sorted=sorted,
        use_bins=use_bins,
        opening_factor=opening_factor,
        current_cost=current_cost,
        findSources=findSources,
    )
    dijkstraState = dijkstra_shortest_paths(TTGraph.graph, src, TTGraph.costMatrix)
    shortestPath = enumerate_paths(dijkstraState, dst)
    removedCost = remove_shortcuts!(shortestPath, TTGraph)
    pathCost = dijkstraState.dists[dst]
    return shortestPath, pathCost - removedCost
end

function parallel_greedy_path(
    solution::Solution,
    TTGraph::TravelTimeGraph,
    TSGraph::TimeSpaceGraph,
    bundle::Bundle,
    src::Int,
    dst::Int,
    CHANNEL::Channel{Vector{Int}},
    use_bins::Bool=true,
    opening_factor::Float64=1.0,
)
    parallel_update_cost_matrix!(
        solution, TTGraph, TSGraph, bundle, CHANNEL, true, use_bins, opening_factor
    )
    dijkstraState = dijkstra_shortest_paths(TTGraph.graph, src, TTGraph.costMatrix)
    shortestPath = enumerate_paths(dijkstraState, dst)
    removedCost = remove_shortcuts!(shortestPath, TTGraph)
    pathCost = dijkstraState.dists[dst]
    return shortestPath, pathCost - removedCost
end

# Compute the path and cost for the greedy insertion of a bundle, handling path admissibility
function parallel_greedy_insertion(
    solution::Solution,
    TTGraph::TravelTimeGraph,
    TSGraph::TimeSpaceGraph,
    bundle::Bundle,
    src::Int,
    dst::Int,
    CHANNEL::Channel{Vector{Int}};
    verbose::Bool=false,
)
    verbose && println("Greedy insertion for bundle $bundle between src $src and dst $dst")
    verbose && println("Bundle arcs : $(TTGraph.bundleArcs[bundle.idx])")
    shortestPath, pathCost = parallel_greedy_path(
        solution, TTGraph, TSGraph, bundle, src, dst, CHANNEL, true, 1.0
    )
    verbose && println("Initial path : $shortestPath for cost $pathCost")
    # If the path is not admissible, re-computing it
    if !is_path_admissible(TTGraph, shortestPath)
        costMatrix = deepcopy(TTGraph.costMatrix)
        shortestPath, pathCost = parallel_greedy_path(
            solution, TTGraph, TSGraph, bundle, src, dst, CHANNEL, true, 0.5
        )
        if !is_path_admissible(TTGraph, shortestPath)
            shortestPath, pathCost = parallel_greedy_path(
                solution, TTGraph, TSGraph, bundle, src, dst, CHANNEL, false, 1.0
            )
            if !is_path_admissible(TTGraph, shortestPath)
                for (aSrc, aDst) in TTGraph.bundleArcs[bundle.idx]
                    TTGraph.costMatrix[aSrc, aDst] =
                        TTGraph.networkArcs[aSrc, aDst].distance
                end
                dijkstraState = dijkstra_shortest_paths(
                    TTGraph.graph, src, TTGraph.costMatrix
                )
                shortestPath = enumerate_paths(dijkstraState, dst)
                removedCost = remove_shortcuts!(shortestPath, TTGraph)
                print("X3")
            else
                print("X2")
            end
        else
            print("X1")
        end
        pathCost = path_cost(shortestPath, costMatrix)
    end
    return shortestPath, pathCost
end

function greedy_insertion(
    solution::Solution,
    TTGraph::TravelTimeGraph,
    TSGraph::TimeSpaceGraph,
    bundle::Bundle,
    src::Int,
    dst::Int,
    CAPACITIES::Vector{Int};
    sorted::Bool=false,
    current_cost::Bool=false,
    findSources::Bool=true,
)
    shortestPath, pathCost = greedy_path(
        solution,
        TTGraph,
        TSGraph,
        bundle,
        src,
        dst,
        CAPACITIES;
        sorted=sorted,
        use_bins=true,
        current_cost=current_cost,
        findSources=findSources,
    )
    # If the path is not admissible, re-computing it
    if !is_path_admissible(TTGraph, shortestPath)
        costMatrix = deepcopy(TTGraph.costMatrix)
        shortestPath, pathCost = greedy_path(
            solution,
            TTGraph,
            TSGraph,
            bundle,
            src,
            dst,
            CAPACITIES;
            sorted=sorted,
            opening_factor=0.5,
            use_bins=true,
            current_cost=current_cost,
            findSources=findSources,
        )
        if !is_path_admissible(TTGraph, shortestPath)
            shortestPath, pathCost = greedy_path(
                solution,
                TTGraph,
                TSGraph,
                bundle,
                src,
                dst,
                CAPACITIES;
                sorted=sorted,
                use_bins=false,
                current_cost=current_cost,
                findSources=findSources,
            )
            if !is_path_admissible(TTGraph, shortestPath)
                for (aSrc, aDst) in TTGraph.bundleArcs[bundle.idx]
                    TTGraph.costMatrix[aSrc, aDst] =
                        TTGraph.networkArcs[aSrc, aDst].distance
                end
                dijkstraState = dijkstra_shortest_paths(
                    TTGraph.graph, src, TTGraph.costMatrix
                )
                shortestPath = enumerate_paths(dijkstraState, dst)
                removedCost = remove_shortcuts!(shortestPath, TTGraph)
                print("X3")
            else
                print("X2")
            end
        else
            print("X1")
        end
        pathCost = path_cost(shortestPath, costMatrix)
    end
    return shortestPath, pathCost
end

function greedy!(solution::Solution, instance::Instance; mode::Int=1)
    TTGraph, TSGraph = instance.travelTimeGraph, instance.timeSpaceGraph
    # Sorting commodities in orders and bundles between them
    sort_order_content!(instance)
    sortedBundleIdxs = sortperm(instance.bundles; by=bun -> bun.maxPackSize, rev=true)
    # Computing the greedy delivery possible for each bundle
    totalCost, totalPathCost = 0.0, 0.0
    print("Greedy introduction progress : ")
    CHANNEL = create_filled_channel()
    percentIdx = ceil(Int, length(sortedBundleIdxs) / 100)
    for (i, bundleIdx) in enumerate(sortedBundleIdxs)
        bundle = instance.bundles[bundleIdx]
        # Retrieving bundle start and end nodes
        suppNode = TTGraph.bundleSrc[bundleIdx]
        custNode = TTGraph.bundleDst[bundleIdx]
        # Computing shortest path
        shortestPath, pathCost = parallel_greedy_insertion(
            solution, TTGraph, TSGraph, bundle, suppNode, custNode, CHANNEL
        )
        totalPathCost += pathCost
        # Adding to solution
        updateCost = update_solution!(solution, instance, bundle, shortestPath; sorted=true)
        # verification
        # if !isapprox(updateCost, pathCost; atol=1.0)
        if updateCost - pathCost > 10
            debug_insertion(instance, solution, bundle, shortestPath, CHANNEL)
        end
        # @assert isapprox(pathCost, updateCost; atol=50 * EPS) "Path cost ($pathCost) and Update cost ($updateCost) don't match \n bundle : $bundle ($suppNode - $custNode) \n shortestPath : $shortestPath \n bundleIdx : $bundleIdx"
        totalCost += updateCost
        i % 10 == 0 && print("|")
        i % percentIdx == 0 && print(" $(round(Int, i/ percentIdx))% ")
    end
    println()
    if !isapprox(totalPathCost, totalCost; atol=10.0)
        @warn "Computed path cost and update cost don't match" totalPathCost totalCost
        throw(ErrorException("Debug"))
    end
    return totalCost
end

function enforce_strict_admissibility!(solution::Solution, instance::Instance)
    TTGraph, TSGraph = instance.travelTimeGraph, instance.timeSpaceGraph
    # Computing the greedy delivery possible for each bundle
    print("Enforcing strict admissibility : ")
    CHANNEL = create_filled_channel()
    percentIdx = ceil(Int, length(instance.bundles) / 100)
    for (i, bundle) in enumerate(instance.bundles)
        bunPath = solution.bundlePaths[bundle.idx]
        # Recomputing if not admissible
        if !is_path_admissible(TTGraph, bunPath)
            # Retrieving bundle start and end nodes
            suppNode = TTGraph.bundleSrc[bundle.idx]
            custNode = TTGraph.bundleDst[bundle.idx]
            # Computing shortest path
            shortestPath, pathCost = parallel_greedy_insertion(
                solution, TTGraph, TSGraph, bundle, suppNode, custNode, CHANNEL
            )
            # Adding to solution
            updateCost = update_solution!(
                solution, instance, bundle, shortestPath; sorted=true
            )
            # verification
            if !isapprox(updateCost, pathCost; atol=1.0)
                debug_insertion(instance, solution, bundle, shortestPath, CHANNEL)
            end
        end
        i % 10 == 0 && print("|")
        i % percentIdx == 0 && print(" $(round(Int, i/ percentIdx))% ")
    end
    return println()
end