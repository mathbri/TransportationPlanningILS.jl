###########################################################################################
#################################   First heuristics   ####################################
###########################################################################################

# Benchmark heuristic where all bundle path are computed as the shortest delivery path on the network graph
function shortest_delivery!(solution::Solution, instance::Instance)
    totCost = 0.0
    # Reconstructing TTGraph to have all entries of the cost matrix to EPS
    TTGraph = TravelTimeGraph(instance.networkGraph, instance.bundles)
    # Sorting commodities
    sort_order_content!(instance)
    # Computing the shortest delivery possible for each bundle
    print("Shortest delivery introduction progress : ")
    percentIdx = ceil(Int, length(instance.bundles) / 100)
    for (i, bundle) in enumerate(instance.bundles)
        # Retrieving bundle start and end nodes
        suppNode = TTGraph.bundleSrc[bundle.idx]
        custNode = TTGraph.bundleDst[bundle.idx]
        # Computing shortest path
        for (aSrc, aDst) in TTGraph.bundleArcs[bundle.idx]
            TTGraph.costMatrix[aSrc, aDst] = TTGraph.networkArcs[aSrc, aDst].distance
        end
        shortestPath = enumerate_paths(
            dijkstra_shortest_paths(TTGraph.graph, suppNode, TTGraph.costMatrix), custNode
        )
        # Adding to solution
        totCost += update_solution!(solution, instance, bundle, shortestPath; sorted=true)
        i % 10 == 0 && print("|")
        i % percentIdx == 0 && print(" $(round(Int, i/ percentIdx))% ")
    end
    println()
    return totCost
end

function direct_delivery!(solution::Solution, instance::Instance)
    totCost = 0.0
    # Reconstructing TTGraph to have all entries of the cost matrix to EPS
    TTGraph = TravelTimeGraph(instance.networkGraph, instance.bundles)
    # Sorting commodities
    sort_order_content!(instance)
    # Computing the shortest delivery possible for each bundle
    print("Shortest delivery introduction progress : ")
    percentIdx = ceil(Int, length(instance.bundles) / 100)
    for (i, bundle) in enumerate(instance.bundles)
        # Retrieving bundle start and end nodes
        suppNode = TTGraph.bundleSrc[bundle.idx]
        custNode = TTGraph.bundleDst[bundle.idx]
        # Computing direct path
        for (aSrc, aDst) in TTGraph.bundleArcs[bundle.idx]
            if TTGraph.networkArcs[aSrc, aDst].type == :shortcut
                TTGraph.costMatrix[aSrc, aDst] = EPS
            elseif TTGraph.networkArcs[aSrc, aDst].type == :direct
                TTGraph.costMatrix[aSrc, aDst] = 1.0
            else
                TTGraph.costMatrix[aSrc, aDst] = 1e6
            end
        end
        shortestPath = enumerate_paths(
            dijkstra_shortest_paths(TTGraph.graph, suppNode, TTGraph.costMatrix), custNode
        )
        # Adding to solution
        totCost += update_solution!(solution, instance, bundle, shortestPath; sorted=true)
        i % 10 == 0 && print("|")
        i % percentIdx == 0 && print(" $(round(Int, i/ percentIdx))% ")
    end
    println()
    return totCost
end

# Benchmark heuristic where all bundle path are computed as the minimum cost average delivery using giant trucks approximation for consolidated arcs
# Can be seen as greedy on the relaxation using averaged bundles
function average_delivery!(solution::Solution, instance::Instance)
    println("Averaging bundles")
    # First step : transforming bundles by averaging bundles orders 
    netGraph, timeHorizon = instance.networkGraph, instance.timeHorizon
    avgBundles = Bundle[
        add_properties(average_bundle(bundle, timeHorizon), netGraph) for
        bundle in instance.bundles
    ]
    for avgBun in avgBundles
        for (o, order) in enumerate(avgBun.orders)
            avgBun.orders[o] = add_properties(order, tentative_first_fit, Int[])
        end
    end
    # Second step : use the already prepared lower bound cost matrix update 
    TTGraph, TSGraph = instance.travelTimeGraph, instance.timeSpaceGraph
    totCost = 0.0
    # Sorting commodities
    sort_order_content!(instance)
    # Computing the average delivery
    print("Average delivery introduction progress : ")
    percentIdx = ceil(Int, length(instance.bundles) / 100)
    for (i, (bundle, avgBundle)) in enumerate(zip(instance.bundles, avgBundles))
        # Retrieving bundle start and end nodes
        bSrc = TTGraph.bundleSrc[bundle.idx]
        bDst = TTGraph.bundleDst[bundle.idx]
        # Computing shortest path
        shortestPath, pathCost = lower_bound_insertion(
            solution, TTGraph, TSGraph, avgBundle, bSrc, bDst; use_bins=true, giant=true
        )
        # Adding to solution
        totCost += update_solution!(solution, instance, bundle, shortestPath; sorted=true)
        i % 10 == 0 && print("|")
        i % percentIdx == 0 && print(" $(round(Int, i/ percentIdx))% ")
    end
    println()
    return totCost
end

# Benchmark heuristic where all bundle path are computed as the minimum cost delivery using random costs on arcs
function random_delivery!(
    solution::Solution, instance::Instance; nSol::Int=5000, check::Bool=false
)
    println("Generating $nSol random deliveries")
    TTGraph = instance.travelTimeGraph
    # Sorting commodities
    sort_order_content!(instance)
    # Computing the best solution in nSol random solutions 
    bestSol = Solution(instance)
    shortest_delivery!(bestSol, instance)
    bestCost = compute_cost(instance, bestSol)
    costs = Float64[]
    percentIdx = ceil(Int, length(instance.bundles) / 100)
    start, timeLimit = time(), 30
    timeRecord, bestCostOverTime = [0.0], [bestCost]
    for s in 1:nSol
        sol = Solution(instance)
        solCost = 0.0
        print("Random delivery introduction progress : ")
        # Updating weights 
        for (i, bundle) in enumerate(instance.bundles)
            # Updating arcs 
            for (aSrc, aDst) in TTGraph.bundleArcs[bundle.idx]
                TTGraph.costMatrix[aSrc, aDst] = rand()
            end
        end
        # Computing the best random deliveries for each bundle
        for (i, bundle) in enumerate(instance.bundles)
            # Retrieving bundle start and end nodes
            suppNode = TTGraph.bundleSrc[bundle.idx]
            custNode = TTGraph.bundleDst[bundle.idx]
            # Computing shortest path
            shortestPath = enumerate_paths(
                dijkstra_shortest_paths(TTGraph.graph, suppNode, TTGraph.costMatrix),
                custNode,
            )
            # Adding to solution
            solCost += update_solution!(sol, instance, bundle, shortestPath; sorted=true)
            # i % 10 == 0 && print("|")
            # i % percentIdx == 0 && print(" $(round(Int, i/ percentIdx))% ")
        end
        println()
        push!(costs, solCost)
        if solCost < bestCost
            bestCost = solCost
            bestSol = sol
        end
        # Checking computations
        if check
            @assert is_feasible(instance, sol)
        end
        # Time recording
        push!(timeRecord, round(time() - start; digits=2))
        push!(bestCostOverTime, bestCost)
        if time() - start > timeLimit
            break
        end
        print(" $s ")
    end
    println("Time record : $timeRecord")
    println("Costs : $costs")
    println("Best cost : $bestCostOverTime")
    # Updating the best has the current solution
    return bestSol
end

###########################################################################################
###############################   MILP based heuristics   #################################
###########################################################################################

function full_perturbation(instance::Instance)
    emptySol = Solution(instance)
    shortest_delivery!(emptySol, instance)
    bundleIdxs = idx(instance.bundles)
    oldPaths = emptySol.bundlePaths
    loads = map(bins -> 0, emptySol.bins)
    return Perturbation(:arc_flow, bundleIdxs, oldPaths, loads)
end

# Construct a lower bound MILP on the full instance
function full_lower_bound_milp(
    instance::Instance; withPacking::Bool=true, withElementarity::Bool=true
)
    # Creating an arc flow perturbation with all bundles 
    perturbation = full_perturbation(instance)
    model = model_with_optimizer(; timeLimit=300.0, verbose=true)
    add_variables!(model, instance, perturbation)
    # Putting the current cost back to default unit costs
    slope_scaling_cost_update!(instance.timeSpaceGraph, Solution(instance))
    add_objective!(model, instance, perturbation; lowerBound=false)
    if !withPacking
        for key in eachindex(model[:tau])
            delete(model, model[:tau][key])
        end
        unregister(model, :tau)
    end
    @info "MILP has $(num_variables(model)) variables ($(count(is_binary, all_variables(model))) binary and $(count(is_integer, all_variables(model))) integer)"
    add_path_constraints!(model, instance, perturbation)
    if withPacking
        add_packing_constraints!(model, instance, perturbation)
        # add_cut_set_inequalities!(model, instance)
        if withElementarity
            add_elementarity_constraints!(model, instance, perturbation)
            @info "MILP has $(num_constraints(model; count_variable_in_set_constraints=false)) constraints ($(length(model[:path])) path, $(length(model[:packing])) packing, $(length(model[:elementarity])) elementarity)"
        else
            @info "MILP has $(num_constraints(model; count_variable_in_set_constraints=false)) constraints ($(length(model[:path])) path, $(length(model[:packing])) packing)"
        end
    else
        @info "MILP has $(num_constraints(model; count_variable_in_set_constraints=false)) constraints (path)"
    end
    return model
end

# Compute a lower bound thanks to the MILP machinery of the lns
function milp_lower_bound!(solution::Solution, instance::Instance; verbose::Bool=false)
    println("MILP lower bound construction")
    # Buidling model
    model = full_lower_bound_milp(instance)
    # Getting bound
    optimize!(model)
    objBound = objective_bound(model)
    println("Lower bound computed = $objBound")
    # Getting the solution paths (if any)
    if has_values(model)
        perturbation = full_perturbation(instance)
        newPaths = get_paths(model, instance, perturbation)
        sort_order_content!(instance)
        update_solution!(solution, instance, instance.bundles, newPaths; sorted=true)
    else
        @warn "No solution found, only bound"
        # If no path computed, returning a shortest delivery solution 
        shortest_delivery!(solution, instance)
    end
    return objBound
end

function milp_lower_bound2!(solution::Solution, instance::Instance; verbose::Bool=false)
    println("MILP lower bound construction")
    # Buidling model
    model = full_lower_bound_milp(instance; withElementarity=true)
    # Getting bound
    optimize!(model)
    objBound = objective_bound(model)
    println("Lower bound computed = $objBound")
    # Getting the solution paths (if any)
    if has_values(model)
        perturbation = full_perturbation(instance)
        newPaths = get_paths(model, instance, perturbation)
        sort_order_content!(instance)
        update_solution!(solution, instance, instance.bundles, newPaths; sorted=true)
    else
        @warn "No solution found, only bound"
        # If no path computed, returning a shortest delivery solution 
        shortest_delivery!(solution, instance)
    end
    return objBound
end

# The perturbations in arc flow formulations like plant, random and supplier are also heuristics to use as benchmark 

function plant_by_plant_milp!(solution::Solution, instance::Instance)
    TTGraph, TSGraph = instance.travelTimeGraph, instance.timeSpaceGraph
    totCost = 0.0
    # emptySol = Solution(instance)
    shortest_delivery!(solution, instance)
    # Sorting commodities
    sort_order_content!(instance)
    # Gathering all plants
    plants = findall(node -> node.type == :plant, TTGraph.networkNodes)
    # Going through all plants in random order to select one
    for (i, plant) in enumerate(shuffle(plants))
        @info "Treating plant : $plant ($i / $(length(plants)))"
        plantIdxs = findall(dst -> dst == plant, TTGraph.bundleDst)
        # If no bundle for this plant, skipping to another directly
        length(plantIdxs) == 0 && continue
        # If too much bundles, seperating into smaller groups
        nCommon, j = length(TSGraph.commonArcs), 1
        totGroups = ceil(
            Int,
            sum(length(TTGraph.bundleArcs[b]) for b in plantIdxs) /
            (MAX_MILP_VAR - nCommon),
        )
        while length(plantIdxs) > 0
            @info "Treating group : $j / $totGroups (plant $i / $(length(plants)))"
            nVars = cumsum([length(TTGraph.bundleArcs[b]) for b in plantIdxs])
            stopIdx = findlast(n -> n <= MAX_MILP_VAR - nCommon, nVars)
            bunGroupIdxs = plantIdxs[1:stopIdx]
            plantIdxs = plantIdxs[(stopIdx+1):end]
            # Computing new paths
            perturbation = arc_flow_perturbation(instance, solution, bunGroupIdxs)
            bunPaths = solve_lns_milp(instance, perturbation; warmStart=false, verbose=true)
            # Adding to solution
            bunGroup = instance.bundles[bunGroupIdxs]
            # TODO : remove old paths before adding them ?
            update_solution!(solution, instance, bunGroup; remove=true)
            totCost += update_solution!(solution, instance, bunGroup, bunPaths; sorted=true)
            println(
                "Bundles added : $(length(bunGroup)) (directs = $(count(x -> length(x) == 2, bunPaths)))",
            )
            j += 1
        end
    end
    return totCost
end

###########################################################################################
###########################   Mixing greedy and Lower Bound   #############################
###########################################################################################

# Construct at the same time the greedy and lower bound solution to allow the construction of the combination of both for free
function mix_greedy_and_lower_bound!(
    solution::Solution, instance::Instance; check::Bool=false
)
    # Initialize the other solution objects
    gSol, lbSol, B = Solution(instance), Solution(instance), length(instance.bundles)
    TTGraph, TSGraph = instance.travelTimeGraph, instance.timeSpaceGraph
    # Run the corresponding heuristics
    sort_order_content!(instance)
    sortedBundleIdxs = sortperm(instance.bundles; by=bun -> bun.maxPackSize, rev=true)
    # Computing the delivery possible for each bundle
    print("All introduction progress : ")
    CAPA, percentIdx = Int[], ceil(Int, B / 100)
    CHANNEL = create_filled_channel()
    lowerBound = 0.0
    greedyCostMatrix = deepcopy(TTGraph.costMatrix)
    for (i, bundleIdx) in enumerate(sortedBundleIdxs)
        bundle = instance.bundles[bundleIdx]
        # Retrieving bundle start and end nodes
        bSrc = TTGraph.bundleSrc[bundleIdx]
        bDst = TTGraph.bundleDst[bundleIdx]
        # Computing greedy shortest path
        gPath, gCost = parallel_greedy_insertion(
            gSol, TTGraph, TSGraph, bundle, bSrc, bDst, CHANNEL
        )
        upCost = update_solution!(gSol, instance, bundle, gPath; sorted=true)
        # verification
        # if !isapprox(upCost, gCost; atol=1.0)
        if upCost - gCost > 10
            println("Update cost : $upCost")
            println("Path cost : $gCost")
            debug_insertion(instance, gSol, bundle, gPath, CHANNEL)
        end
        # Saving cost matrix 
        for (src, dst) in TTGraph.bundleArcs[bundleIdx]
            greedyCostMatrix[src, dst] = TTGraph.costMatrix[src, dst]
        end
        # Computing lower bound shortest path
        lbPath, lbCost = lower_bound_insertion(lbSol, TTGraph, TSGraph, bundle, bSrc, bDst)
        lowerBound += lbCost
        update_solution!(lbSol, instance, bundle, lbPath; sorted=true)
        # Computing mixed shortest path
        mixedCostMatrix = (i / B) .* greedyCostMatrix .+ (B - i / B) .* TTGraph.costMatrix
        dijkstraState = dijkstra_shortest_paths(TTGraph.graph, bSrc, mixedCostMatrix)
        shortestPath = enumerate_paths(dijkstraState, bDst)
        remove_shortcuts!(shortestPath, TTGraph)
        update_solution!(solution, instance, bundle, shortestPath; sorted=true)
        # Record progress
        i % 10 == 0 && print("|")
        i % percentIdx == 0 && print(" $(round(Int, i / percentIdx))% ")
        # Checking computations
        if check
            for (src, dst) in TTGraph.bundleArcs[bundleIdx]
                @assert mixedCostMatrix[src, dst] ≈
                        (i / B) * greedyCostMatrix[src, dst] +
                        (B - i / B) * TTGraph.costMatrix[src, dst]
            end
        end
    end
    println("\nLower bound computed : $lowerBound")
    return gSol, lbSol
end

###########################################################################################
#################################   Fully Outsourcing   ###################################
###########################################################################################

function fully_outsourced!(solution::Solution, instance::Instance)
    # Building the new instance costs
    newInstance = outsource_instance(instance)
    # Solution generation (it is a linear program now)
    lower_bound_filtering!(solution, newInstance)
    println("Fully outsourced cost = $(compute_cost(newInstance, solution))")
    return compute_cost(instance, solution)
end

###########################################################################################
#################################   Load Plan Design   ####################################
###########################################################################################

# The computation as it would be done for the load plan design problem (erera et al. for example)
# Need solving perturbation milps but also retrieving the objective of the milp bound

function load_plan_design_arc_cost(
    TSGraph::TimeSpaceGraph, bins::Vector{Bin}, src::Int, dst::Int
)::Float64
    dstData, arcData = TSGraph.networkNodes[dst], TSGraph.networkArcs[src, dst]
    # Computing useful quantities
    arcVolume = sum(bin.load for bin in bins; init=0)
    stockCost = sum(stock_cost(bin) for bin in bins; init=0.0)
    # Volume and Stock cost 
    cost = dstData.volumeCost * arcVolume / VOLUME_FACTOR
    cost += arcData.carbonCost * arcVolume / arcData.capacity
    cost += arcData.distance * stockCost
    # Transport cost 
    transportUnits = if arcData.isLinear
        (arcVolume / arcData.capacity)
    else
        ceil(arcVolume / arcData.capacity)
    end
    cost += transportUnits * arcData.unitCost
    return cost
end

function load_plan_design_cost(instance::Instance, solution::Solution)
    totalCost, TSGraph = 0.0, instance.timeSpaceGraph
    # Iterate over sparse matrix
    rows = rowvals(solution.bins)
    vals = nonzeros(solution.bins)
    for j in 1:size(solution.bins, 2)
        for idx in nzrange(solution.bins, j)
            i = rows[idx]
            arcBins = vals[idx]
            # Arc cost
            totalCost += load_plan_design_arc_cost(TSGraph, arcBins, i, j)
        end
    end
    return totalCost
end

function load_plan_design_ils!(solution::Solution, instance::Instance; timeLimit::Int=300)
    bestCost = load_plan_design_cost(instance, solution)
    slope_scaling_cost_update!(instance.timeSpaceGraph, Solution(instance))
    @info "Starting Load Plan Design LNS" :start_cost = bestCost
    start, i = time(), 0
    while time() - start < timeLimit
        # Perturbating using single plant
        @info "Starting perturbation $i (single_plant)"
        perturbation = get_perturbation(:single_plant, instance, solution)
        is_perturbation_empty(perturbation; verbose=true) && continue
        println("Bundles : $(length(perturbation.bundleIdxs))")
        # With warm start, guaranteed to get a better solution
        pertPaths = solve_lns_milp(
            instance, perturbation; verbose=true, warmStart=true, lowerBoundObj=true
        )
        !are_new_paths(perturbation.oldPaths, pertPaths; verbose=true) && return 0.0, 0
        # Updating solution
        changedIdxs = get_new_paths_idx(perturbation, perturbation.oldPaths, pertPaths)
        bunIdxs = perturbation.bundleIdxs[changedIdxs]
        println("Changed bundles : $(length(bunIdxs))")
        pertBundles = instance.bundles[bunIdxs]
        oldPaths, pertPaths = perturbation.oldPaths[changedIdxs], pertPaths[changedIdxs]
        update_solution!(solution, instance, pertBundles, oldPaths; remove=true)
        update_solution!(solution, instance, pertBundles, pertPaths; sorted=true)
        # Computing cost change 
        newBestCost = load_plan_design_cost(instance, solution)
        println("New Load Plan Design cost : $newBestCost")
        println("Improvement : $(newBestCost - bestCost)")
        bestCost = newBestCost
        i += 1
    end
end

function perturb_only_ils!(solution::Solution, instance::Instance; timeLimit::Int=300)
    bestCost = load_plan_design_cost(instance, solution)
    realCost = compute_cost(instance, solution)
    slope_scaling_cost_update!(instance.timeSpaceGraph, Solution(instance))
    @info "Starting perturbation only ILS" :start_cost = bestCost :real_cost = realCost
    start, i, noImprov = time(), 0, 0
    while time() - start < timeLimit
        # Perturbating using single plant
        @info "Starting perturbation $i (single_plant)"
        improv, change = perturbate!(solution, instance, :single_plant; verbose=true, objTol=1e-5)
        @info "$neighborhood perturbation(s) applied (without local search)" :improvement =
            improv :change = change
        if change == 0
            noImprov += 1
        end
        i += 1
        if noImprov >= 3
            break
        end
    end
end
