function solve_lns_milp(
    instance::Instance,
    perturbation::Perturbation;
    warmStart::Bool=false,
    verbose::Bool=false,
    optVerbose::Bool=false,
    withCuts::Bool=false,
    lowerBoundObj::Bool=false,
    timeLimit::Float64=180.0
)
    # Buidling MILP
    if length(instance.bundles) < 800
        timeLimit *= 0.5
    end
    model = model_with_optimizer(; verbose=verbose && optVerbose, timeLimit=timeLimit)
    add_variables!(model, instance, perturbation)
    if verbose
        @info "MILP has $(num_variables(model)) variables ($(num_binaries(model)) binary and $(num_integers(model)) integer)"
    end
    add_path_constraints!(model, instance, perturbation)
    add_elementarity_constraints!(model, instance, perturbation)
    add_packing_constraints!(model, instance, perturbation)
    withCuts && add_cut_set_inequalities!(model, instance)
    if verbose
        @info "MILP has $(num_constr(model)) constraints ($(num_path_constr(model)) path, $(num_elem_constr(model)) elementarity, $(num_pack_constr(model)) packing and $(num_cut_constr(model)) cuts)"
    end
    add_objective!(model, instance, perturbation; lowerBound=lowerBoundObj)
    warmStart && warm_start!(model, instance, perturbation)
    # Solving
    start = time()
    optimize!(model)
    # Getting the solution 
    if has_values(model)
        if verbose
            # Can cause conflict with InferOpt objective_value
            value, bound = JuMP.objective_value(model) + 1e-5, objective_bound(model)
            gap = round(min(100, abs(value - bound) / value); digits=3)
            @info "MILP solved in $(round(time() - start; digits=2)) s with an Objective value = $(round(value; digits=2)) (gap = $gap %)"
        end
        return get_paths(model, instance, perturbation)
    else
        @info "MILP found no solution in $(round(time() - start; digits=2)) s"
        return perturbation.oldPaths
    end
end

function perturbate!(
    solution::Solution, instance::Instance, neighborhood::Symbol; verbose::Bool=false, objTol::Float64=0.015
)
    verbose && @info "Starting perturbation with neighborhood $neighborhood"
    # Selecting perturbation based on neighborhood given 
    perturbation = get_perturbation(neighborhood, instance, solution)
    is_perturbation_empty(perturbation; verbose=verbose) && return 0.0, 0
    verbose && println("Bundles : $(length(perturbation.bundleIdxs))")
    # Computing new paths with lns milp
    pertPaths = solve_lns_milp(instance, perturbation; verbose=verbose, optVerbose=false)
    # Filtering bundles to work on
    !are_new_paths(perturbation.oldPaths, pertPaths; verbose=verbose) && return 0.0, 0
    changedIdxs = get_new_paths_idx(perturbation, perturbation.oldPaths, pertPaths)
    bunIdxs = perturbation.bundleIdxs[changedIdxs]
    verbose && println("Changed bundles : $(length(bunIdxs))")
    pertBundles = instance.bundles[bunIdxs]
    # Here we want to have changedIdxs not bunIdxs
    oldPaths, pertPaths = perturbation.oldPaths[changedIdxs], pertPaths[changedIdxs]
    # Applying new paths to the bundles for which it actually changed
    startCost = compute_cost(instance, solution)
    previousBins = save_previous_bins(instance, solution, pertBundles, oldPaths)
    costRemoved = update_solution!(solution, instance, pertBundles, oldPaths; remove=true)
    updateCost = update_solution!(solution, instance, pertBundles, pertPaths; sorted=true)
    improvement = updateCost + costRemoved
    verbose && println(
        "Improvement : $(round(improvement; digits=1)) ($(round((improvement / startCost) * 100; digits=1)) %)) (Cost Removed = $(round(costRemoved; digits=1)), Cost Added = $(round(updateCost; digits=1)))",
    )
    # Reverting if cost augmented by more than 1.5% 
    if updateCost + costRemoved > objTol * startCost
        verbose && println("Update refused")
        revert_solution!(solution, instance, pertBundles, oldPaths, previousBins, pertPaths)
        return 0.0, 0
    else
        verbose && println("Update accpeted")
        return improvement, length(bunIdxs)
    end
end

function ILS!(
    solution::Solution,
    instance::Instance;
    timeLimit::Int=1800,
    perturbTimeLimit::Int=180,
    lsTimeLimit::Int=300,
    lsStepTimeLimit::Int=60,
    resetCost::Bool=false,
    verbose::Bool=true,
    solName::String="world",
    timedelta::Int=0,
)
    startCost = compute_cost(instance, solution)
    start = time()
    changeThreshold = round(Int, 0.15 * length(instance.bundles))
    costThreshold = 0.0175 * startCost
    println("\n")
    @info "Starting ILS" startCost timeLimit timedelta perturbTimeLimit lsTimeLimit solName changeThreshold costThreshold
    bestSol = solution_deepcopy(solution, instance)
    bestCost = startCost
    # Slope scaling cost update
    if resetCost
        slope_scaling_cost_update!(instance.timeSpaceGraph, Solution(instance))
    else
        slope_scaling_cost_update!(instance.timeSpaceGraph, solution)
    end
    # Apply perturbations in random order
    changed, noChange, noImprov, costAdded, i = 0, 0, 0, 0.0, 0
    while time() - start < timeLimit - timedelta
        neighborhood = PERTURBATIONS[i%2+1]
        startMilp = time()
        while time() - startMilp < perturbTimeLimit
            improv, change = perturbate!(solution, instance, neighborhood; verbose=verbose)
            @info "$neighborhood perturbation(s) applied (without local search)" :improvement =
                improv :change = change
            changed += change
            costAdded += improv
            println("Change threshold : $(changed * 100 / changeThreshold)%")
            println("Cost threshold : $(costAdded * 100 / costThreshold)%")
            # If enough path changed
            if changed >= changeThreshold || costAdded >= costThreshold
                break
            end
        end
        # Recording step with no change 
        if changed == 0
            noChange += 1
        else
            noChange = 0
        end
        # If enough path changed, next phase
        if changed >= changeThreshold
            @info "Triggering local search" changed changeThreshold
            # Applying local search 
            large_local_search!(solution, instance; timeLimit=lsTimeLimit)
            changed, costAdded = 0, 0.0
            # If at < 0.25% of best cost, we add some local search
            currentCost = compute_cost(instance, solution)
            if currentCost > bestCost && currentCost < 1.00025 * bestCost
                @info "Extending local search, possible new best solution" changed changeThreshold
                large_local_search!(solution, instance; timeLimit=1800)
            end
            # If new best solution found, store it
            if compute_cost(instance, solution) < bestCost
                bestSol = solution_deepcopy(solution, instance)
                bestCost = compute_cost(instance, solution)
                timeTaken = timedelta + round(Int, time() - start)
                @info "New best solution found" :cost = bestCost :time = timeTaken
                noImprov = 0
                costThreshold = 0.0175 * bestCost
                startSave = time()
                solFile = joinpath(@__DIR__, "solutions", "$(solName)_t$(timeTaken).jld2")
                save(solFile, "solution", bestSol)
                println("Best Solution saved ($(time() - startSave) s)")
            else
                noImprov += 1
            end
            # Applying cost scaling
            slope_scaling_cost_update!(instance, solution)
        end
        # Breaking if multiple times no change
        if noChange >= 3 || noImprov >= 3
            @info "Stopping ILS" noChange noImprov
            break
        end
        i += 1
    end
    println("\n")
    # Final local search
    lastSol = solution_deepcopy(bestSol, instance)
    large_local_search!(
        lastSol, instance; timeLimit=lsTimeLimit
    )
    if compute_cost(instance, lastSol) < bestCost
        bestSol = solution_deepcopy(lastSol, instance)
        bestCost = compute_cost(instance, lastSol)
        timeTaken = timedelta + round(Int, time() - start)
        @info "New best solution found" :cost = bestCost :time = timeTaken
        solFile = joinpath(@__DIR__, "solutions", "$(solName)_t$(timeTaken).jld2")
        save(solFile, "solution", bestSol)
        println("Best Solution saved")
    end
    totImprov = round(Int, bestCost - startCost)
    @info "Full ILS done" :time = round(time() - start + timedelta; digits=2) :improvement = totImprov
    # Reverting if cost augmented by more than 0.75% (heuristic level)
    if bestCost > startCost
        println("Reverting solution because cost degradation")
        revert_solution!(solution, instance, prevSol)
        return 0.0
    else
        return totImprov
    end
end
