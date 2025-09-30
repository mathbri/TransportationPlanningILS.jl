# Functions used to create the MILP model used for the LNS algorithm

# The MILP used for attract / reduce neighborhoods is based on paths as the others are based on arc flows 

###########################################################################################
################################   Adding variables    ####################################
###########################################################################################

# Add to the model the variables for the neighborhood
function add_variables!(model::Model, instance::Instance, perturbation::Perturbation)
    TTGraph, TSGraph = instance.travelTimeGraph, instance.timeSpaceGraph
    if is_attract_reduce(perturbation)
        # Variable z[b] in {0, 1} indicate if bundle b uses old (0) or new (1) path
        @variable(model, z[perturbation.bundleIdxs], Bin)
    else
        # Variable x[b, a] in {0, 1} indicate if bundle b uses arc a in the travel time graph
        @variable(model, x[b in perturbation.bundleIdxs, a in TTGraph.bundleArcs[b]], Bin)
    end
    # Varible tau_a indicate the number of transport units used on arc a in the time space graph
    @variable(model, tau[arc in TSGraph.commonArcs], Int)
end

###########################################################################################
################################   Adding constraints   ###################################
###########################################################################################

# Add path constraints on the travel time graph to the arc flow model
function add_path_constraints!(model::Model, instance::Instance, perturbation::Perturbation)
    # These constraints are implicit in the path model
    if is_attract_reduce(perturbation)
        return nothing
    end
    TTGraph = instance.travelTimeGraph
    bundleIdxs = perturbation.bundleIdxs
    pathConExpr = Dict{Int,Dict{Int,AffExpr}}()
    # Filling expressions bundle by bundle
    for bIdx in bundleIdxs
        pathConExpr[bIdx] = Dict{Int,AffExpr}()
        # Adding +1 to sources and -1 to dests because we want the sum to be 0
        if is_two_shared_node(perturbation)
            pathConExpr[bIdx][perturbation.src] = AffExpr(1)
            pathConExpr[bIdx][perturbation.dst] = AffExpr(-1)
        else
            pathConExpr[bIdx][TTGraph.bundleSrc[bIdx]] = AffExpr(1)
            pathConExpr[bIdx][TTGraph.bundleDst[bIdx]] = AffExpr(-1)
        end
        # For all bundle arcs, adding the corresponding variable to the expressions
        for arc in TTGraph.bundleArcs[bIdx]
            aSrc, aDst = arc
            # With -1 for arc src
            expr = get!(pathConExpr[bIdx], aSrc, AffExpr(0))
            add_to_expression!(expr, -1, model[:x][bIdx, arc])
            # With +1 for arc dst
            expr = get!(pathConExpr[bIdx], aDst, AffExpr(0))
            add_to_expression!(expr, 1, model[:x][bIdx, arc])
        end
    end
    @constraint(
        model, path[b in bundleIdxs, v in keys(pathConExpr[b])], pathConExpr[b][v] == 0
    )
end

# Add all relaxed packing constraints on shared arc of the time space graph
function add_packing_constraints!(
    model::Model, instance::Instance, perturbation::Perturbation
)
    TTGraph, TSGraph = instance.travelTimeGraph, instance.timeSpaceGraph
    tau = model[:tau]
    packConExpr = Dict{Tuple{Int,Int},AffExpr}()
    # Initiating the expressions
    arcDatas = TSGraph.networkArcs
    for (src, dst) in TSGraph.commonArcs
        packConExpr[(src, dst)] = AffExpr(
            perturbation.loads[src, dst], tau[(src, dst)] => -arcDatas[src, dst].capacity
        )
    end
    # Completing the expressions with the path variables
    for bIdx in perturbation.bundleIdxs
        bundle = instance.bundles[bIdx]
        for order in bundle.orders
            # Adding path variable differently is the model is with arc variables or path variables
            if is_attract_reduce(perturbation)
                z = model[:z]
                for (src, dst) in partition(perturbation.oldPaths[bIdx], 2, 1)
                    # Skipping the outsource, direct and shortcut arcs 
                    is_outsource_direct_shortcut(TTGraph, src, dst) && continue
                    # Projecting back on the time space graph
                    tsSrc, tsDst = time_space_projector(TTGraph, TSGraph, src, dst, order)
                    # Getting the corresponding expression
                    packExpr = packConExpr[(tsSrc, tsDst)]
                    # Add to expression order.volume * (1 - z[bIdx])
                    add_to_expression!(packExpr, order.volume)
                    add_to_expression!(packExpr, z[bIdx], -order.volume)
                end
                for (src, dst) in partition(perturbation.newPaths[bIdx], 2, 1)
                    # Skipping the outsource, direct and shortcut arcs 
                    is_outsource_direct_shortcut(TTGraph, src, dst) && continue
                    # Projecting back on the time space graph
                    tsSrc, tsDst = time_space_projector(TTGraph, TSGraph, src, dst, order)
                    # Getting the corresponding expression
                    packExpr = packConExpr[(tsSrc, tsDst)]
                    # Add to expression order.volume * z[bIdx]
                    add_to_expression!(packExpr, z[bIdx], order.volume)
                end
            else
                x = model[:x]
                for (src, dst) in TTGraph.bundleArcs[bIdx]
                    # Skipping the outsource, direct and shortcut arcs 
                    is_outsource_direct_shortcut(TTGraph, src, dst) && continue
                    # Projecting back on the time space graph
                    tsSrc, tsDst = time_space_projector(TTGraph, TSGraph, src, dst, order)
                    # Getting the corresponding expression
                    packExpr = packConExpr[(tsSrc, tsDst)]
                    # Add the corresponding bundle variable and order volume to the expression
                    add_to_expression!(packExpr, x[bIdx, (src, dst)], order.volume)
                end
            end
        end
    end
    @constraint(
        model, packing[(src, dst) in TSGraph.commonArcs], packConExpr[(src, dst)] <= 0
    )
end

###########################################################################################
###########################   Adding elementarity constraints   #############################
###########################################################################################

# Elementarity constraint could be sum of incoming x arc var of the problematic nodes to be less than one

# An interesting question here is that by not enforcing that, the search space is larger and easier to explore, so do we really want to enforce here (an not at the very end) the strict path admissibility ?

function add_elementarity_constraints!(
    model::Model, instance::Instance, perturbation::Perturbation
)
    # These constraints are implicit in the path model
    if is_attract_reduce(perturbation)
        return nothing
    end
    TTGraph = instance.travelTimeGraph
    bundleIdxs = perturbation.bundleIdxs
    elemConExpr = Dict{Int,Dict{String,AffExpr}}()
    # Filling expressions bundle by bundle
    for bIdx in bundleIdxs
        elemConExpr[bIdx] = Dict{String,AffExpr}()
        # For all bundle arcs, adding the corresponding variable to the expressions
        for arc in TTGraph.bundleArcs[bIdx]
            aSrc, aDst = arc
            # If the arc destination is a platform or a port, adding bundle variable to constraint
            dstData = TTGraph.networkNodes[aDst]
            if dstData.type in [:platform, :pol, :pod]
                expr = get!(elemConExpr[bIdx], dstData.account, AffExpr(0))
                add_to_expression!(expr, 1, model[:x][bIdx, arc])
            end
        end
    end
    @constraint(
        model,
        elementarity[b in bundleIdxs, v in keys(elemConExpr[b])],
        elemConExpr[b][v] <= 1
    )
end

###########################################################################################
################################   Adding objective   #####################################
###########################################################################################

function milp_arc_cost(
    instance::Instance, bIdx::Int, src::Int, dst::Int; lowerBound::Bool=false
)
    TTGraph, TSGraph = instance.travelTimeGraph, instance.timeSpaceGraph
    bundle = instance.bundles[bIdx]
    arcData = TTGraph.networkArcs[src, dst]
    # No cost for shortcut arcs
    arcData.type == :shortcut && return EPS
    # Common cost is the sum for orders of volume and stock cost
    cost = sum(volume_stock_cost(TTGraph, src, dst, order) for order in bundle.orders)
    # If it is a oursource or direct arc, computing whole cost
    if is_outsource_direct(TTGraph, src, dst)
        for order in bundle.orders
            tSrc, tDst = time_space_projector(TTGraph, TSGraph, src, dst, order)
            orderTrucks = if lowerBound
                get_lb_transport_units(order, arcData)
            else
                get_transport_units(order, arcData)
            end
            cost += orderTrucks * TSGraph.currentCost[tSrc, tDst]
        end
    end
    return cost
end

# Add the objective function to the model
function add_objective!(
    model::Model, instance::Instance, perturbation::Perturbation; lowerBound::Bool=false
)
    tau = model[:tau]
    TTGraph, TSGraph = instance.travelTimeGraph, instance.timeSpaceGraph
    objExpr = AffExpr()
    # Adding costs of the common arcs in the time space graph 
    for (src, dst) in TSGraph.commonArcs
        add_to_expression!(objExpr, tau[(src, dst)], TSGraph.currentCost[src, dst])
    end
    # Adding costs of arcs in the travel time graph
    for bIdx in perturbation.bundleIdxs
        # Adding a cost per path for attract reduce
        if is_attract_reduce(perturbation)
            z = model[:z]
            oldPathCost = sum(
                milp_arc_cost(instance, bIdx, src, dst; lowerBound=lowerBound) for
                (src, dst) in partition(perturbation.oldPaths[bIdx], 2, 1)
            )
            # Adding to the objective oldPathCost * (1 - z[bIdx])
            add_to_expression!(objExpr, oldPathCost)
            add_to_expression!(objExpr, z[bIdx], -oldPathCost)
            newPathCost = sum(
                milp_arc_cost(instance, bIdx, src, dst; lowerBound=lowerBound) for
                (src, dst) in partition(perturbation.newPaths[bIdx], 2, 1)
            )
            # Adding to the objective newPathCost * z[bIdx]
            add_to_expression!(objExpr, z[bIdx], newPathCost)
            continue
        end
        # Adding a cost per arc for others
        for (src, dst) in TTGraph.bundleArcs[bIdx]
            x = model[:x]
            arcCost = milp_arc_cost(instance, bIdx, src, dst; lowerBound=lowerBound)
            add_to_expression!(objExpr, x[bIdx, (src, dst)], arcCost)
        end
    end
    @objective(model, Min, objExpr)
end

###########################################################################################
###################################   Warm Start   ########################################
###########################################################################################

# Warm start the milp with the current solution
# By default variables, not given are put to 0 which can render the warm start useless
function warm_start!(model::Model, instance::Instance, perturbation::Perturbation)
    TTGraph, TSGraph = instance.travelTimeGraph, instance.timeSpaceGraph
    # Put path variables to 1 for the old path
    for (bIdx, bPath) in zip(perturbation.bundleIdxs, perturbation.oldPaths)
        bundle = instance.bundles[bIdx]
        # For attract reduce, putting the path variable to 0
        if is_attract_reduce(perturbation)
            set_start_value(model[:z][bIdx], 0)
            # Updating the loads along the path
            for (src, dst) in partition(bPath, 2, 1)
                for order in bundle.orders
                    tSrc, tDst = time_space_projector(TTGraph, TSGraph, src, dst, order)
                    perturbation.loads[tSrc, tDst] += order.volume
                end
            end
            continue
        end
        # For two_shared_node, path already between src and dst so no shortcut part needed
        if !is_two_shared_node(perturbation)
            for (src, dst) in partition(get_shortcut_part(TTGraph, bIdx, bPath[1]), 2, 1)
                set_start_value(model[:x][bIdx, (src, dst)], 1)
            end
        end
        # Putting the arc variables to one and adapting loads along the way 
        for (src, dst) in partition(bPath, 2, 1)
            set_start_value(model[:x][bIdx, (src, dst)], 1)
            for order in bundle.orders
                tSrc, tDst = time_space_projector(TTGraph, TSGraph, src, dst, order)
                perturbation.loads[tSrc, tDst] += order.volume
            end
        end
    end
    # Putting packing variables to the current solution
    for (src, dst) in TSGraph.commonArcs
        arcCapacity = TSGraph.networkArcs[src, dst].capacity
        arcLoad = perturbation.loads[src, dst]
        set_start_value(model[:tau][(src, dst)], ceil(arcLoad / arcCapacity))
    end
end

###########################################################################################
################################   Getting solution   #####################################
###########################################################################################

# When getting the solution here, one can check where the objective value differs from the actual update value of the solution path per path 

# Extract the paths from the solution
function get_paths(model::Model, instance::Instance, perturbation::Perturbation)
    TTGraph = instance.travelTimeGraph
    paths = Vector{Vector{Int}}()
    for bIdx in perturbation.bundleIdxs
        if is_attract_reduce(perturbation)
            zVal = value(model[:z][bIdx])
            # If the new path is used, returning it, otherwise returning the old path
            if zVal > 1 - EPS
                push!(paths, perturbation.newPaths[bIdx])
            else
                push!(paths, perturbation.oldPaths[bIdx])
            end
        else
            # Finding all arcs used
            sumX, xAtOne = 0, Tuple{Int,Int}[]
            for (src, dst) in TTGraph.bundleArcs[bIdx]
                xVal = value(model[:x][bIdx, (src, dst)])
                if TTGraph.networkArcs[src, dst].type != :shortcut && xVal > 1 - EPS
                    sumX += 1
                    push!(xAtOne, (src, dst))
                end
                # If the arc is used, it is free, otherwise, he is forbidden
                TTGraph.costMatrix[src, dst] = xVal > 1 - EPS ? EPS : INFINITY
            end
            # Computing a shortest path with all arcs used
            bSrc, bDst = if is_two_shared_node(perturbation)
                perturbation.src, perturbation.dst
            else
                TTGraph.bundleSrc[bIdx], TTGraph.bundleDst[bIdx]
            end
            bundlePath = shortest_path(TTGraph, bSrc, bDst)[1]
            # One check possible is length(bundlePath) - 1 == sum(values(model[:x][bIdx]) without whortcut arcs) (-1 because we count nodes with length and not edges like in x)
            if length(bundlePath) - 1 != sumX
                println("Bundle $bIdx : $(length(bundlePath)) != $(sumX)")
                println(bundlePath)
                println(xAtOne)
            end
            @assert length(bundlePath) - 1 == sumX
            push!(paths, bundlePath)
            # Putting the TTGraph back to normal all free
            for (src, dst) in TTGraph.bundleArcs[bIdx]
                TTGraph.costMatrix[src, dst] = EPS
            end
        end
    end
    println("\nAll paths recovered as modeled and without errors\n")
    return paths
end
