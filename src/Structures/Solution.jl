# Solution structure

# Updating bins and loads in bin packing file

struct Solution
    # Paths used for delivery
    bundlePaths::Vector{Vector{Int}}
    # Bundles on each node of the travel-time common graph + plants
    bundlesOnNode::Dict{Int,Vector{Int}}
    # Transport units completion through time 
    bins::SparseMatrixCSC{Vector{Bin},Int}
end

function Solution(
    travelTimeGraph::TravelTimeGraph,
    timeSpaceGraph::TimeSpaceGraph,
    bundles::Vector{Bundle},
)
    plants = findall(node -> node.type == :plant, travelTimeGraph.networkNodes)
    keys = vcat(travelTimeGraph.commonNodes, plants)
    I, J, V = findnz(timeSpaceGraph.networkArcs)
    bins = [Bin[] for _ in I]
    return Solution(
        [[-1, -1] for _ in 1:length(bundles)],
        Dict{Int,Vector{Int}}(zip(keys, [Int[] for _ in keys])),
        sparse(I, J, bins),
    )
end

function Solution(instance::Instance)
    return Solution(instance.travelTimeGraph, instance.timeSpaceGraph, instance.bundles)
end

# Methods

function update_bundle_path!(
    solution::Solution, bundle::Bundle, path::Vector{Int}; partial::Bool
)
    oldPath = deepcopy(solution.bundlePaths[bundle.idx])
    if partial
        srcIdx = findfirst(node -> node == path[1], oldPath)
        dstIdx = findlast(node -> node == path[end], oldPath)
        if srcIdx === nothing || dstIdx === nothing
            @error "Error : src or dst not found in path" bundle :path = string(path) :oldPath = string(oldPath) partial
        end
        path = vcat(oldPath[1:srcIdx], path[2:(end-1)], oldPath[dstIdx:end])
        oldPath = oldPath[srcIdx:dstIdx]
    end
    solution.bundlePaths[bundle.idx] = path
    return oldPath
end

# As there may be a lot of suppliers given in the paths, change to haskey test instead of get to not create a bundle vector object for nothing
function update_bundle_on_nodes!(
    solution::Solution, bundle::Bundle, path::Vector{Int}; partial::Bool, remove::Bool=false
)
    if partial
        path = path[2:(end-1)]
    end
    for node in path
        bundleVector = get(solution.bundlesOnNode, node, Int[])
        if remove
            # Quick fix for non-admissibility path part removal
            if partial && node in solution.bundlePaths[bundle.idx]
                continue
            end
            filter!(bunIdx -> bunIdx != bundle.idx, bundleVector)
        else
            push!(bundleVector, bundle.idx)
        end
    end
end

# Add / Replace the path of the bundle in the solution with the path given as argument 
# If src and dst are referenced, replace the src-dst part of the current path with the new one
function add_path!(
    solution::Solution, bundle::Bundle, path::Vector{Int}; partial::Bool=false
)
    update_bundle_path!(solution, bundle, path; partial=partial)
    return update_bundle_on_nodes!(solution, bundle, path; partial=partial)
end

# Remove the path of the bundle in the solution
function remove_path!(solution::Solution, bundle::Bundle; src::Int=-1, dst::Int=-1)
    partial = (src != -1) && (dst != -1)
    oldPart = update_bundle_path!(solution, bundle, [src, dst]; partial=partial)
    update_bundle_on_nodes!(solution, bundle, oldPart; partial=partial, remove=true)
    return oldPart
end

# Checking the number of paths
function check_enough_paths(instance::Instance, solution::Solution; verbose::Bool)
    if length(instance.bundles) != length(solution.bundlePaths)
        verbose &&
            @warn "Infeasible solution : $(length(instance.bundles)) bundles and $(length(solution.bundlePaths)) paths"
        return false
    end
    return true
end

# Checking the supplier is correct
function check_supplier(instance::Instance, bundle::Bundle, pathSrc::Int; verbose::Bool)
    pathSupplier = instance.travelTimeGraph.networkNodes[pathSrc]
    if bundle.supplier != pathSupplier
        verbose &&
            @warn "Infeasible solution : bundle $(bundle.idx) has supplier $(bundle.supplier) and its path starts at $(pathSupplier)"
        return false
    end
    return true
end

# Checking the customer is correct
function check_customer(instance::Instance, bundle::Bundle, pathDst::Int; verbose::Bool)
    pathCustomer = instance.travelTimeGraph.networkNodes[pathDst]
    if bundle.customer != pathCustomer
        verbose &&
            @warn "Infeasible solution : bundle $(bundle.idx) has customer $(bundle.customer) and its path ends at $(pathCustomer)"
        return false
    end
    return true
end

# Checking the path is actually feasible on the travel time graph
function check_path_continuity(instance::Instance, path::Vector{Int}; verbose::Bool)
    for (src, dst) in partition(path, 2, 1)
        if !has_edge(instance.travelTimeGraph.graph, src, dst)
            if verbose
                @warn "Infeasible solution : edge $src-$dst doesn't exist in path $path"
            end
            return false
        end
    end
    return true
end

function create_asked_routed_quantities(instance::Instance)
    plantIdxs = findall(node -> node.type == :plant, instance.timeSpaceGraph.networkNodes)
    plantNodes = instance.timeSpaceGraph.networkNodes[plantIdxs]
    plantSteps = instance.timeSpaceGraph.timeStep[plantIdxs]
    askedQuantities = Dict(
        hash(t, p.hash) => Dict{UInt,Int}() for (t, p) in zip(plantSteps, plantNodes)
    )
    routedQuantities = Dict(
        hash(t, p.hash) => Dict{UInt,Int}() for (t, p) in zip(plantSteps, plantNodes)
    )
    return askedQuantities, routedQuantities
end

function update_asked_quantities!(
    askedQuantities::Dict{UInt,Dict{UInt,Int}}, bundle::Bundle
)
    for order in bundle.orders
        timedPlant = hash(order.deliveryDate, bundle.customer.hash)
        for com in order.content
            get!(askedQuantities[timedPlant], com.partNumHash, 0)
            askedQuantities[timedPlant][com.partNumHash] += 1
        end
    end
end

function update_routed_quantities!(
    routedQuantities::Dict{UInt,Dict{UInt,Int}}, instance::Instance, solution::Solution
)
    TSGraph = instance.timeSpaceGraph
    plantIdxs::Vector{Int} = findall(node -> node.type == :plant, TSGraph.networkNodes)
    for pIdx in plantIdxs
        timedPlant = hash(TSGraph.timeStep[pIdx], TSGraph.networkNodes[pIdx].hash)
        for inNeighbor in inneighbors(TSGraph.graph, pIdx)
            for bin in solution.bins[inNeighbor, pIdx]
                for com in bin.content
                    get!(routedQuantities[timedPlant], com.partNumHash, 0)
                    routedQuantities[timedPlant][com.partNumHash] += 1
                end
            end
        end
    end
end

function check_quantities(
    instance::Instance,
    askedQuantities::Dict{UInt,Dict{UInt,Int}},
    routedQuantities::Dict{UInt,Dict{UInt,Int}};
    verbose::Bool=false,
)
    for (timedPlant, quantities) in askedQuantities
        if quantities != routedQuantities[timedPlant]
            allNodes = instance.timeSpaceGraph.networkNodes
            allSteps = instance.timeSpaceGraph.timeStep
            plantIdxs = findall(node -> node.type == :plant, allNodes)
            timedPlantIdx = findfirst(
                i -> hash(allSteps[i], allNodes[i].hash) == timedPlant, plantIdxs
            )
            if verbose
                @warn "Infeasible solution : quantities mismatch between demand and routing"
                println(
                    "timedPlant = $(allNodes[plantIdxs[timedPlantIdx]]) on time step $(allSteps[plantIdxs[timedPlantIdx]])",
                )
                diff = mergewith(-, quantities, routedQuantities[timedPlant])
                filter!(x -> x[2] > 0, diff)
                println("diff = $diff")
            end
            return false
        end
    end
    return true
end

# Check whether a solution is feasible or not 
function is_feasible(
    instance::Instance, solution::Solution; deep::Bool=true, verbose::Bool=false
)
    TTGraph, TSGraph = instance.travelTimeGraph, instance.timeSpaceGraph
    check_enough_paths(instance, solution; verbose=verbose) || return false
    askedQuantities, routedQuantities = create_asked_routed_quantities(instance)
    ALL_COMMODITIES = Commodity[]
    deep && println("Deep feasibility checking...")
    # All paths starts from supplier, end at customer, and are continuous on the graph
    for bundle in instance.bundles
        bundlePath = solution.bundlePaths[bundle.idx]
        pathSrc, pathDst = bundlePath[1], bundlePath[end]
        check_supplier(instance, bundle, pathSrc; verbose=verbose) || return false
        check_customer(instance, bundle, pathDst; verbose=verbose) || return false
        check_path_continuity(instance, bundlePath; verbose=verbose) || return false
        # Updating quantities
        update_asked_quantities!(askedQuantities, bundle)
        if deep
            for order in bundle.orders
                timedPath = time_space_projector(TTGraph, TSGraph, bundlePath, order)
                for (src, dst) in partition(timedPath, 2, 1)
                    # Path continuity
                    if !has_edge(TSGraph.graph, src, dst)
                        @warn "Infeasible solution : edge $src-$dst doesn't exist in path $timedPath"
                        throw(ErrorException("Infeasible"))
                    end
                    # Commodities on timed arcs
                    arcBins = solution.bins[src, dst]
                    allCommodities = get_all_commodities(arcBins, ALL_COMMODITIES)
                    uniqueContent = unique(order.content)
                    for com in uniqueContent
                        orderCount = count(c -> c == com, order.content)
                        arcCount = count(c -> c == com, allCommodities)
                        if !(com in allCommodities)
                            @warn "Infeasible solution : commodity $com doesn't exist in path $timedPath"
                            throw(ErrorException("Infeasible"))
                        end
                        if orderCount < arcCount
                            @warn "$orderCount commodity $com in order and $arcCount in path $timedPath (arc $src->$dst) (bundle $(bundle)) (order $order)"
                        elseif orderCount > arcCount
                            @warn "Infeasible solution : $orderCount commodity $com in order and $arcCount in path $timedPath (arc $src->$dst) (bundle $(bundle)) (order $order)"
                            throw(ErrorException("Infeasible"))
                        end
                    end
                end
            end
        end
    end
    # Getting routed quantities to plants
    update_routed_quantities!(routedQuantities, instance, solution)
    return check_quantities(instance, askedQuantities, routedQuantities; verbose=verbose)
end

# Compute arc cost with respect to the bins on it
function compute_arc_cost(
    TSGraph::TimeSpaceGraph, bins::Vector{Bin}, src::Int, dst::Int; current_cost::Bool
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
    transportUnits = arcData.isLinear ? (arcVolume / arcData.capacity) : 1.0 * length(bins)
    transportCost = current_cost ? TSGraph.currentCost[src, dst] : arcData.unitCost
    cost += transportUnits * transportCost
    return cost
end

# Compute the cost of a solution : node cost + arc cost + commodity cost
function compute_cost(
    instance::Instance, solution::Solution; current_cost::Bool=false
)
    totalCost = 0.0
    # Iterate over sparse matrix
    rows = rowvals(solution.bins)
    vals = nonzeros(solution.bins)
    for j in 1:size(solution.bins, 2)
        for idx in nzrange(solution.bins, j)
            i = rows[idx]
            arcBins = vals[idx]
            # Arc cost
            arcCost = compute_arc_cost(
                instance.timeSpaceGraph, arcBins, i, j; current_cost=current_cost
            )
            totalCost += arcCost
        end
    end
    return totalCost
end

# Compute the extracted cost of a solution
function compute_extracted_cost(instance::Instance, solution::Solution)
    totalCost = 0.0
    for arc in edges(instance.timeSpaceGraph.graph)
        arcBins = solution.bins[src(arc), dst(arc)]
        # If the arc is not direct, it is not extracted
        arcData = instance.timeSpaceGraph.networkArcs[src(arc), dst(arc)]
        srcData = instance.timeSpaceGraph.networkNodes[src(arc)]
        dstData = instance.timeSpaceGraph.networkNodes[dst(arc)]
        trueDirect =
            arcData.type == :direct && srcData.type == :supplier && dstData.type == :plant
        if !trueDirect
            continue
        end
        # If there is no bins, skipping arc
        length(arcBins) == 0 && continue
        # Arc cost
        totalCost += compute_arc_cost(
            instance.timeSpaceGraph, arcBins, src(arc), dst(arc); current_cost=false
        )
    end
    return totalCost
end

# Project a path on the sub instance
function project_on_sub_instance(
    path::Vector{Int}, instance::Instance, subInstance::Instance
)
    TTGraph, subTTGraph = instance.travelTimeGraph, subInstance.travelTimeGraph
    newPath = [new_node_index(subTTGraph, TTGraph, node) for node in path]
    all(node -> node != -1, newPath) || return Int[]
    return newPath
end

# Repair paths by putting direct paths for bundles with errors
function repair_paths!(
    paths::Vector{Vector{Int}}, instance::Instance; directRepair::Bool=false
)
    TTGraph, TSGraph = instance.travelTimeGraph, instance.timeSpaceGraph
    repaired, fullrepair = 0, 0
    for (idx, path) in enumerate(paths)
        startNode = TTGraph.networkNodes[path[1]]
        bundleStart = TTGraph.networkNodes[TTGraph.bundleSrc[idx]]
        if startNode != bundleStart
            repaired += 1
            bundle = instance.bundles[idx]
            for (src, dst) in TTGraph.bundleArcs[bundle.idx]
                # If the arc doesn't need an update, skipping
                is_update_candidate(TTGraph, src, dst, bundle) || continue
                # Otherwise, computing the new cost
                if directRepair
                    TTGraph.costMatrix[src, dst] = TTGraph.networkArcs[src, dst].distance
                else
                    TTGraph.costMatrix[src, dst] = arc_lb_filtering_update_cost(
                        TTGraph, TSGraph, bundle, src, dst
                    )
                end
            end
            # Starting with partial repair
            suppNode, custNode = TTGraph.bundleSrc[idx], path[1]
            shortestPath, pathCost = shortest_path(TTGraph, suppNode, custNode)
            if length(shortestPath) >= 2
                # Partial repair worked
                paths[idx] = vcat(shortestPath[1:(end-1)], path)
            else
                # Partial repair didn't worked : doing full repair 
                custNode = TTGraph.bundleDst[idx]
                shortestPath, pathCost = shortest_path(TTGraph, suppNode, custNode)
                fullrepair += 1
                paths[idx] = shortestPath
            end
        end
    end
    @info "Repaired $(repaired) paths, with full repair for $(fullrepair) paths"
    return count
end

# Extract a sub solution according to the instance given
function extract_sub_solution(solution::Solution, instance::Instance, subInstance::Instance)
    # Extracting corresponding paths
    bundleIdxs = findall(bun -> bun in subInstance.bundles, instance.bundles)
    subPaths = solution.bundlePaths[bundleIdxs]
    # Travel time graph indexing may have changed due to the extraction
    subPaths = [project_on_sub_instance(path, instance, subInstance) for path in subPaths]
    repaired = repair_paths!(subPaths, subInstance)
    # Computing a new solution for the sub instance
    subSolution = Solution(subInstance)
    update_solution!(subSolution, subInstance, subInstance.bundles, subPaths)
    # Checking feasibility and cost
    feasible = is_feasible(subInstance, subSolution; verbose=true)
    totalCost = compute_cost(subInstance, subSolution)
    @info "Extracted solution properties" :repaired = repaired :feasible = feasible :total_cost =
        totalCost
    return subSolution
end

# Extract a filtered instance from the instance and solution given
function extract_filtered_instance(instance::Instance, solution::Solution)
    # Redifining network and bundles
    bundleIdxsToKeep = findall(x -> length(x) > 2, solution.bundlePaths)
    newBundles = instance.bundles[bundleIdxsToKeep]
    length(newBundles) == 0 && @warn "No bundles in the sub instance"
    extractedCost = compute_extracted_cost(instance, solution)
    newVertices = filter(
        n -> !is_node_filterable(instance.networkGraph, n, newBundles),
        vertices(instance.networkGraph.graph),
    )
    # Filtering bundle and orders
    newBundles = [change_idx(bundle, idx) for (idx, bundle) in enumerate(newBundles)]
    newNetGraph, _ = induced_subgraph(instance.networkGraph.graph, newVertices)
    newNetwork = NetworkGraph(newNetGraph)
    nNode, nLeg, nBun = nv(newNetGraph), ne(newNetGraph), length(newBundles)
    nOrd = sum(length(bundle.orders) for bundle in newBundles; init=0)
    nCom = sum(
        sum(length(order.content) for order in bundle.orders) for bundle in newBundles;
        init=0,
    )
    @info "Filtered instance has $nNode nodes, $nLeg legs, $nBun bundles, $nOrd orders and $nCom commodities"
    @info "Extracted cost is $extractedCost"
    return Instance(
        newNetwork,
        TravelTimeGraph(newNetwork, newBundles),
        TimeSpaceGraph(newNetwork, instance.timeHorizon),
        newBundles,
        instance.timeHorizon,
        instance.dates,
        instance.partNumbers,
    )
end

function merge_solutions(
    subSolution::Solution, fullSolution::Solution, instance::Instance, subInstance::Instance
)
    @info "Fusing solutions"
    fusedPaths = [Int[] for _ in 1:length(instance.bundles)]
    # For each bundle of the instance, look for the sub instance index correxponding
    for bundle in instance.bundles
        subIdx = findfirst(x -> x == bundle, subInstance.bundles)
        # If it exist, take the sub instance index path in the sub solution
        if subIdx !== nothing
            subPath = subSolution.bundlePaths[subIdx]
            fusedPaths[bundle.idx] = project_on_sub_instance(subPath, subInstance, instance)
        else
            # otherwise, take the full solution path
            fusedPaths[bundle.idx] = fullSolution.bundlePaths[bundle.idx]
        end
    end
    # Travel time graph indexing may have changed due to the extraction
    repaired = repair_paths!(fusedPaths, instance)
    # Compute the new solution with update solution
    fusedSolution = Solution(instance)
    update_solution!(fusedSolution, instance, instance.bundles, fusedPaths; sorted=true)
    # Checking feasibility and cost
    feasible = is_feasible(instance, fusedSolution; verbose=true)
    totalCost = compute_cost(instance, fusedSolution)
    @info "Fused solution properties" :repaired = repaired :feasible = feasible :total_cost =
        totalCost
    return fusedSolution
end

function project_on_split_instance(solution::Solution, instance::Instance, instanceSplit::Instance)
    @info "Splitting solution by fusing with full shortest solution"
    shortestSol = Solution(instanceSplit)
    shortest_delivery!(shortestSol, instanceSplit)
    sort_order_content!(instanceSplit)
    # Fusing solutions being careful that now there is multiple bundle per OD
    @info "Fusing solutions"
    fusedPaths = [Int[] for _ in 1:length(instanceSplit.bundles)]
    # For each bundle of the instance, look for the sub instance index correxponding
    for bundle in instanceSplit.bundles
        suppCustH = hash(bundle)
        subIdx = findfirst(x -> hash(x) == suppCustH, instance.bundles)
        # If it exist, take the sub instance index path in the sub solution
        if subIdx !== nothing
            subPath = solution.bundlePaths[subIdx]
            fusedPaths[bundle.idx] = project_on_sub_instance(subPath, instance, instanceSplit)
        else
            # otherwise, take the full solution path (only last node to trigger the repair)
            fusedPaths[bundle.idx] = shortestSol.bundlePaths[bundle.idx][end:end]
        end
    end
    # Travel time graph indexing may have changed due to the extraction
    repaired = repair_paths!(fusedPaths, instanceSplit)
    # Compute the new solution with update solution
    fusedSolution = Solution(instanceSplit)
    update_solution!(fusedSolution, instanceSplit, instanceSplit.bundles, fusedPaths; sorted=true)
    # Checking feasibility and cost
    feasible = is_feasible(instanceSplit, fusedSolution; verbose=true)
    totalCost = compute_cost(instanceSplit, fusedSolution)
    @info "Fused solution properties" :repaired = repaired :feasible = feasible :total_cost =
        totalCost
    return fusedSolution
end

# Own deepcopy function to remove runtime dispatch with deepcopy
function solution_deepcopy(solution::Solution, instance::Instance)
    newSol = Solution(instance)
    for bundle in instance.bundles
        empty!(newSol.bundlePaths[bundle.idx])
        append!(newSol.bundlePaths[bundle.idx], solution.bundlePaths[bundle.idx])
    end
    for node in keys(solution.bundlesOnNode)
        append!(newSol.bundlesOnNode[node], solution.bundlesOnNode[node])
    end
    # Efficient iteration over sparse matrix
    rows = rowvals(solution.bins)
    bins = nonzeros(solution.bins)
    for j in 1:size(solution.bins, 2)
        for i in nzrange(solution.bins, j)
            append!(newSol.bins[rows[i], j], my_deepcopy(bins[i]))
        end
    end
    return newSol
end