# Instance structure to store problem metadata

struct Instance
    # Network 
    networkGraph::NetworkGraph
    # Travel time graph
    travelTimeGraph::TravelTimeGraph
    # Time space graph
    timeSpaceGraph::TimeSpaceGraph
    # Commodities ordered in bundles
    bundles::Vector{Bundle}
    # Time Horizon 
    timeHorizon::Int
    # Fields needed for writing the solution
    dates::Vector{String}
    partNumbers::Dict{UInt,String}
end

# Methods

# Computing all objects properties
function add_properties(instance::Instance, bin_packing::Function, CAPACITIES::Vector{Int})
    @info "Adding properties to instance"
    start = time()
    newBundles = Bundle[
        add_properties(bundle, instance.networkGraph) for bundle in instance.bundles
    ]
    newTTGraph = TravelTimeGraph(instance.networkGraph, newBundles; reachability=false)
    # Checking a path exists for every bundle in the travel time graph
    checkedBundles = Int[]
    for bundle in newBundles
        if has_path(
            newTTGraph.graph,
            newTTGraph.bundleSrc[bundle.idx],
            newTTGraph.bundleDst[bundle.idx],
        )
            push!(checkedBundles, bundle.idx)
        else
            suppNode = code_for(instance.networkGraph.graph, bundle.supplier.hash)
            custNode = code_for(instance.networkGraph.graph, bundle.customer.hash)
            supplier = bundle.supplier
            customer = bundle.customer
            has_path_in_flat_network = has_path(
                instance.networkGraph.graph, suppNode, custNode
            )
            @error "Bundle $(bundle.idx) doesn't have path in the travel time graph" supplier customer has_path_in_flat_network
        end
    end
    @warn "$(length(newBundles) - length(checkedBundles)) bundles have no path in the travel time graph"
    newBundles = newBundles[checkedBundles]
    newBundles = [change_idx(bundle, idx) for (idx, bundle) in enumerate(newBundles)]
    println("Properties added to bundles")
    for bundle in newBundles
        newOrders = [
            add_properties(order, bin_packing, CAPACITIES) for order in bundle.orders
        ]
        empty!(bundle.orders)
        append!(bundle.orders, newOrders)
    end
    println("Properties added to orders")
    newTTGraph = TravelTimeGraph(instance.networkGraph, newBundles)
    @info "Travel-time graph has $(nv(newTTGraph.graph)) nodes and $(ne(newTTGraph.graph)) arcs"
    newTSGraph = TimeSpaceGraph(instance.networkGraph, instance.timeHorizon)
    @info "Time-space graph has $(nv(newTSGraph.graph)) nodes and $(ne(newTSGraph.graph)) arcs"
    timeTaken = round(time() - start; digits=1)
    @info "Full properties added" :time = timeTaken
    return Instance(
        instance.networkGraph,
        newTTGraph,
        newTSGraph,
        newBundles,
        instance.timeHorizon,
        instance.dates,
        instance.partNumbers,
    )
end

# Sorting once and for all order contents
function sort_order_content!(instance::Instance)
    for bundle in instance.bundles
        for order in bundle.orders
            sort!(order.content; rev=true)
        end
    end
end

# Extract a sub instance (either by country or continent) from the instance given
function extract_sub_instance(
    instance::Instance; country::String="", continent::String="", timeHorizon::Int=3
)
    noExtraction = (country == "" && continent == "") || (country != "" && continent != "")
    noExtraction && return instance
    # Redifining network and bundles
    newNetwork = deepcopy(instance.networkGraph)
    newBundles = Bundle[]
    newVertices = Int[]
    # If country arg is not empty, filtering with it
    if country != ""
        newBundles = filter(bun -> is_bundle_in_country(bun, country), instance.bundles)
        newVertices = filter(
            n -> is_node_in_country(newNetwork, n, country), vertices(newNetwork.graph)
        )
    else
        # Otherwise, the contient arg is not empty
        newBundles = filter(bun -> is_bundle_in_continent(bun, continent), instance.bundles)
        newVertices = filter(
            n -> is_node_in_continent(newNetwork, n, continent), vertices(newNetwork.graph)
        )
    end
    length(newBundles) == 0 && @warn "No bundles in the sub instance"
    # Filtering bundle and orders
    newBundles = [
        remove_orders_outside_horizon(bundle, timeHorizon) for bundle in newBundles
    ]
    newBundles = [bundle for bundle in newBundles if length(bundle.orders) > 0]
    newBundles = [change_idx(bundle, idx) for (idx, bundle) in enumerate(newBundles)]
    newNetGraph, _ = induced_subgraph(instance.networkGraph.graph, newVertices)
    newNetwork = NetworkGraph(newNetGraph)
    nNode, nLeg, nBun = nv(newNetGraph), ne(newNetGraph), length(newBundles)
    nOrd = sum(length(bundle.orders) for bundle in newBundles; init=0)
    nCom = sum(
        sum(length(order.content) for order in bundle.orders) for bundle in newBundles;
        init=0,
    )
    nComU = sum(
        sum(length(unique(order.content)) for order in bundle.orders) for
        bundle in newBundles;
        init=0,
    )
    @info "Extracted instance has $nNode nodes, $nLeg legs, $nBun bundles, $nOrd orders and $nCom commodities ($nComU unique) on a $timeHorizon steps time horizon"
    return Instance(
        newNetwork,
        TravelTimeGraph(newNetwork, newBundles),
        TimeSpaceGraph(newNetwork, timeHorizon),
        newBundles,
        timeHorizon,
        instance.dates[1:timeHorizon],
        instance.partNumbers,
    )
end

function extract_sub_instance2(
    instance::Instance; continents::Vector{String}=[""], timeHorizon::Int=3
)
    # Redifining network and bundles
    newNetwork = deepcopy(instance.networkGraph)
    newBundles = filter(bun -> is_bundle_in_continents(bun, continents), instance.bundles)
    filter!(bun -> bun.maxDelTime <= timeHorizon, newBundles)
    newVertices = filter(
        n -> is_node_in_continents(newNetwork, n, continents), vertices(newNetwork.graph)
    )
    length(newBundles) == 0 && @warn "No bundles in the sub instance"
    # Filtering bundle and orders
    newBundles = [
        remove_orders_outside_horizon(bundle, timeHorizon) for bundle in newBundles
    ]
    newBundles = [bundle for bundle in newBundles if length(bundle.orders) > 0]
    newBundles = [change_idx(bundle, idx) for (idx, bundle) in enumerate(newBundles)]
    newNetGraph, _ = induced_subgraph(instance.networkGraph.graph, newVertices)
    newNetwork = NetworkGraph(newNetGraph)
    nNode, nLeg, nBun = nv(newNetGraph), ne(newNetGraph), length(newBundles)
    nOrd = sum(length(bundle.orders) for bundle in newBundles; init=0)
    nCom = sum(
        sum(length(order.content) for order in bundle.orders) for bundle in newBundles;
        init=0,
    )
    nComU = sum(
        sum(length(unique(order.content)) for order in bundle.orders) for
        bundle in newBundles;
        init=0,
    )
    @info "Extracted instance has $nNode nodes, $nLeg legs, $nBun bundles, $nOrd orders and $nCom commodities ($nComU unique) on a $timeHorizon steps time horizon"
    return Instance(
        newNetwork,
        TravelTimeGraph(newNetwork, newBundles),
        TimeSpaceGraph(newNetwork, timeHorizon),
        newBundles,
        timeHorizon,
        instance.dates[1:timeHorizon],
        instance.partNumbers,
    )
end

# Split the instance into multiple instances based on the common time horizon given, 13 weeks by default
function split_instance(instance::Instance, newHorizon::Int=13)
    nFullInstances = div(instance.timeHorizon, newHorizon)
    network = instance.networkGraph
    splitInstances = Instance[]
    for i in 1:(nFullInstances+1)
        # Computing the corresponding time frame
        tStart, tEnd = if i <= nFullInstances
            (i - 1) * newHorizon + 1, i * newHorizon
        else
            nFullInstances * newHorizon + 1, instance.timeHorizon
        end
        # Computing the new bundles involved
        newBundles = [
            remove_orders_outside_frame(bundle, tStart, tEnd) for bundle in instance.bundles
        ]
        filter!(bun -> length(bun.orders) > 0, newBundles)
        # If there is no bundle left, skipping
        length(newBundles) == 0 && continue
        newBundles = [change_idx(bundle, idx) for (idx, bundle) in enumerate(newBundles)]
        # Creating the corresponding instance
        newInstance = Instance(
            network,
            TravelTimeGraph(network, newBundles),
            TimeSpaceGraph(network, newHorizon),
            newBundles,
            newHorizon,
            instance.dates[tStart:tEnd],
            instance.partNumbers,
        )
        push!(splitInstances, newInstance)
    end
    return splitInstances
end

# Split the instance into multiple instances based on the common time horizon given, 12 weeks by default
function split_all_bundles_by_part(instance::Instance)
    # Splitting all bundles
    newBundles = vcat([split_bundle_by_part(bundle, 1) for bundle in instance.bundles]...)
    nBun = length(newBundles)
    println("Split $(length(instance.bundles)) bundles into $nBun bundles")
    newBundles = [change_idx(bundle, idx) for (idx, bundle) in enumerate(newBundles)]
    nBun = length(newBundles)
    nOrd = sum(length(bundle.orders) for bundle in newBundles; init=0)
    nCom = sum(
        sum(length(order.content) for order in bundle.orders) for bundle in newBundles;
        init=0,
    )
    nComU = sum(
        sum(length(unique(order.content)) for order in bundle.orders) for
        bundle in newBundles;
        init=0,
    )
    timeHorizon = instance.timeHorizon
    @info "Split instance (by parts) $nBun bundles, $nOrd orders and $nCom commodities ($nComU unique) on a $timeHorizon steps time horizon"
    return Instance(
        instance.networkGraph,
        TravelTimeGraph(instance.networkGraph, newBundles),
        TimeSpaceGraph(instance.networkGraph, timeHorizon),
        newBundles,
        timeHorizon,
        instance.dates[1:timeHorizon],
        instance.partNumbers,
    )
end

function split_all_bundles_by_time(instance::Instance, newHorizon::Int=13)
    # Splitting all bundles
    newBundles = vcat(
        [split_bundle_by_time(bundle, 1, newHorizon) for bundle in instance.bundles]...
    )
    newBundles = [bundle for bundle in newBundles if length(bundle.orders) > 0]
    nBun = length(newBundles)
    println("Split $(length(instance.bundles)) bundles into $nBun bundles")
    newBundles = [change_idx(bundle, idx) for (idx, bundle) in enumerate(newBundles)]
    nOrd = sum(length(bundle.orders) for bundle in newBundles; init=0)
    nCom = sum(
        sum(length(order.content) for order in bundle.orders) for bundle in newBundles;
        init=0,
    )
    nComU = sum(
        sum(length(unique(order.content)) for order in bundle.orders) for
        bundle in newBundles;
        init=0,
    )
    timeHorizon = instance.timeHorizon
    @info "Split instance (by time) $nBun bundles, $nOrd orders and $nCom commodities ($nComU unique) on a $timeHorizon steps time horizon"
    return Instance(
        instance.networkGraph,
        TravelTimeGraph(instance.networkGraph, newBundles),
        TimeSpaceGraph(instance.networkGraph, timeHorizon),
        newBundles,
        timeHorizon,
        instance.dates[1:timeHorizon],
        instance.partNumbers,
    )
end

function split_all_bundles_into_n(instance::Instance, n::Int=2)
    # Splitting all bundles
    newBundles = vcat(
        [split_bundle_into_n(bundle, 1, n) for bundle in instance.bundles]...
    )
    newBundles = [bundle for bundle in newBundles if length(bundle.orders) > 0]
    nBun = length(newBundles)
    println("Split $(length(instance.bundles)) bundles into $nBun bundles")
    newBundles = [change_idx(bundle, idx) for (idx, bundle) in enumerate(newBundles)]
    nOrd = sum(length(bundle.orders) for bundle in newBundles; init=0)
    nCom = sum(
        sum(length(order.content) for order in bundle.orders) for bundle in newBundles;
        init=0,
    )
    nComU = sum(
        sum(length(unique(order.content)) for order in bundle.orders) for
        bundle in newBundles;
        init=0,
    )
    timeHorizon = instance.timeHorizon
    @info "Split instance (into $n) $nBun bundles, $nOrd orders and $nCom commodities ($nComU unique) on a $timeHorizon steps time horizon"
    return Instance(
        instance.networkGraph,
        TravelTimeGraph(instance.networkGraph, newBundles),
        TimeSpaceGraph(instance.networkGraph, timeHorizon),
        newBundles,
        timeHorizon,
        instance.dates[1:timeHorizon],
        instance.partNumbers,
    )
end

function outsource_arc(arc::NetworkArc, oursourceCost::Float64)
    # If already outsource or is a direct, returning the arc
    arc.type in [:outsource, :direct, :shortcut] && return arc
    # Otherwise, modifying the arc
    newCost = arc.distance * (LAND_CAPACITY / VOLUME_FACTOR) * oursourceCost
    D, T, C, cap = arc.distance, arc.travelTime, arc.carbonCost, arc.capacity
    # Putting isLinear to true
    return NetworkArc(arc.type, D, T, arc.isCommon, newCost, true, C, cap)
end

# Construct an instance where all arcs are linear with oursource costs 
function outsource_instance(instance::Instance)
    # Getting the outsource km m3 cost thanks to the mean of the oursource arcs 
    nOutsourceArcs = count(a -> a.type == :outsource, instance.travelTimeGraph.networkArcs)
    sumOutsourceCost = sum(
        a -> if a.type == :outsource
            a.unitCost / (a.distance * (LAND_CAPACITY / VOLUME_FACTOR))
        else
            0
        end,
        instance.travelTimeGraph.networkArcs,
    )
    oursourceCost = sumOutsourceCost / nOutsourceArcs
    println("Mean outsourcing cost : $(round(oursourceCost; digits=3)) € / m3 / km")
    # Building the new travel time graph with outsource costs
    newInstance = deepcopy(instance)
    I, J, V = findnz(instance.travelTimeGraph.networkArcs)
    newArcs = map(a -> outsource_arc(a, oursourceCost), V)
    for (i, j, newArc) in zip(I, J, newArcs)
        newInstance.travelTimeGraph.networkArcs[i, j] = newArc
    end
    # Buidling the new time space graph
    I, J, V = findnz(instance.timeSpaceGraph.networkArcs)
    newArcs = map(a -> outsource_arc(a, oursourceCost), V)
    for (i, j, newArc) in zip(I, J, newArcs)
        newInstance.timeSpaceGraph.networkArcs[i, j] = newArc
    end
    return newInstance
end
