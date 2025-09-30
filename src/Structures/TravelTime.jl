# Graph structure used to compute bundle paths 

# Travel Time Graph
struct TravelTimeGraph
    # Core fields
    graph::SimpleDiGraph # graph structure
    networkNodes::Vector{NetworkNode} # nodes data
    networkArcs::SparseMatrixCSC{NetworkArc,Int} # arcs data 
    stepToDel::Vector{Int} # nodes time step until delivery step
    # Properties
    costMatrix::SparseMatrixCSC{Float64,Int} # cost matrix used for path computation
    commonNodes::Vector{Int} # common nodes of the graph
    bundleSrc::Vector{Int} # start nodes for bundles
    bundleDst::Vector{Int} # end nodes for bundles
    hashToIdx::Dict{UInt,Int} # dict to easily recover nodes from time space to travel time
    bundleArcs::Vector{Vector{Tuple{Int,Int}}} # arcs usable for each bundle to ease looping through them
end

function TravelTimeGraph()
    return TravelTimeGraph(
        SimpleDiGraph(),
        NetworkNode[],
        sparse([], [], NetworkArc[]),
        Int[],
        sparse(zeros(Float64, 0, 0)),
        Int[],
        Int[],
        Int[],
        Dict{UInt,Int}(),
        Vector{Tuple{Int,Int}}[],
    )
end

function TravelTimeGraph(bundles::Vector{Bundle})
    return TravelTimeGraph(
        SimpleDiGraph(),
        NetworkNode[],
        sparse([], [], NetworkArc[]),
        Int[],
        sparse(zeros(Float64, 0, 0)),
        Int[],
        fill(-1, length(bundles)),
        fill(-1, length(bundles)),
        Dict{UInt,Int}(),
        [Tuple{Int,Int}[] for _ in bundles],
    )
end

function TravelTimeGraph(
    travelTimeGraph::TravelTimeGraph,
    I::Vector{Int},
    J::Vector{Int},
    arcs::Vector{NetworkArc},
    costs::Vector{Float64},
)
    return TravelTimeGraph(
        travelTimeGraph.graph,
        travelTimeGraph.networkNodes,
        sparse(I, J, arcs),
        travelTimeGraph.stepToDel,
        sparse(I, J, costs),
        travelTimeGraph.commonNodes,
        travelTimeGraph.bundleSrc,
        travelTimeGraph.bundleDst,
        travelTimeGraph.hashToIdx,
        travelTimeGraph.bundleArcs,
    )
end

# Methods

function get_bundle_on_nodes(bundles::Vector{Bundle})
    bundlesOnNodes = Dict{UInt,Vector{Bundle}}()
    for bundle in bundles
        supplierBundles = get!(bundlesOnNodes, bundle.supplier.hash, Bundle[])
        push!(supplierBundles, bundle)
        customerBundles = get!(bundlesOnNodes, bundle.customer.hash, Bundle[])
        push!(customerBundles, bundle)
    end
    return bundlesOnNodes
end

function add_timed_node!(
    travelTimeGraph::TravelTimeGraph, nodeData::NetworkNode, stepToDel::Int
)
    add_vertex!(travelTimeGraph.graph)
    push!(travelTimeGraph.networkNodes, nodeData)
    push!(travelTimeGraph.stepToDel, stepToDel)
    nodeIdx = nv(travelTimeGraph.graph)
    travelTimeGraph.hashToIdx[hash(stepToDel, nodeData.hash)] = nodeIdx
    return nodeIdx
end

function add_timed_supplier!(
    travelTimeGraph::TravelTimeGraph,
    nodeData::NetworkNode,
    stepToDel::Int,
    bundlesOnNodes::Dict{UInt,Vector{Bundle}},
)
    nodeIdx = add_timed_node!(travelTimeGraph, nodeData, stepToDel)
    bundleOnSupplier = get(bundlesOnNodes, nodeData.hash, Bundle[])
    startOnNode = filter(bundle -> bundle.maxDelTime == stepToDel, bundleOnSupplier)
    for bundle in startOnNode
        travelTimeGraph.bundleSrc[bundle.idx] = nodeIdx
    end
end

function add_timed_customer!(
    travelTimeGraph::TravelTimeGraph,
    nodeData::NetworkNode,
    bundlesOnNodes::Dict{UInt,Vector{Bundle}},
)
    nodeIdx = add_timed_node!(travelTimeGraph, nodeData, 0)
    bundleOnCustomer = get(bundlesOnNodes, nodeData.hash, Bundle[])
    for bundle in bundleOnCustomer
        travelTimeGraph.bundleDst[bundle.idx] = nodeIdx
    end
end

function add_timed_platform!(
    travelTimeGraph::TravelTimeGraph, nodeData::NetworkNode, stepToDel::Int
)
    nodeIdx = add_timed_node!(travelTimeGraph, nodeData, stepToDel)
    return push!(travelTimeGraph.commonNodes, nodeIdx)
end

function add_network_node!(
    travelTimeGraph::TravelTimeGraph,
    nodeData::NetworkNode,
    bundlesOnNodes::Dict{UInt,Vector{Bundle}},
    maxTime::Int,
)
    # Number of times we have to add a timed copy of the node 
    # plant = 0, suppliers = max of bundle del time, platforms = overall max del time
    if nodeData.type == :supplier
        bundleOnSupplier = get(bundlesOnNodes, nodeData.hash, Bundle[])
        extraCopies = maximum(bundle -> bundle.maxDelTime, bundleOnSupplier; init=0)
        for stepToDel in 0:extraCopies
            add_timed_supplier!(travelTimeGraph, nodeData, stepToDel, bundlesOnNodes)
        end
    elseif nodeData.type == :plant
        add_timed_customer!(travelTimeGraph, nodeData, bundlesOnNodes)
    else
        for stepToDel in 0:maxTime
            add_timed_platform!(travelTimeGraph, nodeData, stepToDel)
        end
    end
end

function add_network_arc!(
    travelTimeGraph::TravelTimeGraph,
    srcData::NetworkNode,
    dstData::NetworkNode,
    arcData::NetworkArc,
)
    maxTime = maximum(travelTimeGraph.stepToDel)
    srcs, dsts = Int[], Int[]
    for stepToDel in 0:maxTime
        # Adding timed copy of network arc
        src = get(travelTimeGraph.hashToIdx, hash(stepToDel, srcData.hash), nothing)
        dst = get(
            travelTimeGraph.hashToIdx,
            hash(stepToDel - arcData.travelTime, dstData.hash),
            nothing,
        )
        if src !== nothing && dst !== nothing
            add_edge!(travelTimeGraph.graph, src, dst)
            push!(srcs, src)
            push!(dsts, dst)
        end
    end
    return srcs, dsts
end

function add_arc_to_vectors!(
    vectors::Tuple{Vector{Int},Vector{Int},Vector{NetworkArc},Vector{Float64}},
    srcs::Vector{Int},
    dsts::Vector{Int},
    arcData::NetworkArc,
)
    I, J, arcs, costs = vectors
    append!(I, srcs)
    append!(J, dsts)
    append!(arcs, fill(arcData, length(srcs)))
    return append!(costs, fill(EPS, length(srcs)))
end

function add_bundle_arcs!(
    travelTimeGraph::TravelTimeGraph,
    bundle::Bundle,
    commonNetworkReach::Dict{Tuple{Int,Int},Bool},
)
    bunArcs = Tuple{Int,Int}[]
    # Discovering all nodes reachable from the bundle source 
    reachableNodes = bfs_parents(
        travelTimeGraph.graph, travelTimeGraph.bundleSrc[bundle.idx]
    )
    # If a node cannot be reached, it as 0 as parent
    reachableNodes = [i for (i, n) in enumerate(reachableNodes) if n != 0]
    # Constructing bundle arcs by adding all outgoing arcs of all reachable nodes
    bunDst = travelTimeGraph.bundleDst[bundle.idx]
    for arcSrc in reachableNodes
        srcData = travelTimeGraph.networkNodes[arcSrc]
        # If the source is in the common network and it can't reach the pant
        if srcData.type != :supplier && !haskey(commonNetworkReach, (arcSrc, bunDst))
            # We can skip it directly 
            continue
        end
        # Otherwise, it is either a supplier (a plant can't be a source) or a common node that can reach the plant
        # We will need to check its outneighbors to see which arcs we want to keep
        for outNeighbor in outneighbors(travelTimeGraph.graph, arcSrc)
            # For every outneighbor 
            outData = travelTimeGraph.networkNodes[outNeighbor]
            # Either it is a common node and we need to check if it can reach the plant
            if outData.type != :plant && outData.type != :supplier
                if haskey(commonNetworkReach, (outNeighbor, bunDst))
                    # If it can reach, adding it 
                    push!(bunArcs, (arcSrc, outNeighbor))
                end
            elseif outData.type == :supplier
                # Either it is a supplier because shortcut and we add it if has path
                if has_path(travelTimeGraph.graph, outNeighbor, bunDst)
                    push!(bunArcs, (arcSrc, outNeighbor))
                end
            else
                # Either it is a plant and we need to check if it is the right plant
                if outNeighbor == bunDst
                    push!(bunArcs, (arcSrc, outNeighbor))
                end
            end
        end
    end
    # Adding the complete list to the Travel Time Graph
    return travelTimeGraph.bundleArcs[bundle.idx] = bunArcs
end

function TravelTimeGraph(
    network::NetworkGraph, bundles::Vector{Bundle}; reachability::Bool=true
)
    # Computing for each node which bundles starts and which bundles end at this node 
    bundlesOnNodes = get_bundle_on_nodes(bundles)
    maxTime = maximum(bundle -> bundle.maxDelTime, bundles)
    # Initializing structure
    travelTimeGraph = TravelTimeGraph(bundles)
    # Adding all nodes from the network graph
    for nodeHash in labels(network.graph)
        add_network_node!(travelTimeGraph, network.graph[nodeHash], bundlesOnNodes, maxTime)
    end
    # Initializing vectors for sparse matrices
    I, J, arcs, costs = Int[], Int[], NetworkArc[], Float64[]
    # Adding all arcs form the network graph
    for (srcHash, dstHash) in edge_labels(network.graph)
        srcData, dstData = network.graph[srcHash], network.graph[dstHash]
        arcData = network.graph[srcHash, dstHash]
        srcs, dsts = add_network_arc!(travelTimeGraph, srcData, dstData, arcData)
        add_arc_to_vectors!((I, J, arcs, costs), srcs, dsts, arcData)
    end
    println(
        "Common arcs : $((count(arc -> arc.type in COMMON_ARC_TYPES, arcs))) (travel time graph)",
    )
    # Computing common network reach
    if reachability
        println("Reachbility properties : ")
        plantNodes = findall(node -> node.type == :plant, travelTimeGraph.networkNodes)
        commonNetworkReach = Dict{Tuple{Int,Int},Bool}()
        for node in travelTimeGraph.commonNodes
            for plant in plantNodes
                if has_path(travelTimeGraph.graph, node, plant)
                    commonNetworkReach[(node, plant)] = true
                end
            end
        end
        print("Common network done, bundles ")
        # Computing bundle arcs
        allUsableArcs = Set{Tuple{Int,Int}}()
        for bundle in bundles
            add_bundle_arcs!(travelTimeGraph, bundle, commonNetworkReach)
            union!(allUsableArcs, travelTimeGraph.bundleArcs[bundle.idx])
            bundle.idx % 100 == 0 && print("|")
        end
        println()
        meanArcs = mean(length.(travelTimeGraph.bundleArcs))
        println(
            "Mean bundle arcs : $(round(meanArcs; digits=1)) ($(minimum(length.(travelTimeGraph.bundleArcs))) -> $(maximum(length.(travelTimeGraph.bundleArcs))))",
        )
        println("All usable arcs : $(length(allUsableArcs))")
    end
    # Creating final structures (because of sparse matrices)
    return TravelTimeGraph(travelTimeGraph, I, J, arcs, costs)
end

function is_path_elementary(travelTimeGraph::TravelTimeGraph, path::Vector{Int})
    return is_path_elementary(hash.(travelTimeGraph.networkNodes[path]))
end

function is_port(travelTimeGraph::TravelTimeGraph, node::Int)
    return travelTimeGraph.networkNodes[node].type in [:pol, :pod]
end

function is_platform(travelTimeGraph::TravelTimeGraph, node::Int)
    return travelTimeGraph.networkNodes[node].type in [:xdock, :iln]
end

function remove_shortcuts!(path::Vector{Int}, travelTimeGraph::TravelTimeGraph)
    firstNode = 1
    for (src, dst) in partition(path, 2, 1)
        if travelTimeGraph.networkArcs[src, dst].type == :shortcut
            firstNode += 1
        else
            break
        end
    end
    deleteat!(path, 1:(firstNode-1))
    return (firstNode - 1) * EPS
end

# Shortcut for computing shortest paths
function shortest_path(travelTimeGraph::TravelTimeGraph, src::Int, dst::Int)
    dijkstraState = dijkstra_shortest_paths(
        travelTimeGraph.graph, src, travelTimeGraph.costMatrix; maxdist=INFINITY
    )
    shortestPath = enumerate_paths(dijkstraState, dst)
    removedCost = remove_shortcuts!(shortestPath, travelTimeGraph)
    pathCost = dijkstraState.dists[dst]
    return shortestPath, pathCost - removedCost
end

# get new index of a node according to a new travel time graph
function new_node_index(
    newTTGraph::TravelTimeGraph, oldTTGraph::TravelTimeGraph, oldTTNode::Int
)
    oldHash = hash(oldTTGraph.stepToDel[oldTTNode], oldTTGraph.networkNodes[oldTTNode].hash)
    !haskey(newTTGraph.hashToIdx, oldHash) && return -1
    return newTTGraph.hashToIdx[oldHash]
end