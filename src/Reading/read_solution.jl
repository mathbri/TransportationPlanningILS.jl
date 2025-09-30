# Check that the bundle exists in the instance and return it if found
function check_bundle(instance::Instance, row::CSV.Row)
    suppNode = NetworkNode(row.supplier_account, :supplier, "", "", true, 0.0, 0, 0.0)
    custNode = NetworkNode(row.customer_account, :plant, "", "", true, 0.0, 0, 0.0)
    bundleHash = hash(suppNode, hash(custNode))
    bundleIdx = findfirst(b -> b.hash == bundleHash, instance.bundles)
    if bundleIdx === nothing
        # @warn "Bundle unknown in the instance" :supplier = suppNode :customer = custNode :row =
        #     row
    else
        return instance.bundles[bundleIdx]
    end
end

# Check that the node exists ins the instance and return it if found
function check_node(instance::Instance, row::CSV.Row)
    nodeHash = hash(row.point_account, hash(Symbol(row.point_type)))
    if haskey(instance.networkGraph.graph, nodeHash)
        return instance.networkGraph.graph[nodeHash]
    else
        # @warn "Node unknown in the network" :account = row.point_account :type =
        #     row.point_type :row = row
    end
end

# Add the node to the path, updating path length if needed
function add_node_to_path!(path::Vector{NetworkNode}, node::NetworkNode, idx::Int)
    if idx > length(path)
        for _ in 1:(idx-length(path))
            push!(path, zero(NetworkNode))
        end
    end
    return path[idx] = node
end

# Detect errors in paths
function check_paths(paths::Vector{Vector{NetworkNode}})
    emptyPaths = findall(x -> length(x) == 0, paths)
    if length(emptyPaths) > 0
        emptyPathBundles = join(emptyPaths, ", ")
        @warn "Found $(length(emptyPaths)) empty paths for bundles $(emptyPathBundles)"
    end
    missingPointPaths = findall(
        x -> findfirst(y -> y == zero(NetworkNode), x) !== nothing, paths
    )
    if length(missingPointPaths) > 0
        if length(missingPointPaths) > 10
            missingPointBundles = join(missingPointPaths[1:10], ", ")
            @warn "Missing points in $(length(missingPointPaths)) paths for bundles $missingPointBundles ..."
        else
            missingPointBundles = join(missingPointPaths, ", ")
            @warn "Missing points in $(length(missingPointPaths)) paths for bundles $missingPointBundles"
        end
    end
end

function filter_missing_nodes!(paths::Vector{Vector{NetworkNode}})
    for (i, path) in enumerate(paths)
        paths[i] = filter(x -> x != zero(NetworkNode), path)
    end
end

# Detect whether the path already has an error or not
function is_path_projectable(path::Vector{NetworkNode})
    missingPoints = count(x -> x == zero(NetworkNode), path)
    return (length(path) > 0) && (missingPoints == 0)
end

# Find the next node in the projected path if it exists
function find_next_node(TTGraph::TravelTimeGraph, ttNode::Int, node::NetworkNode)
    inNodes = inneighbors(TTGraph.graph, ttNode)
    nextNodeIdx = findfirst(idx -> TTGraph.networkNodes[idx] == node, inNodes)
    nextNodeIdx === nothing && return nothing
    return inNodes[nextNodeIdx]
end

function project_path(path::Vector{NetworkNode}, TTGraph::TravelTimeGraph, idx::Int)
    # The travel time path is re-created backward by searching for corresponding nodes
    ttPath = [TTGraph.bundleDst[idx]]
    if length(path) <= 1
        return ttPath, true
    end
    # Paths in data files are already backwards
    for node in path[2:end]
        # For each node of the path, we search its inneighbor having the same information
        nextNode = find_next_node(TTGraph, ttPath[end], node)
        if nextNode === nothing
            pathStr = join(string.(path), ", ")
            prev_node = TTGraph.networkNodes[ttPath[end]]
            # @warn "Full path not projectable for bundle $(idx)" :node = node :prev_node =
            #     prev_node :path = pathStr
            break
        end
        push!(ttPath, nextNode)
    end
    error = length(ttPath) < length(path)
    return reverse(ttPath), error
end

# Paths read on the network needs to be projected on the travel-time graph
function project_all_paths(paths::Vector{Vector{NetworkNode}}, TTGraph::TravelTimeGraph)
    ttPaths = [Int[] for _ in 1:length(paths)]
    errors = 0
    for (idx, path) in enumerate(paths)
        # Paths with errors are skipped (only empty paths now)
        # is_path_projectable(path) || continue
        ttPath, error = project_path(path, TTGraph, idx)
        ttPaths[idx] = ttPath
        # If not projected completly, leaving it empty, adding it otherwise 
        errors += error
    end
    @info "Projected $(length(ttPaths)) paths, encoutered errors for $(errors) paths"
    return ttPaths
end

function read_solution(instance::Instance, solution_file::String)
    start = time()
    # Reading .csv file
    csv_reader = CSV.File(
        solution_file;
        types=Dict(
            "supplier_account" => String,
            "customer_account" => String,
            "point_account" => String,
            "point_number" => Int,
        ),
    )
    @info "Reading solution from CSV file $(basename(solution_file)) ($(length(csv_reader)) lines)"
    paths = [NetworkNode[] for _ in 1:length(instance.bundles)]
    ignored = Dict(:unknown_bundle => 0, :unknown_node => 0)
    # Reading paths
    for row in csv_reader
        # Check bundle and node existence (warnings if not found)
        bundle = check_bundle(instance, row)
        if bundle === nothing
            ignored[:unknown_bundle] += 1
            continue
        end
        node = check_node(instance, row)
        if node === nothing
            ignored[:unknown_node] += 1
            continue
        end
        # Add node to path 
        add_node_to_path!(paths[bundle.idx], node, row.point_number)
    end
    check_paths(paths)
    ignoredStr = join(pairs(ignored), ", ")
    timeTaken = round(time() - start; digits=1)
    @info "Read $(sum(length.(paths))) nodes in $(length(paths)) paths" :ignored =
        ignoredStr :time = timeTaken
    # Filtering missing nodes 
    filter_missing_nodes!(paths)
    allPaths = project_all_paths(paths, instance.travelTimeGraph)
    repair_paths!(allPaths, instance)
    # Creating and updating solution
    print("Constructing current solution : ")
    solution = Solution(instance)
    for (path, bundle) in zip(allPaths, instance.bundles)
        update_solution!(solution, instance, bundle, path; sorted=true)
        bundle.idx % 100 == 0 && print("|")
    end
    println()
    feasible = is_feasible(instance, solution; verbose=true)
    totalCost = compute_cost(instance, solution)
    timeTaken = round(time() - start; digits=1)
    @info "Current solution properties" :feasible = feasible :total_cost = totalCost :time =
        timeTaken
    return solution
end