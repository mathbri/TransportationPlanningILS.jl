# File containing all functions to read an instance

# TODO : data validation
# Verify capacity of arcs : either land or sea
# Verify travel times are in time steps
# Verify stock cost are not aberrant 
# Verify there is no false oversea or cross plat
# Verify there is max delivery time for bundles and add it to bundle reading 
# Verify there is capacity for node and add it to node reading

function are_node_data_missing(row::CSV.Row)
    return ismissing(row.point_account) ||
           ismissing(row.point_type) ||
           ismissing(row.point_m3_cost) ||
           ismissing(row.point_m3_capacity)
end

function read_node!(counts::Dict{Symbol,Int}, row::CSV.Row)
    nodeType = Symbol(row.point_type)
    haskey(counts, nodeType) && (counts[nodeType] += 1)
    account, country, continent = promote(
        row.point_account, row.point_country, row.point_continent
    )
    return NetworkNode(
        account,
        nodeType,
        country,
        continent,
        nodeType in COMMON_NODE_TYPES,
        row.point_m3_cost,
        round(Int, row.point_m3_capacity),
        row.point_m3_cost
    )
end

function read_and_add_nodes!(network::NetworkGraph, node_file::String; verbose::Bool=false)
    start = time()
    counts = Dict([(nodeType, 0) for nodeType in NODE_TYPES])
    # Reading .csv file
    csv_reader = CSV.File(
        node_file; types=Dict("point_account" => String, "point_type" => String)
    )
    @info "Reading nodes from CSV file $(basename(node_file)) ($(length(csv_reader)) lines)"
    ignored = Dict(:same_node => 0, :unknown_type => 0, :missing_data => 0)
    for row in csv_reader
        if are_node_data_missing(row)
            ignored[:missing_data] += 1
            continue
        end
        node = read_node!(counts, row)
        added, ignore_type = add_node!(network, node; verbose=verbose)
        added || (ignored[ignore_type] += 1)
    end
    ignoredStr = join(pairs(ignored), ", ")
    timeTaken = round(time() - start; digits=1)
    @info "Read $(nv(network.graph)) nodes : $counts" :ignored = ignoredStr :time =
        timeTaken
end

function are_leg_data_missing(row::CSV.Row)
    return ismissing(row.src_account) ||
           ismissing(row.dst_account) ||
           ismissing(row.src_type) ||
           ismissing(row.dst_type) ||
           ismissing(row.leg_type) ||
           ismissing(row.distance) ||
           ismissing(row.travel_time) ||
           ismissing(row.shipment_cost) ||
           ismissing(row.capacity) ||
           ismissing(row.carbon_cost)
end

function src_dst_hash(row::CSV.Row)
    return hash(row.src_account, hash(Symbol(row.src_type))),
    hash(row.dst_account, hash(Symbol(row.dst_type)))
end

function is_common_arc(row::CSV.Row)
    return Symbol(row.leg_type) in COMMON_ARC_TYPES
end

function read_leg!(counts::Dict{Symbol,Int}, row::CSV.Row, isCommon::Bool)
    arcType = Symbol(row.leg_type)
    haskey(counts, arcType) && (counts[arcType] += 1)
    return NetworkArc(
        arcType,
        row.distance,
        row.travel_time,
        isCommon,
        row.shipment_cost,
        row.is_linear,
        row.carbon_cost,
        round(Int, row.capacity * VOLUME_FACTOR),
    )
end

function read_and_add_legs!(network::NetworkGraph, leg_file::String; verbose::Bool=false)
    start = time()
    neShortcut = ne(network.graph)
    counts = Dict([(arcType, 0) for arcType in ARC_TYPES])
    counts[:shortcut] = neShortcut
    # Reading .csv file
    columns = ["src_account", "dst_account", "src_type", "dst_type", "leg_type"]
    csv_reader = CSV.File(leg_file; types=Dict([(column, String) for column in columns]))
    @info "Reading legs from CSV file $(basename(leg_file)) ($(length(csv_reader)) lines)"
    ignored = Dict(
        :same_arc => 0,
        :unknown_type => 0,
        :unknown_source => 0,
        :unknown_dest => 0,
        :missing_data => 0,
    )
    minCost, maxCost, meanCost = INFINITY, 0.0, 0.0
    for row in csv_reader
        if are_leg_data_missing(row)
            ignored[:missing_data] += 1
            continue
        end
        src, dst = src_dst_hash(row)
        arc = read_leg!(counts, row, is_common_arc(row))
        added, ignore_type = add_arc!(network, src, dst, arc; verbose=verbose)
        added || (ignored[ignore_type] += 1)
        if added
            minCost = min(minCost, row.shipment_cost)
            maxCost = max(maxCost, row.shipment_cost)
            meanCost += row.shipment_cost
        end
    end
    ignoredStr = join(pairs(ignored), ", ")
    meanCost = meanCost / ne(network.graph)
    timeTaken = round(time() - start; digits=1)
    @info "Read $(ne(network.graph)) legs : $counts" :ignored = ignoredStr :time = timeTaken :min_cost =
        minCost :max_cost = maxCost :mean_cost = meanCost
end

function are_commodity_data_missing(row::CSV.Row)
    return ismissing(row.supplier_account) ||
           ismissing(row.customer_account) ||
           ismissing(row.delivery_time_step) ||
           ismissing(row.size) ||
           ismissing(row.delivery_date) ||
           ismissing(row.part_number) ||
           ismissing(row.quantity) ||
           ismissing(row.lead_time_cost) ||
           ismissing(row.max_delivery_time)
end

function bundle_hash(row::CSV.Row)
    return hash(row.supplier_account, hash(row.customer_account))
end

function order_hash(row::CSV.Row)
    return hash(row.delivery_time_step + 1, bundle_hash(row))
end

function com_size(row::CSV.Row)
    baseSize = min(round(Int, max(1, row.size * 100)), SEA_CAPACITY)
    return baseSize
end

function get_bundle!(bundles::Dict{UInt,Bundle}, row::CSV.Row, network::NetworkGraph)
    # Get supplier and customer nodes
    if !haskey(network.graph, hash(row.supplier_account, hash(:supplier)))
        @warn "Supplier unknown in the network" :supplier = row.supplier_account :row = row
    elseif !haskey(network.graph, hash(row.customer_account, hash(:plant)))
        @warn "Customer unknown in the network" :customer = row.customer_account :row = row
    else
        supplierNode = network.graph[hash(row.supplier_account, hash(:supplier))]
        customerNode = network.graph[hash(row.customer_account, hash(:plant))]
        maxDelTime = row.max_delivery_time
        return get!(
            bundles,
            bundle_hash(row),
            Bundle(supplierNode, customerNode, length(bundles) + 1, maxDelTime),
        )
    end
end

# Vectors starts at index 0 in Python so adding 1 to get the right index in Julia
function get_order!(orders::Dict{UInt,Order}, row::CSV.Row, bundle::Bundle)
    return get!(orders, order_hash(row), Order(bundle, row.delivery_time_step + 1))
end

function add_date!(dateHorizon::Vector{String}, row::CSV.Row)
    dateIdx = row.delivery_time_step + 1
    if length(dateHorizon) < dateIdx
        append!(dateHorizon, ["" for _ in 1:(dateIdx-length(dateHorizon))])
    end
    if dateHorizon[dateIdx] == ""
        dateHorizon[dateIdx] = row.delivery_date
    end
end

function read_commodities(networkGraph::NetworkGraph, commodities_file::String)
    start = time()
    orders, bundles = Dict{UInt,Order}(), Dict{UInt,Bundle}()
    dates, partNums = Vector{String}(), Dict{UInt,String}()
    comCount, comUnique = 0, 0
    # Reading .csv file
    csv_reader = CSV.File(
        commodities_file;
        types=Dict(
            "supplier_account" => String,
            "customer_account" => String,
            "delivery_date" => String,
        ),
    )
    @info "Reading commodity orders from CSV file $(basename(commodities_file)) ($(length(csv_reader)) lines)"
    # Creating objects : each line is a commodity order
    ignored = Dict(:unknown_bundle => 0, :missing_data => 0)
    for row in csv_reader
        if are_commodity_data_missing(row)
            ignored[:missing_data] += 1
            continue
        end
        # Getting bundle, order and commodity data
        bundle = get_bundle!(bundles, row, networkGraph)
        if bundle === nothing
            ignored[:unknown_bundle] += 1
            continue
        end
        order = get_order!(orders, row, bundle)
        # If the order is new (no commodities) we have to add it to the bundle
        length(order.content) == 0 && push!(bundle.orders, order)
        # Creating (and Duplicating) commodity
        partNumHash = hash(row.part_number)
        partNums[partNumHash] = string(row.part_number)
        commodity = Commodity(order.hash, partNumHash, com_size(row), row.lead_time_cost)
        rowQuantity = round(Int, row.quantity)
        append!(order.content, [commodity for _ in 1:(rowQuantity)])
        comCount += rowQuantity
        comUnique += 1
        # Is it a new time step ?
        add_date!(dates, row)
    end
    # Transforming dictionnaries into vectors (sorting the vector so that the idx field correspond to the actual idx in the vector)
    bundleVector = sort(collect(values(bundles)); by=bundle -> bundle.idx)
    ignoreStr = join(pairs(ignored), ", ")
    timeTaken = round(time() - start; digits=1)
    @info "Read $(length(bundles)) bundles, $(length(orders)) orders and $comCount commodities ($comUnique without quantities) on a $(length(dates)) steps time horizon" :ignored =
        ignoreStr :time = timeTaken
    return bundleVector, dates, partNums
end

function read_instance(node_file::String, leg_file::String, commodities_file::String)
    networkGraph = NetworkGraph()
    read_and_add_nodes!(networkGraph, node_file)
    read_and_add_legs!(networkGraph, leg_file)
    # Adding general properties 
    netGraph = networkGraph.graph
    maxLeg = maximum(netGraph[srcHash, dstHash].travelTime for (srcHash, dstHash) in edge_labels(netGraph))
    bundles, dates, partNums = read_commodities(networkGraph, commodities_file)
    horizon = length(dates)
    if maxLeg > length(dates)
        @warn "The longest leg is $(maxLeg) steps long, which is longer than the time horizon ($(length(dates)) steps), the new horizon will be $(maxLeg + 1) steps"
        horizon = maxLeg + 1
        lastWeek = Date(DateTime(dates[end][1:10]))
        for i in 1:(horizon-length(dates))
            push!(dates, string(lastWeek + Dates.Week(i)))
        end
    end
    return Instance(
        networkGraph, TravelTimeGraph(), TimeSpaceGraph(), bundles, horizon, dates, partNums
    )
end