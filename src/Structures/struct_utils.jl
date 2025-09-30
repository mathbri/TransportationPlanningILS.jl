# Useful functions using / connecting multiple structures

function Order(bundle::Bundle, deliveryDate::Int)
    return Order(bundle.hash, deliveryDate)
end

function add_properties(bundle::Bundle, network::NetworkGraph)
    maxPackSize = maximum(order -> maximum(com -> com.size, order.content), bundle.orders)
    supp, cust, idx, h = bundle.supplier, bundle.customer, bundle.idx, bundle.hash
    return Bundle(supp, cust, bundle.orders, idx, h, maxPackSize, bundle.maxDelTime)
end

function get_lb_transport_units(order::Order, arcData::NetworkArc)
    # If the arc is shared or already linear
    arcData.type != :direct && return (order.volume / arcData.capacity)
    # If the arc is direct
    return ceil(order.volume / arcData.capacity)
end

function get_transport_units(order::Order, arcData::NetworkArc)
    # If the arc has linear cost
    arcData.isLinear && return (order.volume / arcData.capacity)
    # If the arc is consolidated
    return get(order.bpUnits, arcData.type, 0)
end

function is_node_filterable(networkGraph::NetworkGraph, node::Int, bundles::Vector{Bundle})
    nodeData = networkGraph.graph[label_for(networkGraph.graph, node)]
    !(nodeData.type in [:supplier, :plant]) && return false
    if nodeData.type == :supplier
        return findfirst(b -> b.supplier == nodeData, bundles) === nothing
    else
        return findfirst(b -> b.customer == nodeData, bundles) === nothing
    end
end