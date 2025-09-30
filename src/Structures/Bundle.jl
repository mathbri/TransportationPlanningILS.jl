# Bundle structure (group of orders with the same origin and destination)

struct Bundle
    # Core fields
    supplier::NetworkNode  # supplier node
    customer::NetworkNode  # customer node
    orders::Vector{Order}  # vector of order
    # Properties
    idx::Int               # index in the instance vector (fast and easy retrieval of informations)
    hash::UInt             # hash of the bundle (fast and easy retrieval of informations)
    maxPackSize::Int       # size of the largest commodity in the bundle
    maxDelTime::Int        # maximum number of steps for delivery authorized
end

function Bundle(supplier::NetworkNode, customer::NetworkNode, idx::Int)
    return Bundle(supplier, customer, Order[], idx, hash(supplier, hash(customer)), 0, 0)
end

function Bundle(supplier::NetworkNode, customer::NetworkNode, idx::Int, maxDelTime::Int)
    return Bundle(supplier, customer, Order[], idx, hash(supplier, hash(customer)), 0, maxDelTime)
end

function Base.hash(bundle::Bundle)
    return hash(bundle.supplier, hash(bundle.customer))
end

function Base.:(==)(bun1::Bundle, bun2::Bundle)
    return bun1.hash == bun2.hash
end

function idx(bundles::Vector{Bundle})
    return map(bundle -> bundle.idx, bundles)
end

function Base.show(io::IO, bundle::Bundle)
    return print(io, "Bundle($(bundle.supplier), $(bundle.customer), idx=$(bundle.idx))")
end

function is_bundle_in_country(bundle::Bundle, country::String)
    return bundle.supplier.country == country && bundle.customer.country == country
end

function is_bundle_in_continent(bundle::Bundle, continent::String)
    return bundle.supplier.continent == continent && bundle.customer.continent == continent
end

function is_bundle_in_continents(bundle::Bundle, continents::Vector{String})
    return bundle.supplier.continent in continents &&
           bundle.customer.continent in continents
end

function change_idx(bundle::Bundle, idx::Int)
    return Bundle(
        bundle.supplier,
        bundle.customer,
        bundle.orders,
        idx,
        bundle.hash,
        bundle.maxPackSize,
        bundle.maxDelTime,
    )
end

function remove_orders_outside_horizon(bundle::Bundle, timeHorizon::Int)
    return Bundle(
        bundle.supplier,
        bundle.customer,
        [order for order in bundle.orders if order.deliveryDate <= timeHorizon],
        bundle.idx,
        bundle.hash,
        bundle.maxPackSize,
        bundle.maxDelTime,
    )
end

function remove_orders_outside_frame(bundle::Bundle, timeStart::Int, timeEnd::Int)
    newOrders = filter(order -> timeStart <= order.deliveryDate <= timeEnd, bundle.orders)
    return Bundle(
        bundle.supplier,
        bundle.customer,
        newOrders,
        bundle.idx,
        bundle.hash,
        bundle.maxPackSize,
        bundle.maxDelTime,
    )
end

# Split the bundle according to part numbers 
# The properties are not recomputed here, no its needs recomputation
function split_bundle_by_part(bundle::Bundle, startIdx::Int)
    # Gather al part numbers in the bundle
    partNums = Set{UInt}()
    for order in bundle.orders
        partNums = union(partNums, map(com -> com.partNumHash, order.content))
    end
    # Creating one bundle per part number
    newBundles = Bundle[]
    for (i, partNum) in enumerate(partNums)
        newHash = hash(partNum, bundle.hash)
        newOrders = Order[]
        # Creating orders only composed of the partNum 
        for order in bundle.orders
            if count(com -> com.partNumHash == partNum, order.content) > 0
                newContent = filter(com -> com.partNumHash == partNum, order.content)
                newBunOrder = Order(order.hash, order.deliveryDate, newContent)
                push!(newOrders, newBunOrder)
            end
        end
        newIdx = startIdx + i - 1
        newBundle = Bundle(
            bundle.supplier, bundle.customer, newOrders, newIdx, newHash, 0, 0
        )
        push!(newBundles, newBundle)
    end
    return newBundles
end

# Split the bundle according to time frames
# The properties are not recomputed here, no its needs recomputation
function split_bundle_by_time(bundle::Bundle, startIdx::Int, newHorizon::Int)
    maxDelDate = maximum(order -> order.deliveryDate, bundle.orders)
    nFullHorizon = div(maxDelDate, newHorizon)
    # Creating one bundle per time frame
    newBundles = Bundle[]
    for i in 1:(nFullHorizon+1)
        # Computing the corresponding time frame (accounting for the last time frame)
        tStart, tEnd = if i <= nFullHorizon
            (i - 1) * newHorizon + 1, i * newHorizon
        else
            nFullHorizon * newHorizon + 1, maxDelDate
        end
        newHash = hash(i, bundle.hash)
        # Computing the new orders involved
        newOrders = filter(order -> tStart <= order.deliveryDate <= tEnd, bundle.orders)
        # If there is orders, adding it to the new bundles
        if length(newOrders) > 0
            newOrders = [Order(newHash, order.deliveryDate, order.content) for order in newOrders]
            newIdx = startIdx + i - 1
            newBundle = Bundle(
                bundle.supplier, bundle.customer, newOrders, newIdx, newHash, 0, 0
            )
            push!(newBundles, newBundle)
        end
    end
    return newBundles
end

# Split the bundle according to time frames
# The properties are not recomputed here, no its needs recomputation
function split_bundle_into_n(bundle::Bundle, startIdx::Int, n::Int)
    # No splitting if all orders have 1 unique commodity00
    if all(o -> length(o.content) == 1, bundle.orders) || n <= 1
        return [bundle]
    end
    # Creating one bundle per part number
    newBundles = Bundle[]
    for i in 1:n
        newHash = hash(i, bundle.hash)
        newOrders = Order[]
        # Creating orders only composed of the partNum 
        for order in bundle.orders
            quantInOrder = div(length(order.content), n)
            if length(order.content) < n
                if i == 1
                    newBunOrder = Order(order.hash, order.deliveryDate, order.content)
                    push!(newOrders, newBunOrder)
                end
            else
                cStart, cEnd = if i < n
                    (i - 1) * quantInOrder + 1, i * quantInOrder
                else
                    (n - 1) * quantInOrder + 1, length(order.content)
                end
                # Also need to update order and commodity hash
                newOrderH = hash(order.deliveryDate, newHash)
                newContent = [Commodity(newOrderH, com.partNumHash, com.size, com.stockCost) for com in order.content[cStart:cEnd]]
                newBunOrder = Order(newOrderH, order.deliveryDate, newContent)
                push!(newOrders, newBunOrder)
            end
        end
        newIdx = startIdx + i - 1
        newBundle = Bundle(
            bundle.supplier, bundle.customer, newOrders, newIdx, newHash, 0, 0
        )
        push!(newBundles, newBundle)
    end
    return newBundles
end

# Average the bundle orders on the whole horizon
function average_bundle(bundle::Bundle, timeHorizon::Int)
    # Computing totals (volume, nb of com, stock cost) on the time horizon and the commodity mean volume
    totVolume = sum(order.volume for order in bundle.orders)
    totCom = sum(length(order.content) for order in bundle.orders)
    totStockCost = sum(order.stockCost for order in bundle.orders)
    # Computing news
    newOrderVolume = totVolume / timeHorizon
    meanComSize = totVolume / totCom
    nCom = ceil(newOrderVolume / meanComSize)
    # Rounding won't have puch impact as we are with m3 / 100
    newComSize = round(newOrderVolume / nCom)
    newComStockCost = totStockCost / nCom
    # Creating 1 new order with commodities of the mean volume  
    newDelDate = bundle.orders[1].deliveryDate
    newContent = [Commodity(0, 0, newComSize, newComStockCost) for _ in 1:nCom]
    newOrder = Order(bundle.hash, newDelDate, newContent)
    return Bundle(
        bundle.supplier, bundle.customer, [newOrder], bundle.idx, bundle.hash, 0, 0
    )
end

function split_bundle_in_half(bundle::Bundle, startIdx::Int)
    # No splitting if only one commodity
    if all(o -> length(o.content) == 1, bundle.orders)
        return [bundle]
    end
    # Dividing bundle in two
    supp, cust, bIdx = bundle.supplier, bundle.customer, startIdx
    newBundles = [
        Bundle(supp, cust, Order[], bIdx, hash(1, bundle.hash), 0, 0),
        Bundle(supp, cust, Order[], bIdx + 1, hash(2, bundle.hash), 0, 0),
    ]
    # Dividing orders 
    for order in bundle.orders
        quantity = max(1, round(Int, length(order.content) / n))
        order1 = Order(order.hash, order.deliveryDate, order.content[1:quantity])
        push!(newBundles[1].orders, order1)
        if length(order.content) > 1
            order2 = Order(
                order.hash, order.deliveryDate, order.content[(quantity+1):end]
            )
            push!(newBundles[2].orders, order2)
        end
    end
    return newBundles
end