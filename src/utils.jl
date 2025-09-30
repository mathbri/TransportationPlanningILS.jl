# Utility functions

# These functions are made to be used across the whole project and are not specifically made for a stucture or specific algorithm

function get_path_nodes(path::Vector{T}) where {T<:AbstractEdge}
    return vcat([src(e) for e in path], [dst(path[end])])
end

function is_path_elementary(path::Vector{UInt})
    if length(path) >= 4
        for (nodeIdx, nodeHash) in enumerate(path)
            if nodeHash in path[(nodeIdx+1):end]
                # println("Non elementary path found : $path")
                return false
            end
        end
    end
    return true
end

function init_counters(labels::Vector{String})
    return Dict{String,Int}(labels .=> 0)
end

function print_counters(counters::Dict{String,Int})
    for (key, value) in pairs(counters)
        println("$key : $value")
    end
end

function Base.zero(::Type{Vector{Int}})
    return Int[]
end

function get_elapsed_time(startTime::Float64)
    return round((time() - startTime) * 1000) / 1000
end

function num_binaries(model::Model)
    return count(is_binary, all_variables(model))
end

function num_integers(model::Model)
    return count(is_integer, all_variables(model))
end

function num_constr(model::Model)
    return num_constraints(model; count_variable_in_set_constraints=false)
end

function num_path_constr(model::Model)
    try
        n = length(model[:path])
        return n
    catch e
        return 0
    end
end

function num_elem_constr(model::Model)
    try
        n = length(model[:elementarity])
        return n
    catch e
        return 0
    end
end

function num_pack_constr(model::Model)
    try
        n = length(model[:packing])
        return n
    catch e
        return 0
    end
end

function num_cut_constr(model::Model)
    try
        n = length(model[:cutSet])
        return n
    catch e
        return 0
    end
end