# File used to launch all kinds of scripts using OFOND package 

function julia_main_test(;
    instanceName::String="small",
    timeLimit::Int=60,
)
    println("\n######################################\n")
    println("Launching TransportationPlanningILS tests")
    println("\n######################################\n")

    @info "Test parameters" instanceName timeLimit

    #####################################################################
    # 1. Read instance and solution
    #####################################################################

    directory = joinpath(Base.dirname(@__DIR__), "data")
    println("Reading data from $directory")
    println("Reading instance $instanceName")
    node_file = joinpath(directory, "$(instanceName)_nodes.csv")
    leg_file = joinpath(directory, "$(instanceName)_legs.csv")
    com_file = joinpath(directory, "$(instanceName)_commodities.csv")
    # read instance 
    instance = read_instance(node_file, leg_file, com_file)
    # adding properties to the instance
    CAPACITIES = Int[]
    instance = add_properties(instance, tentative_first_fit, CAPACITIES)

    totVol = sum(sum(o.volume for o in b.orders) for b in instance.bundles)
    println("Instance volume : $(round(Int, totVol / VOLUME_FACTOR)) m3")

    # read solution
    sol_file = joinpath(directory, "$(instanceName)_routes.csv")
    solution = read_solution(instance, joinpath(directory, sol_file))

    println("\n######################################\n")

    #####################################################################
    # 2. Run all heuristics
    #####################################################################

    # Linear lower bound 
    println("\n####   Lower Bound    #####\n")
    solLB = run_simple_heursitic(instance, lower_bound!)
    println("\n######################################\n")

    # Load plan design lower bound 
    if lower_bound_milp_variables(instance) < MAX_MILP_VAR
        run_simple_heursitic(instance, milp_lower_bound!)
        println("\n######################################\n")
    end

    # Shortest delivery
    run_simple_heursitic(instance, shortest_delivery!)
    println("\n######################################\n")

    # Fully outsourced 
    run_simple_heursitic(instance, fully_outsourced!)
    println("\n######################################\n")

    # Local search on current
    println("\n###   Trying local search on current  ###\n")
    run_local_search(instance, local_search!, solution, timeLimit)
    println("\n######################################\n")

    # Load plan design ils on current 
    println("\n###   Trying load plan design LNS on current  ###\n")
    run_local_search(instance, load_plan_design_ils!, solution, timeLimit)
    println("\n######################################\n")

    # # Plan by plant milp 
    # run_simple_heursitic(instance, plant_by_plant_milp!)
    # println("\n######################################\n")

    # Average, Greedy and ILS
    @info "Filtering with standard procedure"
    solution_LBF = run_simple_heursitic(instance, lower_bound_filtering!)
    println("Bundles filtered : $(count(x -> length(x) == 2, solution_LBF.bundlePaths))")

    instanceSub = extract_filtered_instance(instance, solution_LBF)
    instanceSub = add_properties(instanceSub, tentative_first_fit, CAPACITIES)

    # Average delivery 
    run_simple_heursitic(instanceSub, average_delivery!)
    println("\n######################################\n")

    @info "Constructing greedy, lower bound and mixed solution"
    solution_Mix = Solution(instanceSub)
    solution_G, solution_LB = mix_greedy_and_lower_bound!(solution_Mix, instanceSub)
    feasibles = [
        is_feasible(instanceSub, sol) for sol in [solution_Mix, solution_G, solution_LB]
    ]
    mixCost = compute_cost(instanceSub, solution_Mix)
    gCost = compute_cost(instanceSub, solution_G)
    lbCost = compute_cost(instanceSub, solution_LB)
    @info "Mixed heuristic results" :feasible = string(feasibles) :mixed_cost = mixCost :greedy_cost =
        gCost :lower_bound_cost = lbCost

    # Choosing the best solution as the starting solution
    solutionSub = solution_G
    choiceSolution = argmin([mixCost, gCost, lbCost])
    if choiceSolution == 3
        solutionSub = solution_LB
        @info "Lower bound solution chosen"
    elseif choiceSolution == 1
        solutionSub = solution_Mix
        @info "Mixed solution chosen"
    else
        @info "Greedy solution chosen"
    end

    # Full ILS 
    pertLimit = round(Int, min(180, timeLimit / 5))
    lsLimit = min(2700, timeLimit)

    println("\n###   Initial local search   ###\n")
    local_search!(solutionSub, instanceSub; timeLimit=lsLimit)

    ILS!(
        solutionSub,
        instanceSub;
        timeLimit=min(21600, 2 * timeLimit),
        perturbTimeLimit=pertLimit,
        lsTimeLimit=lsLimit,
        solName="$(instanceName)",
        timedelta=lsLimit,
    )

    return 0 # if things finished successfully
end

function julia_main_logged(
    instanceNames::Vector{String}=["tiny"], timeLimits::Vector{Int}=[300]
)
    for (instanceName, timeLimit) in zip(instanceNames, timeLimits)
        logFileName = joinpath("logs", "log2_$(instanceName).txt")
        open(logFileName, "w") do file
            logger = SimpleLogger(file)
            global_logger(logger)
            redirect_stdout(file) do
                julia_main_test(instanceName=instanceName, timeLimit=timeLimit)
            end
        end
    end
end
