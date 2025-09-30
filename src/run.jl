# Common part of all the heuristic solving process

function run_heuristic(
    instance::Instance,
    heuristic::Function;
    timeLimit::Int=-1,
    preSolve::Bool=false,
    startSol::Solution=Solution(instance),
)
    @info "Running heuristic $heuristic"
    # Initialize start time
    startTime = time()
    # Complete Instance object with all properties needed
    if preSolve
        CAPACITIES = Int[]
        instance = add_properties(instance, tentative_first_fit, CAPACITIES)
        @info "Pre-solve done" :pre_solve_time = get_elapsed_time(startTime)
    end
    # Initialize solution object
    solution = if Base.size(startSol.bins) != (0, 0)
        deepcopy(startSol)
    else
        Solution(instance)
    end

    # Run the corresponding heuristic
    heuristic(solution, instance)
    println("Cost after initial heuristic: $(compute_cost(instance, solution))")
    improvement = 1.0
    while get_elapsed_time(startTime) < timeLimit && improvement > 1e-3
        improvement = heuristic(solution, instance)
    end

    solveTime = get_elapsed_time(startTime)
    feasible = is_feasible(instance, solution)
    totalCost = compute_cost(instance, solution)
    # detect_infeasibility(instance, solution)
    @info "$heuristic heuristic results" :solve_time = solveTime :feasible = feasible :total_cost =
        totalCost
    return instance, solution
end

function run_simple_heursitic(
    instance::Instance, heuristic::Function, solution::Solution=Solution(instance)
)
    @info "Running heuristic $heuristic"
    startTime = time()
    heuristic(solution, instance)
    # println("Cost after heuristic: $(compute_cost(instance, solution))")
    solveTime = get_elapsed_time(startTime)
    feasible = is_feasible(instance, solution)
    totalCost = compute_cost(instance, solution)
    @info "Results" :solve_time = solveTime :feasible = feasible :total_cost = totalCost
    return solution
end

function run_local_search(
    instance::Instance, heuristic::Function, solution::Solution, timeLimit::Int
)
    @info "Running heuristic $heuristic (on current)"
    startTime = time()
    solutionLS = deepcopy(solution)
    heuristic(solutionLS, instance; timeLimit=timeLimit)
    # println("Cost after heuristic: $(compute_cost(instance, solutionLS))")
    solveTime = get_elapsed_time(startTime)
    feasible = is_feasible(instance, solutionLS)
    totalCost = compute_cost(instance, solutionLS)
    @info "Results" :solve_time = solveTime :feasible = feasible :total_cost = totalCost
end