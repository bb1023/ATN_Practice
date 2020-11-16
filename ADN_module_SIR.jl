module ADN_module_SIR

using Statistics

function infected_c(
    agent_state_1,
    agent_state_2,
    activity_good,
    activity_bad,
    m,
    mu,
    lambda,
    alpha,
    N_good,
    N_bad,
    # C_1, 
    # C_2
)

    N = N_good + N_bad
    AgentStateNext_1 = zeros(Int, N_good)
    AgentStateNext_2 = zeros(Int, N_bad)
    for k = 1:N_good
        AgentStateNext_1[k] = agent_state_1[k]
    end
    for k = 1:N_bad
        AgentStateNext_2[k] = agent_state_2[k]
    end

    for k = 1:N_good
        if agent_state_1[k] == 1
            if rand() < mu
                AgentStateNext_1[k] = 2
            end
        end
    end

    for k = 1:N_bad
        if agent_state_2[k] == 1
            if rand() < mu
                AgentStateNext_2[k] = 2
            end
        end
    end

    # loop over good agents
    for i = 1:N
        # if active
        # choose a random index
        if i <= N_good
            if rand() < activity_good[i]
                for j = 1:m
                    rInt = rand(1:N)

                    if rInt <= N_good

                        # the good active node gets the good neigbor sick
                        if agent_state_1[i] == 1
                            if (
                                (agent_state_1[rInt] == 0) &&
                                (rand() < alpha^2 * lambda)
                            )
                                AgentStateNext_1[rInt] = 1
                            end
                        else
                            # agent_state_1[i] == 0
                            # the good active node gets sick from the good neighbor node
                            if (
                                (agent_state_1[rInt] == 1) &&
                                (rand() < alpha^2 * lambda)
                            )
                                AgentStateNext_1[i] = 1
                            end
                        end
                    else
                        rInt = rInt - N_good
                        if agent_state_1[i] == 1
                            if (
                                (agent_state_2[rInt] == 0) &&
                                (rand() < alpha * lambda)
                            )
                                AgentStateNext_2[rInt] = 1
                            end
                        else
                            # agent_state_1[i] == 1
                            # the good active node gets sick from the good neighbor node
                            if (
                                (agent_state_2[rInt] == 1) &&
                                (rand() < alpha * lambda)
                            )
                                AgentStateNext_1[i] = 1
                            end
                        end


                    end
                end
            end
        else
            i = i - N_good
            if rand() < activity_bad[i]
                for j = 1:m
                    rInt = rand(1:N)

                    if rInt <= N_good

                        # the good active node gets the good neigbor sick
                        if agent_state_2[i] == 1
                            if (
                                (agent_state_1[rInt] == 0) &&
                                (rand() < alpha * lambda)
                            )
                                AgentStateNext_1[rInt] = 1
                            end
                        else
                            # agent_state_2[i] == 0
                            # the good active node gets sick from the good neighbor node
                            if (
                                (agent_state_1[rInt] == 1) &&
                                (rand() < alpha * lambda)
                            )
                                AgentStateNext_2[i] = 1
                            end
                        end
                    else
                        rInt = rInt - N_good
                        if agent_state_2[i] == 1
                            if ((agent_state_2[rInt] == 0) && (rand() < lambda))
                                AgentStateNext_2[rInt] = 1
                            end
                        else
                            # agent_state_2[i] == 1
                            # the good active node gets sick from the good neighbor node
                            if ((agent_state_2[rInt] == 1) && (rand() < lambda))
                                AgentStateNext_2[i] = 1
                            end
                        end


                    end
                end
            end
        end
    end


    # if ((agent_state_2[rInt] == 0) && (rand() < alpha * lambda))
    #     AgentStateNext_2[rInt] = 1
    # end


    # C_1 = C_1 .| .!iszero.(AgentStateNext_1)
    # C_2 = C_2 .| .!iszero.(AgentStateNext_2)
    return AgentStateNext_1, AgentStateNext_2#, C_1, C_2
end

    function infected_totals(
    tmax,
    activity_good,
    activity_bad,
    m,
    mu,
    lambda,
    alpha,
    N_good,
    N_bad,
)
    agent_state_1 = zeros(Int, N_good)
    agent_state_2 = zeros(Int, N_bad)
    C_1 = .!iszero.(agent_state_1)
    C_2 = .!iszero.(agent_state_2)
    # Iend_1 = zeros(tmax)
    for k = 1:100
        agent_state_1[k] = 1
        agent_state_2[k] = 1
    end


    for i_time = 1:tmax
        agent_state_1, agent_state_2 = infected_c(
            agent_state_1,
            agent_state_2,
            activity_good,
            activity_bad,
            m,
            mu,
            lambda,
            alpha,
            N_good,
            N_bad,
            # C_1,
            # C_2
        )

    R=(sum(.!iszero.(agent_state_1)) + sum(.!iszero.(agent_state_2)))/(N_good+N_bad)
    T = sum([a==1 for a in agent_state_1]) + sum([a==1 for a in  agent_state_2])
    if R>.95 || T<4
        break
    end
    end

    return  (sum(.!iszero.(agent_state_1)) + sum(.!iszero.(agent_state_2)))/(N_good+N_bad)
end

   function parameter_sweep(n_trials, activity, grid)
    N = length(activity)
    proportion = grid[1]
    proportion = 1 - proportion
    # proportion =.1
    lambda = 0.992681
    alpha =grid[2]
    N_bads = Int(ceil((proportion * N)))
    N_good = N - N_bads
    activity_bad = zeros(N_bads)
    activity_good = zeros(N_good)

    for k = 1:N_good
        activity_good[k] = activity[k]
    end
    for k = 1:N_bads
        activity_bad[k] = activity[k + N_good]
    end


    eta_ratio = 1

    # recovery rate
    mu = 0.05
    # infection rate
# *(1-.77 *proportion)

    m = 5
    tmax = 2000
    total_1 = 0.0
    total_2 = 0.0

    for k = 1:n_trials
        infected_1 = infected_totals(
            tmax,
            activity_good,
            activity_bad,
            m,
            mu,
            lambda,
            alpha,
            N_good,
            N_bads,
        )
        total_1 = total_1 + infected_1
        # total_2 = total_2 + infected_2 / (N)

    end
    return total_1 / n_trials
end

end
