module ADN_module_SIR

using Statistics

function infected_c(
    agent_state_1,
    agent_state_2,
    activity_good,
    activity_asshole,
    m,
    mu,
    lambda,
    alpha,
    N_good,
    N_asshole,
    C_1, 
    C_2
)
    latency=3*mu
    N = N_good + N_asshole
    AgentStateNext_1 = zeros(Int, N_good)
    AgentStateNext_2 = zeros(Int, N_asshole)
    for k = 1:N_good
        AgentStateNext_1[k] = agent_state_1[k]
    end
    for k = 1:N_asshole
        AgentStateNext_2[k] = agent_state_2[k]
    end

    for k = 1:N_good
        if agent_state_1[k] == 1
            if rand() < mu
                AgentStateNext_1[k] = 2
            end
        end
    end

    for k = 1:N_asshole
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
                                AgentStateNext_1[rInt] = -1
                            end
                        else
                            # agent_state_1[i] == 0
                            # the good active node gets sick from the good neighbor node
                            if (
                                (agent_state_1[rInt] == 1) &&
                                (rand() < alpha^2 * lambda)
                            )
                                AgentStateNext_1[i] = -1
                            end
                        end
                    else
                        rInt = rInt - N_good
                        if agent_state_1[i] == 1
                            if (
                                (agent_state_2[rInt] == 0) &&
                                (rand() < alpha * lambda)
                            )
                                AgentStateNext_2[rInt] = -1
                            end
                        else
                            # agent_state_1[i] == 1
                            # the good active node gets sick from the good neighbor node
                            if (
                                (agent_state_2[rInt] == 1) &&
                                (rand() < alpha * lambda)
                            )
                                AgentStateNext_1[i] = -1
                            end
                        end


                    end
                end
            end
        else
            i = i - N_good
            if rand() < activity_asshole[i]
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


    C_1 = C_1 .| .!iszero.(AgentStateNext_1)
    C_2 = C_2 .| .!iszero.(AgentStateNext_2)
    return AgentStateNext_1, AgentStateNext_2, C_1, C_2
end

    function infected_totals(
    tmax,
    activity_good,
    activity_asshole,
    m,
    mu,
    lambda,
    alpha,
    N_good,
    N_asshole,
)
    agent_state_1 = zeros(Int, N_good)
    agent_state_2 = zeros(Int, N_asshole)
    C_1 = .!iszero.(agent_state_1)
    C_2 = .!iszero.(agent_state_2)
    # Iend_1 = zeros(tmax)
    for k = 1:100
        agent_state_1[k] = 1
        agent_state_2[k] = 1
    end


    for i_time = 1:tmax
        agent_state_1, agent_state_2, C_1, C_2 = infected_c(
            agent_state_1,
            agent_state_2,
            activity_good,
            activity_asshole,
            m,
            mu,
            lambda,
            alpha,
            N_good,
            N_asshole,
            C_1,
            C_2
        )
  
    end
    return sum(C_1) + sum(C_2)
end

   function parameter_sweep(n_trials, activity, grid)
    N = length(activity)
    proportion = grid[1]
    proportion = 1 - proportion
    # proportion =.1
    lambda = .1
    alpha = .7
    N_assholes = Int(ceil((proportion * N)))
    N_good = N - N_assholes
    activity_asshole = zeros(N_assholes)
    activity_good = zeros(N_good)

    for k = 1:N_good
        activity_good[k] = activity[k]
    end
    for k = 1:N_assholes
        activity_asshole[k] = activity[k + N_good]
    end


    eta_ratio = 1

    # recovery rate
    mu = 0.1
    # infection rate


    m = 5
    tmax = 30
    total_1 = 0.0
    total_2 = 0.0

    for k = 1:n_trials
        infected_1 = infected_totals(
            tmax,
            (2.5 + eta_ratio * 2.5) * activity_good,
            5 * activity_asshole,
            m,
            mu,
            lambda,
            alpha,
            N_good,
            N_assholes,
        )
        total_1 = total_1 + infected_1 / (N)
        # total_2 = total_2 + infected_2 / (N)

    end
    return total_1 / n_trials
end

end
