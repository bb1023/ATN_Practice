module ATN_Coupled_SIRD
export time_series_par
export time_series_eta
export time_series
using Statistics

function infected_c(
    AgentState_1,
    AgentState_2,
    AgentAct_1,
    AgentAct_2,
    m,
    mu_1,
    mu_2,
    lambda_1,
    lambda_2,
    beta_1,
    beta_2,
    gamma_1,
    gamma_2,
)

    N = length(AgentAct_1)
    AgentStateNext_1 = zeros(Int, N)
    AgentStateNext_2 = zeros(Int, N)
    for k = 1:N
        AgentStateNext_1[k] = AgentState_1[k]
        AgentStateNext_2[k] = AgentState_2[k]
    end
    for i = 1:N
        #Change someone sick in group 1 recovers, mu_1
        if AgentState_1[i] == 1
            if rand() < mu_1
                AgentStateNext_1[i] = 2
            end
        end
        #Change someone sick in group 2 recovers, mu_2
        if AgentState_2[i] == 1
            if rand() < mu_2
                AgentStateNext_2[i] = 2
            end
        end

        # if AgentState_2[i] == 2
        #     if rand() < 0.01
        #         AgentStateNext_2[i] = 0
        #     end
        # end
        if AgentState_1[i] == 1
            if rand() < gamma_1
                AgentStateNext_1[i] = 3
            end
        end
        #Change someone sick in group 2 recovers, mu_2
        if AgentState_2[i] == 1
            if rand() < gamma_2
                AgentStateNext_2[i] = 3
            end
        end



        #If node is active
        if rand() < 5 * AgentAct_1[i]
            #it makes m connections

            for j = 1:m
                #choose a random index
                rInt = rand(1:N)
                #probability of staying is beta_1

                if AgentState_1[i] == 1
                    #the active node gets the new neigbor sick in POP_1
                    if ((AgentState_1[rInt] == 0) && (rand() < lambda_1))
                        AgentStateNext_1[rInt] = 1
                    end
                else
                    #the active node gets sick from the new neighbor node
                    if ((AgentState_1[rInt] == 1) && (rand() < lambda_1))
                        AgentStateNext_1[i] = 1
                    end
                end

                if rand() < beta_1
                    #the active node gets the new neigbor sick in POP_2
                    if AgentState_1[i] == 1
                        if ((AgentState_1[rInt] == 0) && (rand() < lambda_2))
                            AgentStateNext_2[rInt] = 1
                        end
                    else
                        if ((AgentState_2[rInt] == 1) && (rand() < lambda_1))
                            AgentStateNext_1[i] = 1
                        end
                    end
                end
            end
        end

        if rand() < 5 * AgentAct_2[i]
            #it makes m connections
            for j = 1:m
                #choose a random index
                rInt = rand(1:N)
                #probability of staying is beta_2

                if AgentState_2[i] == 1
                    #the active node gets the new neigbor sick in POP_1
                    if ((AgentState_2[rInt] == 0) && (rand() < lambda_2))
                        AgentStateNext_2[rInt] = 1
                    end
                else
                    #the active node gets sick from the new neighbor node
                    if ((AgentState_2[rInt] == 1) && (rand() < lambda_2))
                        AgentStateNext_2[i] = 1
                    end
                end

                if rand() < beta_2

                #the active node gets the new neigbor sick in POP_2
                    if AgentState_2[i] == 1
                        if ((AgentState_2[rInt] == 0) && (rand() < lambda_1))
                            AgentStateNext_1[rInt] = 1
                        end
                    else
                        if ((AgentState_1[rInt] == 1) && (rand() < lambda_2))
                            AgentStateNext_2[i] = 1
                        end
                    end

                end
            end


        end
    end
    return AgentStateNext_1, AgentStateNext_2

end

function infected_totals(
    tmax,
    agent_act_1,
    agent_act_2,
    m,
    mu_1,
    mu_2,
    lambda_1,
    lambda_2,
    beta_1,
    beta_2,
    gamma_1,
    gamma_2,
)
    N = length(agent_act_1)
    agent_state_1 = zeros(Int, N)
    agent_state_2 = zeros(Int, N)

    # Iend_1 = zeros(tmax)
    Iend_2 = zeros(tmax)
    for k = 1:100
        # agent_state_1[k] = 1
        agent_state_2[k] = 1
    end
    t_window = 500

    for i_time = 1:tmax
        agent_state_1, agent_state_2 = infected_c(
            agent_state_1,
            agent_state_2,
            agent_act_1,
            agent_act_2,
            m,
            mu_1,
            mu_2,
            lambda_1,
            lambda_2,
            beta_1,
            beta_2,
            gamma_1,
            gamma_2,
        )
    end
    dead_1 = 0.0
    dead_2 = 0.0
    for k = 1:N
        if agent_state_1[k] == 3
            dead_1 = dead_1 + 1
        end
        # if agent_state_2[k] == 3
        #     dead_2 = dead_2 + 1
        # end
    end
    return dead_1
end

function parameter_sweep(n_trials, agent_act_1, agent_act_2, grid)
    mu_1 = .1
    mu_2 = .1
    N = length(agent_act_1)
    beta_2 = grid[1]
    beta_1 = beta_2
    lambda_2 =  grid[2]
    lambda_1 = lambda_2
    m = 3
    tmax = 1000
    total_1 = 0.0
    total_2 = 0.0
    gamma_1 = .1
    gamma_2 = .01
    for k = 1:n_trials
        infected_1 = infected_totals(
            tmax,
            agent_act_1,
            agent_act_2,
            m,
            mu_1,
            mu_2,
            lambda_1,
            lambda_2,
            beta_1,
            beta_2,
            gamma_1,
            gamma_2,
        )
        total_1 = total_1 + infected_1 / (N)
        # total_2 = total_2 + infected_2 / (N)

    end
    return total_1 / n_trials
end
end
