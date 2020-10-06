module atn_coupled_module
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
)

    N = length(AgentAct_1)
    AgentStateNext_1 = zeros(Bool, N)
    AgentStateNext_2 = zeros(Bool, N)
    for k = 1:N
        AgentStateNext_1[k] = AgentState_1[k]
        AgentStateNext_2[k] = AgentState_2[k]
    end
    for i = 1:N
        #Change someone sick in group 1 recovers, mu_1
        if AgentState_1[i] == 1
            if rand() < mu_1
                AgentStateNext_1[i] = 0
            end
        end
        #Change someone sick in group 2 recovers, mu_2
        if AgentState_2[i] == 1
            if rand() < mu_2
                AgentStateNext_2[i] = 0
            end
        end


        #If node is active
        if rand() < AgentAct_1[i]
            #it makes m connections
            for j = 1:m
                #choose a random index
                rInt = rand(1:N)
                #probability of staying is beta_1

                if rand() < beta_1
                    if AgentState_2[i] == 1
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
                else
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

        if rand() < AgentAct_2[i]
            #it makes m connections
            for j = 1:m
                #choose a random index
                rInt = rand(1:N)
                #probability of staying is beta_2

                if rand() < beta_2
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
                else
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

function avg_infected(
    N,
    NTrials,
    m,
    mu,
    eta,
    lambda,
    etaS,
    AgentAct,
    tmax,
    TWindow,
)
    InfectedOverTrials = 0.0
    for nt = 1:NTrials
        AgentState = zeros(Bool, N)
        AgentState[1:Int(0.01 * N)] .= 1
        InfectedOverATrial = 0.0
        for i_time = 1:tmax
            NumInf, AgentState = infected_c(
                AgentState,
                AgentAct,
                m,
                mu,
                lambda,
                etaS,
                eta,
                InfectedOverATrial,
                tmax,
                TWindow,
            )
            if i_time > (tmax - TWindow)
                NumInf = sum(AgentState)
                InfectedOverATrial = InfectedOverATrial + NumInf
            end
        end
        InfectedOverATrial = InfectedOverATrial / TWindow
        InfectedOverTrials = InfectedOverTrials + InfectedOverATrial
    end
    return (InfectedOverTrials / NTrials) / N
end

function avg_infected_par(
    N,
    NTrials,
    m,
    mu,
    grid,
    etaS,
    AgentAct,
    tmax,
    TWindow,
)
    eta = grid[1]
    lambda = grid[2]
    InfectedOverTrials = 0.0
    for nt = 1:NTrials
        AgentState = zeros(Bool, N)
        for k = 1:100
            AgentState[k] = 1
        end
        InfectedOverATrial = 0.0
        for i_time = 1:tmax
            AgentState = infected_c(
                AgentState,
                AgentAct,
                m,
                mu,
                lambda,
                etaS,
                eta,
                InfectedOverATrial,
                tmax,
                TWindow,
            )
            if i_time > (tmax - TWindow)
                NumInf = sum(AgentState)
                InfectedOverATrial = InfectedOverATrial + NumInf
            end
        end
        InfectedOverATrial = InfectedOverATrial / TWindow
        InfectedOverTrials = InfectedOverTrials + InfectedOverATrial
    end
    return (InfectedOverTrials / NTrials) / N
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
)
    N = length(agent_act_1)
    agent_state_1 = zeros(Bool, N)
    agent_state_2 = zeros(Bool, N)

    Iend_1 = zeros(tmax)
    Iend_2 = zeros(tmax)
    for k = 1:100
        agent_state_1[k] = 1
        # agent_state_2[k] = 1
    end


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
        )
        Iend_1[i_time] = sum(agent_state_1)
        Iend_2[i_time] = sum(agent_state_2)
    end
    return Iend_1, Iend_2
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
)
    N = length(agent_act_1)
    agent_state_1 = zeros(Bool, N)
    agent_state_2 = zeros(Bool, N)

    # # Iend_1 = zeros(tmax)
    # Iend_2 = zeros(tmax)
    for k = 1:100
        # agent_state_1[k] = 1
        agent_state_2[k] = 1
    end
    t_window = 500
    infected = 0.0
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
        )
        # Iend_1[i_time] = sum(agent_state_1)
        if i_time > (tmax - t_window)
            infected = infected + sum(agent_state_1)
        end
    end
    #
    return infected / t_window
end

function parameter_sweep(
    n_trials,
    agent_act_1,
    agent_act_2,
    mu_1,
    mu_2,
    lambda_1,
    # beta_1,
    grid,
)
    N = length(agent_act_1)
    beta_2 = 1 - grid[1]
    beta_1 = beta_2
    lambda_2 = grid[2]
    m = 5
    tmax = 2000
    total = 0.0
    for k = 1:n_trials
        infected = infected_totals(
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
        )
        total = total + infected / N
    end
    return total / n_trials
end
end
