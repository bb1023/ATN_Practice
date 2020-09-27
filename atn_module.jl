module atn_module
export time_series_par
export time_series_eta
export time_series
using Statistics

function time_series_par(matrix, Ni, N, m, EtaMax, mu, tmax)
    eta = matrix[1]
    lambda = matrix[2]
    ITrend = zeros(Int, 0)

    Ni = Int(0.01 * N)  # initial number of infects
    gamma = -2.22       # exponent of the power law
    tolerance = 1e-10   # termination tolerance
    TWindow = 500       # time window for evaluating steady state
    TMinConvergence = 500
    ITrend = zeros(Int, 0)
    ActMin = 0.001      # minimum activation probability
    ActMax = 1    
    AgentAct = ((ActMax.^(gamma + 1) - ActMin.^(gamma + 1)) .* rand(N) .+
            ActMin.^(gamma + 1)).^(1.0 / (gamma + 1))
    # AvgDegree = 2 * m * mean(AgentAct)
    AgentState = zeros(Bool, N)
    AgentState[1:Ni] .= 1

    for t = 1:tmax
        AgentStat, ITrend, termination = number_infected_eta(N, AgentState, AgentAct, lambda, eta, EtaMax,  mu, m,  tolerance, t, ITrend, TWindow, TMinConvergence)
        if termination
            break
        end
    end
    return ITrend[end]
end

function time_series(Ni, N, m, gamma, lambda, mu, tolerance, TWindow, TMinConvergence, tmax)
    ITrend = zeros(0)
    ActMin = 0.001      # minimum activation probability
    ActMax = 1    
    AgentAct = ((ActMax.^(gamma + 1) - ActMin.^(gamma + 1)) .* rand(N) .+
            ActMin.^(gamma + 1)).^(1.0 / (gamma + 1))
    # AvgDegree = 2 * m * mean(AgentAct)
    AgentState = zeros(Bool, N)
    AgentState[1:Ni] .= 1

    for t = 1:tmax
        AgentStat, ITrend, termination = number_infected(N, AgentState, AgentAct, lambda,  mu, m,  tolerance, t, ITrend, TWindow, TMinConvergence)
        if termination
            break
        end

    end
    return ITrend
end

function time_series_eta(Ni, N, m, eta, EtaMax, gamma, lambda, mu, tolerance, TWindow, TMinConvergence, tmax)
    ITrend = zeros(0)
    ActMin = 0.001      # minimum activation probability
    ActMax = 1    
    AgentAct = ((ActMax.^(gamma + 1) - ActMin.^(gamma + 1)) .* rand(N) .+
            ActMin.^(gamma + 1)).^(1.0 / (gamma + 1))
    # AvgDegree = 2 * m * mean(AgentAct)
    AgentState = zeros(Bool, N)
    AgentState[1:Ni] .= 1

    for t = 1:tmax
        AgentStat, ITrend, termination = number_infected_eta(N, AgentState, AgentAct, lambda, eta, EtaMax,  mu, m,  tolerance, t, ITrend, TWindow, TMinConvergence)
        if termination
            break
        end

    end
    return ITrend
end

function number_infected(N, AgentState, AgentAct, lambda,  mu, m,  tolerance, t, ITrend, TWindow, TMinConvergence)
    AdjMatrix = zeros(Bool, N, N)
    ActiveIndicesBool = [rand() < agentAct for agentAct in AgentAct]
    indx = [k[1] for k in findall(ActiveIndicesBool)]

    for i in indx
        pick_from = setdiff(1:N, i)
        ContactedAgents = rand(pick_from, m)
        AdjMatrix[i, ContactedAgents] .= 1
        AdjMatrix[ContactedAgents, i] .= 1 # connections are bilateral!
    end

    AgentStateNext = AgentState

    InfectedIndicesBool = [agentState == 1 for agentState in AgentState]
    InfectedIndices = [k[1] for k in findall(InfectedIndicesBool)]

    for j in InfectedIndices
        ConnectedSusceptibles = AdjMatrix[j, :] .& (.~AgentState) # find individuals connected to an infected that are susceptible
        ConnSuscIndices = findall(ConnectedSusceptibles) # find their indices
        AgentStateNext[ConnSuscIndices] = [rand() < lambda for v in ConnSuscIndices] # infection process S - > I on the susceptible and connected;
    end

    InfectedIndices = [k[1] for k in findall(AgentStateNext)]# find Infected again and
    AgentStateNext[InfectedIndices] = [~(rand() < mu) for v in InfectedIndices] # process I - > S (it is independent from the network connections)
    push!(ITrend, Int(sum(AgentState)))
    if t > (TMinConvergence + TWindow)
        termination = (abs(mean(ITrend[t - TWindow:t - 1]) - mean(ITrend[t - TWindow + 1:t])) < tolerance);
    else
        termination = false
    end

    return AgentStateNext, ITrend, termination

end
function number_infected_eta(N, AgentState, AgentAct, lambda, eta, EtaMax,  mu, m,  tolerance, t, ITrend, TWindow, TMinConvergence)
    AdjMatrix = zeros(Bool, N, N)

    AgentAct = eta * (ones(N) - AgentState) .* AgentAct + AgentState .* AgentAct * EtaMax
    ActiveIndicesBool = [rand() < agentAct for agentAct in AgentAct]
    indx = [k[1] for k in findall(ActiveIndicesBool)]

    for i in indx
        pick_from = setdiff(1:N, i)
        ContactedAgents = rand(pick_from, m)
        for contact in ContactedAgents
            AdjMatrix[i, contact] = 1
            AdjMatrix[contact, i] = 1 # connections are bilateral!
        end
    end

    AgentStateNext = AgentState

    InfectedIndicesBool = [agentState == 1 for agentState in AgentState]
    InfectedIndices = [k[1] for k in findall(InfectedIndicesBool)]

    for j in InfectedIndices
        ConnectedSusceptibles = zeros(Bool, N)
        for k in 1:N
            ConnectedSusceptibles[k] = AdjMatrix[j, k] & (~AgentState[k]) 
        end # find individuals connected to an infected that are susceptible
        ConnSuscIndices = findall(ConnectedSusceptibles) # find their indices
        for k in ConnSuscIndices
            AgentStateNext[k] = rand() < lambda 
        end# infection process S - > I on the susceptible and connected;
    end

    InfectedIndices = [k[1] for k in findall(AgentStateNext)]# find Infected again and
    for k in InfectedIndices
        AgentStateNext[k] = (rand() < mu) 
    end
    # process I - > S (it is independent from the network connections)
    push!(ITrend, Int(sum(AgentStateNext)))
    if t > (TMinConvergence + TWindow)
        termination = (abs(mean(ITrend[t - TWindow:t - 1]) - mean(ITrend[t - TWindow + 1:t])) < tolerance);
    else
        termination = false
    end

    return AgentStateNext, ITrend, termination

end



function infected_c(i_time, AgentState, AgentAct, m, mu, lambda, etaS, eta, InfectedOverATrial, tmax, TWindow)
    # AgentState = zeros(Bool, N)
    N = length(AgentAct)
    # AgentState = zeros(Bool, N)
    # AgentState[1:Ni] .= 1
    AgentStateNext = AgentState
    for i in 1:N
        if AgentState[i] == 1
            if rand() < mu
                AgentStateNext[i] = 0
            end
        end
        Activity = (1 - AgentState[i]) * (etaS * AgentAct[i]) + AgentState[i] * (eta * AgentAct[i])
        if Activity > rand()
            for j in 1:m
                rInt = Int(round(rand() * (N - 1) + 1))
                if AgentState[i] == 1
                    if ((AgentState[rInt] == 0) && (rand() < lambda) )
                        AgentStateNext[rInt] = 1
                    end
                else
                    if ((AgentState[rInt] == 1) && (rand() < lambda) )
                        AgentStateNext[i] = 1
                    end
                end
            end
        end
    end
    NumInf = 0
    for i = 1:N
        AgentState[i] = AgentStateNext[i]
        if AgentState[i] == 1
            NumInf = NumInf + 1
        end
    end
    if i_time > (tmax - TWindow)
        InfectedOverATrial = InfectedOverATrial + NumInf
    end
    return InfectedOverATrial, AgentState
 
end
    
function avg_infected(N, NTrials, m, mu, eta, lambda, etaS,  AgentAct, tmax, TWindow)
    InfectedOverTrials = 0
    for nt in 1:NTrials
        AgentState = zeros(Bool, N)
        AgentState[1:Int(.01 * N)] .= 1
        InfectedOverATrial = 0.0
        for i_time in 1:tmax
            InfectedOverATrial, AgentState = infected_c(i_time, AgentState, AgentAct, m, mu, lambda, etaS, eta, InfectedOverATrial, tmax, TWindow)
        end
        InfectedOverATrial = InfectedOverATrial / TWindow
        InfectedOverTrials = InfectedOverTrials + InfectedOverATrial
    end
    return InfectedOverTrials / NTrials
end

end