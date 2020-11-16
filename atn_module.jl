module atn_module

using Statistics

function infected_c(
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
    N = length(AgentAct)
    AgentStateNext = zeros(Bool, N)
    # AgentState = zeros(Bool, N)

    # AgentState = zeros(Bool, N)
    # AgentState[1:Ni] .= 1
    for k = 1:N
        AgentStateNext[k] = AgentState[k]
    end
    for i = 1:N
        if AgentState[i] == 1
            if rand() < mu
                AgentStateNext[i] = 0
            end
        end
        # Activity =
        #     (1 - AgentState[i]) * (etaS * AgentAct[i]) +
        #     AgentState[i] * (eta * AgentAct[i])
        Activity =
            (1 - AgentState[i]) * (eta * AgentAct[i]) +
            AgentState[i] * (eta * AgentAct[i])

        if Activity > rand()
            for j = 1:m
                rInt = rand(1:N)
                if AgentState[i] == 1
                    r = rand()
                    if ((AgentState[rInt] == 0) && (rand() < lambda))
                        AgentStateNext[rInt] = 1
                    end
                else
                    if ((AgentState[rInt] == 1) && (rand() < lambda))
                        AgentStateNext[i] = 1
                    end
                end
            end
        end
    end
    for k = 1:N
        AgentState[k] = AgentStateNext[k]
        # if AgentState[k] == 1
        #     NumInf = NumInf + 1
        # end
    end

    return AgentState

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

end
