using Statistics
using ProgressMeter
using MATLAB
include("atn_module.jl")

N = 10_000          # number of nodes
m = 5               # number of active links per active node
mu = 0.1            # probability I -> S
tmax = 1000         # maximum duration of the simulation
Ni = Int(0.01 * N)  # initial number of infects
gamma = -2.1      # exponent of the power law
NTrials = 1       # number of trials for each value of lambda
TWindow = 500.0      # time window for evaluating steady state
TMinConvergence = 500

LambdaMin = 0.0
LambdaMax = 1.0
LambdaStep = 0.05
Lambda = LambdaMin:LambdaStep:LambdaMax
LambdaSize = length(Lambda)

EtaMin = 0.0
EtaMax = 15
EtaStep = 0.5
Eta = EtaMin:EtaStep:EtaMax
EtaSize = length(Eta)
AvgInfectedLambda = zeros(LambdaSize)
Ni = Int(0.01 * N)  # initial number of infects
gamma = -2.1       # exponent of the power law
TWindow = 500       # time window for evaluating steady state
TMinConvergence = 500
ITrend = zeros(Int, 0)
ActMin = 0.001      # minimum activation probability
ActMax = 1
AgentAct =
    (
        (ActMax .^ (gamma + 1) - ActMin .^ (gamma + 1)) .* rand(N) .+
        ones(N) * ActMin .^ (gamma + 1)
    ) .^ (1.0 / (gamma + 1))
# AvgDegree = 2 * m * mean(AgentAct)
# AgentState = zeros(Bool, N)
# AgentState[1:Ni] .= 1
etaS = 1
eta = 1

# grid=Reshape([(eta, lambda) for eta in Eta, lambda in Lambda],:)
# @showprogress for (ietaI, eta) in enumerate(Eta)
# @showprogress for (ilambda, lambda) in enumerate(Lambda)
lambda = 1
InfectedOverTrials = zeros(0)
InfectedOverATrial = zeros(0)

AgentState = zeros(Int, N)
AgentStateNext = zeros(Int, N)
for k = 1:100
    AgentState[k] = 1
end
# InfectedOverATrial = 0.0
for i_time = 1:tmax
    for l = 1:N
        AgentStateNext[l] = AgentState[l]
    end
    for i = 1:N
        if (AgentState[i] == 1)
            if (rand() < mu)
                AgentStateNext[i] = 0
            end
        end
        # println(i)
        if (AgentAct[i] > rand())
            for j = 1:m
                rInt = rand(1:N)
                # println(rInt)
                if (AgentState[i] == 1)
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
    for l = 1:N
        AgentState[l] = AgentStateNext[l]
    end
    NumInf = sum(AgentState)
    push!(InfectedOverTrials, NumInf)
end
# InfectedOverATrial = InfectedOverATrial / TWindow
# push!(InfectedOverTrials, InfectedOverATrial)
# end


mat"figure();set(gcf, 'Position',  [0, 0, 2500, 1500]); hold on;"
# # mat"axis([ -0,1,-0,10 ])"
# # mat"imagesc([0,1],[0,10],$AvgInfectedEtaLambda)"
# # mat"colormap jet"
#
mat"plot(1:3000,$InfectedOverTrials,'.-','MarkerSize',20)"
