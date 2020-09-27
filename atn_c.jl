using Statistics
using ProgressMeter
using MATLAB
include("atn_module.jl")


N = 10_000          # number of nodes
m = 5               # number of active links per active node
mu = 0.1            # probability I -> S
tmax = 3000         # maximum duration of the simulation
Ni = Int(0.01 * N)  # initial number of infects
gamma = -2.1      # exponent of the power law
NTrials = 50        # number of trials for each value of lambda
tolerance = 1e-10   # termination tolerance
TWindow = 500       # time window for evaluating steady state
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
AvgInfectedEtaLambda = zeros(EtaSize, LambdaSize)
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
# AgentState = zeros(Bool, N)
# AgentState[1:Ni] .= 1
etaS = 15

# grid=Reshape([(eta, lambda) for eta in Eta, lambda in Lambda],:)
@showprogress for (ietaI, eta) in enumerate(Eta)
    for (ilambda, lambda) in enumerate(Lambda)

        AvgInfectedEtaLambda[ietaI,ilambda] = atn_module.avg_infected(N, NTrials, m, mu, eta, lambda, etaS,  AgentAct, tmax, TWindow) / N
    end
end                                                                                                                                      