using Statistics
using ProgressMeter
using MATLAB
include("atn_module.jl")
N = 10_000          # number of nodes
m = 5               # number of active links per active node
mu = 0.1            # probability I -> S
tmax = 3000         # maximum duration of the simulation
Ni = Int(0.01 * N)  # initial number of infects
gamma = -2.22       # exponent of the power law
NTrials = 5         # number of trials for each value of lambda
tolerance = 1e-10   # termination tolerance
TWindow = 500      # time window for evaluating steady state
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
IT = zeros(EtaSize, LambdaSize)
@showprogress for (i, lambda) in enumerate(Lambda)
    for (j, eta) in enumerate(Eta)
        ITrend = time_series_eta(Ni, N, m, eta, EtaMax, gamma, lambda, mu, tolerance,
                                 TWindow, TMinConvergence, tmax)
        IT[j, i] = ITrend[end]
    end
end
# mat"plot($ITrend)"