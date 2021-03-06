using Distributed
using Statistics
using ProgressMeter
using MATLAB
include("atn_module.jl")

@everywhere begin
    using Statistics
    using ProgressMeter
    using MATLAB
    include("atn_module.jl")
    # using .time_series_par
    N = 10_000          # number of nodes
    m = 5               # number of active links per active node
    mu = 0.1            # probability I -> S
    tmax = 10000         # maximum duration of the simulation
    Ni = Int(0.01 * N)  # initial number of infects
    gamma = -2.1   # exponent of the power law
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
    EtaMax = 5
    EtaStep = 0.5
    Eta = 5
    EtaSize = length(Eta)
end

infected = @showprogress pmap(
    atn_module.time_series_par,
    [(eta, lambda) for eta in Eta, lambda in Lambda],
    Ni * ones(Int, EtaSize, LambdaSize),
    N * ones(Int, EtaSize, LambdaSize),
    m * ones(Int, EtaSize, LambdaSize),
    EtaMax * ones(EtaSize, LambdaSize),
    mu * ones(EtaSize, LambdaSize),
    tmax * ones(Int, EtaSize, LambdaSize),
)

infected = float64.(transpose(infected))
mat"imagesc($infected)"
