using Statistics
using ProgressMeter
using MATLAB
using MAT
using Dates
include("atn_module.jl")

@everywhere begin
    using Statistics
    using ProgressMeter
    using MATLAB
    using MAT
    using Dates
    include("atn_module.jl")
    N = 10_000          # number of nodes
    m = 5               # number of active links per active node
    mu = 0.1            # probability I -> S
    tmax = 1_000         # maximum duration of the simulation
    Ni = Int(0.01 * N)  # initial number of infects
    gamma = -2.1        # exponent of the power law
    NTrials = 5        # number of trials for each value of lambda
    TWindow = 200       # time window for evaluating steady state
    TMinConvergence = 100
    LambdaSize = 1000
    LambdaMin = 0.0
    LambdaMax = 1.0
    LambdaStep = 0.05
    # Lambda = LambdaMin:LambdaStep:LambdaMax
    # LambdaSize = length(Lambda)
    Lambda = range(LambdaMin, stop = LambdaMax, length = LambdaSize)
    EtaSize = 1000
    EtaMin = 0.0
    EtaMax = 15
    EtaStep = 0.5
    # Eta = EtaMin:EtaStep:EtaMax
    Eta = range(EtaMin, stop = EtaMax, length = EtaSize)
    ActMin = 0.001      # minimum activation probability
    ActMax = 1
    AgentAct =
        (
            (ActMax .^ (gamma + 1) - ActMin .^ (gamma + 1)) .* rand(N) .+
            ActMin .^ (gamma + 1)
        ) .^ (1.0 / (gamma + 1))
    Activity = [AgentAct]
    N_iters = EtaSize * LambdaSize
    for k = 1:(N_iters-1)
        push!(Activity, AgentAct)
    end
    Activity = reshape(Activity, LambdaSize, EtaSize)
    # AvgDegree = 2 * m * mean(AgentAct)
    # AgentState = zeros(Bool, N)
    # AgentState[1:Ni] .= 1
    etaS = 15
end
#  EtaSize,LambdaSize
# grid = reshape([(eta, lambda) for lambda in Lambda, eta in Eta], :)
# grid = reshape([(eta, lambda) for eta in Eta, lambda in Lambda], :)
grid = [(eta, lambda) for lambda in Lambda, eta in Eta]
AvgInfectedEtaLambda = @showprogress pmap(
    atn_module.avg_infected_par,
    N * ones(Int, LambdaSize, EtaSize),
    NTrials * ones(Int, LambdaSize, EtaSize),
    m * ones(Int, LambdaSize, EtaSize),
    mu * ones(LambdaSize, EtaSize),
    grid,
    etaS * ones(LambdaSize, EtaSize),
    Activity,
    tmax * ones(LambdaSize, EtaSize),
    TWindow * ones(LambdaSize, EtaSize),
)

# AvgInfectedEtaLambda =Float64.(transpose(reshape(AvgInfectedEtaLambda, EtaSize, LambdaSize)))
# AvgInfectedEtaLambda = Float64.(transpose(reshape(AvgInfectedEtaLambda, LambdaSize,EtaSize)))
# AvgInfectedEtaLambda = Float64.(transpose(reshape(AvgInfectedEtaLambda, EtaSize,LambdaSize)))
# AvgInfectedEtaLambda = Float64.(AvgInfectedEtaLambda)


# new_1=reshape(grid, 5, 3)
location_mat = "/mnt/storage/MAT/"
right_now = replace(
    replace(replace(string(Dates.now()), "." => "_"), ":" => "_"),
    "-" => "_",
)

file_name = "practice" * right_now
mat_name = location_mat * file_name * ".mat"

matwrite(
    mat_name,
    Dict("AvgInfectedEtaLambda" => AvgInfectedEtaLambda);
    compress = true,
)
# AvgInfectedEtaLambda = Float64.(transpose(AvgInfectedEtaLambda))
mat"figure();set(gcf, 'Position',  [0, 0, 2500, 1500]); hold on;"
mat"axis([ -0,1,-0,10 ])"
mat"imagesc([0,1],[0,10],$AvgInfectedEtaLambda)"
mat"colormap jet"
