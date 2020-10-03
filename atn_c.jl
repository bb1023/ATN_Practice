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
AgentAct = ((ActMax .^ (gamma + 1) - ActMin .^ (gamma + 1)) .* rand(N) .+
            ActMin .^ (gamma + 1)) .^ (1.0 / (gamma + 1))
# AvgDegree = 2 * m * mean(AgentAct)
# AgentState = zeros(Bool, N)
# AgentState[1:Ni] .= 1
etaS = 15
eta = 1

lambda = 0.5
InfectedOverTrials = zeros(0)
InfectedOverATrial = zeros(0)
Activitylist = zeros(0)
for nt in 1:NTrials
    AgentState = zeros(Bool, N)
    AgentState[1:Int(0.01 * N)] .= 1
    InfectedOverATrial = 0.0
    for i_time in 1:tmax
        AgentStateNext = AgentState
        for i in 1:N
            if AgentState[i] == 1
                r = rand()
                if (r < mu)
                    AgentStateNext[i] = 0
                r = rand()
                if (r < mu)
            Activity = AgentState[i]
            push!(Activitylist, Activity)
            r = rand()
            if (Activity > r)
                for j in 1:m
                    rInt = rand(1:N)
                    if (AgentState[i] == 1)
                        r = rand()
                        if ((AgentState[rInt] == 0) && (r < lambda))
                    if (AgentState[i] == 1)
                        r = rand()
                        if ((AgentState[rInt] == 0) && (r < lambda))
                            # println(rInt)
                            AgentStateNext[rInt] = 1
                        end
                    else
                        r = rand()
                        if ((AgentState[rInt] == 1) && (r < lambda))
                            # println(rInt)
                end
            end
        end
        NumInf = 0
        for k in 1:N
            AgentState[k] = AgentStateNext[k]
            if (AgentState[k] == 1)
                NumInf = NumInf + 1
            endntState[k] = AgentStateNext[k]
            if (AgentState[k] == 1)
                NumInf = NumInf + 1
            end
        end
        push!(InfectedOverTrials, NumInf)
    end
    # InfectedOverATrial = InfectedOverATrial / TWindow
    # push!(InfectedOverTrials, InfectedOverATrial)
end

# println(sum(InfectedOverTrials) / NTrials)
# AvgInfectedLambda[ilambda] = InfectedOverTrials / NTrials
# end
# end
# AvgInfectedEtaLambda = Float64.(transpose(AvgInfectedEtaLambda))
# mat"figure();set(gcf, 'Position',  [0, 0, 2500, 1500]); hold on;"
# mat"axis([ -0,1,-0,10 ])"
# mat"imagesc([0,1],[0,10],$AvgInfectedEtaLambda)"
# mat"colormap jet"
# location_mat = "/mnt/storage/MAT/"
# right_now = replace(
#     replace(replace(string(Dates.now()), "." => "_"), ":" => "_"),
#     "-" => "_",
# )

# file_name = "practice" * right_now
# mat_name = location_mat * file_name * ".mat"
#
# matwrite(
#     mat_name,
#     Dict("AvgInfectedLambda" => AvgInfectedLambda);
#     compress = true,
# )
# AvgInfectedEtaLambda = Float64.(transpose(AvgInfectedEtaLambda))
mat"figure();set(gcf, 'Position',  [0, 0, 2500, 1500]); hold on;"
# # mat"axis([ -0,1,-0,10 ])"
# # mat"imagesc([0,1],[0,10],$AvgInfectedEtaLambda)"
# # mat"colormap jet"
#
mat"plot(1:3000,$InfectedOverTrials,'.-','MarkerSize',20)"
