using Statistics
using ProgressMeter
using MATLAB
using MAT
using Dates
include("atn_coupled_module.jl")

N = 10_000          # number of nodes
m = 5               # number of active links per active node
tmax = 1000         # maximum duration of the simulation
Ni = Int(0.01 * N)  # initial number of infects
gamma = -2.1        # exponent of the power law

ActMin = 0.001      # minimum activation probability
ActMax = 1
agent_act_1 =
    (
        (ActMax .^ (gamma + 1) - ActMin .^ (gamma + 1)) .* rand(N) .+
        ActMin .^ (gamma + 1)
    ) .^ (1.0 / (gamma + 1))

agent_act_2 =
    (
        (ActMax .^ (gamma + 1) - ActMin .^ (gamma + 1)) .* rand(N) .+
        ActMin .^ (gamma + 1)
    ) .^ (1.0 / (gamma + 1))




mu_1 = 0.1
mu_2 = 0.01

lambda_1 = 0.2
lambda_2 = 0.3

beta_1 = 0.9
beta_2 = beta_1

Iend_1, Iend_2 = atn_coupled_module.infected_totals(
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


location_mat = "/mnt/storage/MAT/"
right_now = replace(
    replace(replace(string(Dates.now()), "." => "_"), ":" => "_"),
    "-" => "_",
)

file_name = "practice" * right_now
mat_name = location_mat * file_name * ".mat"

# matwrite(
#     mat_name,
#     Dict("AvgInfectedEtaLambda" => AvgInfectedEtaLambda);
#     compress = true,
# )
time = Float64.(1:tmax)
mat"figure(); hold on; plot($time,$Iend_1,'.-');plot($time,$Iend_2,'.-')"
