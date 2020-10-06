using Statistics
using ProgressMeter
using MATLAB
using MAT
using Dates
include("atn_coupled_module.jl")
@everywhere begin
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
    lambda_size = 200
    lambda_min = 0.0
    lambda_max = 1.0
    n_trials = 5
    Lambda = range(lambda_min, stop = lambda_max, length = lambda_size)
    beta_size = 200
    beta_min = 0.0
    beta_max = 1.0
    Beta = range(beta_min, stop = beta_max, length = beta_size)
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
    activity_1 = [agent_act_1]
    activity_2 = [agent_act_2]

    N_iters = beta_size * lambda_size
    for k = 1:(N_iters-1)
        push!(activity_1, agent_act_1)
        push!(activity_2, agent_act_2)
    end
    activity_1 = reshape(activity_1, lambda_size, beta_size)
    activity_2 = reshape(activity_2, lambda_size, beta_size)


    mu_1 = 0.1
    mu_2 = 0.1

    lambda_1 = 1
    # lambda_2 = 0.3

    beta_1 = 0.9
    # beta_2 = beta_1
end

grid = [(beta, lambda) for lambda in Lambda, beta in Beta]

Iend_2 = @showprogress pmap(
    atn_coupled_module.parameter_sweep,
    n_trials * ones(Int, lambda_size, beta_size),
    activity_1,
    activity_2,
    mu_1 * ones(Int, lambda_size, beta_size),
    mu_2 * ones(Int, lambda_size, beta_size),
    lambda_1 * ones(Int, lambda_size, beta_size),
    # beta_1 * ones(Int, lambda_size, beta_size),
    grid,
)


location_mat = "/mnt/storage/MAT/"
right_now = replace(
    replace(replace(string(Dates.now()), "." => "_"), ":" => "_"),
    "-" => "_",
)

file_name = "ATN_SIS_COUPLED" * right_now
mat_name = location_mat * file_name * ".mat"

matwrite(mat_name, Dict("Iend_2" => Iend_2); compress = true)
Iend_2 = Float64.(Iend_2)
mat"figure();set(gcf, 'Position',  [0, 0, 2500, 1500]); hold on;"
mat"axis([ -0,1,-0,10 ])"
mat"imagesc([0,1],[0,10],$Iend_2)"
mat"colormap jet"
