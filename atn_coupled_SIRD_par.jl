using Statistics
using ProgressMeter
using MATLAB
using MAT
using Dates
include("ATN_Coupled_SIRD.jl")
@everywhere begin
    using Statistics
    using ProgressMeter
    using MATLAB
    using MAT
    using Dates
    include("ATN_Coupled_SIRD.jl")
    N = 100_000          # number of nodes
    m = 5               # number of active links per active node
    tmax = 1000         # maximum duration of the simulation
    Ni = Int(0.01 * N)  # initial number of infects
    gamma = -2.1        # exponent of the power law
    lambda_size = 20
    lambda_min = 0.0
    lambda_max =1
    n_trials = 1
    Lambda = range(lambda_min, stop = lambda_max, length = lambda_size)
    beta_size =20
    beta_min = 0.0
    beta_max = 1
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
end

grid = [(beta, lambda) for lambda in Lambda, beta in Beta]

Dead= @showprogress pmap(
    ATN_Coupled_SIRD.parameter_sweep,
    n_trials * ones(Int, lambda_size, beta_size),
    activity_1,
    activity_2,
    grid,
)


location_mat = "/mnt/storage/MAT/"
right_now = replace(
    replace(replace(string(Dates.now()), "." => "_"), ":" => "_"),
    "-" => "_",
)

file_name = "ATN_SIRD_COUPLED" * right_now
mat_name = location_mat * file_name * ".mat"
Dead = Float64.(Dead)
matwrite(mat_name, Dict("Dead" => Dead); compress = true)

mat"figure();set(gcf, 'Position',  [0, 0, 2500, 1500]); hold on;"
mat"axis([ -0,1,-0,10 ])"
mat"imagesc([0,1],[0,10],$Dead,[0,.5])"
mat"colorbar;colormap('inferno')"
