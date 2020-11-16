using Statistics
using ProgressMeter
using MATLAB
using MAT
using Dates
include("ADN_Coupled_assholes_threshold.jl")


@everywhere begin
    using Statistics
    using ProgressMeter
    using MATLAB
    using MAT
    using Dates
    include("ADN_Coupled_assholes_threshold.jl")
    N = 10_000          # number of nodes
    m = 5               # number of active links per active node
    tmax = 1000         # maximum duration of the simulation
    Ni = Int(0.01 * N)  # initial number of infects
    gamma = -2.1        # exponent of the power law
    n_trials = 100

    proportion_size = 120
    proportion_min = 0.01
    proportion_max = 0.99
    Proportion = range(proportion_min, stop = proportion_max, length = proportion_size)
    ActMin = 0.001      # minimum activation probability
    ActMax = 1
    agent_act = ((ActMax.^(gamma + 1) - ActMin.^(gamma + 1)) .* rand(N) .+
                 ActMin.^(gamma + 1)).^(1.0 / (gamma + 1))
    activity = [agent_act]
    # N_iters = proportion_size * alpha_size
    for k in 1:proportion_size - 1
        push!(activity, agent_act)
    end
    # activity = reshape(activity, alpha_size, proportion_size)
    activity = reshape(activity, proportion_size)

end


extinct = @showprogress pmap(ADN_Coupled_assholes_threshold.parameter_sweep,
                          n_trials * ones(Int, proportion_size), activity, Proportion)

location_mat = "/mnt/storage/MAT/"
right_now = replace(replace(replace(string(Dates.now()), "." => "_"), ":" => "_"),
                    "-" => "_")

file_name = "ADN_extinct_" * right_now
mat_name = location_mat * file_name * ".mat"
matwrite(mat_name, Dict("extinct" => extinct); compress = true)
mat"figure();set(gcf, 'Position',  [0, 0, 2500, 1500]); hold on;"
mat"plot($Proportion, $extinct,'b--')"