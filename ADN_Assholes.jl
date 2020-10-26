using Statistics
using ProgressMeter
using MATLAB
using MAT
using Dates
include("ATN_Coupled_assholes.jl")
@everywhere begin
    using Statistics
    using ProgressMeter
    using MATLAB
    using MAT
    using Dates
    include("ATN_Coupled_assholes.jl")
    N = 10_000          # number of nodes
    m = 5               # number of active links per active node
    tmax = 1000         # maximum duration of the simulation
    Ni = Int(0.01 * N)  # initial number of infects
    gamma = -2.1        # exponent of the power law
    alpha_size = 100
    alpha_min = 0.0
    alpha_max = 1
    n_trials = 5
    Alpha = range(alpha_min, stop = alpha_max, length = alpha_size)
    proportion_size = 100
    proportion_min = 0.1
    proportion_max = 0.9
    Proportion =
        range(proportion_min, stop = proportion_max, length = proportion_size)
    ActMin = 0.001      # minimum activation probability
    ActMax = 1
    agent_act =
        (
            (ActMax .^ (gamma + 1) - ActMin .^ (gamma + 1)) .* rand(N) .+
            ActMin .^ (gamma + 1)
        ) .^ (1.0 / (gamma + 1))
    activity = [agent_act]
    N_iters = proportion_size * alpha_size
    for k = 1:(N_iters-1)
        push!(activity, agent_act)

    end
    activity = reshape(activity, alpha_size, proportion_size)

end

grid = [(proportion, alpha) for alpha in Alpha, proportion in Proportion]

Dead = @showprogress pmap(
    ATN_Coupled_SIRD.parameter_sweep,
    n_trials * ones(Int, alpha_size, proportion_size),
    activity,
    grid,
)


location_mat = "/mnt/storage/MAT/"
right_now = replace(
    replace(replace(string(Dates.now()), "." => "_"), ":" => "_"),
    "-" => "_",
)

file_name = "ATN_AssHats" * right_now
mat_name = location_mat * file_name * ".mat"
Dead = Float64.(Dead)
matwrite(mat_name, Dict("Dead" => Dead); compress = true)

mat"figure();set(gcf, 'Position',  [0, 0, 2500, 1500]); hold on;"
# mat"axis([ 0,15,0,10 ])"
# mat"axis([ 0,1,0,1 ])"
mat"[X,Y]=meshgrid($Proportion,$Alpha);contourf(X,Y,$Dead,15)"

mat"colorbar;colormap('inferno')"
