using Statistics
using MATLAB
include("atn_module.jl")
# using .time_series

N = 10_000          # number of nodes
m = 30               # number of active links per active node
mu = 0.1            # probability I -> S
tmax = 1000          # maximum duration of the simulation
Ni = Int(0.1 * N)   # initial number of infects
gamma = -2.22       # exponent of the power law
NTrials = 5         # number of trials for each value of lambda
tolerance = 1e-15   # termination tolerance
TWindow = 1000      # time window for evaluating steady state
TMinConvergence = 500

# health state of agents S=0 I=1

# diffferent lambda
# lambda=0.1:0.1:1;
# IAvgLambda=zeros(length(lambda))
# beta=zeros(length(lambda))

lambda = 0.2

ITrend = atn_module.time_series(Ni, N, m, gamma, lambda, mu, tolerance, TWindow,
                                TMinConvergence, tmax);

mat"plot($ITrend)"