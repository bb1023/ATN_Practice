using Statistics
N = 10000         #number of nodes
m = 3               #number of active links per active node
mu = 0.1            #probability I -> S
tmax = 6000         #maximum duration of the simulation
Ni = Int(0.1*N)     #initial number of infects
ActMin = 0.001      #minimum activation probability
ActMax = 1          #maximum activation probability
gamma = -2.22       #exponent of the power law

AgentAct=rand(N,1);
AgentAct = ((ActMax.^(gamma+1) - ActMin.^(gamma+1)).*AgentAct.+ ActMin.^(gamma+1)).^(1. /(gamma+1))
AvgDegree = 2*m*mean(AgentAct)

NTrials = 5         #number of trials for each value of lambda

tolerance = 1e-15   #termination tolerance
TWindow = 1000      #time window for evaluating steady state
TMinConvergence=500


#health state of agents S=0 I=1
AgentState=zeros(N)
AgentStateNext=zeros(N)
AgentState[1:Ni].=1.0

#diffferent lambda
# lambda=0.1:0.1:1;
# IAvgLambda=zeros(length(lambda))
# beta=zeros(length(lambda))

lambda=.2
termination=1.0
t=1;


termination = false 
t=2
AdjMatrix=zeros(N,N)
ActiveIndicesBool=[ rand()<agentAct for agentAct in AgentAct]

println(sum(ActiveIndicesBool))
indx=[ k[1] for k in findall(ActiveIndicesBool)]
for i in indx
    pick_from=setdiff(1:N, i)
    ContactedAgents=rand(pick_from,m)
    AdjMatrix[i, ContactedAgents] .= 1
    AdjMatrix[ContactedAgents,i] .= 1 # connections are bilateral!
end

# connections are bilateral!
AgentStateNext=AgentState

InfectedIndicesBool = [agentState==1 for agentState in AgentState]
InfectedIndices=[ k[1] for k in findall(InfectedIndicesBool)]

for j in InfectedIndices
    ConnectedSusceptibles  = double(AdjMatrix(InfectedIndices(j),:) & (~AgentState')); #find individuals connected to an infected that are susceptible
    ConnSuscIndices = find(ConnectedSusceptibles == 1); #find their indices
    AgentStateNext(ConnSuscIndices)=(rand(length(ConnSuscIndices),1)<lambda); #infection process S->I on the susceptible and connected;
end
