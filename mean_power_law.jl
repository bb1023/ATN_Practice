using Statistics
N=100000;
val=zeros(N)
for i=1:N
    gamma=-2.1 # minimum activation probability
    ActMax = 1
    agent_act =
        (
            (ActMax .^ (gamma + 1) - ActMin .^ (gamma + 1)) .* rand(N) .+
            ActMin .^ (gamma + 1)
        ) .^ (1.0 / (gamma + 1))
    val[i]=mean(agent_act)+sqrt(mean(agent_act .* agent_act))
end
# mmmean=   
# s=1.5104924655094145

# s1=-0.019221243887172085



println(mean(val))
