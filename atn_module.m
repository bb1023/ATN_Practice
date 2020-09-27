function number_infected()

    if t > (TMinConvergence + TWindow)
        termination = (abs(mean(ITrend[t - TWindow:t - 1]) - mean(ITrend[t - TWindow + 1:t])) < tolerance);

        if termination
            break
        end

    end

    return AgentState, ITrend, terminantion

end
