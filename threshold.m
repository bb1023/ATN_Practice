function [alpha] = threshold(beta, gamma,r)
    delta=1/gamma;
%     beta=( (1-delta) ./(1-alpha.^2));
    alpha=sqrt((beta+delta-1) ./(r .* beta));
end

