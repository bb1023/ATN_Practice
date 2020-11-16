% uiopen('/home/brandonbehring/Dropbox/github/ATN_Practice/adn4.fig',1)
close all
alpha_size = 100;
alpha_min = 0.0;
alpha_max = 1;

proportion_size = 100;
proportion_min = 0.01;
proportion_max = 0.99;
Alpha = linspace(alpha_min,  alpha_max,alpha_size);
Proportion =linspace(proportion_min, proportion_max, proportion_size);
figure();set(gcf, 'Position',  [0, 0, 2500, 1500]); hold on;
axis([proportion_min, proportion_max, alpha_min, alpha_max])
axis equal
[X,Y]=meshgrid(Proportion,Alpha);contourf(X,Y,Dead,15)

colorbar;colormap('inferno')
result=compute_threshold(Dead,Alpha);
plot(Proportion,result, 'w--', 'LineWidth',5)
Alpha = linspace(alpha_min,  alpha_max,1000);
lambda_mu=fa(.5,.8,Alpha)/2.2655;
hold on
% result=real(threshold( .6, 0.8828/2, Proportion));
plot(Alpha,.1*lambda_mu,  'LineWidth',10)
% plot(Proportion,result,  'LineWidth',5)
% plot(Proportion,result, 'w--', 'LineWidth',5)

% title("Meta-population: Good, $\alpha=0.5$,$\beta=.1$",'FontSize',40)

title(" Proportion Infected: GOOD ",'FontSize',40)
% title("Meta-population: Good, $\alpha=0.25$",'FontSize',40)
% xlabel("$\eta^{(g)}/\eta^{(b)}$",'FontSize',40)
ylabel("$\lambda$",'FontSize',40)

% xlabel("$\alpha$",'FontSize',40)
xlabel("$\beta$",'FontSize',40)

% ylabel("$\alpha$",'FontSize',40)

colormap('inferno')
axis([proportion_min, proportion_max, alpha_min, alpha_max])
