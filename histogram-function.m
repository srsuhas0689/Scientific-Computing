%%  Copyright Erik Dali, GPL-3.0 License
%%  https://github.com/erik-dali

%%  1.1 Histogram validation
clear all; format compact; format long e;

%   Declar variables
rng(198127);
x_min = 0;
x_max = 10;
number_bins = 50;
number_samples = 2500;

%   Create the histogram
[mybins bin_frq f_x ci_p mymean, var_dist, var_mean]=hst(@exponential_sampler, number_bins, number_samples, x_min,x_max,'y');


%   Plot the actual distribution
hold on
fplot(@exponential_pdf,[x_min x_max],'g', 'LineWidth',1);
hold off
%   Figure Options
set(gca,'ycolor','k')
ylabel({'Probability'},'Interpreter','Latex', 'FontSize',12)
legend({'Histogram','95\% CI','$e^{-t}$'},'Interpreter','Latex', 'FontSize',12, 'Location','North East')

%   Save figures as EPSC
%   Naming convention: "Problem_#_Tittle_Erik_Dali"
saveas(1,'Problem_1_1_Histogram_Erik_Dali','epsc')


%   Define the histogram function
function [bins, frq, fx, cip, mu, var_dist, var_mean] = hst(sampler, num_bins, num_samp, a, b, plotYN)
    bin_bounds = linspace(a,b,num_bins+1)';
    samples = sampler(num_samp);
    bins = bin_bounds(1:(end-1));
    frq = zeros(num_bins,1); %  Pre-allocate for the frequency
    fx = zeros(num_bins,1); %  Pre-allocate for the f(x)

    
    %   Compute mean
    mu = sum(samples)/num_samp;
    
    %   Compute variance
    var_dist = (1/num_samp)*sum( (samples-mu).^2); %    Variance of the sample
    var_mean = (1/num_samp)*var_dist; % Variance of the mean estimator
    
    %   Compute the frequencies
    for i = 1:(num_bins-1)
        bin_idx = ( (bin_bounds(i) < samples) & (samples < bin_bounds(i+1)) );
        frq(i) = sum( bin_idx );
    end
    bin_idx = (bin_bounds(end) < samples);
    frq(end) = sum( bin_idx );
    
    %   Compute the empirical probabilities
    dx = bin_bounds(2)-bin_bounds(1);
    fx = (1/dx)*(1/num_samp)*frq; %    Probability of x_i
    
    %   Compute the confidence interval
    cip = (2/dx)*sqrt(frq)/num_samp;
    
    if plotYN == 'y'
        %   Plot histogram
        figure('Name','Histogram', 'WindowStyle','docked')
        bar(bins,fx, 'histc');

        %   Figure Options
        title('Histogram','FontSize',12, 'FontWeight','normal')
        xlabel({'x'},'Interpreter','Latex','FontSize',12)
        ylabel({'Frequency'},'Interpreter','Latex', 'FontSize',12)
        xlim([a b])
        pbaspect([1.75 1 1]) % Aspect Ratio


        %   Plot errors
        hold on
        errorbar(bins+(dx/2),fx,cip,'r','LineStyle','none');
        hold off

        %   Figure Options
        ylim([0 inf])
        legend({'Histogram','95% CI'},'FontSize',12, 'Location','North East')
    end
end

%   The sampler function for the exponential distribution with lambda =1
function V = exponential_sampler(num_samp)
    U = rand(num_samp,1);
    V = -log(U);
end

%   The pdf for the exponential distribution with lambda =1
function v = exponential_pdf(u)
    v = exp(-u);
end
