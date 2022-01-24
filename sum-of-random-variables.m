%%  1.2 Simple sampler
clear all; format compact; format long e;

%   Declar variables
rng(198127);
x_min = 0;
x_max = 10;
number_bins = 50;
number_samples = 1e4;
number_bins2 = 100;
number_samples2 = 1e5;

%   For 1e4 observations
[bins1, frq1, fx1, cip1, mean1, var_dist1, var_mean1]= hst(@new_exponential_sampler, number_bins, number_samples, x_min,x_max,'n');
mean1
var_dist1
mean_CI = [mean1-(2*sqrt(var_mean1)) mean1+(2*sqrt(var_mean1))] % Confidence Interval of the mean


%   For 1e5 observations and 100 bins
[bins2, frq2, fx2, cip2, mean2, var_dist2, var_mean2] = hst(@new_exponential_sampler, number_bins2, number_samples2, x_min,x_max,'y');
mean2
var_dist2
mean2_CI = [mean2-(2*sqrt(var_mean2)) mean2+(2*sqrt(var_mean2))] % Confidence Interval of the mean


%   Plot the actual distribution
hold on
fplot(@new_exponential_pdf,[x_min x_max],'g', 'LineWidth',1);
hold off
%   Figure Options
set(gca,'ycolor','k')
ylabel({'Probability'},'Interpreter','Latex', 'FontSize',12)
legend({'Histogram','95\% CI','$te^{-t}$'},'Interpreter','Latex', 'FontSize',12, 'Location','North East')

%   Save figures as EPSC
%   Naming convention: "Problem_#_Tittle_Erik_Dali"
saveas(1,'Problem_1_2_Histogram_Erik_Dali','epsc')


%   Define the histogram function
%   For this part we ignore the values larger than the x_max
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
        bin_idx = ( (bin_bounds(i) <= samples) & (samples < bin_bounds(i+1)) );
        frq(i) = sum( bin_idx );
    end
    %bin_idx = (bin_bounds(end) < samples);
    %frq(end) = sum( bin_idx );
    
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

%   The new sampler function using the trick
%   It draws two sets of independent observations and adds them Y=X_1+X_2
function V = new_exponential_sampler(num_samp)
    U = rand(num_samp,2);
    V = -log(U(:,1))-log(U(:,2));
end

%   The pdf for the new exponential distribution
function v = new_exponential_pdf(u)
    v = u.*exp(-u);
end
