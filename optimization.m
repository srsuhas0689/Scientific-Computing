%%  Copyright Erik Dali, GPL3 Licesnse

%%  2.1 Synthetic Data
clear all; format compact; format long e;

%   Fix the seed number
rng(719271)


m = 100; %  Number of points
epsilon = 1e-2; %   Magnitude of errors
x = sort(rand(m,1)*10); %    Generate random numbers between 0 and 10
y0 = zeros(m,1); %    Pre-allocate for y=f(x,c)
y = zeros(m,1); %    Pre-allocate for y=f(x,c)+error


%   Generate the synthetic data
y0 = exp(-x/2).*sin(2*x);
y = y0 + epsilon*randn(m,1);


%   Plot y
figure('Name','Function', 'WindowStyle','docked')
plot(x,y0,'b-', x,y,'ro--', 'LineWidth',1);

%   Figure Options
legend({'f(x;c)','f(x;c)+\delta'},'FontSize',12, 'Location','North East')
title('Function Fitting','FontSize',12, 'FontWeight','normal')
xlabel('x','FontSize',12)
ylabel('y','FontSize',12,'Rotation',0)
set(gca, 'YGrid','on','XGrid','off','YMinorTick','off','XMinorTick','off','YMinorGrid','off');
%yticks([3 pi (pi-3)+pi 2*(pi-3)+pi 3*(pi-3)+pi]);
%xticks(1:n)
set(gca,'linewidth',1)
pbaspect([1.75 1 1]) % Aspect Ratio

%   Save figures as EPSC
%   Naming convention: "Problem_#_Tittle_Erik_Dali"
saveas(1,'Problem_2_1_Synthetic_Erik_Dali','epsc')


%%  2.2 Gauss-Newton Method
clear all; format compact; format long e;

%   Fix the seed number
rng(719271)

nrm = 2; %  Norm used in calculations
m = 100; %  Number of points
n = 4; %   Dimension of unknown parameters
k = 20; %  Number of iterations
%epsilon = 1e-2; %   Magnitude of errors
epsilon = 0; %   Magnitude of errors
guess = [0.95 0.45 1.95 0.05]'; %   Initial guess
%guess = [1 1 1 1]'; %   Initial guess
%guess = [0.5 -0.1 1 0.5]'; %   Initial guess
x = sort(rand(m,1)*10); %    Generate random numbers between 0 and 10
y0 = zeros(m,1); %    Pre-allocate for y=f(x,c)
y = zeros(m,1); %    Pre-allocate for y=f(x,c)+error
c0 = [1 1/2 2 0]'; %   The actual parameters
c = zeros(n,k); % Pre-allocate for the c_k (by row) amplitude c1, decay c2, period c3,  phase c4
J = zeros(m,n); % Pre-allocate for the Jacobian (by column) amplitude c1, decay c2, period c3,  phase c4
c(:,1) = guess; % Set the initial guess
K = zeros(k,1); % Pre-allocate for the conditioning number of the (J.')*J matrix
err = zeros(1,k); % Pre-allocate for the error
C = 0; %   Pre-allocate for C

%   Define the function
f = @(d) d(1)*exp(-d(2)*x).*sin((d(3)*x) + d(4));

%   Define the derivates
fc1 = @(d) exp(-d(2)*x).*sin(d(3)*x + d(4));
fc2 = @(d) -x*d(1).*exp(-d(2)*x).*sin(d(3)*x + d(4));
fc3 = @(d) x*d(1).*exp(-d(2)*x).*cos(d(3)*x + d(4));
fc4 = @(d) d(1).*exp(-d(2)*x).*cos(d(3)*x + d(4));

%   Generate the synthetic data
y0 = f(c0);
y = y0 + (epsilon*randn(m,1));


%   Apply the Gauss-Newton method
for i = 1:k
    J(:,1) = fc1( c(:,i) );
    J(:,2) = fc2( c(:,i) );
    J(:,3) = fc3( c(:,i) );
    J(:,4) = fc4( c(:,i) );
    
    %   Applying the normal method
    LHS = (J')*J;
    RHS = (J')*(y0 - f(c(:,i)) );
    R = chol(LHS);
    z = R'\RHS;
    delta_c_k = R\z;
    c(:,i+1) = c(:,i) + delta_c_k;
    K(i) = cond(LHS,2);
end


%   Convergence
for i = 1:(k-1)
    err(:,i) = norm(c(:,i+1)-c0,nrm)/(norm(c(:,i)-c0,nrm)^1.1);
end


%%  2.3 Levenberg-Marquardt Algorithm
clear all; format compact; format long e;

%   Fix the seed number
rng(719271)

nrm = 2; %  Norm used in calculations
m = 100; %  Number of points
n = 4; %   Dimension of unknown parameters
k = 20; %  Number of iterations
epsilon = 1e-2; %   Magnitude of errors
guess = [1 1 1 1]'; %   Initial guess
x = sort(rand(m,1)*10); %    Generate random numbers between 0 and 10
y0 = zeros(m,1); %    Pre-allocate for y=f(x,c)
y = zeros(m,1); %    Pre-allocate for y=f(x,c)+error
c0 = [1 1/2 2 0]; %   The actual parameters
c = zeros(n,k); % Pre-allocate for the c_k (by column) amplitude c1, decay c2, period c3,  phase c4
J = zeros(m,n); % Pre-allocate for the Jacobian (by column) amplitude c1, decay c2, period c3,  phase c4
c(:,1) = guess; % Set the initial guess
K = zeros(k,1); % Pre-allocate for the conditioning number of the (J.')*J matrix
lambda = zeros(k,1); %   Pre-allocate for lambda
lambda(1) = 0.00511; %   Initial lambda


%   Define the function
f = @(d) d(1).*exp(-d(2).*x).*sin(d(3).*x + d(4));

%   Define the derivates
fc1 = @(d) exp(-d(2).*x).*sin(d(3).*x + d(4));
fc2 = @(d) -x.*d(1).*exp(-d(2).*x).*sin(d(3).*x + d(4));
fc3 = @(d) x.*d(1).*exp(-d(2).*x).*cos(d(3).*x + d(4));
fc4 = @(d) d(1).*exp(-d(2).*x).*cos(d(3).*x + d(4));

%   Generate the synthetic data
y0 = f(c0);
y = y0 + (epsilon*randn(m,1));


%   Apply the Gauss-Newton method
for i = 1:k
    J(:,1) = fc1( c(:,i) );
    J(:,2) = fc2( c(:,i) );
    J(:,3) = fc3( c(:,i) );
    J(:,4) = fc4( c(:,i) );
    
    %   Applying the normal method
    LHS = ((J')*J) + lambda(i)*diag(diag((J')*J));
    RHS = (J')*(y0 - f(c(:,i)) );
    R = chol(LHS);
    z = R'\RHS;
    delta_c_k = R\z;
    c(:,i+1) = c(:,i) + delta_c_k;
    K(i) = cond(LHS,2);
    
    lambda(i+1) = lambda(i)/2;
    
    %   Check if overall fit is better to half the lambda, or to double it
    %   if the fit is worse
    %err1 = norm(y - f(c(:,i)),nrm);
    %err2 = norm(y - f(c(:,i+1)),nrm);
    %if err2 < err1
    %    lambda(i+1) = lambda(i)/2;
    %else
    %    lambda(i+1) = 2*lambda(i);
    %end
end


%   Plot y
x=(1:k);
figure('Name','Convergence', 'WindowStyle','docked')
plot(x,c(1,x),'bo-', x,c(2,x),'rs-', x,c(3,x),'g^-', x,c(4,x),'cd-', 'LineWidth',1);

%   Figure Options
legend({'c_{1}','c_{2}','c_{3}','c_{4}'},'FontSize',12, 'Location','North West')
title('Convergence','FontSize',12, 'FontWeight','normal')
xlabel('iteration','FontSize',12)
ylabel('estimate','FontSize',12,'Rotation',90)
set(gca, 'YGrid','on','XGrid','off','YMinorTick','off','XMinorTick','off','YMinorGrid','off');
yticks([-2 -1 0 1/2 1 2 3 4]);
%xticks(1:n)
set(gca,'linewidth',1)
pbaspect([1.75 1 1]) % Aspect Ratio

%   Save figures as EPSC
%   Naming convention: "Problem_#_Tittle_Erik_Dali"
saveas(1,'Problem_2_3_Convergence_Erik_Dali','epsc')
