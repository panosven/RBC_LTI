%-----------------------------------------------------------------------

% Basic RBC (McKay's website)
% Solution Method: Linear Time Iteration (Rendahl 2017)
% Panagiotis Veneris, University of Liverpool
% Date: 26/10/2020

%------------------------------------------------------------------------


clear all;
close all;
clc;


%Define Parameters
alppha = 0.36;         % capital share
gamma =  1;            % Risk aversion
delta = 0.025;         % Depreciation rate
rho_z = 0.95;          % shock persistence
beta = 0.97;           % discount factor
it = 0;


% Define Steady State  
R = 1/beta;


% Define Coefficients 
coeff0 = (R-1+delta)/R;
coeff1 = (1/gamma) * coeff0;
coeff2 = ((alppha-1)/gamma) * coeff0;
coeff3 = R;
coeff4 = (R-1+delta)/alppha;
coeff5 = (R-1+(delta*(1-alppha)))/alppha;


%---------------------------------------------------
% Matrix form of Linear RE model
%---------------------------------------------------


A = [ 1 0 -coeff1;
    0 0 0;
    0 0 0];

B = [ -1 -coeff2 0;
    coeff5 1 -coeff4;
    0 0 1];

C = [ 0 0 0;
    0 -coeff3 0;
    0 0 -rho_z];

W = [ 0; 0; -1];    % for negative prod.shock remove (-)


%-----------------------------------------------------------------------------------------------
% Main Loop
%-----------------------------------------------------------------------------------------------

Tol = 1e-13;        % tolerance level for time iteration algorithm, default: 1e-13

metric = 1; 
P = 0;              % initial guess for P


tic
while metric>Tol
    it = it+1;
    
    P = -(A*P+B)\C; 
                    
                    
    metric = max(max(abs([A*P*P+B*P+C])));   
                                            
                                         
Q = -inv(A*P+B)*W;   
                   
  fprintf('iteration = %d, metric = %e\n', it, metric);                
end
toc


% IRFs
nresp = 80; % horizon of responses (IRFs)

 
x = zeros(3, nresp);             
                      
u = zeros(1, nresp);             


% Initialize u   
u(:,1)= [1];  
                

for t=2:nresp                         
    x(:,t) = P*x(:,t-1)+ Q*u(:,t-1);                                      
end                                                                      

% Remaining IRFs

k_irf = x(1,2:end);                       
z_irf = x(3,2:end);                      
                                      
y_irf = z_irf + alppha*k_irf;             


R_irf = coeff0*(y_irf - k_irf);           



% Plotting
figure(1);
subplot(3,1,1); plot(x(1,2:end), 'r','Linewidth',2);    
                                                       
title('Consumption');
grid on;

subplot(3,1,2); plot(y_irf, 'm', 'Linewidth',2);
title('Output');
grid on;

subplot(3,1,3); plot(x(3,2:end), 'g', 'Linewidth',2);
title('Productivity');
grid on;

figure(2);
subplot(2,1,1); plot(x(2,2:end), 'b', 'Linewidth',2);
title('Capital');
grid on;

subplot(2,1,2); plot(R_irf, 'c', 'Linewidth',2);
title('Real Interest Rate');
grid on;



