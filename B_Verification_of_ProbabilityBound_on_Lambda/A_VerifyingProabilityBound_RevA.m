clear all;
close all;
clc;

% Define the STD shifts
k = (0:0.1:10);

% Generate Lambda
n = (1:400)';
alpha = 25;
amp = 1;
lambda = n.^2.*exp(-n.^2/(2*alpha^2))/alpha^3;
lambda = (amp/max(lambda)).*lambda;

% Precalculated Right Side of the inequality
H1 = sum(2*lambda.^2+lambda);
H2 = sum(lambda);

% Experiments
Re = 100;                           % number of experiments
Rt = 100;                          % number of captured spectra per experiment
lambda_matrix = lambda*ones(1,Rt); % Generate the Lambda matrix (saves computation time)
y = zeros(length(lambda),Rt);      % Poisson realizations of lambda

% Run simulation
Prob_k = zeros(length(k),1);
for ki = 1:length(k)
    
    Prob_lambda = zeros(Re,1);     % Probability Vector
    
    J = (H2 + k(ki)*sqrt(H1));     % Precalculate right side of inequality
    
    for j = 1:Re
        % Generate Poisson Realizations
        y = poissrnd(lambda_matrix);
        
        % Calculate T, L2-norm
        T = var(sum((y-lambda_matrix).^2));
        
        % Add up the probability
        if (T >= J)
            Prob_lambda(j) = 1;
        end
        
    end
    
    Prob_k(ki) = nnz(Prob_lambda)/length(Prob_lambda);
    
    % Debug
%     disp(' ');
%     disp(['k = ',num2str(k(ki))]);
%     disp(['J = ',num2str(J)]);
%     disp(['P(k) = ',num2str(Prob_k(ki))]);
%     pause;
    
end

% Plot the Results
figure(1);
plot(k,Prob_k,'b');hold on;plot(k,(1./(1+k.^2)),'r');
title('Probability Bound Plot');
xlabel('k');

