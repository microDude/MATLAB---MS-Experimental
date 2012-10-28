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
%H1 = sum(2*lambda.^2+lambda);
%H2 = sum(lambda);

% Experiments
Re = 300;                           % number of experiments
Rt = 100;                          % number of captured spectra per experiment
alpha = 15;
lambda_matrix = lambda*ones(1,Rt); % Generate the Lambda matrix (saves computation time)
y = zeros(length(lambda),Rt);      % Poisson realizations of lambda

% Run simulation
Prob_k1 = zeros(length(k),1);
for ki = 1:length(k)
    
    %J = (H2 + k(ki)*sqrt(H1));     % Precalculate right side of inequality
    
    Prob_k2 = zeros(alpha,1);
    for pq = 1:alpha
        Prob_lambda = zeros(Re,1);     % Probability Vector
 
        for j = 1:Re
            y = poissrnd(lambda_matrix);        % Generate Poisson Realizations 
            
            Sy = var(sum(y,1));
            Gamma = Sy + k(ki)^2/2 + sqrt((Sy*k(ki)^2+(k(ki)^4/4)));
            J = Gamma + k(ki)*sqrt(2*Gamma^2+Gamma);
            
            T = var(sum((y-lambda_matrix).^2)); % Calculate T, L2-norm
            
            %pause;
            
            % Add up the probability
            if (T >= J)
                Prob_lambda(j) = 1;
            end
        end
        Prob_k2(pq) = nnz(Prob_lambda)/length(Prob_lambda);
    end
    Prob_k1(ki) = sum(Prob_k2)/alpha;
end

% Plot the Results
figure(1);
plot(k,Prob_k1,'b');hold on;plot(k,(1./(1+k.^2)),'--r');
title('Probability Bound Plot');
xlabel('k');

