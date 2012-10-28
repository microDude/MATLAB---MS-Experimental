clear all;
close all;

% Generate a set of Lambda Amplitudes
amp_space = linspace(0.1,1000,100);
STDP = zeros(length(amp_space),1);

% Test different Amplitudes
for q = 1:length(amp_space)
    % Generate Lambda
    n = (1:400)';
    alpha = 25;
    amp = amp_space(q);
    lambda = n.^2.*exp(-n.^2/(2*alpha^2))/alpha^3;
    lambda = (amp/max(lambda)).*lambda;
    
    % Precalculated Right Side of the equality
    H = sum(2*lambda.^2+lambda);
    
    % Experiments
    Re = 1000;                         % number of experiments
    Rt = 100;                          % number of captured spectra per experiment
    lambda_matrix = lambda*ones(1,Rt); % Generate the Lambda matrix (saves computation time)
    y = zeros(length(lambda),Rt);      % Poisson realizations of lambda
    T = zeros(Rt,1);                   % L2-norm
    Psi = zeros(1,Re);                 % error of Chebyshev Bound
    for j = 1:Re
        % Generate Poisson Realizations
        y = poissrnd(lambda_matrix);
        
        % Calculate T, L2-norm
        T = sum((y-lambda_matrix).^2);
        
        % Calculate the error of the Chebyshev Bound
        Psi(j) = var(T)-H;
    end
    STDP(q) = std(Psi);
end

% Plot the Results
% figure(1);
% histfit(Psi,50);title(['\Psi (error in Chebyshev Bound) for ',num2str(Re),...
%     ' Experiments with ',num2str(Rt),' spectra per Experiment and \lambda normalized to ',num2str(amp)]);
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','b','EdgeColor','w');

% Evaluate Equality
%[mean(T),sum(lambda)]
%[std(T),sqrt(sum(2*lambda.^2+lambda))] % Standard Deviation
