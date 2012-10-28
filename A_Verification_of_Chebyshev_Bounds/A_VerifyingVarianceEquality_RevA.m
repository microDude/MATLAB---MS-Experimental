% Generate Lambda
n = (1:400)';
alpha = 25;
lambda = 10000.*n.^2.*exp(-n.^2/(2*alpha^2))/alpha^3;

% Generate Poisson Realizations
R = 100000; %number of expierments
y = zeros(length(lambda),R);
for i = 1:R
   y(:,i) = poissrnd(lambda);
end

% Calculate T
% L2-norm
T = sum((y-lambda*ones(1,R)).^2); 

% Evaluate Equality
%[mean(T),sum(lambda)]
%[std(T),sqrt(sum(2*lambda.^2+lambda))] % Standard Deviation
j = [var(T),(sum(2*lambda.^2+lambda))]      % Variance
diff = j(1) - j(2) 