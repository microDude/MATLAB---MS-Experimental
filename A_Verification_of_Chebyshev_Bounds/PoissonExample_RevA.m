clear;
clc;

alpha = 0.5; % Max Amplitude of Lambda
n = (1:1000); % Number of samples in Lambda

N = 100;      % Number of experiments

lambda = alpha/length(n)*n;                  % Lambda
lambda_matrix = lambda'*ones(1,N);  % Lambda matrix

%Generate y
y = poissrnd(lambda_matrix);

%L2-norm
T = sum((lambda_matrix - y).^2);

%Var[T] vs. sum of Lambdas
LeftSide = var(T);
RightSide_single = sum(2*lambda.^2+lambda);
RightSide_matrix = sum(sum(2*lambda_matrix.^2+lambda_matrix));

% Psi
Psi = abs(LeftSide-RightSide_single);


