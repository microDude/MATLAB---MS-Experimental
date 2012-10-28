clear all;
close all;
clc;

% Run the MMV algorithm
%     ---------------------
%     | Calculate Epsilon |
%     ---------------------
%              |
%              v
%     ---------------------
%     | Run with Update   |
%     |  (Adjust Lambda)  |
%     ---------------------
%              |
%              v
%     ---------------------
%     | return L0 Norm    |
%     | of X_estimate     |
%     ---------------------

% Define Constants
R = 100;            % Number of Statistical Observations
m = 200;            % Length of each spectra
n = 5;             % Number of captured spectra per experiment
niter = 100;          % Number of regression iterations
ki = 0.121;         % The standard deviation that was found to yield 91% = P(T<E) < 1/(K^2+1) 

% Section:01 ----------------------------------------------------
% Generate Test Data

% Generate the Convolution Kernel, Maxwell-Boltzman
b = (1:m);
alpha = 15;  % Controls the fatness of the distribution
amp = 50;    % Controls the amplitude
a = b.^2.*exp(-b.^2/(2*alpha^2))/alpha^3;
a = (amp/max(a)).*a;
A = toeplitz(zeros(1,m),a)';

% Generate Sparse x vector (same non-zero indexes)
x = zeros(m,n);
x([10,12,20,28,52],:) = 1 + (5-1).*rand(5,n); %Uniform random values assigned to fixed x values

% Generate Lambda
lambda = A*x;

% Generate y
y = poissrnd(lambda);

% Section:02 ----------------------------------------------------
% Calculate epsilon
Sy = sum(sum(y,1));
Gamma = Sy + ki^2/2 + sqrt((Sy*ki^2+(ki^4/4)));
epsilon = Gamma + ki*sqrt(2*Gamma^2+Gamma);
fprintf('Calculated Epsilson : %d \n',epsilon);

% Section:03 ----------------------------------------------------
% Run Iterations
%  Update Lambda on succesive iterations

niter_langrangian = 100; % number of iterations to find the langrangian

lagrangianMin = 1e-8; % Starting Lagrangian min
lagrangianMax = 100;  % Starting Lagrangian max

Xsol = rand(m,n);

debug = true;

% Preallocations
constraint = zeros(niter_langrangian,1);
lagrangian = zeros(niter_langrangian,1);
StateMachine.now        = 0;
StateMachine.previous   = 0;
   % StateMachine = 0 , intial
   % StateMachine = 1 , Adjust Lambda Min
   % StateMachine = 2 , Adjust Lambda Max

for iterL = 1:niter_langrangian
    
    % Define/Update lambda
    lagrangian(iterL) = (lagrangianMax+lagrangianMin)/2;
    
    % Solve for X
    if (iterL>20)
        Xsol = fSolveForX_revB(A,y,lagrangian(iterL),niter,Xsol);
    else
        Xsol = fSolveForX_revB(A,y,lagrangian(iterL),niter,rand(m,n));
    end 
    
    % Check against Constraint
    %constraint(iterL) = constraint(iterL)+(norm(A*Xsol-y,2)^2);
    constraint(iterL) = constraint(iterL)+sum(sum((A*Xsol-y).^2));
    
    switch StateMachine.now
        case 0 % intial
            if (constraint(iterL) > epsilon)
                StateMachine.now = 2;
            else
                StateMachine.now = 1;
            end
            str = sprintf('Intial State on iter: %d, Move to State: %d, Intial Constraint: %d  epsilon: %d ',iterL,StateMachine.now, constraint(iterL), epsilon);
        case 1 % Adjust Lambda Min
            if (StateMachine.now == StateMachine.previous)
                if ((constraint(iterL) < epsilon) && constraint(iterL) >= constraint(iterL-1))
                    lagrangianMin = lagrangian(iterL);
                    str = sprintf('  Adjust lagrangianMin on iter: %d, Constraint: %d < epsilon: %d ',iterL, constraint(iterL), epsilon);
                else
                    lagrangianMin = lagrangian(iterL - 1);
                    StateMachine.now = 2;
                    str = sprintf('Switch to State2 on iter: %d, Adjust lagrangianMin back, Constraint: %d  epsilon: %d ',iterL, constraint(iterL), epsilon);
                end
            else
                StateMachine.previous = StateMachine.now;
                str = 'State-Transition';
                if ((constraint(iterL) < epsilon))
                    lagrangianMin = lagrangian(iterL);
                    str = sprintf('State-Transitioned on iter: %d, Adjust lagrangianMin, Constraint: %d < epsilon: %d ',iterL, constraint(iterL), epsilon);
                end
            end
        case 2 % Adjust Lambda Max
            if (StateMachine.now == StateMachine.previous)
                if ((constraint(iterL) > epsilon) && (constraint(iterL) <= constraint(iterL-1)))
                    lagrangianMax = lagrangian(iterL);
                    str = sprintf('  Adjust lagrangianMax on iter: %d, Constraint: %d > epsilon: %d ', iterL,constraint(iterL), epsilon);
                else
                    lagrangianMax = lagrangian(iterL - 1);
                    StateMachine.now = 1;
                    str = sprintf('Switch to State1 on iter: %d, Adjust lagrangianMax back, Constraint: %d  epsilon: %d ',iterL, constraint(iterL), epsilon);
                end
            else
                StateMachine.previous = StateMachine.now;
                str = 'State-Transition';
                if ((constraint(iterL) > epsilon))
                    lagrangianMax = lagrangian(iterL);
                    str = sprintf('State-Transitioned on iter: %d, Adjust lagrangianMax, Constraint: %d > epsilon: %d ', iterL,constraint(iterL), epsilon);
                end
            end
        otherwise
            disp('StateMachine Error');
    end
    
    if (debug == true)
        disp(str);
        fprintf('     Constraint-Epsilson error: %d \n',abs(constraint(iterL) - epsilon));
        x_est = sum(Xsol,2);
        true_x = sum(x,2);
        figure(902);
                stem(true_x*(max(x_est)/max(true_x)),'b');hold on;stem(x_est,'r');hold off;
                title(['MMV Lagrangian Iternation Number = ',num2str(iterL)]);
                xlabel('eV sample');ylabel('counts');
    end
    
end

% Calculate the L0 norm (Sum along rows) of Xsol matrix
x_est = sum(Xsol,2);
true_x = sum(x,2);

% Plot the results
if(debug == true)
    close(902);
    figure('Name','MMV results for estimate of x','NumberTitle','off');
    stem(true_x*(max(x_est)/max(true_x)),'b');hold on;stem(x_est,'r');hold off;
    title('MMV Results');
    xlabel('eV sample');ylabel('counts');
end
