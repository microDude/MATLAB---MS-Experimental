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
Riter = 5;            % Number of Statistical Observations
m = 100;            % Length of each spectra
n = 5;             % Number of captured spectra per experiment
niter = 100;          % Number of regression iterations
ki = 0.121;         % The standard deviation that was found to yield 91% = P(T<E) < 1/(K^2+1)

L0_norm_diff = zeros(Riter,1);
epsilon_track = zeros(Riter,1);

for R = 1:Riter
    % Section:01 ----------------------------------------------------
    % Generate Test Data
    
    % Generate the Convolution Kernel, Maxwell-Boltzman
    b = (1:m);
    alpha = 5;  % Controls the fatness of the distribution
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
    
    langrangian.niter = 100; % number of iterations to find the langrangian
    
    langrangian.min = 1e-8; % Starting Lagrangian min
    langrangian.max = 50;  % Starting Lagrangian max
    
    Xsol = rand(m,n);        %intialize the Xsolution with random entries
    
    % Preallocations
    constraint = zeros(langrangian.niter,1);
    lagrangian = zeros(langrangian.niter,1);
    StateMachine.now        = 0;
    StateMachine.previous   = 0;
    % StateMachine = 0 , intial
    % StateMachine = 1 , Adjust Lambda Min
    % StateMachine = 2 , Adjust Lambda Max
    
    for iterL = 1:langrangian.niter
        
        % Define/Update lambda
        lagrangian(iterL) = (langrangian.max+langrangian.min)/2;
        
        % Solve for X
        if (iterL > 20)
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
            case 1 % Adjust Lambda Min
                if (StateMachine.now == StateMachine.previous)
                    if ((constraint(iterL) < epsilon) && constraint(iterL) >= constraint(iterL-1))
                        langrangian.min = lagrangian(iterL);
                    else
                        langrangian.min = lagrangian(iterL - 1);
                        StateMachine.now = 2;
                    end
                else
                    StateMachine.previous = StateMachine.now;
                    if ((constraint(iterL) < epsilon))
                        langrangian.min = lagrangian(iterL);
                    end
                end
            case 2 % Adjust Lambda Max
                if (StateMachine.now == StateMachine.previous)
                    if ((constraint(iterL) > epsilon) && (constraint(iterL) <= constraint(iterL-1)))
                        langrangian.max = lagrangian(iterL);
                    else
                        langrangian.max = lagrangian(iterL - 1);
                        StateMachine.now = 1;
                    end
                else
                    StateMachine.previous = StateMachine.now;
                    if ((constraint(iterL) > epsilon))
                        langrangian.max = lagrangian(iterL);
                    end
                end
            otherwise
                disp('StateMachine Error');
        end
        
    end
    
    % Section:05 ----------------------------------------------------
    % Calculate L0-norm
    
    % Form a single vector for the result
    x_est = sum(Xsol,2);
    true_x = sum(x,2);
    
    % Truncate off any elements that are below the mean
    x_trunc_index = (x_est < mean(x_est));
    x_est(x_trunc_index) = 0;
    
    % Calculate the L0-norm
    L0.x_est = nnz(x_est);
    L0.true_x = nnz(true_x);

    % Save the L0-norm distance and the epsilon used to generate it
    L0_norm_diff(R) = (L0.x_est - L0.true_x);
    epsilon_track(R) = epsilon;
end


% Plot the results
% if(debug == true)
%     figure('Name','MMV results for estimate of x','NumberTitle','off');
%     stem(true_x*(max(x_est)/max(true_x)),'b');hold on;stem(x_est,'r');hold off;
%     title('MMV Results');
%     xlabel('eV sample');ylabel('counts');
% end
