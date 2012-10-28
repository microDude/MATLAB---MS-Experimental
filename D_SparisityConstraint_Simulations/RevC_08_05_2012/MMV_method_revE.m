clear all;
close all;
clc;

tic;
% Turn off warnings about Toeplitz matrix
warning('off','all');

% Run the MMV algorithm
%     ---------------------
%     | Generate Test Data |
%     ---------------------
%              |
%              v
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
Riter = 1;                 % Number of Statistical Expirements
m = 100;                   % Length of each spectra
n = 50;                    % Number of captured spectra per experiment
niter = 150;               % Number of regression iterations for langrangian Hunting
niter_final = 800;         % Number of regression iterations for the final solution
ki = 3;                    % The standard deviation that was found to yield 91% = P(T>E) < 1/(K^2+1)
langrangian.niter = 30;    % Number of hunting steps to find the langrangian
langrangian.min = 1e-8; % Starting Lagrangian min
langrangian.max = 1e-2;     % Starting Lagrangian max

L0_norm_diff = zeros(Riter,1); % Number of differnt non-zero indexes between true_x and x_est
x_est_saved = zeros(m,Riter);  % x_est for Statistical Expirement 'R'

for R = 1:Riter
    fprintf('Current Statistical Expirement Iteration : %d \n',R);
    
    % Section:01 ----------------------------------------------------
    % Generate Test Data
    
    % Generate the Convolution Kernel, Maxwell-Boltzman
    b = (1:m);
    A = cell(1,n);
    
    for t=1:n
        alpha = 4 + (20-4).*(t-1)/(n-1);                 % Controls the fatness of the distribution
        amp = 90 + (110-90).*(t-1)/(n-1);                % Controls the amplitude
        a = b.^2.*exp(-b.^2/(2*alpha^2))/alpha^3; % Generate Maxwell-Boltzman Curve
        a = (amp/max(a)).*a;                      % Normalize to defined amplitude
        A{1,t} = tril(toeplitz(a));         % Cell array which holds A_i
    end
    
    % Generate Sparse x vector (same non-zero indexes)
    x = zeros(m,n);
  %  x([10,12,20,28,52],:) = 4 + (5-4).*rand(5,n); %Uniform random values assigned to fixed x values
    x([10,12],:) = 2 + (10-2).*rand(2,n); %Uniform random values assigned to fixed x values
    
    % Populate the observation space
    y = zeros(m,n);
    for t=1:n
        lambda = A{1,t}*x(:,t);      % Generate lambda_i
        y(:,t) = poissrnd(lambda);   % Generate y_i
    end
       
    % Section:02 ----------------------------------------------------
    % Calculate epsilon
    Sy = sum(sum(y,1));
    Gamma = Sy + ki^2/2 + sqrt((Sy*ki^2+(ki^4/4)));
    epsilon = Gamma + ki*sqrt(2*Gamma^2+Gamma);
    epsilon = 1*epsilon;
    fprintf('Calculated Epsilson : %d \n',epsilon);
    
    % Section:03 ----------------------------------------------------
    % Run Iterations
    % Update langrangian on succesive iterations
    
    Xsol = rand(m,n);        %intialize the Xsolution with random entries
    
    % Preallocations for lagrangian hunting
    constraint = zeros(langrangian.niter,1);
    lagrangian = zeros(langrangian.niter,1);
    StateMachine.now        = 0;
    StateMachine.previous   = 0;
        % --------------------------------------
        % StateMachine = 0 , intial
        % StateMachine = 1 , Adjust Lambda Min
        % StateMachine = 2 , Adjust Lambda Max
        % --------------------------------------
    
    for iterL = 1:langrangian.niter
        iterL
        % Define/Update lambda
        lagrangian(iterL) = (langrangian.max+langrangian.min)/2;
        lagrangian(iterL)
        % Solve for X
        %if (iterL > 25)
        %    Xsol = fSolveForX_revB(A,y,lagrangian(iterL),niter,(spones(Xsol).*(0+(1-0).*rand(m,n))));
        %else
            Xsol = fSolveForX_revC(A,y,lagrangian(iterL),niter,(1+(5-1).*rand(m,n)),1);
        %end
        
        % Check against Constraint
        %constraint(iterL) = constraint(iterL)+(norm(A*Xsol-y,2)^2);
        %constraint(iterL) = constraint(iterL)+sum(sum((A*Xsol-y).^2));
        sqr_error = zeros(1,n);
        for t=1:n
            sqr_error(t) = sum(sum((A{1,t}*Xsol(:,t)-y(:,t)).^2));
        end
        constraint(iterL) = constraint(iterL) + sum(sqr_error);
        
        % Adjust langrangian
%         switch StateMachine.now
%             case 0 % intial
%                 if (constraint(iterL) > epsilon)
%                     StateMachine.now = 2;
%                 else
%                     StateMachine.now = 1;
%                 end
%             case 1 % Adjust Lambda Min
%                 if (StateMachine.now == StateMachine.previous)
%                     if ((constraint(iterL) < epsilon) && constraint(iterL) >= constraint(iterL-1))
%                         langrangian.min = lagrangian(iterL);
%                     else
%                         langrangian.min = lagrangian(iterL - 1);
%                         StateMachine.now = 2;
%                     end
%                 else
%                     StateMachine.previous = StateMachine.now;
%                     if ((constraint(iterL) < epsilon))
%                         langrangian.min = lagrangian(iterL);
%                     end
%                 end
%             case 2 % Adjust Lambda Max
%                 if (StateMachine.now == StateMachine.previous)
%                     if ((constraint(iterL) > epsilon) && (constraint(iterL) <= constraint(iterL-1)))
%                         langrangian.max = lagrangian(iterL);
%                     else
%                         langrangian.max = lagrangian(iterL - 1);
%                         StateMachine.now = 1;
%                     end
%                 else
%                     StateMachine.previous = StateMachine.now;
%                     if ((constraint(iterL) > epsilon))
%                         langrangian.max = lagrangian(iterL);
%                     end
%                 end
%             otherwise
%                 disp('StateMachine Error');
%         end

if (constraint(iterL) > epsilon)
    langrangian.min=lagrangian(iterL);
else
    langrangian.max=lagrangian(iterL);
end
    end
    
    % After Correct Lagrangian is found, Calculate the Xsol with many iterations
    %Xsol = fSolveForX_revB(A,y,lagrangian(iterL),niter_final,(1+(5-1).*rand(m,n)));
    Xsol = fSolveForX_revC(A,y,lagrangian(iterL),niter_final,(1+(5-1).*rand(m,n)),1);
    
    % Section:05 ----------------------------------------------------
    % Calculate L0-norm
    
    % Form a single vector for the result
    eps3 = 1e-2;
    x_est = sum(Xsol.^2,2) > eps3;
    true_x = sum(x.^2,2) > eps3;
    
%     % Truncate off any elements that are below the mean
%     index.x_est_truncate = (x_est < mean(x_est)/2);
%     x_est(index.x_est_truncate) = 0;
%         
%     % Match the index locations for x and x_est
%     % If index match, reduce L0-norm penality
%     % If no index match, add a L0-norm penality
%     index.x = find(true_x);
%     index.penality = 0;
%     for i=1:length(index.x)
%         if (x_est(index.x(i)) < 2*eps)
%             index.penality = index.penality + 1; % No Match
%         else
%             index.penality = index.penality - 1; % Match
%         end
%     end
    
    % Save the L0-norm distance
    L0_norm_diff(R) = nnz(x_est-true_x);
    x_est_saved(:,R) = x_est;
end

%Save Workspace
save 'xyz';
toc;

warning('on','all');

% end
