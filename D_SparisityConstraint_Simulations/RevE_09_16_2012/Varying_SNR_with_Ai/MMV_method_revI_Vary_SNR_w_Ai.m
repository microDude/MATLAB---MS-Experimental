% ***************************************
% GRID-SEARCH: MMV Method across parameter space
% 09/28/2012
% revI uses Evgenia's updated methods
% ***************************************
clear all;
close all;
clc;
mainTIC = tic;

% Define Parameter Constants
iter.statistical = 1;                       % Number of Statistical Expirements
iter.final = 200;                           % Number of regression iterations for the final solution
LagrangianHunting.numSteps = 35;            % Number of steps allowed to hone in on the true Lagrangain
LagrangianHunting.iter = 10000;              % Number of regression iterations per step for lagrangian Hunting
LagrangianHunting.min = 1e-8;               % Starting Lagrangian min
LagrangianHunting.max = 1;                  % Starting Lagrangian max
m = 60;                                    % Length of a_i, x_i, y_i, i.e. Length of each spectra
ki = 3;                                     % The standard deviation that was found to yield 91% = P(T>E) < 1/(K^2+1)
Cardx = 2;                                  % The cardinality of the x vector
ampAlpha = 75;                              % Amplitude of convolution kernel 'a'
eps3 = 1e-3;                                % Error threshold for x_est
alpha = 2:0.1:14;                           % Sample space used to generate 'a' and likewise A.

% Grid Search Parameters
n = [20];                                   % Column space of X and y, i.e. number of captured spectra
delta = [5];                                % Card(x) pattern spacing
SNR = 0:30/(n-1):30;                        % SNR (dB) of the signal

% Choose Appropriate Solver
% 1 = fSolveForX_revD   --> Orginal solver
% 2 = fLSSP3_revB       --> New LS solver
% 3 = fMLSP_revC        --> New ML solver
solver.type = 3;

% **************** Begin GRID SEARCH ***********************
% iN            = iterator for Column space of X and y
% iDELTA        = iterator for Card(x) pattern spacing
% iSNR          = iterator for SNR (dB) of the signal
% iVARIETY      = iterator for variety in sample space of alpha
% iSTATISTICAL  = iterator for statistical expirements
nTIC = tic;

% Sample alpha uniformly over n samples
alphaINDEX = round(linspace(1,length(alpha),n));
aVec = alpha(alphaINDEX); % create a vector to loop through to grab alpha values

% --------------------------------------------------
% Generate a linear of SNR to match Ai
% Generate Sparse x vector (same non-zero indexes)
minAMPx = exp(SNR/20); % Figure out the min for the x amplitude
X = zeros(m,n);

for iSNR=1:length(SNR)
    X([5,5+delta],iSNR) = minAMPx(iSNR) + minAMPx(iSNR).*rand(2,1); %Uniform random values assigned to fixed x values
    %XTrueGrid(n,:) = mean(X,2);
end
% --------------------------------------------------

% Generate the Convolution Kernel, Maxwell-Boltzman
b = 1:m;                            % vector used to generate 'a' and then A.
A = cell(1,n);                  % preallocate for A

for t=1:length(A)
    a = b.^2.*exp(-b.^2/(2*aVec(t)^2))/aVec(t)^3; % Generate Maxwell-Boltzman Curve
    a = (ampAlpha/max(a)).*a;                     % Normalize to defined amplitude
    A{t} = tril(toeplitz(a));                     % Cell array which holds A_i
end

XsolSTAT = zeros(m,iter.statistical);             % Preallocate matrix for statistical iterations
for iSTATISTICAL=1:iter.statistical
    timeiSTAT = tic;
    display(['    iSTATISTICAL = ',num2str(iSTATISTICAL),'/',num2str(iter.statistical)]);
    % Section:01 ----------------------------------------------------
    % Populate the observation space
    y = zeros(m,n);
    for t=1:n
        lambda = A{t}*X(:,t);           % Generate lambda_i
        y(:,t) = poissrnd(lambda);      % Generate y_i
    end
    
    % Preallocations for lagrangian hunting
    Xsol = zeros(m,n);        %intialize the Xsolution
    constraint = zeros(LagrangianHunting.numSteps,1);
    lagrangian = zeros(LagrangianHunting.numSteps,1);
    
    % Section:02 ----------------------------------------------------
    % Solve for X
    switch solver.type
        case 1
            % Section:03 ----------------------------------------------------
            % calculate epsilon for LSSP method
            % probability s.t. P(constraint<epsilonLSSP)>1-p
            fp = 0.1;
            Sy = sum(sum(y));
            % epsilon without simplification (look in the notes method 2 epsilon 2)
            fhlp = sqrt((2/fp)-1);
            epsilon = Sy+(fhlp^2)/2+sqrt(fhlp^2*Sy+fhlp^4/4)+fhlp*sqrt(2*Sy^2+Sy*(4*fhlp^2+1)+fhlp^4+fhlp^2/2+(4*Sy+2*fhlp^2+1)*sqrt(fhlp^2*Sy+fhlp^4/4));
            
            % Peform Lagrangian Hunting
            for iterL = 1:LagrangianHunting.numSteps
                
                % Define/Update lambda
                lagrangian(iterL) = (LagrangianHunting.max + LagrangianHunting.min)/2;
                
                % Solve for X
                Xsol = fSolveForX_revD(A,y,lagrangian(iterL),LagrangianHunting.iter,(1+(5-1).*rand(m,n)),0);
                
                % Evaluate the Constraint, sum of squared error (least squares)
                sqr_error = zeros(1,n);
                for t=1:n
                    sqr_error(t) = sum(sum((A{t}*Xsol(:,t)-y(:,t)).^2));
                end
                constraint(iterL) = sum(sqr_error);
                
                % Apply Lagrangian update rule
                if (constraint(iterL) > epsilon)
                    LagrangianHunting.min = lagrangian(iterL);
                else
                    LagrangianHunting.max = lagrangian(iterL);
                end
            end
            
            % After Correct Lagrangian is found, Calculate the Xsol with many iterations
            Xsol = fSolveForX_revD(A,y,lagrangian(iterL),iter.final,(1+(5-1).*rand(m,n)),1);
        case 2
            % Section:03 ----------------------------------------------------
            % calculate epsilon for LSSP method
            % probability s.t. P(constraint<epsilonLSSP)>1-p
            fp = 0.1;
            Sy = sum(sum(y));
            % epsilon without simplification (look in the notes method 2 epsilon 2)
            fhlp = sqrt((2/fp)-1);
            epsilon = Sy+(fhlp^2)/2+sqrt(fhlp^2*Sy+fhlp^4/4)+fhlp*sqrt(2*Sy^2+Sy*(4*fhlp^2+1)+fhlp^4+fhlp^2/2+(4*Sy+2*fhlp^2+1)*sqrt(fhlp^2*Sy+fhlp^4/4));
            
            [Xsol,~,~,constraint,lagrangian] = fLSSP3_revB(A,y,(1+(5-1).*rand(m,n)),epsilon,LagrangianHunting.numSteps,LagrangianHunting.iter,...
                LagrangianHunting.min,LagrangianHunting.max,0);
        case 3
            % Section:03 ----------------------------------------------------
            % calculate epsilon for MLSP method
            % probability of being in the confidence region
            confidProbabMLSP = 0.99;
            % parameter for one-sided Chebyshev inequality
            chebParam = sqrt(1/(1-confidProbabMLSP)-1);
            % bound on expectation of D_ij
            Ce = 0.6;
            % bound on variance of D_ij
            Cv = 0.55;
            % radius of the confidence constrain set
            epsilon = m*n*Ce+chebParam*sqrt(m*n*Cv);
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Calculate matrix Arowsum
            %~~~~~~~~~~~~~~~~~~~~~~~~~~
            % row j of matrix Arowsum contains sum of rows of matrix A{j}
            fArowsum = zeros(n,m);
            for t = 1:length(A)
                fArowsum(t,:) = sum(A{t},1);
            end;
            
            [Xsol,~,~,constraint,lagrangian] = fMLSP_revC(A,fArowsum,y,(1+(5-1).*rand(m,n)),epsilon,LagrangianHunting.iter,...
                LagrangianHunting.numSteps,LagrangianHunting.min,LagrangianHunting.max,0);
        otherwise
            warning(' !!! Unexpected solver type');
    end
    
    % Section:05 ----------------------------------------------------
    x_est = sqrt(sum(Xsol.^2,2));            % Form a single vector for the result
    x_est(x_est < (max(x_est)*eps3)) = 0;    % Zero-out any element less then 3-order of magnitudes of the peak
    
    % Store the result
    XsolSTAT(:,iSTATISTICAL) = x_est;
    display(['Time to complete one solution : ',num2str(toc(timeiSTAT)),' secs']);
end

% ---------------------------------------------------------
% Store the Final Solution for this Grid Point
XSolFinal = mean(XsolSTAT,2);

display(['Time to complete with column space (',num2str(n),'): ',num2str(toc(nTIC)),' secs']);

%Save Workspace
save(['SNR[',num2str(min(SNR)),'_',num2str(max(SNR)),']dB across i[',num2str(n),'] with iter = ',num2str(LagrangianHunting.iter)]);
toc(mainTIC);

% end
