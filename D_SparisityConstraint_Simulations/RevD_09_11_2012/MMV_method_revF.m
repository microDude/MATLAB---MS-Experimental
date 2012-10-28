% ***************************************
% GRID-SEARCH: MMV Method across parameter space
% 09/06/2012
% ***************************************
clear all;
close all;
clc;
mainTIC = tic;

% Define Parameter Constants
iter.statistical = 4;                       % Number of Statistical Expirements
iter.final = 200;                           % Number of regression iterations for the final solution
LagrangianHunting.numSteps = 30;            % Number of steps allowed to hone in on the true Lagrangain
LagrangianHunting.iter = 80;                % Number of regression iterations per step for lagrangian Hunting
LagrangianHunting.min = 1e-8;               % Starting Lagrangian min
LagrangianHunting.max = 1;                  % Starting Lagrangian max
m = 100;                                    % Length of a_i, x_i, y_i, i.e. Length of each spectra
ki = 3;                                     % The standard deviation that was found to yield 91% = P(T>E) < 1/(K^2+1)
Cardx = 2;                                  % The cardinality of the x vector
ampAlpha = 100;                             % Amplitude of convolution kernel 'a'
eps6 = 1e-6;                                % Error threshold for x_est
alpha = 5:0.1:24.9;                         % Sample space used to generate 'a' and likewise A.

% Grid Search Parameters
n = 10:10:100;                                 % Column space of X and y, i.e. number of captured spectra
delta = 1:1:10;                             % Card(x) pattern spacing
SNR = 5:5:30;                              % SNR (dB) of the signal
AlphaSampleVariety = {'Low','Half','Full'}; % Low Variety, Half Variety, Full Variety

% Define the solution space
XSolGrid  = cell(length(n),length(delta),length(SNR),length(AlphaSampleVariety));
XTrueGrid = cell(length(n),length(delta),length(SNR));

% **************** Begin GRID SEARCH ***********************
% iN            = iterator for Column space of X and y
% iDELTA        = iterator for Card(x) pattern spacing
% iSNR          = iterator for SNR (dB) of the signal
% iVARIETY      = iterator for variety in sample space of alpha
% iSTATISTICAL  = iterator for statistical expirements
for iN=1:length(n)
    display(['iN = ',num2str(iN),'/',num2str(length(n))]);
    nTIC = tic;
    % Choose which convolution kernel widths will be used to generate the variety of A.
    alphaINDEX = cell(1,length(AlphaSampleVariety));
    alphaINDEX{1} = 1:n(iN);
    alphaINDEX{2} = round(linspace(1,length(alpha)/2,n(iN)));
    alphaINDEX{3} = round(linspace(1,length(alpha),n(iN)));
    
    for iDELTA=1:length(delta)
        display([' iDELTA = ',num2str(iDELTA),'/',num2str(length(delta))]);
        
        for iSNR=1:length(SNR)
            display(['  iSNR = ',num2str(iSNR),'/',num2str(length(SNR))]);
            % Generate Sparse x vector (same non-zero indexes)
            minAMPx = exp(SNR(iSNR)/20); % Figure out the min for the x amplitude
            X = zeros(m,n(iN));
            %X([10,10+delta(iDELTA)],:) = (minAMPx + ((2*minAMPx)- minAMPx).*rand(Cardx,1))*ones(1,n(iN)); %Similar Uniform random values assigned to fixed x values
            X([10,10+delta(iDELTA)],:) = minAMPx + ((2*minAMPx)- minAMPx).*rand(Cardx,n(iN)); %Uniform random values assigned to fixed x values
            XTrueGrid{iN,iDELTA,iSNR} = mean(X,2);
    
            for iVARIETY=1:length(AlphaSampleVariety)
                display(['   iVARIETY = ',num2str(iVARIETY),'/',num2str(length(AlphaSampleVariety))]);
                % Generate the Convolution Kernel, Maxwell-Boltzman
                b = 1:m;                            % vector used to generate 'a' and then A.
                A = cell(1,n(iN));                  % preallocate for A
                aVec = alpha(alphaINDEX{iVARIETY}); % create a vector to loop through to grab alpha values
                for t=1:length(A)
                    a = b.^2.*exp(-b.^2/(2*aVec(t)^2))/aVec(t)^3; % Generate Maxwell-Boltzman Curve
                    a = (ampAlpha/max(a)).*a;                     % Normalize to defined amplitude
                    A{t} = tril(toeplitz(a));                     % Cell array which holds A_i
                end
                
                XsolSTAT = zeros(m,iter.statistical);             % Preallocate matrix for statistical iterations
                for iSTATISTICAL=1:iter.statistical
                    % Section:01 ----------------------------------------------------
                    % Populate the observation space
                    y = zeros(m,n(iN));
                    for t=1:n(iN)
                        lambda = A{t}*X(:,t);           % Generate lambda_i
                        y(:,t) = poissrnd(lambda);      % Generate y_i
                    end
                    
                    % Section:02 ----------------------------------------------------
                    % Calculate epsilon
                    Sy = sum(sum(y,1));
                    Gamma = Sy + ki^2/2 + sqrt((Sy*ki^2+(ki^4/4)));
                    epsilon = Gamma + ki*sqrt(2*Gamma^2+Gamma);
                    %epsilon = 1*epsilon;
                    %fprintf('Calculated Epsilson : %d \n',epsilon);
                    
                    % Section:03 ----------------------------------------------------
                    % Perform Lagrangian Hunting
                    
                    % Preallocations for lagrangian hunting
                    Xsol = zeros(m,n(iN));        %intialize the Xsolution
                    constraint = zeros(LagrangianHunting.numSteps,1);
                    lagrangian = zeros(LagrangianHunting.numSteps,1);

                    % Peform Lagrangian Hunting
                    for iterL = 1:LagrangianHunting.numSteps

                        % Define/Update lambda
                        lagrangian(iterL) = (LagrangianHunting.max + LagrangianHunting.min)/2;
                        
                        % Solve for X
                        Xsol = fSolveForX_revD(A,y,lagrangian(iterL),LagrangianHunting.iter,(1+(5-1).*rand(m,n(iN))),0);
                        
                        % Evaluate the Constraint, sum of squared error (least squares)
                        sqr_error = zeros(1,n(iN));
                        for t=1:n(iN)
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
                    
                    % Section:04 ----------------------------------------------------
                    % After Correct Lagrangian is found, Calculate the Xsol with many iterations
                    Xsol = fSolveForX_revD(A,y,lagrangian(iterL),iter.final,(1+(5-1).*rand(m,n(iN))),1);
                    
                    x_est = sum(Xsol.^2,2);                  % Form a single vector for the result
                    x_est(x_est < (max(x_est)*eps6)) = 0;    % Zero-out any element less then 6-order of magnitudes of the peak
                    
                    % Store the result
                    XsolSTAT(:,iSTATISTICAL) = x_est;
                end
                
                % ---------------------------------------------------------
                % Store the Final Solution for this Grid Point
                XSolGrid{iN,iDELTA,iSNR,iVARIETY} = mean(XsolSTAT,2);
            end
        end
    end
    display(['Time to complete with column space (',num2str(n(iN)),'): ',num2str(toc(nTIC)),' secs']);
end
 
%Save Workspace
save 'xyz';
toc(mainTIC);

% end
