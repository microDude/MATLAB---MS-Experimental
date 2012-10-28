
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                            DESCRIPTION                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code calculates:
%
%  1) CRLB and MSE's for different estimators of matrix X, namely:
%        - MLNOSP (Maximum Likelihood NO sparsity enforsed),
%        - MLSP (Maximum Likelihood with sparsity enforsed),
%        - MLNOSPonSupport (Maximum Likelihood NO sparsity enforsed, but
%          updates only on support),
%        - LSNOSP (Least Squares NO sparsity enforsed),
%        - LSSP (Least Squares with enforsed sparsity),
%        - LSNOSPonSupport (Least Squares NO sparsity enforsed, updates
%          only on support),
%         NB!!! support can be either found by other algorithms or known
%         true sparsity pattren of X can be taken as support
%         (we call this case oracle).
%
%  2) probability of correct row-sparsity recovery
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc
format long
clk = clock;
f1=fullfile(strcat('exp'));
if exist(f1,'dir') == 0
   mkdir (f1);
end;
f2=fullfile(strcat('CRLB_',num2str(date),'_',num2str(clk(4:5))));
if exist(f2,'dir') == 0
   mkdir (f2);
end;
tic
diary on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Initialize parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Define this general parameters
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% path for saving all plots
saveToPath = 'C:\Users\Evgenia\Documents\Work\Research\Project Poisson\Codes\Poisson Case ML and LS frameworks\September 25\exp\';

disp('General parameters of the problem:');
disp('**********************************');
% N<<M<<K
% A's are MxK, X is KxN, Y is MxN
% number of measurements
N = 9;
disp(strcat('number of measurements N = ',num2str(N)));
% dimension of vector y
M = 11;
disp(strcat('dimension of vector y M = ',num2str(M)));
% dimension of vector x - column of X
K = 12;
disp(strcat('dimension of vector x K = ',num2str(K)));
% row-sparsity of matrix X
rowSparsity = 2;
disp(strcat('row-sparsity of X = ',num2str(rowSparsity)));
% initialize Xsolution for all algorithms
XsolInitial = rand(K,N);

% epsilonThresholdFixed - epsilon with fixed value for thresholding
% solution matrix X
epsilonThresholdFixed = 1e-5;

% epsilonThresholdPercentage - this epsilon is defined later
% (in fGenerateData.m) as 0.1% of the smallest value of the non-zero 2 norm
% of X. (this value is defined after X is defined)
% epsilonThresholdPercentage =

% epsilonThreshold is defined after data is generated. Can have 3 values:
% epsilonThresholdFixed, epsilonThresholdPercentage, 0 (zero, i.e. no
% thresholding only projection on the domain R^{N*K}_{+})
% epsilonThreshold =

% number of runs for every multiple of matrix X
runNo = 1;
disp(strcat('number of runs for every multiple of matrix X = ',num2str(runNo)));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% initialize parameters for MLNOSP method
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('  ');
disp('Parameters of the MLNOSP method:');
disp('**********************************');
NiterMLNOSP = 15;
disp(strcat('NiterMLNOSP = ',num2str(NiterMLNOSP)));


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% initialize parameters for MLSP method
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('  ');
disp('Parameters of the MLSP method:');
disp('**********************************');
% probability of being in the confidence region
confidProbabMLSP = 0.99;
% parameter for one-sided Chebyshev inequality
chebParam = sqrt(1/(1-confidProbabMLSP)-1);
% bound on expectation of D_ij
Ce = 0.6;
% bound on variance of D_ij
Cv = 0.55;
% radius of the confidence constrain set
epsilonMLSP = M*N*Ce+chebParam*sqrt(M*N*Cv);
disp(strcat('epsilonMLSP = ',num2str(epsilonMLSP)));
% number of iterations for finding optimal gamma
NiterGamma = 2;
disp(strcat('NiterGamma = ',num2str(NiterGamma)));
% min values of gamma range
gammaMin = 1e-1;
% max value of gamma range
gammaMax = 1e1;
% number of iteration for update of solution matrix
NiterMLSP = 10;
disp(strcat('NiterMLSP = ',num2str(NiterMLSP)));


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% initialize parameters for MLNOSPonSupport (MLNOSP on support) method,
% support (non-zero rows) is either found by MLSP method or true support of
% X (oracle)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('  ');
disp('Parameters of the MLNOSPonSUpport method:');
disp('**********************************');
NiterMLNOSPsup = 15;
disp(strcat('NiterMLNOSPsup = ',num2str(NiterMLNOSPsup)));


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% initialize parameters for LSNOSP method
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('  ');
disp('Parameters of the LSNOSP method:');
disp('**********************************');
NiterLSNOSP = 15;
disp(strcat('NiterLSNOSP = ',num2str(NiterLSNOSP)));


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% initialize parameters for LSSP method
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('  ');
disp('Parameters of the LSSP method:');
disp('**********************************');
% number of iterations to find optimal delta
NiterDelta = 3;
disp(strcat('NiterDelta = ',num2str(NiterDelta)));
% min value for delta range
deltaMinG = 1e-1;
% max value for delta range
deltaMaxG = 1e1;
% number of iterations for LSSP algorithm
NiterLSSP = 15;
disp(strcat('NiterLSSP = ',num2str(NiterLSSP)));
% probability for calculation of epsilonLSSP (epsilon is calculated later)
pLSSP = 0.1;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% initialize parameters for LSNOSPonSupport (LSNOSP on support) method,
% support (non-zero rows) is either found by LSSP method or true support of
% X (oracle)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('  ');
disp('Parameters of the LSNOSPonSupport method:');
disp('**********************************');
NiterLSNOSPsup = 15;
disp(strcat('NiterLSNOSPsup = ',num2str(NiterLSNOSPsup)));


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                FLAGS
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% thresholding flag, if 0 then no thresholding, if 1 then thresholding with
% epsilonThresholdFixed, if 2 then thresholding solution matrix with
% epsilonThresholdPercentage
thresholdFlag = 2;

% if 1 then displays X and all solution matrices
displayMatr = 0;
% if 1 then displays all auxiliary plots (e.g. objective functions, etc.)
displayAuxiliary = 1;
% if 1 then displays 12norms of X and solution matrix for each method
display12norm = 1;
%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Main part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                          Generate data
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[A,Arowsum,X,X12norm,sppatX,Y,ripconst,epsilonThresholdPercentage,epsilonLSSP] = fGenerateData(K,N,M,rowSparsity);

% MLSP: lower bound for min value of 2 row norm that guarantees correct sparsity
% recovery
lbMLSP = 2*sqrt(epsilonMLSP)/(1-ripconst);

% vectorize true X
xtrue = reshape(X,N*K,1);

if displayMatr==1
    disp('   ')
    disp('initial row-sparse matrix X');
    disp('***************************');
    disp(X);
end;

% epsilonThreshold is defined
if thresholdFlag == 1
    epsilonThreshold = epsilonThresholdFixed;
else
    if thresholdFlag == 2
        epsilonThreshold = epsilonThresholdPercentage;
    else
        epsilonThreshold = 0;
    end;
end;
disp('   ')
disp(strcat('epsilonThreshold = ',num2str(epsilonThreshold)));
disp('***************************');
pause
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% CRLB and MSEs for different methods, Probability of correct sparsity
% recovery
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% multiples of matrix X
multParam = [1e-3,1,1e3];

disp('   ');
disp('multiples of matrix X')
disp('************************');
disp(num2str(multParam'));

% number of multiples of X
ParamLength = length(multParam);

% Probability of correct sparsity recovery for different methods
ProbabilityMLSP = zeros(ParamLength,1);
ProbabilityLSSP = zeros(ParamLength,1);

% total MSE for different methods
MSEtotalErrorMLNOSP = zeros(ParamLength,1);
MSEtotalErrorMLSP = zeros(ParamLength,1);
MSEtotalErrorLSNOSP = zeros(ParamLength,1);
MSEtotalErrorLSSP = zeros(ParamLength,1);
% MSE of MLNOSP on support found by MLSP method
MSEtotalErrorMLNOSPsupMLSP = zeros(ParamLength,1);
% MSE of MLNOSP on support of true matrix X
MSEtotalErrorMLNOSPsupTrue = zeros(ParamLength,1);
% MSE of LSNOSP on support found by LSSP method
MSEtotalErrorLSNOSPsupLSSP = zeros(ParamLength,1);
% MSE of LSNOSP on support of true matrix X
MSEtotalErrorLSNOSPsupTrue = zeros(ParamLength,1);

% initialize CRLBX - array that contains CRLBX for each mult parameter
CRLBX = zeros(ParamLength,1);

for imultParam = 1:ParamLength
    disp('   ');
    disp('                        LOOP OVER MULTIPLES OF X')
    disp('*****************************************************************************');
    disp('   ');
    disp(strcat('current multiple parameter is_',num2str(imultParam),'_out of_',num2str(ParamLength)));
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
    
    % correct sparsity pattern counter for each method
    sppatCntrMLSP = zeros(1,runNo);
    sppatCntrLSSP = zeros(1,runNo);
    
    % initialize errors for each method
    errorMLNOSP = zeros(N*K,runNo);
    errorMLSP = zeros(N*K,runNo);
    errorLSNOSP = zeros(N*K,runNo);
    errorLSSP = zeros(N*K,runNo);
    errorMLNOSPsupportMLSP = zeros(N*K,runNo);
    errorMLNOSPsupportTrue = zeros(N*K,runNo);
    errorLSNOSPsupportLSSP = zeros(N*K,runNo);
    errorLSNOSPsupportTrue = zeros(N*K,runNo);
    
    % define new X and its 12norm
    newX = X*multParam(imultParam);
    newX12norm = sqrt(sum(newX.*newX,2));
      
    % initialize Lambda - matrix that has A{iN}*X(:,iN) values as columns
    Lambda = zeros(M,N);
    % define Lambda
    for iN = 1:N
        Lambda(:,iN) = A{iN}*newX(:,iN);
    end;
    
    % define CRLBX for current parameter
    for iN=1:N
        CRLBX(imultParam) =  CRLBX(imultParam) + trace(inv( A{iN}(:,sppatX==1)'*diag(1./Lambda(:,iN))*A{iN}(:,sppatX==1)));
    end;
    
    % for current mult parameter run algorithms runNo times
    for irun = 1:runNo
        
        disp('  ');
        disp(strcat('current run is_',num2str(irun),'_out of_',num2str(runNo)));
        disp('..................................');
        
        % generate Y(:,*) = Poisson(A{i}*X(:,i)) for each new run
        Y = zeros(M,N);
        for iN = 1:N
            Y(:,iN) = poissrnd(Lambda(:,iN));
        end;
        disp('   ');
        disp('matrix of Poisson counts Y for current multiple parameter');
        disp('*******************************');
        disp(Y);
        
        fsumY = sum(sum(Y));
        % epsilon calculated without simplification (look in the notes method 2 epsilon 2)
        hlpEpsLSSP = sqrt((2/pLSSP)-1);
        epsilonLSSP = fsumY+(hlpEpsLSSP^2)/2+sqrt(hlpEpsLSSP^2*fsumY+hlpEpsLSSP^4/4)+hlpEpsLSSP*sqrt(2*fsumY^2+fsumY*(4*hlpEpsLSSP^2+1)+hlpEpsLSSP^4+hlpEpsLSSP^2/2+(4*fsumY+2*hlpEpsLSSP^2+1)*sqrt(hlpEpsLSSP^2*fsumY+hlpEpsLSSP^4/4));
        disp('   ');
        disp(strcat('current epsilonLSSP = ', num2str(epsilonLSSP)));
               
        % check the inequality for min 2 norm of rows of X
        if min(newX12norm(sppatX==1))>2*sqrt(epsilonLSSP)/(1-ripconst)
            disp('LSSP: Inequality for min 2 norm is satisfied');
            disp('****************************');
        else
            disp('LSSP: Inequality for min 2 norm is NOT satisfied');
            disp('****************************');
        end;
        
        % initialize Xsolution for all algorithms
        XsolInitial = rand(K,N);
        
        
        %*********************************
        %             MLNOSP
        %*********************************
        disp('   ');
        disp('----------------');
        disp('MLNOSP');
        disp('----------------');
        
        [XsolMLNOSP,XsolMLNOSP12norm,objFuncMLNOSP] = fMLNOSP(A,Arowsum,Y,XsolInitial,NiterMLNOSP,epsilonThreshold);
        
        if displayMatr == 1
            disp('   ');
            disp('solution of MLNOSP matrix XsolMLNOSP');
            disp('************************************');
            disp(XsolMLNOSP);
        end;
        if displayAuxiliary == 1
            % display objective function for MLNOSP
            fig = figure('Name','Objective function MLNOSP');
            plot(1:NiterMLNOSP,objFuncMLNOSP);
            grid on
            xlabel('Iteration')
            ylabel('Objective function MLNOSP')
            title('Objective function MLNOSP')
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'Objective function MLNOSP') '.fig']);
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'Objective function MLNOSP') '.jpg']);
            close(fig)
        end;
        if display12norm == 1
            % display numerically and compare 12norm of true X and a solution to MLNOSP
            disp('   ');
            disp('12norm of trueX*multParam and XsolMLNOSP ');
            disp('********************************');
            disp([newX12norm,XsolMLNOSP12norm]);
        end;
        % vectorize current solution
        xMLNOSP = reshape(XsolMLNOSP,K*N,1);
        % add errors over runs
        errorMLNOSP(:,irun) = (xMLNOSP - xtrue*multParam(imultParam)).^2;
        
        
        %*********************************
        %              MLSP
        %*********************************
        disp('   ');
        disp('----------------');
        disp('MLSP');
        disp('----------------');
        
        [XsolMLSP,XsolMLSP12norm,objFuncMLSP,constrainMLSP,gamma] = fMLSP(A,Arowsum,Y,XsolInitial,epsilonMLSP,NiterMLSP,NiterGamma,gammaMin,gammaMax,epsilonThreshold);
        
        if displayMatr==1
            disp('   ');
            disp('solution of MLSP matrix XsolMLSP');
            disp('********************************');
            disp(XsolMLSP);
        end;
        if displayAuxiliary == 1
            % display objective function for MLSP for first iterGamma
            fig = figure('Name','Objective function MLSP');
            subplot(2,1,1)
            plot(1:NiterMLSP,objFuncMLSP(1,:));
            grid on
            xlabel('Iteration')
            ylabel('Objective function MLSP')
            title('objective function MLSP for on the first gamma iteration')
            % display objective function for MLSP for last iterGamma
            subplot(2,1,2)
            plot(1:NiterMLSP,objFuncMLSP(end,:));
            grid on
            xlabel('Iteration')
            ylabel('Objective function MLSP')
            title('objective function MLSP for on the last gamma iteration')
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'Objective function MLSP') '.fig']);
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'Objective function MLSP') '.jpg']);
            close(fig)
                       
            disp('   ');
            disp('Lagrangian parameter gamma (for MLSP)');
            disp('*************************************');
            disp(num2str(gamma(end)));
            
            % display constrainMLSP and true epsilonMLSP
            fig = figure('Name','constrainMLSP and true epsilonMLSP');
            plot(1:NiterGamma,constrainMLSP,'-b*',1:NiterGamma,epsilonMLSP*ones(1,NiterGamma),'r');
            hlpLegendconstrMLSP = legend('constrainMLSP','epsilonMLSP');
            set(hlpLegendconstrMLSP,'Location','SouthOutside')
            grid on
            xlabel('Iteration for gamma')
            ylabel('constrains and epsilon')
            title('constrainMLSP and true epsilonMLSP');
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'constrainMLSP and true epsilonMLSP') '.fig']);
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'constrainMLSP and true epsilonMLSP') '.jpg']);
            close(fig);
            
            % display values of gamma vs iterGamma
            fig = figure('Name','values of gamma vs iterGamma');
            plot(1:NiterGamma,gamma);
            xlabel('iterGamma')
            ylabel('values of gamma ')
            grid on
            title('values of gamma vs iterGamma');
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'values of gamma vs iterGamma') '.fig']);
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'values of gamma vs iterGamma') '.jpg']);
            close(fig);
        end;
        if display12norm == 1
            % display numerically and compare 12norm of true X and a solution to MLSP
            disp('   ');
            disp('12norm of trueX*multParam and XsolMLSP');
            disp('*****************************');
            disp([newX12norm,XsolMLSP12norm(end,:)']);
        end;
        % vectorize current solution
        xMLSP = reshape(XsolMLSP,K*N,1);
        % add errors over runs
        errorMLSP(:,irun) = (xMLSP - xtrue*multParam(imultParam)).^2;
        
        % find sparsity pattern of the MLSP solution
        sparsPatternMLSP = (XsolMLSP12norm(end,:)>0); 
        % check if the sparsity pattern was found correctly
        if sum(abs(sparsPatternMLSP-sppatX'))==0
            % then sparsity pattern was found correctly 
            sppatCntrMLSP(irun) = 1;
        else
            % sparsity pattren was found incorrectly
            sppatCntrMLSP(irun) = 0;
        end;
                
        %*********************************
        % MLNOSP on support found by MLSP
        %*********************************
        disp('   ');
        disp('----------------');
        disp('MLNOSP on support found by MLSP');
        disp('----------------');
                
        [XsolMLNOSPsupFull,XsolMLNOSPfull12norm,objFuncMLNOSPsup] = fMLNOSPonSupport(A,Arowsum,Y,NiterMLNOSPsup,epsilonThreshold,sparsPatternMLSP);
        
        if displayMatr==1
            disp('   ');
            disp('solution of MLNOSP on MLSP support ');
            disp('******************************');
            disp(XsolMLNOSPsupFull);
        end;
        if displayAuxiliary == 1
            % display objective function for MLNOSPonSupport
            fig = figure('Name','Objective function MLNOSPsup');
            plot(1:NiterMLNOSPsup,objFuncMLNOSPsup);
            grid on
            xlabel('Iteration')
            ylabel('Objective function MLNOSPsupMLSP')
            title('Objective function MLNOSPsupMLSP')
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'Objective function MLNOSPsupMLSP') '.fig']);
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'Objective function MLNOSPsupMLSP') '.jpg']);
            close(fig)
        end;
        if display12norm == 1
            % display numerically and compare 12norm of true X and a solution to
            % MLNOSP on support
            disp('   ');
            disp('12norm of trueX*multParam and XsolMLNOSPsupportMLSP');
            disp('********************************');
            disp([newX12norm,XsolMLNOSPfull12norm]);
        end;
        % vectorize current solution
        xMLNOSPsupportMLSP = reshape(XsolMLNOSPsupFull,K*N,1);
        % add errors over runs
        errorMLNOSPsupportMLSP(:,irun) = (xMLNOSPsupportMLSP - xtrue*multParam(imultParam)).^2;
       
        
        %*********************************
        % MLNOSP on true support ORACLE
        %*********************************
        disp('   ');
        disp('----------------');
        disp('MLNOSP on true support ORACLE');
        disp('----------------');
               
        [XsolMLNOSPfullO,XsolMLNOSPfull12normO,objFuncMLNOSPsupO] = fMLNOSPonSupport(A,Arowsum,Y,NiterMLNOSPsup,epsilonThreshold,sppatX);
        
        if displayMatr==1
            disp('   ');
            disp('solution of MLNOSP on True support ');
            disp('******************************');
            disp(XsolMLNOSPfullO);
        end;
        if displayAuxiliary == 1
            % display objective function for MLNOSPonSupport
            fig = figure('Name','Objective function MLNOSPsup True');
            plot(1:NiterMLNOSPsup,objFuncMLNOSPsupO);
            grid on
            xlabel('Iteration')
            ylabel('Objective function MLNOSPsupTrue')
            title('Objective function MLNOSPsupTrue')
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'Objective function MLNOSPsupTrue') '.fig']);
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'Objective function MLNOSPsupTrue') '.jpg']);
            close(fig)
        end;
        if display12norm == 1
            % display numerically and compare 12norm of true X and a solution to
            % MLNOSP on support
            disp('   ');
            disp('12norm of trueX*multParam and XsolMLNOSPsupportTrue');
            disp('********************************');
            disp([newX12norm,XsolMLNOSPfull12normO]);
        end;
        % vectorize current solution
        xMLNOSPsupportTrue = reshape(XsolMLNOSPfullO,K*N,1);
        % add errors over runs
        errorMLNOSPsupportTrue(:,irun) = (xMLNOSPsupportTrue - xtrue*multParam(imultParam)).^2;
       
        
        %*********************************
        %              LSNOSP
        %*********************************
        disp('   ');
        disp('----------------');
        disp('LSNOSP');
        disp('----------------');
        
        [XsolLSNOSP,XsolLSNOSP12norm,objFuncLSNOSP] = fLSNOSP(A,Y,XsolInitial,NiterLSNOSP,epsilonThreshold);
        
        if displayMatr==1
            disp('   ');
            disp('solution of LSNOSP matrix XsolLSNOSP');
            disp('********************************');
            disp(XsolLSNOSP);
        end;
        if displayAuxiliary == 1
            % display objective function for LSNOSP
            fig = figure('Name','Objective function LSNOSP');
            plot(1:NiterLSNOSP,objFuncLSNOSP);
            grid on
            xlabel('Iteration')
            ylabel('Objective function LSNOSP')
            title('Objective function LSNOSP')
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'Objective function LSNOSP') '.fig']);
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'Objective function LSNOSP') '.jpg']);
            close(fig);
        end;
        if display12norm == 1
            % display numerically and compare 12norm of true X and a solution to LSNOSP
            disp('   ');
            disp('12norm of trueX*multParam and XsolLSNOSP');
            disp('*******************************');
            disp([newX12norm,XsolLSNOSP12norm]);
        end;
        % vectorize current solution
        xLSNOSP = reshape(XsolLSNOSP,K*N,1);
        % add errors over runs
        errorLSNOSP(:,irun) = (xLSNOSP - xtrue*multParam(imultParam)).^2;
        
        
        %*********************************
        %               LSSP
        %*********************************
        disp('   ');
        disp('----------------');
        disp('LSSP');
        disp('----------------');
        
        [XsolLSSP,XsolLSSP12norm,objFuncLSSP,constrainLSSP,delta] = fLSSP3(A,Y,XsolInitial,epsilonLSSP,NiterLSSP,NiterDelta,deltaMinG,deltaMaxG,epsilonThreshold);
        
        if displayMatr==1
            disp('   ');
            disp('solution of LSSP matrix XsolLSSP');
            disp('********************************');
            disp(XsolLSSP);
        end;
        if displayAuxiliary == 1
            % display objective function for LSSP for first iterDelta
            fig = figure('Name','Objective function LSSP');
            subplot(2,1,1)
            plot(1:NiterLSSP,objFuncLSSP(1,:));
            grid on
            xlabel('Iteration')
            ylabel('Objective function LSSP')
            title('objective function LSSP for on the first delta iteration')
            % display objective function for LSSP for last iterDelta
            subplot(2,1,2)
            plot(1:NiterLSSP,objFuncLSSP(end,:));
            grid on
            xlabel('Iteration')
            ylabel('Objective function LSSP')
            title('objective function LSSP for on the last delta iteration')
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'Objective function LSSP') '.fig']);
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'Objective function LSSP') '.jpg']);
            close(fig);
                       
            % display constrainLSSP and true epsilonLSSP
            fig = figure('Name','constrainLSSP and true epsilonLSSP');
            plot(1:NiterDelta,constrainLSSP,'-b*',1:NiterDelta,epsilonLSSP*ones(1,NiterDelta),'r');
            hlpLegendconstrLSSP = legend('constrainLSSP','epsilonLSSP');
            set(hlpLegendconstrLSSP,'Location','SouthOutside')
            grid on
            xlabel('Iteration for delta')
            ylabel('constrains and epsilon')
            title('constrainLSSP and true epsilonLSSP');
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'constrainLSSP and true epsilonLSSP') '.fig']);
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'constrainLSSP and true epsilonLSSP') '.jpg']);
            close(fig);
            
            % display values of delta vs iterDelta
            fig = figure('Name','values of delta vs iterDelta');
            plot(1:NiterDelta,delta);
            xlabel('iterDelta')
            ylabel('values of delta ')
            grid on
            title('values of delta vs iterDelta');
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'values of delta vs iterDelta') '.fig']);
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'values of delta vs iterDelta') '.jpg']);
            close(fig);
        end;
        if display12norm == 1
            % display numerically and compare 12norm of true X and a solution to LSSP
            disp('   ');
            disp('12norm of trueX*multParam and XsolLSSP');
            disp('*****************************');
            disp([newX12norm,XsolLSSP12norm(end,:)']);
        end;
        % vectorize current solution
        xLSSP = reshape(XsolLSSP,K*N,1);
        % add errors over runs
        errorLSSP(:,irun) = (xLSSP - xtrue*multParam(imultParam)).^2;
        
        % find sparsity pattern of the LSSP solution
        sparsPatternLSSP = (XsolLSSP12norm(end,:)>0); 
        % check if the sparsity pattern was found correctly
        if sum(abs(sparsPatternLSSP-sppatX'))==0
            % then sparsity pattern was found correctly 
            sppatCntrLSSP(irun) = 1;
        else
            % sparsity pattren was found incorrectly
            sppatCntrLSSP(irun) = 0;
        end;
        
                
        %*********************************
        % LSNOSP on support found by LSSP
        %*********************************
        disp('   ');
        disp('----------------');
        disp('LSNOSP on support found by LSSP');
        disp('----------------');
               
        [XsolLSNOSPsupFull,XsolLSNOSPsupFull12norm,objFuncLSNOSPsup] = fLSNOSPonSupport(A,Y,NiterLSNOSPsup,epsilonThreshold,sparsPatternLSSP);
        
        if displayMatr==1
            disp('   ');
            disp('solution of LSNOSP on LSSP Support');
            disp('********************************');
            disp(XsolLSNOSPsupFull);
        end;
        if displayAuxiliary == 1
            % display objective function for LSNOSPonSupport
            fig = figure('Name','Objective function LSNOSPsupLSSP');
            plot(1:NiterLSNOSPsup,objFuncLSNOSPsup);
            grid on
            xlabel('Iteration')
            ylabel('Objective function LSNOSPsupLSSP')
            title('Objective function LSNOSPsupLSSP')
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'Objective function LSNOSPsupLSSP') '.fig']);
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'Objective function LSNOSPsupLSSP') '.jpg']);
            close(fig);
        end;
        if display12norm == 1
            % display numerically and compare 12norm of true X and a solution
            % to LSNOSPonSupport
            disp('   ');
            disp('12norm of trueX*multParam and XsolLSNOSPonLSSPsupport');
            disp('*****************************');
            disp([newX12norm,XsolLSNOSPsupFull12norm]);
        end;
        % vectorize current solution
        xLSNOSPsupportLSSP = reshape(XsolLSNOSPsupFull,K*N,1);
        % add errors over runs
        errorLSNOSPsupportLSSP(:,irun) = (xLSNOSPsupportLSSP - xtrue*multParam(imultParam)).^2;
        
        
        %*********************************
        % LSNOSP on true support ORACLE
        %*********************************
        disp('   ');
        disp('----------------');
        disp('LSNOSP on true support ORACLE');
        disp('----------------');
        
        [XsolLSNOSPsupFullO,XsolLSNOSPsupFull12normO,objFuncLSNOSPsupO] = fLSNOSPonSupport(A,Y,NiterLSNOSPsup,epsilonThreshold,sppatX);
        
        if displayMatr==1
            disp('   ');
            disp('solution of LSNOSP on True Support');
            disp('********************************');
            disp(XsolLSNOSPsupFullO);
        end;
        if displayAuxiliary == 1
            % display objective function for LSNOSPonSupport
            fig = figure('Name','Objective function LSNOSPsupTrue');
            plot(1:NiterLSNOSPsup,objFuncLSNOSPsupO);
            grid on
            xlabel('Iteration')
            ylabel('Objective function LSNOSPsupTrue')
            title('Objective function LSNOSPsupTrue')
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'Objective function LSNOSPsupTrue') '.fig']);
            saveas(fig,[saveToPath strcat(num2str(imultParam),'+',num2str(irun),'Objective function LSNOSPsupTrue') '.jpg']);
            close(fig);
        end;
        if display12norm == 1
            % display numerically and compare 12norm of true X and a solution
            % to LSNOSPonSupport
            disp('   ');
            disp('12norm of trueX*multParam and XsolLSNOSPonTruesupport');
            disp('*****************************');
            disp([newX12norm,XsolLSNOSPsupFull12normO]);
        end;
        xLSNOSPsupportTrue = reshape(XsolLSNOSPsupFullO,K*N,1);
        % add errors over runs
        errorLSNOSPsupportTrue(:,irun) = (xLSNOSPsupportTrue - xtrue*multParam(imultParam)).^2;
        
                
        % display all 12 norms of all 8 solutions
        %****************************************
        disp('  ');
        disp('12norms of the solutions: true*multParam, MLNOSP, MLSP, MLNOSPsupMLSP, MLNOSPsupTRUE, LSNOSP, LSSP, LSNOSPsupLSSP, LSNOSPsupTRUE: ');
        disp('----------------------------------------------------------------------------------------------------');
        disp([newX12norm, XsolMLNOSP12norm, XsolMLSP12norm(end,:)', XsolMLNOSPfull12norm, XsolMLNOSPfull12normO, XsolLSNOSP12norm, XsolLSSP12norm(end,:)',XsolLSNOSPsupFull12norm,XsolLSNOSPsupFull12normO]);
        disp('  ');
        
    end; % irun
    
    % Probability of correct sparsity recovery for MLSP 
    ProbabilityMLSP(imultParam) = sum(sppatCntrMLSP)/runNo;
    % Probability of correct sparsity recovery for LSSP 
    ProbabilityLSSP(imultParam) = sum(sppatCntrLSSP)/runNo;
    
    % average MSE error for MLNOSP over number of runs
    MSEtotalErrorMLNOSP(imultParam) = sum(sum(errorMLNOSP,2)/runNo);
    % average MSE error for MLSP over number of runs
    MSEtotalErrorMLSP(imultParam) = sum(sum(errorMLSP,2)/runNo);
    % average MSE error for MLNOSP on support found by MLSP over number of runs
    MSEtotalErrorMLNOSPsupMLSP(imultParam) = sum(sum(errorMLNOSPsupportMLSP,2)/runNo);
    % average MSE error for MLNOSP on true support over number of runs
    MSEtotalErrorMLNOSPsupTrue(imultParam) = sum(sum(errorMLNOSPsupportTrue,2)/runNo);
    % average MSE error for LSNOSP over number of runs
    MSEtotalErrorLSNOSP(imultParam) = sum(sum(errorLSNOSP,2)/runNo);
    % average MSE error for LSSP over number of runs
    MSEtotalErrorLSSP(imultParam) = sum(sum(errorLSSP,2)/runNo);
    % average MSE error for LSNOSP on support found by LSSP over number of runs
    MSEtotalErrorLSNOSPsupLSSP(imultParam) = sum(sum(errorLSNOSPsupportLSSP,2)/runNo);
    % average MSE error for LSNOSP on true support over number of runs
    MSEtotalErrorLSNOSPsupTrue(imultParam) = sum(sum(errorLSNOSPsupportTrue,2)/runNo);
end; % multParam

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Display results
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% display all 8 MSE and CRLB 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig = figure('Name','MSE for 8 methods');
loglog(multParam,MSEtotalErrorMLNOSP'./(multParam.^2),'-y*','LineWidth',1)
hold on
grid on
loglog(multParam,MSEtotalErrorMLSP'./(multParam.^2),'-ms','LineWidth',1);
loglog(multParam,MSEtotalErrorMLNOSPsupMLSP'./(multParam.^2),'--go','LineWidth',1);
loglog(multParam,MSEtotalErrorMLNOSPsupTrue'./(multParam.^2),'-.mh','LineWidth',1);
loglog(multParam, MSEtotalErrorLSNOSP'./(multParam.^2),'-ko','LineWidth',1);
loglog(multParam,MSEtotalErrorLSSP'./(multParam.^2),'-bp','LineWidth',1);
loglog(multParam, MSEtotalErrorLSNOSPsupLSSP'./(multParam.^2),'--co','LineWidth',1);
loglog(multParam,MSEtotalErrorLSNOSPsupTrue'./(multParam.^2),'-.bh','LineWidth',1);
loglog(multParam,CRLBX'./(multParam.^2),':rd','LineWidth',1);
hold off
hlpLegendconstr = legend('MSE MLNOSP','MSE MLSP','MSE MLNOSPsupMLSP','MSE MLNOSPsupTRUE','MSE LSNOSP','MSE LSSP','MSE LSNOSPsupLSSP','MSE LSNOSPsupTRUE','CRLBx');
set(hlpLegendconstr,'Location','SouthOutside')
xlabel('number of parameters')
ylabel('MSEs ')
title('MSEs for X');
saveas(fig,[saveToPath 'MSEs and CRLB for X' '.fig']);
saveas(fig,[saveToPath 'MSEs and CRLB for X' '.jpg']);
saveas(fig,[saveToPath 'MSEs and CRLB for X' '.eps']);
close(fig);

% probability of correct sparsity recovery for LSSP and MLSP methods
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% together
fig = figure('Name','Probability of correct sparsity recovery for MLSP and LSSP');
plot(1:ParamLength,ProbabilityMLSP,'-c*','LineWidth',2);
hold on
grid on
plot(1:ParamLength,ProbabilityLSSP,'-b*','LineWidth',2);
%plot(lbMLSP,0:0.001:1,'g','LineWidth',4);
hold off
hlpLegendconstr = legend('probability for MLSP','probability for LSSP');
set(hlpLegendconstr,'Location','SouthOutside')
xlabel('multiplication parameters')
ylabel('probability ')
title('Probability of correct sparsity recovery for MLSP and LSSP');
saveas(fig,[saveToPath 'Probability of correct sparsity recovery for MLSP and LSSP' '.fig']);
saveas(fig,[saveToPath 'Probability of correct sparsity recovery for MLSP and LSSP' '.jpg']);
saveas(fig,[saveToPath 'Probability of correct sparsity recovery for MLSP and LSSP' '.eps']);
close(fig);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% probability of correct sparsity recovery for LSSP and MLSP methods
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% separately
% MLSP
fig = figure('Name','Probability of correct sparsity recovery for MLSP');
plot(1:ParamLength,ProbabilityMLSP,'-b*','LineWidth',2);
grid on
%hold on
%plot(lbMLSP,0:0.001:1,'g','LineWidth',4);
%hold off
xlabel('Multiplication parameters')
ylabel('Probability ')
title('Probability of correct sparsity recovery for MLSP');
saveas(fig,[saveToPath 'Probability of correct sparsity recovery for MLSP' '.fig']);
saveas(fig,[saveToPath 'Probability of correct sparsity recovery for MLSP' '.jpg']);
saveas(fig,[saveToPath 'Probability of correct sparsity recovery for MLSP' '.eps']);
close(fig);

% LSSP
fig = figure('Name','Probability of correct sparsity recovery for LSSP');
plot(1:ParamLength,ProbabilityLSSP,'-b*','LineWidth',2);
grid on
xlabel('multiplication parameters')
ylabel('probability ')
title('Probability of correct sparsity recovery for LSSP');
saveas(fig,[saveToPath 'Probability of correct sparsity recovery for LSSP' '.fig']);
saveas(fig,[saveToPath 'Probability of correct sparsity recovery for LSSP' '.jpg']);
saveas(fig,[saveToPath 'Probability of correct sparsity recovery for LSSP' '.eps']);
close(fig);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tm = toc;
disp(strcat('running time = ', num2str(tm/60)));
save all
diary off