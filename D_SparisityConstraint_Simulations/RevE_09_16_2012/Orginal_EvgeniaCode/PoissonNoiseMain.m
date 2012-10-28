
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                            DESCRIPTION                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code finds solution to main problem using different methods:
% MLNOSP (Maximum Likelihood NO sparsity enforsed),
% MLSP (Maximum Likelihood with sparsity enforsed),
% MLNOSPonSupport (Maximum Likelihood NO sparsity enforsed, but updates
% only on support),
% LSNOSP (Least Squares NO sparsity enforsed),
% LSSP (Least Squares with enforsed sparsity),
% LSNOSPonSupport (Least Squares NO sparsity enforsed, updates only on
% support),
% where support can be either found by other algorithms or known true
% sparsity pattren of X can be taken as support (we call this case oracle).
%
% Moreover this code allows to prodice row-sparsity plots vs epsilon for
% LSSP and MLSP methods.
% To choose which method to run go to the FLAG section of this code.
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
f2=fullfile(strcat('MAIN_',num2str(date),'_',num2str(clk(4:5))));
if exist(f2,'dir') == 0
   mkdir (f2);
end;
tic
diary on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      PART 1:  Initialize parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Define this general parameters
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% path for saving all plots
saveToPath = 'C:\Users\Evgenia\Documents\Work\Research\Project Poisson\Codes\Poisson Case ML and LS frameworks\September 25\exp\';

% N<<M<<K
% A's are MxK, X is KxN, Y is MxN
% number of measurements
N = 20;
% dimension of vector y
M = 23;
% dimension of vector x - column of X
K = 25;
% row-sparsity of matrix X
rowSparsity = 2;
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


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% initialize parameters for MLNOSP method
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NiterMLNOSP = 30;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% initialize parameters for MLSP method
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
% number of iterations for finding optimal gamma
NiterGamma = 2;
% min values of gamma range
gammaMin = 1e-1;
% max value of gamma range
gammaMax = 1e1;
% number of iteration for update of solution matrix
NiterMLSP = 15;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% initialize parameters for MLNOSPonSupport (MLNOSP on support) method,
% support (non-zero rows) is either found by MLSP method or true support of
% X (oracle)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NiterMLNOSPsup = 10;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% initialize parameters for LSNOSP method
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NiterLSNOSP = 30;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% initialize parameters for LSSP method
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% number of iterations to find optimal delta
NiterDelta = 25;
% min value for delta range
deltaMinG = 1e-1;
% max value for delta range
deltaMaxG = 1e1;
% number of iterations for LSSP algorithm
NiterLSSP = 100;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% initialize parameters for LSNOSPonSupport (LSNOSP on support) method,
% support (non-zero rows) is either found by LSSP method or true support of
% X (oracle)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NiterLSNOSPsup = 10;



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                              FLAGS
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% thresholding flag, if 0 then no thresholding, if 1 then thresholding with
% epsilonThresholdFixed, if 2 then thresholding solution matrix with
% epsilonThresholdPercentage
thresholdFlag = 2;

% if 1 then displays X and all solution matrices
displayMatr = 0;
% if 1 then displays all auxiliary plots (e.g. objective functions, etc.)
displayAuxiliary = 1;
% if 1 then displays 12norms of X and all solution matrices
display12norm = 1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% choose which algorithms to run
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% if 1 run MLNOSP
runMLNOSPflag = 0;
% if 1 run MLSP
runMLSPflag = 0;
% if 1 run MLNOSPonSupport
runMLNOSPonSupportFlag = 0;
% choose support: 1 - support found by MLSP method (then turn runMLSPflag
% on), 0 - true support of X (OracleMLNOSP)
SupportMLNOSPflag = 0;
if (SupportMLNOSPflag == 1) && (runMLNOSPonSupportFlag == 1)
    % turn the MLSP method on to get support
    runMLSPflag = 1;
end;
% run MLSP on the range of epsilons and displays row-sparsity vs epsilons
runRowSparsityVSepsilonMLSP = 0;

% if 1 run LSNOSP
runLSNOSPflag = 0;
% if 1 run LSSP
runLSSPflag = 1;
% if 1 run LSNOSPonSupport
runLSNOSPonSupportFlag = 0;
% choose support: 1 - support found by LSSP method (then turn runLSSPflag
% on), 0 - true support of X (OracleLSNOSP)
SupportLSNOSPflag = 0;
if (SupportLSNOSPflag == 1) && (runLSNOSPonSupportFlag == 1)
    % turn LSSP method on to get support
    runLSSPflag = 1;
end;
% run LSSP on the range of epsilons and displays row-sparsity vs epsilons
runRowSparsityVSepsilonLSSP = 0;

%%%%%%%%%%%%%%%%%%%% PART 1:  Initialize parameters %%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  PART 2: Main functional part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                          Generate data
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[A,Arowsum,X,X12norm,sppatX,Y,ripconst,epsilonThresholdPercentage,epsilonLSSP] = fGenerateData(K,N,M,rowSparsity);

% matrix Lambda = A*X, to check if true solution is inso=ide epsilon ball
% for LSSP and MLSP algorithms
Lambda = zeros(M,N);
for iN = 1:N
    Lambda(:,iN) = A{iN}*X(:,iN);
end;

if displayMatr==1
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
pause
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!! ML FRAMEWORKS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% MLNOSP framework: ML, NO SPARSITY enforced
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if runMLNOSPflag == 1
    disp('*****************************');
    disp('MLNOSP framework');
    disp('*****************************');
    
    % run fMLNOSP function
    [XsolMLNOSP,XsolMLNOSP12norm,objFuncMLNOSP] = fMLNOSP(A,Arowsum,Y,XsolInitial,NiterMLNOSP,epsilonThreshold);
    
    if displayMatr == 1
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
        saveas(fig,[saveToPath 'Objective function MLNOSP' '.fig']);
        saveas(fig,[saveToPath 'Objective function MLNOSP' '.jpg']);
        close(fig)
    end;
    
    if display12norm == 1
        % display numerically and compare 12norm of true X and a solution to MLNOSP
        disp('12norm of true X and XsolMLNOSP ');
        disp('********************************');
        disp([X12norm,XsolMLNOSP12norm]);
    end;
end; % MLNOSP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% MLSP framewotk SPARSITY enforced
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if runMLSPflag == 1
    disp('*****************************');
    disp('MLSP framework');
    disp('*****************************');
    
    % check that true solution is inside epsilonMLSP ball
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    constrTrueX = 0;
    for iN = 1:N
        % indices where y_i(j)~=0
        ihlp = find(Y(:,iN)~=0);
        Hlp1 = zeros(M,1);
        Hlp1(ihlp) = Y(ihlp,iN)./(Lambda(ihlp,iN));
        Hlp2 = sum(Y(ihlp,iN).*log(Hlp1(ihlp)));
        constrTrueX = constrTrueX + Hlp2;
    end;
    % value of the constrain for true X
    constrTrueX = constrTrueX + sum(sum(Lambda)) - sum(sum(Y));
    % check if Xtrue is inside epsilon ball
    if constrTrueX < epsilonMLSP
        disp('   ');
        disp('True X is inside epsilonMLSP ball');
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
    else
        disp('   ');
        disp('True X is NOT inside epsilonMLSP ball');
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
    end;
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    % run fMLSP function
    [XsolMLSP,XsolMLSP12norm,objFuncMLSP,constrainMLSP,gamma] = fMLSP(A,Arowsum,Y,XsolInitial,epsilonMLSP,NiterMLSP,NiterGamma,gammaMin,gammaMax,epsilonThreshold);
    
    % find support/row-sparsity pattern of the XsolMLSP
    sppatXsolMLSP = XsolMLSP12norm(end,:);
    indsppatXsolMLSP = XsolMLSP12norm(end,:)~=0;
    sppatXsolMLSP(indsppatXsolMLSP) = 1;
    
    if displayMatr==1
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
        saveas(fig,[saveToPath 'Objective function MLSP' '.fig']);
        saveas(fig,[saveToPath 'Objective function MLSP' '.jpg']);
        close(fig);
        
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
        saveas(fig,[saveToPath 'ConstrainMLSP and epsilon' '.fig']);
        saveas(fig,[saveToPath 'ConstrainMLSP and epsilon' '.jpg']);
        close(fig);
        
        % display values of gamma vs iterGamma
        fig = figure('Name','values of gamma vs iterGamma');
        plot(1:NiterGamma,gamma);
        xlabel('iterGamma')
        ylabel('values of gamma ')
        grid on
        title('values of gamma vs iterGamma');
        saveas(fig,[saveToPath 'Values of gamma vs iterGamma' '.fig']);
        saveas(fig,[saveToPath 'Values of gamma vs iterGamma' '.jpg']);
        close(fig);
    end;
    
    if display12norm == 1
        % display numerically and compare 12norm of true X and a solution to MLSP
        disp('12norm of true X and XsolMLSP');
        disp('*****************************');
        disp([X12norm,XsolMLSP12norm(end,:)']);
    end;
end; % MLSP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%        MLNOSP on SUPPORT
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if runMLNOSPonSupportFlag==1
    disp('*****************************');
    disp('MLNOSP on Support framework');
    disp('*****************************');
    disp('');
    
    % set the support that will be used for this code: either found by MLSP
    % or true support of X (i.e. OracleMLNOSP)
    if SupportMLNOSPflag == 1
        % then use support found by MLSP method
        SupportMLNOSPonSup = sppatXsolMLSP;
        disp('support found by MLSP');
    else
        % then use support of true X
        SupportMLNOSPonSup = sppatX;
        disp('Oracle MLNOSP');
    end;
    
    % run MLNOSPonSupport
    [XsolMLNOSPfull,XsolMLNOSPfull12norm,objFuncMLNOSPsup] = fMLNOSPonSupport(A,Arowsum,Y,NiterMLNOSPsup,epsilonThreshold,SupportMLNOSPonSup);
    
    if displayMatr==1
        disp('');
        disp('solution of MLNOSP on support ');
        disp('******************************');
        disp(XsolMLNOSPfull);
    end;
    
    if displayAuxiliary == 1
        % display objective function for MLNOSPonSupport
        fig =  figure('Name','Objective function MLNOSPsup');
        plot(1:NiterMLNOSPsup,objFuncMLNOSPsup);
        grid on
        xlabel('Iteration')
        ylabel('Objective function MLNOSPsup')
        title('Objective function MLNOSPsup')
        saveas(fig,[saveToPath 'Objective function MLNOSPsup' '.fig']);
        saveas(fig,[saveToPath 'Objective function MLNOSPsup' '.jpg']);
        close(fig);
    end;
    
    if display12norm == 1
        % display numerically and compare 12norm of true X and a solution to
        % MLNOSP on support
        disp('12norm of true X and XsolMLNOSPsupport ');
        disp('********************************');
        disp([X12norm,XsolMLNOSPfull12norm]);
    end;
    
end; % run MLNOSP on support
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Row-sparsity VS epsilon analysis for MLSP framework
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% if 1 then produce the plot depicting row sparsity of MLSP solutions for
% different values of epsilon
if runRowSparsityVSepsilonMLSP == 1
    disp('****************************************************');
    disp('Row-sparsity VS epsilon analysis for MLSP framework ');
    disp('****************************************************');
    
    % define range for epsilons
    % epsilonMLSPrange = [epsilonMLSP/64,epsilonMLSP/32,epsilonMLSP/16,epsilonMLSP/8,epsilonMLSP/4,epsilonMLSP/2,epsilonMLSP,2*epsilonMLSP,4*epsilonMLSP,8*epsilonMLSP,16*epsilonMLSP,32*epsilonMLSP,64*epsilonMLSP];
    epsilonMLSPrange = [epsilonMLSP/3,epsilonMLSP,epsilonMLSP*3];
    
    % number of differents epsilons
    epsilonMLSPrangeLength = length(epsilonMLSPrange);
    % array that contains row sparsity of the solution for current epsilon
    rowSparsityMLSPVsEpsilon = zeros(epsilonMLSPrangeLength,1);
    % initialize counter
    epsCounter = 0;
    % save 12 norms of the solutions for different epsilon
    XsolMLSP12normForEachEps = zeros(K,epsilonMLSPrangeLength);
    % save optimal gamma found for each epsilon
    gammaEpsilonMLSP = zeros(1,epsilonMLSPrangeLength);
    
    % for every epsilon in the range find solution and its 12 norm
    for iterEpsilonMLSP = epsilonMLSPrange
        % number of the epsilon in the range of epsilons
        epsCounter = epsCounter + 1;
        disp(strcat('current epsilon is epsilon #',num2str(epsCounter),' out of ',num2str(epsilonMLSPrangeLength)));
        % for a given epsilonMLSP solve the MLSP problem
        [XsolMLSP,XsolMLSP12norm,objFuncMLSP,constrainMLSP,gamma] = fMLSP(A,Arowsum,Y,XsolInitial,iterEpsilonMLSP,NiterMLSP,NiterGamma,gammaMin,gammaMax,epsilonThreshold);
                
        % save 12norm that we got for current epsilon
        XsolMLSP12normForEachEps(:,epsCounter) = XsolMLSP12norm(end,:)';
        
        if displayAuxiliary == 1
            %display objective function for MLSP for first iterGamma
            fig = figure('Name',strcat('Objective function MLSP for epsilon #',num2str(epsCounter)));
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
            saveas(fig,[saveToPath strcat('Objective function MLSP for epsilon #',num2str(epsCounter)) '.fig']);
            saveas(fig,[saveToPath strcat('Objective function MLSP for epsilon #',num2str(epsCounter)) '.jpg']);
            close(fig);
            
            % display constrainMLSP and current epsilonMLSP, i.e. iterEpsilonMLSP
            fig = figure('Name',strcat('constrainMLSP, true epsilonMLSP and current epsilon #',num2str(epsCounter)));
            plot(1:NiterGamma,constrainMLSP,'-b*',1:NiterGamma,iterEpsilonMLSP*ones(1,NiterGamma),'r',1:NiterGamma,epsilonMLSP*ones(1,NiterGamma),'m');
            hlpLegendconstrMLSP = legend('constrainMLSP','iterEpsilonMLSP','true epsilon MLSP');
            set(hlpLegendconstrMLSP,'Location','SouthOutside')
            grid on
            xlabel('Iteration for gamma')
            ylabel(strcat('constrainMLSP, true epsilonMLSP and current epsilon #',num2str(epsCounter)));
            title(strcat('constrainMLSP, true epsilonMLSP and current epsilon #',num2str(epsCounter)));
            saveas(fig,[saveToPath strcat('constrainMLSP, true epsilonMLSP and current epsilon #',num2str(epsCounter)) '.fig']);
            saveas(fig,[saveToPath strcat('constrainMLSP, true epsilonMLSP and current epsilon #',num2str(epsCounter)) '.jpg']);
            close(fig);
            
            % display values of gamma vs iterGamma
            fig = figure('Name',strcat('values of gamma vs iterGamma for current epsilon #',num2str(epsCounter)));
            plot(1:NiterGamma,gamma);
            xlabel('iterGamma')
            ylabel('values of gamma ')
            grid on
            title(strcat('values of gamma vs iterGamma for current epsilon #',num2str(epsCounter)));
            saveas(fig,[saveToPath strcat('values of gamma vs iterGamma for current epsilon #',num2str(epsCounter)) '.fig']);
            saveas(fig,[saveToPath strcat('values of gamma vs iterGamma for current epsilon #',num2str(epsCounter)) '.jpg']);
            close(fig);
        end;
        
        % save optimal gamma found for current epsilon
        gammaEpsilonMLSP(epsCounter) = gamma(end);
        
        % find row sparsity of the tresholded solution
        rowSparsityMLSPVsEpsilon(epsCounter) = sum(XsolMLSP12normForEachEps(:,epsCounter)~=0);
        
    end; % epsilonMLSPRange
    
    % display the row sparsity of the solution vs different values of
    % epsilon
    fig = figure('Name','Row sparsity of MLSP solution vs different values of epsilon');
    epsilonMLSPpositionIndex = find(epsilonMLSPrange == epsilonMLSP);
    plot(1:epsilonMLSPrangeLength,rowSparsityMLSPVsEpsilon,'-*',1:epsilonMLSPrangeLength,rowSparsity*ones(1,epsilonMLSPrangeLength),'-ro',epsilonMLSPpositionIndex, 0:0.02:max(rowSparsityMLSPVsEpsilon)+1,'g');
    hlpLegendRowSPMLSP = legend('rowSparsityMLSPVsEpsilon','true row-sparsity','true value of epsilonMLSP found by the algorithm');
    set(hlpLegendRowSPMLSP,'Location','SouthOutside')
    grid on
    xlabel('Different values of epsilon')
    ylabel('Row sparsity of MLSP solution')
    title('Row sparsity of MLSP solution vs different values of epsilon')
    saveas(fig,[saveToPath 'Row sparsity of MLSP solution vs different values of epsilon''.fig']);
    saveas(fig,[saveToPath 'Row sparsity of MLSP solution vs different values of epsilon' '.jpg']);
    saveas(fig,[saveToPath 'Row sparsity of MLSP solution vs different values of epsilon' '.eps']);
    close(fig);
    
    if display12norm == 1
        disp('XsolMLSP12norms for different epsilons and 12norm of true X ');
        disp('**************************************************************');
        disp(XsolMLSP12normForEachEps);
        disp(X12norm);
    end;
    
    if displayAuxiliary == 1
        % display optimal gammas found for each epsilon in the range
        disp('optimal gammas found for each epsilon in the range ');
        disp('***************************************************');
        disp(gammaEpsilonMLSP);
    end;
    
end; % if rowSparsityVSepsilonMLSP == 1
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ML FRAMEWORK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!! LS FRAMEWORKS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% LS framewotk NO SPARSITY enforced
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if runLSNOSPflag == 1
    disp('********************');
    disp('LSNOSP framework ');
    disp('********************');
    
    [XsolLSNOSP,XsolLSNOSP12norm,objFuncLSNOSP] = fLSNOSP(A,Y,XsolInitial,NiterLSNOSP,epsilonThreshold);
    
    if displayMatr==1
        disp('solution of LSNOSP matrix XsolLSNOSP');
        disp('********************************');
        disp(XsolLSNOSP);
    end;
    
    if displayAuxiliary == 1
        % display objective function for LSNOSP
        fig = figure('Name','Objective function LSNOSP');
        semilogy(1:NiterLSNOSP,objFuncLSNOSP);
        grid on
        xlabel('Iteration')
        ylabel('Objective function LSNOSP')
        title('Objective function LSNOSP')
        saveas(fig,[saveToPath 'Objective function LSNOSP' '.fig']);
        saveas(fig,[saveToPath 'Objective function LSNOSP' '.jpg']);
        close(fig);
    end;
    
    if display12norm == 1
        % display numerically and compare 12norm of true X and a solution to LSNOSP
        disp('12norm of true X and XsolLSNOSP');
        disp('*******************************');
        disp([X12norm,XsolLSNOSP12norm]);
    end;
    
end; % LSNOSP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% LS framewotk SPARSITY enforced
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if runLSSPflag == 1
    disp('********************');
    disp('LSSP framework ');
    disp('********************');
    
    % check that true solution is inside epsilonLSSP ball
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    constrTrueX = 0;
    for iN = 1:N
        % value of the constrain for true X
        constrTrueX = constrTrueX + norm(Lambda(:,iN)-Y(:,iN),2)^2;
    end;
    % check if Xtrue is inside epsilon ball
    if constrTrueX < epsilonLSSP
        disp('   ');
        disp('True X is inside epsilonLSSP ball');
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
    else
        disp('   ');
        disp('True X is NOT inside epsilonLSSP ball');
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
    end;
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
    [XsolLSSP,XsolLSSP12norm,objFuncLSSP,constrainLSSP,delta] = fLSSP3(A,Y,XsolInitial,epsilonLSSP,NiterLSSP,NiterDelta,deltaMinG,deltaMaxG,epsilonThreshold);
    
    % find support/row-sparsity pattern of the XsolLSSP
    sppatXsolLSSP = XsolLSSP12norm(end,:);
    indsppatXsolLSSP = XsolLSSP12norm(end,:)~=0;
    sppatXsolLSSP(indsppatXsolLSSP) = 1;
    
    if displayMatr==1
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
        saveas(fig,[saveToPath 'Objective function LSSP' '.fig']);
        saveas(fig,[saveToPath 'Objective function LSSP' '.jpg']);
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
        saveas(fig,[saveToPath 'ConstrainLSSP and true epsilonLSSP' '.fig']);
        saveas(fig,[saveToPath 'ConstrainLSSP and true epsilonLSSP' '.jpg']);
        close(fig);
        
        % display values of delta vs iterDelta
        fig = figure('Name','values of delta vs iterDelta');
        plot(1:NiterDelta,delta);
        xlabel('iterDelta')
        ylabel('values of delta ')
        grid on
        title('values of delta vs iterDelta');
        saveas(fig,[saveToPath 'values of delta vs iterDelta' '.fig']);
        saveas(fig,[saveToPath 'values of delta vs iterDelta' '.jpg']);
        close(fig);
    end;
    
    if display12norm == 1
        % display numerically and compare 12norm of true X and a solution to LSSP
        disp('12norm of true X and XsolLSSP');
        disp('*****************************');
        disp([X12norm,XsolLSSP12norm(end,:)']);
    end;
end; % LSSP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% LSNOSP on Support framework
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if runLSNOSPonSupportFlag == 1
    disp('********************');
    disp('LSNOSP on support framework ');
    disp('********************');
    disp('');
    
    % set the support that will be used for this code: either found by LSSP
    % or true support of X (i.e. OracleLSNOSP)
    if SupportLSNOSPflag == 1
        % then use support found by LSSP method
        SupportLSNOSPonSup = sppatXsolLSSP;
        disp('support found by LSSP');
    else
        % then use support of true X
        SupportLSNOSPonSup = sppatX;
        disp('Oracle LSNOSP');
    end;
    
    % run LSNOSP on support
    [XsolLSNOSPsupFull,XsolLSNOSPsupFull12norm,objFuncLSNOSPsup] = fLSNOSPonSupport(A,Y,NiterLSNOSPsup,epsilonThreshold,SupportLSNOSPonSup);
    
    if displayMatr==1
        disp('solution of LSNOSP on Support');
        disp('********************************');
        disp(XsolLSNOSPsupFull);
    end;
    
    if displayAuxiliary == 1
        % display objective function for LSNOSPonSupport
        fig = figure('Name','Objective function LSNOSPsup');
        plot(1:NiterLSNOSPsup,objFuncLSNOSPsup);
        grid on
        xlabel('Iteration')
        ylabel('Objective function LSNOSPsup')
        title('Objective function LSNOSPsup')
        saveas(fig,[saveToPath 'Objective function LSNOSPsup' '.fig']);
        saveas(fig,[saveToPath 'Objective function LSNOSPsup' '.jpg']);
        close(fig);
    end;
    
    if display12norm == 1
        % display numerically and compare 12norm of true X and a solution
        % to LSNOSPonSupport
        disp('12norm of true X and XsolLSNOSPonSupport');
        disp('*****************************');
        disp([X12norm,XsolLSNOSPsupFull12norm]);
    end;
end; % LSNOSPonSupport
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Row-sparsity VS epsilon analysis for LSSP framework
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% if 1 then produce the plot depicting row sparsity of LSSP solutions for
% different values of epsilon
if runRowSparsityVSepsilonLSSP == 1
    disp('****************************************************');
    disp('Row-sparsity VS epsilon analysis for LSSP framework ');
    disp('****************************************************');
    
    % define range for epsilons
    % epsilonLSSPrange = [epsilonLSSP/64,epsilonLSSP/32,epsilonLSSP/16,epsilonLSSP/8,epsilonLSSP/4,epsilonLSSP/2,epsilonLSSP,2*epsilonLSSP,4*epsilonLSSP,8*epsilonLSSP,16*epsilonLSSP,32*epsilonLSSP,64*epsilonLSSP];
    epsilonLSSPrange = [epsilonLSSP/2,epsilonLSSP,epsilonLSSP*2];
    % number of differents epsilons
    epsilonLSSPrangeLength = length(epsilonLSSPrange);
    % array that contains row sparsity of the solution for current epsilon
    rowSparsityLSSPVsEpsilon = zeros(epsilonLSSPrangeLength,1);
    % initialize counter
    epsCounter = 0;
    % save 12 norms of the solutions for different epsilon
    XsolLSSP12normForEachEps = zeros(K,epsilonLSSPrangeLength);
    % save optimal delta found for each epsilon
    deltaEpsilonLSSP = zeros(1,epsilonLSSPrangeLength);
    % for every epsilon in the range find solution and its 12 norm
    for iterEpsilonLSSP = epsilonLSSPrange
        % number of the epsilon in the range of epsilons
        epsCounter = epsCounter + 1;
        disp(strcat('current epsilon is epsilon #',num2str(epsCounter),' out of ',num2str(epsilonLSSPrangeLength)));
        
        % check if the inequality is satisfied
        indsppatX = X12norm~=0;
        disp('   ');
        if min(X12norm(indsppatX))>2*sqrt(iterEpsilonLSSP)/(1-ripconst)
            disp('Inequality for gamma is satisfied');
            disp('****************************');
        else
            disp('Inequality for gamma is NOT satisfied');
            disp('****************************');
        end;
        disp('   ');
        
        % for a given epsilonLSSP solve the LSSP problem
        [XsolLSSP,XsolLSSP12norm,objFuncLSSP,constrainLSSP,delta] = fLSSP3(A,Y,XsolInitial,iterEpsilonLSSP,NiterLSSP,NiterDelta,deltaMinG,deltaMaxG,epsilonThreshold);
        
        % save 12norm that we got for current epsilon
        XsolLSSP12normForEachEps(:,epsCounter) = XsolLSSP12norm(end,:)';
        
        if displayAuxiliary == 1
            % display objective function for LSSP for first iterDelta
            fig = figure('Name',strcat('Objective function LSSP for epsilon #',num2str(epsCounter)));
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
            saveas(fig,[saveToPath strcat('Objective function LSSP for epsilon #',num2str(epsCounter)) '.fig']);
            saveas(fig,[saveToPath strcat('Objective function LSSP for epsilon #',num2str(epsCounter)) '.jpg']);
            close(fig);
            
            % display constrainLSSP and current epsilonLSSP, i.e. iterEpsilonLSSP
            fig = figure('Name',strcat('constrainLSSP, true epsilonLSSP and current epsilon #',num2str(epsCounter)));
            plot(1:NiterDelta,constrainLSSP,'-b*',1:NiterDelta,iterEpsilonLSSP*ones(1,NiterDelta),'r',1:NiterDelta,epsilonLSSP*ones(1,NiterDelta),'m');
            hlpLegendconstrLSSP = legend('constrainLSSP','iterEpsilonLSSP','true epsilon LSSP');
            set(hlpLegendconstrLSSP,'Location','SouthOutside')
            grid on
            xlabel('Iteration for delta')
            ylabel(strcat('constrainLSSP, true epsilonLSSP and current epsilon #',num2str(epsCounter)));
            title(strcat('constrainLSSP, true epsilonLSSP and current epsilon #',num2str(epsCounter)));
            saveas(fig,[saveToPath strcat('constrainLSSP, true epsilonLSSP and current epsilon #',num2str(epsCounter)) '.fig']);
            saveas(fig,[saveToPath strcat('constrainLSSP, true epsilonLSSP and current epsilon #',num2str(epsCounter)) '.jpg']);
            close(fig);
                        
            % display values of delta vs iterDelta
            fig = figure('Name',strcat('values of delta vs iterDelta for current epsilon #',num2str(epsCounter)));
            plot(1:NiterDelta,delta);
            xlabel('iterDelta')
            ylabel('values of delta ')
            grid on
            title(strcat('values of delta vs iterDelta for current epsilon #',num2str(epsCounter)));
            saveas(fig,[saveToPath strcat('values of delta vs iterDelta for current epsilon #',num2str(epsCounter)) '.fig']);
            saveas(fig,[saveToPath strcat('values of delta vs iterDelta for current epsilon #',num2str(epsCounter)) '.jpg']);
            close(fig);
        end;
        
        % save optimal delta found for current epsilon
        deltaEpsilonLSSP(epsCounter) = delta(end);
        
        % find row sparsity of the tresholded solution
        rowSparsityLSSPVsEpsilon(epsCounter) = sum(XsolLSSP12normForEachEps(:,epsCounter)~=0);
        
    end; % epsilonLSSPRange
    
    % display the row sparsity of the solution vs different values of
    % epsilon
    fig = figure('Name','Row sparsity of LSSP solution vs different values of epsilon');
    epsilonLSSPpositionIndex = find(epsilonLSSPrange == epsilonLSSP);
    plot(1:epsilonLSSPrangeLength,rowSparsityLSSPVsEpsilon,'-*',1:epsilonLSSPrangeLength,rowSparsity*ones(1,epsilonLSSPrangeLength),'-ro',epsilonLSSPpositionIndex, 0:0.02:max(rowSparsityLSSPVsEpsilon)+1,'g');
    hlpLegendRowSPLSSP = legend('rowSparsityLSSPVsEpsilon','true row-sparsity','true value of epsilonLSSP found by the algorithm');
    set(hlpLegendRowSPLSSP,'Location','SouthOutside')
    grid on
    xlabel('Different values of epsilon')
    ylabel('Row sparsity of LSSP solution')
    title('Row sparsity of LSSP solution vs different values of epsilon')
    saveas(fig,[saveToPath 'Row sparsity of LSSP solution vs different values of epsilon' '.fig']);
    saveas(fig,[saveToPath 'Row sparsity of LSSP solution vs different values of epsilon' '.jpg']);
    saveas(fig,[saveToPath 'Row sparsity of LSSP solution vs different values of epsilon' '.eps']);
    close(fig);
    
    if display12norm == 1
        disp('XsolLSSP12norms for different epsilons and 12norm of true X ');
        disp('**************************************************************');
        disp(XsolLSSP12normForEachEps);
        disp(X12norm);
    end;
    
    if displayAuxiliary == 1
        % display optimal deltas found for each epsilon in the range
        disp('optimal deltas found for each epsilon in the range ');
        disp('***************************************************');
        disp(deltaEpsilonLSSP);
    end;
    
end; % if rowSparsityVSepsilonLSSP == 1
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LS FRAMEWORK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%%%%%%%%%%%%%%%%%%%% PART 2:  MAIN FUNCTIONAL PART %%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      PART 3: DISPLAY RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%       Display the results
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% display numerically and compare 12norm of true X and all solutions
if runMLSPflag && runMLNOSPflag && runLSNOSPflag && runLSSPflag && ~runMLNOSPonSupportFlag && ~runLSNOSPonSupportFlag && ~runRowSparsityVSepsilonMLSP && ~runRowSparsityVSepsilonLSSP
    disp('12norm of true X and 4 solutions: MLNOSP, MLSP, LSNOSP, LSSP');
    disp('**********************************');
    disp([X12norm,XsolMLNOSP12norm,XsolMLSP12norm(end,:)',XsolLSNOSP12norm,XsolLSSP12norm(end,:)']);
end;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 3:  RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tm = toc;
disp(strcat('running time = ', num2str(tm/60)));
save all
diary off
