function [fA,fArowsum,fX,fX12norm,fsppatX,fY,fripconst,fepsilonThresholdPercentage,fepsilonLSSP] = fGenerateData(fK,fN,fM,frowSparsity)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates mixing matrices A's, initial true row-sparse
% matrix X and matrix of Poisson counts Y.
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% NB!!!: choose one of two options of calculating RIP of A{i}'s:
%        either RIP of concatinated A (set fconcatinatedARIPflag to 1) or
%        choose max of RIP's for each A{i} (set fseparatelyARIPflag to 1)
%
fconcatenatedARIPflag = 0;
fconcatIter = fK*fN*10;      % number of random samplings from column space
fseparatelyARIPflag = 1;
fseparatIter = fK*10;        % number of random samplings from column space
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% INPUT: K,N,M - dimensions of matrices A's,X,Y such that N<<M<<K and
%        A_i is MxK, i = 1..N, X is KxN, Y is MxN.
%        rowSparsity - desired row-sparsity of matrix X, i.e. number of
%        nonzero rows of X
% OUTPUT: A - cell(1,N) that contains mixing matrices A's
%         Arowsum - row j of Arowsum contains sum of rows of matrix A{j}
%         X - true row-sparse matrix
%         X12norm - 12norm of matrix X, such that min 2 norm of rows of X
%         (called gamma in the paper) is greater then the
%         (2*sqrt(epsilon)/1-delta_2s_o)
%         sppatX - sparsity pattern of X, i.e. vector of 1's and 0's (1's
%         are located on the positions of the non-zero rows)
%         Y - matrix of Poisson counts
%         ripconst - RIP constant of the matrices A{i} (two options :
%         concatinated or separate)
%         epsilonThresholdPercentage - second epsilon (first was
%         epsilonThresholdFixed) to threshold the solution matrices X, is
%         defined as 0.1% of the smallest value of non-zero elements of X
%         epsilonLSSP - epsilon for method LSSP, is used for LSSP function
%         is used to check that inequality for the min 2 norm of rows of X
%         is satisfied, i.e. (min 2 norm of rows is called gamma in the
%         paper) : gamma > (2*sqrt(epsilonLSSP)/1-delta_2s_o)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% NB!!! These parameters are defined/changed in this code:
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% for defining matrices A's
% fdefinitionParameterA = 0.7;
% for defining matrix X large enought (to satisf inequality for min 2 norm)
 fmultParam = 1e4;
% for thresholding: percentage of min 2 norm 
% fPercentThre = 0.001;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

disp('   ');
disp('******************** GENERATED DATA ********************')
disp('   ');

% A_i satisfy 2so-RIP
fripSp = 2*frowSparsity;

 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Analize the dependence of fripconst on fdefinitionParameterA:
% % what fdefinitionParameterA should be in order to generate matrices A such
% % that they will satisfy RIP.
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % first generate matrices A for fdefinitionParameterA = [0.1,0.2,...1],
% % then calculate fripconst for A's and take the max one. Plot fripconst for
% % every 1-fdefinitionParameterA (or fdefinitionParameterA) and chose such
% % fdefinitionParameterA for which ripconst will be less then 1, i.e. for
% % which A's will satisfy RIP
% % define array that contains ripconst for different fdefinitionParameterA
% fripcontsArray = zeros(1,10);
% 
% for fidp = 1:10
%     
%     %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     %           Define mixing matrices A's
%     %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     % for every i = 1,N: element A_i(*,*)>fdefinitionParameterA
%     % then A_i(*,*) = 1 otherwise A_i(*,*) = 0
%     % define fdefinitionParameterA
%     fdefinitionParameterA = 0.1*fidp;
%     
%     fA = cell(1,fN);
%     for fiA = 1:fN
%         fA{fiA} = double(rand(fM,fK)>fdefinitionParameterA);
%         for fiColumn = 1:fK
%             % for every column check if this column consist only of zeros then
%             % add 1 on random position of the column
%             if sum(fA{fiA}(:,fiColumn))==0
%                 fA{fiA}(randi(fM,1,1),fiColumn) = 1;
%             end;
%             % normalize all columns
%             fA{fiA}(:,fiColumn) = fA{fiA}(:,fiColumn)/norm(fA{fiA}(:,fiColumn),2);
%         end;
%         
%         % NB!!! need this part to make every element of A_i nonzero, to get
%         % nonzero lambdas for CRLB(X). I general can comment this line
%         fA{fiA} = fA{fiA}+1e-7;
%     end;
%     
%     %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     %    Calculate RIPS of A's
%     %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     
%     % OPTION 1: calculate RIP of concatinated A's
%     %----------
%     if fconcatenatedARIPflag == 1
%         fconcatenatedA = zeros(fM,fN*fK);
%         for fiN = 1:fN
%             % concatinate all matrices A into one matrix MxKN
%             fconcatenatedA(:,1+(fiN-1)*fK:fK*fiN) = fA{fiN};
%         end;
%         % iteratively choose random columns of concatinatedA and find RIP that
%         % corresponds to submatrices consisting of those columns
%         fhlp2 = zeros(fconcatIter,1);
%         for iconcat=1:fconcatIter
%             fhlp1=randperm(fN*fK);
%             fhlp1=fhlp1(1:fripSp);
%             fxtmp=zeros(fN*fK,1);
%             fxtmp(fhlp1)=randn(fripSp,1);
%             fhlp2(iconcat)=norm(fconcatenatedA*fxtmp,2)/norm(fxtmp,2);
%         end
%         fhlp3=sort(fhlp2,'ascend');
%         fripconstLowConcat=1-fhlp3(round(fconcatIter*0.05));
%         fripconstHighConcat=fhlp3(round(fconcatIter*0.95))-1;
%         % choose the RIP constant
%         fripconst = max(fripconstLowConcat,fripconstHighConcat);
%     end; % OPTION 1 (if fconcatinatedARIPflag ==1)
%     
%     % OPTION 2: calculate RIP for each Ai separately and choose max of them
%     %----------
%     if fseparatelyARIPflag == 1
%         % find RIP for each A{i} separately
%         fRIPA = zeros(fN,1);
%         for fiN = 1:fN
%             fhlp1 = zeros(1,fseparatIter);
%             for fiter=1:fseparatIter;
%                 fhlp2=randperm(fK);
%                 fhlp2=fhlp2(1:fripSp);
%                 fxtmp=zeros(fK,1);
%                 fxtmp(fhlp2)=randn(fripSp,1);
%                 fhlp1(fiter)=norm(fA{fiN}*fxtmp,2)^2/norm(fxtmp,2)^2;
%             end;
%             fhlp3 = sort(fhlp1,'ascend');
%             fripconstLowAi = 1-fhlp3(round(fseparatIter*0.05));
%             fripconstHighAi = fhlp3(round(fseparatIter*0.95))-1;
%             % save RIP for current A
%             fRIPA(fiN) = max(fripconstLowAi,fripconstHighAi);
%         end;
%         % display RIPs for all A{i}'s
%         disp('   ');
%         disp('RIPs for A{i}s');
%         disp('**********************');
%         disp(fRIPA);
%         disp('   ');
%         % as RIP for all A's choose min of individual RIPs
%         fripconst = max(fRIPA);
%     end; % OPTION 2 (if fseparatelyARIPflag == 1)
%     
%     % save fripconst for current fdefinitionParameterA
%     fripcontsArray(fidp) = fripconst;
%     
%     % display the RIP
%     if fseparatelyARIPflag == 1
%         disp(strcat('RIP const = ',num2str(fripconst),' {separately}'));
%         disp('****************************');
%     else
%         disp(strcat('RIP const = ',num2str(fripconst),' {concatenated}'));
%         disp('****************************');
%     end;
%     
% end; % fidp
% 
% % Plot ripconts vs fdefinitionParameterA
% fhlpprob = 1*ones(10,1)-0.1*(1:10)';
% plot(fhlpprob,fripcontsArray)
% title('ripcont vs 1-definitionParameterA')
% grid on
% disp('ripcont vs 1-definitionParameterA');
% disp('**********************************');
% disp([fhlpprob,fripcontsArray']);
% 
% pause
% 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Analize ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%           Define mixing matrices A's
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% for every i = 1,N: element A_i(*,*)>fdefinitionParameterA
% then A_i(*,*) = 1 otherwise A_i(*,*) = 0
% define fdefinitionParameterA
fdefinitionParameterA = 0.7;

fA = cell(1,fN);
for fiA = 1:fN
    fA{fiA} = double(rand(fM,fK)>fdefinitionParameterA);
    for fiColumn = 1:fK
        % for every column check if this column consist only of zeros then
        % add 1 on random position of the column
        if sum(fA{fiA}(:,fiColumn))==0
            fA{fiA}(randi(fM,1,1),fiColumn) = 1;
        end;
        % normalize all columns
        fA{fiA}(:,fiColumn) = fA{fiA}(:,fiColumn)/norm(fA{fiA}(:,fiColumn),2);
    end;
    
    % NB!!! need this part to make every element of A_i nonzero, to get
    % nonzero lambdas for CRLB(X). I general can comment this line
    fA{fiA} = fA{fiA}+1e-7;
end;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%    Calculate RIPS of A's
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% OPTION 1: calculate RIP of concatinated A's
%----------
if fconcatenatedARIPflag == 1
    fconcatenatedA = zeros(fM,fN*fK);
    for fiN = 1:fN
        % concatinate all matrices A into one matrix MxKN
        fconcatenatedA(:,1+(fiN-1)*fK:fK*fiN) = fA{fiN};
    end;
    % iteratively choose random columns of concatinatedA and find RIP that
    % corresponds to submatrices consisting of those columns
    fhlp2 = zeros(fconcatIter,1);
    for iconcat=1:fconcatIter
        fhlp1=randperm(fN*fK);
        fhlp1=fhlp1(1:fripSp);
        fxtmp=zeros(fN*fK,1);
        fxtmp(fhlp1)=randn(fripSp,1);
        fhlp2(iconcat)=norm(fconcatenatedA*fxtmp,2)/norm(fxtmp,2);
    end
    fhlp3=sort(fhlp2,'ascend');
    fripconstLowConcat=1-fhlp3(round(fconcatIter*0.05));
    fripconstHighConcat=fhlp3(round(fconcatIter*0.95))-1;
    % choose the RIP constant
    fripconst = max(fripconstLowConcat,fripconstHighConcat);
end; % OPTION 1 (if fconcatinatedARIPflag ==1)

% OPTION 2: calculate RIP for each Ai separately and choose max of them
%----------
if fseparatelyARIPflag == 1
    % find RIP for each A{i} separately
    fRIPA = zeros(fN,1);
    for fiN = 1:fN
        fhlp1 = zeros(1,fseparatIter);
        for fiter=1:fseparatIter;
            fhlp2=randperm(fK);
            fhlp2=fhlp2(1:fripSp);
            fxtmp=zeros(fK,1);
            fxtmp(fhlp2)=randn(fripSp,1);
            fhlp1(fiter)=norm(fA{fiN}*fxtmp,2)^2/norm(fxtmp,2)^2;
        end;
        fhlp3 = sort(fhlp1,'ascend');
        fripconstLowAi = 1-fhlp3(round(fseparatIter*0.05));
        fripconstHighAi = fhlp3(round(fseparatIter*0.95))-1;
        % save RIP for current A
        fRIPA(fiN) = max(fripconstLowAi,fripconstHighAi);
    end;
    % display RIPs for all A{i}'s
    disp('   ');
    disp('RIPs for A{i}s');
    disp('**********************');
    disp(fRIPA);
    disp('   ');
    % as RIP for all A's choose min of individual RIPs
    fripconst = max(fRIPA);
end; % OPTION 2 (if fseparatelyARIPflag == 1)

% display the RIP
if fseparatelyARIPflag == 1
    disp(strcat('RIP const = ',num2str(fripconst),' {separately}'));
    disp('****************************');
else
    disp(strcat('RIP const = ',num2str(fripconst),' {concatenated}'));
    disp('****************************');
end;


%~~~~~~~~~~~~~~~~~~~~~~~~~~
% Calculate matrix Arowsum
%~~~~~~~~~~~~~~~~~~~~~~~~~~
% row j of matrix Arowsum contains sum of rows of matrix A{j}
fArowsum = zeros(fN,fK);
for fiN = 1:fN
    fArowsum(fiN,:) = sum(fA{fiN},1);
end;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     Define X - original row-sparse matrix
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%fmultParam = 1e2; % is defined at the beginning of this code
fX = zeros(fK,fN);
fp = randperm(fK);
for fiColumn = 1:fN
    fx = zeros(fK,1);
    fx(fp(1:frowSparsity)) = fmultParam*randn(frowSparsity,1);
    fX(:,fiColumn) = fx;
end;
% make X nonnegative
fX = abs(fX);
% row-norm (12norm) of the true matrix X
fX12norm = sqrt(sum(fX.*fX,2));
disp('   ');
disp('2 norm of rows of true X');
disp('****************************');
disp(fX12norm);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% find row-sparsity pattern
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fsppatX = zeros(fK,1);
findsppatX = fX12norm~=0;
fsppatX(findsppatX) = 1;
disp('   ');
disp('sparsity pattern of true X');
disp('****************************');
disp(fsppatX);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%           Define Y
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% generate Y(:,*) = Poisson(A{i}*X(:,i))
fY = zeros(fM,fN);
for fiN = 1:fN
    fparametr = fA{fiN}*fX(:,fiN);
    fY(:,fiN) = poissrnd(fparametr);
end;
disp('   ');
disp('matrix of Poisson counts Y');
disp('****************************');
disp(fY);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% calculate epsilon for LSSP method
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% probability s.t. P(constraint<epsilonLSSP)>1-p
fp = 0.1;
fsumY = sum(sum(fY));
% epsilon without simplification (look in the notes method 2 epsilon 2)
fhlp = sqrt((2/fp)-1);
fepsilonLSSP = fsumY+(fhlp^2)/2+sqrt(fhlp^2*fsumY+fhlp^4/4)+fhlp*sqrt(2*fsumY^2+fsumY*(4*fhlp^2+1)+fhlp^4+fhlp^2/2+(4*fsumY+2*fhlp^2+1)*sqrt(fhlp^2*fsumY+fhlp^4/4));

% % epsilon using simplification
% fepsilonLSSP = ( fsumY + (fp-1/2) + sqrt( (2/fp-1)*fsumY + ((1/fp)-(1/2))^2 ) )*(1 + sqrt(6/fp-3));
% % epsilon that uses Chebyshev's and Markov's inequalities
% fbeta = 1.01;
% fepsilonLSSP31 = (2/fp)*max(fbeta*fsumY/(fbeta-1),2*fbeta^2/fp);
% fbeta = 1.99;
% fepsilonLSSP32 = (2/fp)*max(fbeta*fsumY/(fbeta-1),2*fbeta^2/fp);
% % display all epsilons
% disp('   epsilon1,    epsilon2,    epsilon31,    epsilon32');
% fepsArray = [fepsilonLSSP,fepsilonLSSP,fepsilonLSSP31,fepsilonLSSP32];
% disp(fepsArray)


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% check inequality for gamma
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('   ');
disp(strcat('value of the epsilonLSSP = ',num2str(fepsilonLSSP)));
% disp(strcat('gamma = ',num2str(min(fX12norm(findsppatX)))));
% disp(strcat('2*sqrt(fepsilonLSSP2)/(1-fripconst) = ',num2str(2*sqrt(fepsilonLSSP)/(1-fripconst))));
if min(fX12norm(findsppatX))>2*sqrt(fepsilonLSSP)/(1-fripconst)
    disp('Inequality for gamma is satisfied');
    disp('****************************');
else
    disp('Inequality for gamma is NOT satisfied');
    disp('****************************');
end;
disp('   ');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% set epsilonThresholdPercentage as fPercentThre*100% of of the smallest
% 2 norm value of non-zero row of X
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fPercentThre = 0.001;
fepsilonThresholdPercentage  = fPercentThre*min(fX12norm(findsppatX));
