function [fXsolMLSP,fXsolMLSP12norm,fobjFuncMLSP,fconstrainMLSP,fgamma] = fMLSP(fA,fArowsum,fY,fXsolInitial,fepsilonMLSP,fNiterMLSP,fNiterGamma,fgammaMin,fgammaMax,fepsilonThreshold)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function finds ML solution with enforced sparsity and thresholds 12 norm
% of the solution with epsilonThreshold, then nulls those rows of X that
% were thresholded in 12norm
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% INPUT: A - cell(1,N) that contains mixing matrices A's
%        Arowsum - matrix, ith row of which is sum of rows of A{i}
%        Y - matrix of Poisson counts
%        XsolInitial - initial matrix X solution: one initial X solution
%        for all algorithms
%        epsilonMLSP - epsilon for MLSP method
%        NiterMLSP - number of iterations for iterative algorithm MLSP
%        NiterGamma - number of iteration to find optimal gamma
%        gammaMin - minimal value of gamma range
%        gammaMax - maximal value of gamma range
%        epsilonThreshold - epsilon for thresholding solution mtrix
%        XsolMLSP
% OUTPUT: XsolMLSP - solution matrix
%         XsolMLSP12norm - 12norm of the solution matrix XsolMLSP, saved
%         for every iterGamma
%         objFuncMLSP - objective function that is minimized, save for every
%         iterGamma and every iretation for update X (fiter), when solution
%         matrix gets updated
%         constrainMLSP - values of the constains saved for every iterGamma
%         gamma - Lagrangian parameter saved for each iterGamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% NB!!! These parameters are defined/changed in this code:
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% max number of times gamma range can be increased/decreased
% fstepCntrMax = 20;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[~,fK] = size(fA{1});
[fM,fN] = size(fY);

% need to check that if zero matrix is inside the constrains then solution
% must be automatically defined as zero matrix
% if sum_{i=1,N,j=1,M}(y_ij)*(log(y_ij)-1) <= epsilonMLSP
fYindhlp = find(fY~=0);
fconstrHlp = sum(fY(fYindhlp).*log(fY(fYindhlp))-fY(fYindhlp));

if fconstrHlp <= fepsilonMLSP
    % it means zero point inside espilon ball, hence solution is zero
    % matrix
    fXsolMLSP = zeros(fK,fN);
    % 12 norm is a zero vector
    fXsolMLSP12norm = zeros(fNiterGamma,fK);
    % null objective function
    fobjFuncMLSP = zeros(1,fNiterMLSP);
    % constraints are const sum_{i=1,N,j=1,M}(y_ij)*(log(y_ij)-1)
    fconstrainMLSP = fconstrHlp*ones(fNiterGamma,1);
    % gamma can be any but here I make it zeros
    fgamma = zeros(fNiterGamma,1);
    disp('zero point is inside the constrained region - therefore solution is zero matrix!');
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
    
else % if epsilon ball doesn't contain zero point then run algorithm
    
    % initialize objective function calculated for each gamma, each iteration
    fobjFuncMLSP = zeros(fNiterGamma,fNiterMLSP);
    % initialize array that contains constrains for each gamma
    fconstrainMLSP = zeros(fNiterGamma,1);
    % initialize array that contains values of 12 norm of XsolMLSP for each
    % iterGamma
    fXsolMLSP12norm = zeros(fNiterGamma,fK);
    % initialize array that contains values of gamma for each iterGamma
    fgamma = zeros(fNiterGamma,1);
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %                     Find range for Gamma
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    disp('**********************************');
    disp('      Find range for Gamma');
    disp('**********************************');
    % flag that indicates if correct gamma range was found or not
    fgammaRangeFoundFlag = 0;
    
    % counter of times gamma range was increased/decreased
    fstepCntr = 0;
    % max number of times gamma range can be increased/decreased
    fstepCntrMax = 20;
    
    % check that for max gamma constraints are less then epsilon, for min gamma
    % constraints are greater then epsilon, if not decrease/increase gamma
    % range, do it no more then fstepCntrMax times
    
    while fgammaRangeFoundFlag == 0 && (fstepCntr<=fstepCntrMax)
        
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
        fstepCntr = fstepCntr+1;
        disp(strcat('current step = ',num2str(fstepCntr)));
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
        disp(strcat('current fGammaMax = ',num2str(fgammaMax)))
        disp(strcat('current fGammaMin = ',num2str(fgammaMin)))
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Check 1: check that for max delta constraints are less then epsilon
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % initialize Xsol
        fXsolMLSP = rand(fK,fN);
        % set current gamma to fgammaMax
        fgammaHlp = fgammaMax;
        
        % find solution corresponding for fgammaMax
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % find solution matrix XsolMLSP iteratively
        for fiter = 1:fNiterMLSP
            
            % vector of Eucledian norms of the rows of XsolMLSP
            fHlp0 = sqrt(sum(fXsolMLSP.^2,2));
            fHlp2 = repmat(fHlp0',fN,1);
            fHlp3 = -fgammaHlp*(fArowsum.*fHlp2);
            fbeta = zeros(fK,fN);
            for fiN = 1:fN % calculate matrix beta
                % A_i*x_i term
                fHlp6 = fA{fiN}*fXsolMLSP(:,fiN);
                % indices where y_i(j)~=0
                fiHlp5 = find(fY(:,fiN)~=0);
                fHlp1 = zeros(fM,1);
                fHlp1(fiHlp5) = fY(fiHlp5,fiN)./(fHlp6(fiHlp5));
                % column of auxiliary matrix beta
                fbeta(:,fiN) = (fHlp1'*fA{fiN}).*(fXsolMLSP(:,fiN)');
            end;
            fHlp4 = 4*fgammaHlp*(fbeta.*fHlp2');
            
            % update solution matrix XsolMLSP
            fXsolMLSP = ((fHlp3+sqrt(fHlp3.^2+fHlp4'))/2)';
            
            % projection on positive semi-axis, since Problem Domain is R^{K*N}_{+}
            fXsolMLSPind = fXsolMLSP<0;
            fXsolMLSP(fXsolMLSPind) = 0;
            
            % calculate constrains
            fHlp7 = 0;
            fHlp8 = 0;
            % calculate auxiliary terms for constrains
            for fiN = 1:fN
                fiHlp5 = find(fY(:,fiN)~=0);
                fHlp9 = fA{fiN}*fXsolMLSP(:,fiN);
                fHlp10 = fY(fiHlp5,fiN)'*log(fHlp9(fiHlp5));
                fHlp7 = fHlp7 + sum(fHlp10);
                fHlp8 = fHlp8 + fY(fiHlp5,fiN)'*log(fY(fiHlp5,fiN));
            end;
            fconstrainHlpMax = sum(diag(fArowsum*fXsolMLSP)) - fHlp7 - sum(sum(fY)) + fHlp8;
            
        end; % for iter
        
        % check that constraints for gammaMax are less then epsilon
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if fconstrainHlpMax > fepsilonMLSP
            disp('constr for max gamma > eps');
            fconstrGammaMaxFlag = 0;
        else
            fconstrGammaMaxFlag = 1;
            disp('constr for max gamma < eps');
        end;
        %~~~~~~~~~~~~~~~~~~~~~~% MAX %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Check 2: check that for min gamma constraints are greater then epsilon
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % initialize Xsol
        fXsolMLSP = rand(fK,fN);
        % set current gamma to min
        fgammaHlp = fgammaMin;
        
        % find solution corresponding for fgammaMin
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % find solution matrix XsolMLSP iteratively
        for fiter = 1:fNiterMLSP
            
            % vector of Eucledian norms of the rows of XsolMLSP
            fHlp0 = sqrt(sum(fXsolMLSP.^2,2));
            fHlp2 = repmat(fHlp0',fN,1);
            fHlp3 = -fgammaHlp*(fArowsum.*fHlp2);
            fbeta = zeros(fK,fN);
            for fiN = 1:fN % calculate matrix beta
                % A_i*x_i term
                fHlp6 = fA{fiN}*fXsolMLSP(:,fiN);
                % indices where y_i(j)~=0
                fiHlp5 = find(fY(:,fiN)~=0);
                fHlp1 = zeros(fM,1);
                fHlp1(fiHlp5) = fY(fiHlp5,fiN)./(fHlp6(fiHlp5));
                % column of auxiliary matrix beta
                fbeta(:,fiN) = (fHlp1'*fA{fiN}).*(fXsolMLSP(:,fiN)');
            end;
            fHlp4 = 4*fgammaHlp*(fbeta.*fHlp2');
            
            % update solution matrix XsolMLSP
            fXsolMLSP = ((fHlp3+sqrt(fHlp3.^2+fHlp4'))/2)';
            
            % projection on positive semi-axis, since Problem Domain is R^{K*N}_{+}
            fXsolMLSPind = fXsolMLSP<0;
            fXsolMLSP(fXsolMLSPind) = 0;
            
            % calculate constrains
            fHlp7 = 0;
            fHlp8 = 0;
            % calculate auxiliary terms for constrains
            for fiN = 1:fN
                fiHlp5 = find(fY(:,fiN)~=0);
                fHlp9 = fA{fiN}*fXsolMLSP(:,fiN);
                fHlp10 = fY(fiHlp5,fiN)'*log(fHlp9(fiHlp5));
                fHlp7 = fHlp7 + sum(fHlp10);
                fHlp8 = fHlp8 + fY(fiHlp5,fiN)'*log(fY(fiHlp5,fiN));
            end;
            fconstrainHlpMin = sum(diag(fArowsum*fXsolMLSP)) - fHlp7 - sum(sum(fY)) + fHlp8;
     
        end; % for iter
        
        % check that constraints for gammaMin are greater then epsilon
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if fconstrainHlpMin  > fepsilonMLSP
            disp('constr for min gamma > eps');
            fconstrGammaMinFlag = 1;
        else
            disp('constr for min gamma < eps');
            fconstrGammaMinFlag = 0;
        end;
        %~~~~~~~~~~~~~~~~~~~~~~% MIN %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        % check that constraints for gammaMax are less then epsilon and
        % constraints for gammaMin are greater then epsilon
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % if not - then change the constraints range and check again
        if fconstrGammaMaxFlag && fconstrGammaMinFlag
            % the correct range for gamma is found - exit while loop
            fgammaRangeFoundFlag = 1;
        else
            % if gammaMax is too small - increase gammaMax
            if fconstrGammaMaxFlag == 0
                % if constraints for gammaMin are greater then epsilon and
                % simultaneously constraints for gammaMax are greater then
                % epsilon then move gammaMin to gammaMax position
                if fconstrGammaMinFlag 
                    fgammaMin = fgammaMax;
                end;
                fgammaMax = fgammaMax*1e1;
            end;
            % if gammaMin is too big - decrease deltaMin
            if fconstrGammaMinFlag == 0
                % if constraints for gammaMin is less then epsilon and
                % simultaneously constraints for gammaMax is less then
                % epsilon then move gammaMax to gammaMin position
                if fconstrGammaMaxFlag 
                    fgammaMax = fgammaMin;
                end;
                fgammaMin = fgammaMin*1e-1;
            end;
        end;
    end; % while
    
    % display DeltaMin and DeltaMax that were found
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    disp('**********************************');
    disp(strcat('fGammaMax = ',num2str(fgammaMax)))
    disp(strcat('fGammaMin = ',num2str(fgammaMin)))
    if fgammaRangeFoundFlag
        disp('**********************************');
        disp('      Range for gamma is found ');
        disp('**********************************');
    else
        disp('**********************************');
        disp('Exit gamma range search loop (max number of steps)');
        disp('**********************************');
    end;
    %~~~~~~~~~~~~~~~~~~~~~~~ Find range for Gamma ~~~~~~~~~~~~~~~~~~~~~~~~~
    
       
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %                     Find solution matrix
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % initialize Xsol
    fXsolMLSP = fXsolInitial;
       
    % gamma iterations
    for fiGamma = 1:fNiterGamma
        
        % set Lagrangian parameter gamma
        fgamma(fiGamma)=(fgammaMax+fgammaMin)/2;
        
        % find solution matrix XsolMLSP iteratively
        for fiter = 1:fNiterMLSP
            
            % vector of Eucledian norms of the rows of XsolMLSP
            fHlp0 = sqrt(sum(fXsolMLSP.^2,2));
            fHlp2 = repmat(fHlp0',fN,1);
            fHlp3 = -fgamma(fiGamma)*(fArowsum.*fHlp2);
            fbeta = zeros(fK,fN);
            for fiN = 1:fN % calculate matrix beta
                % A_i*x_i term
                fHlp6 = fA{fiN}*fXsolMLSP(:,fiN);
                
                % indices where y_i(j)~=0
                fiHlp5 = find(fY(:,fiN)~=0);
                fHlp1 = zeros(fM,1);
                fHlp1(fiHlp5) = fY(fiHlp5,fiN)./(fHlp6(fiHlp5));
                                               
                % column of auxiliary matrix beta
                fbeta(:,fiN) = (fHlp1'*fA{fiN}).*(fXsolMLSP(:,fiN)');
               
            end;
           
            fHlp4 = 4*fgamma(fiGamma)*(fbeta.*fHlp2');
                       
            % update solution matrix XsolMLSP
            fXsolMLSP = ((fHlp3+sqrt(fHlp3.^2+fHlp4'))/2)';
                                               
            % projection on positive semi-axis, since Problem Domain is R^{K*N}_{+}
            fXsolMLSPind = fXsolMLSP<0;
            fXsolMLSP(fXsolMLSPind) = 0;
            
            % calculate constrains
            fHlp7 = 0;
            fHlp8 = 0;
            % calculate auxiliary terms for constrains
            for fiN = 1:fN
                fiHlp5 = find(fY(:,fiN)~=0);
                fHlp9 = fA{fiN}*fXsolMLSP(:,fiN);
                fHlp10 = fY(fiHlp5,fiN)'*log(fHlp9(fiHlp5));
                fHlp7 = fHlp7 + sum(fHlp10);
                fHlp8 = fHlp8 + fY(fiHlp5,fiN)'*log(fY(fiHlp5,fiN));
            end;
            fconstrainMLSPhlp = sum(diag(fArowsum*fXsolMLSP)) - fHlp7 - sum(sum(fY)) + fHlp8;
            
            % calculate objective function on current iteration
            fobjFuncMLSP(fiGamma,fiter) = sum(sqrt(sum(fXsolMLSP.^2,2))) + fgamma(fiGamma)*(fconstrainMLSPhlp-fepsilonMLSP);
            %fobjFuncMLSP(fiGamma,fiter) = sum(sqrt(sum(fXsolMLSP.^2,2)));
        end; % for iter
        
        % save constrain for current iteration gamma
        fconstrainMLSP(fiGamma) = fconstrainMLSPhlp;
        
        % calculate new range for gamma
        if fconstrainMLSP(fiGamma) > fepsilonMLSP
            % then increase gamma
            fgammaMin = fgamma(fiGamma);
        else
            % then decrease gamma
            fgammaMax = fgamma(fiGamma);
        end;
        
        % save value of 12 norm of the solution for current value of Gamma
        fXsolMLSP12norm(fiGamma,:) = sqrt(sum(fXsolMLSP.*fXsolMLSP,2));
        
    end; % for iGamma
    
    % threshold the 12 norm of the rows of the solution matrix XsolMLSP
    findThr = fXsolMLSP12norm(fiGamma,:)<=fepsilonThreshold;
    fXsolMLSP12norm(fiGamma,findThr) = 0;
    
    % null those rows of the solution matrix XsolMLSP which were thresholded by
    % the epsilonThreshold in XsolMLSP12norm (we need this for MSE analysis)
    fXsolMLSP(findThr==1,:) = 0;
    
end; %{if zero point is in epsilon ball}
