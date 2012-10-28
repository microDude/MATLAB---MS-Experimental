function [fXsolLSSP,fXsolLSSP12norm,fobjFuncLSSP,fconstrainLSSP,fdelta] = fLSSP3_revA(fA,fY,fXsolInitial,fepsilonLSSP,fNiterE,fNiterDelta,fdeltaMin,fdeltaMax,fepsilonThreshold)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function finds LS (sparsity enforces) solution, thresholds solution's 12
% norm with epsilonThreshold, and then nulls those rows of the solution X
% which were thresholded in 12norm
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% INPUT: A - cell(1,N) that contains mixing matrices A's
%        Y - matrix of Poisson counts
%        XsolInitial - initial matrix X solution: one initial X solution
%        for all algorithms
%        epsilonLSSP - epsilon for LSSP method
%        Niter - number of iterations for iterative algorithm
%        NiterDelta - number of iteration to find optimal Delta
%        deltaMin - minimal value of delta range
%        deltaMax - maximal value of delta range
%        epsilonThreshold - epsilon for thresholding 12 norm of the
%        solution matrix XsolLSSP
% OUTPUT: XsolLSSP - solution matrix
%         XsolLSSP12norm - 12norm of the solution matrix XsolLSSP, saved
%         for every iterDelta
%         objFuncLSSP - objective function that is minimized, save for every
%         iterDelta and every iretation iter, when solution matrix gets
%         updated
%         constrainLSSP - values of the constains saved for every iterDelta
%         delta - Lagrangian parameter saved for each iterDelta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% NB!!! These parameters are defined/changed in this code:
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% number of internal iterations
fNiterI = 1000;
% max number of internal iterations for while loop for internal iterations
fNiterIMax = 100;
% max number of times delta range can be increased/decreased
% fstepCntrMax = 20;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[~,fK] = size(fA{1});
[~,fN] = size(fY);

% need to check that if zero matrix is inside the constrains then solution
% must be automatically defined as zero matrix
% if sum_{i=1,N}||y_i||_2^2 <= epsilonLSSP
fconstrHlp = sum(sum(fY.^2));
if fconstrHlp <= fepsilonLSSP
    % it means zero point inside espilon ball, hence solution is zero
    % matrix
    fXsolLSSP = zeros(fK,fN);
    % 12 norm is a zero vector
    fXsolLSSP12norm = zeros(fNiterDelta,fK);
    % null objective function
    fobjFuncLSSP = zeros(1,fNiterE);
    % constraints are const sum_{i=1,N}||y_i||_2^2
    fconstrainLSSP = fconstrHlp*ones(fNiterDelta,1);
    % delta can be any but here I make it zeros
    fdelta = zeros(fNiterDelta,1);
%     disp('zero point is inside the constrained region - therefore solution is zero matrix!');
%     disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
    
else % if epsilon ball doesn't contain zero point then run algorithm
    
    % number of internal iterations
    % fNiterI = 20;
    % initialize objective function calculated for each delta, each iteration
    fobjFuncLSSP = zeros(fNiterDelta,fNiterE);
    % initialize array that contains constrains for each delta
    fconstrainLSSP = zeros(fNiterDelta,1);
    % initialize array that contains values of 12 norm of XsolLSSP for each
    % iterDelta
    fXsolLSSP12norm = zeros(fNiterDelta,fK);
    % initialize array that contains values of delta for each iterDelta
    fdelta = zeros(fNiterDelta,1);
    % matrix of N*1 basis vectors
    En = eye(fN,fN);
    
    %calculate auxiliaty matrices M_l, are not changing with iterations
    fHlp1 = zeros(fK,fN); % auxiliary matrix for finding M_l
    for fiN = 1:fN
        % each column of this matrix contains sum of the columns of
        % A'A+diag(A'A)
        fHlp1(:,fiN) = sum(fA{fiN}'*fA{fiN},2) + diag(fA{fiN}'*fA{fiN});
    end;
    % cell contains matrices M_l
    fM = cell(1,fK);
    for fiK = 1:fK
        fM{fiK} = diag(fHlp1(fiK,:));
    end;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %                     Find range for Delta
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     disp('**********************************');
%     disp('      Find range for Delta');
%     disp('**********************************');
    % flag that indicates if correct delta range was found or not
    fdeltaRangeFoundFlag = 0;
    
    % counter of times delta range was increased/decreased
    fstepCntr = 0;
    % max number of times delta range can be increased/decreased
    fstepCntrMax = 20;
    
    % check that for max delta constraints are less then epsilon, for min delta
    % constraints are greater then epsilon, if not decrease/increase delta
    % range, do it no more then fstepCntrMax times
    
    while fdeltaRangeFoundFlag == 0 && (fstepCntr<=fstepCntrMax)
        
%         disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
         fstepCntr = fstepCntr+1;
         disp(strcat('current Delta Ranging step = ',num2str(fstepCntr)));
%         disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
%         disp(strcat('current fDeltaMax = ',num2str(fdeltaMax)))
%         disp(strcat('current fDeltaMin = ',num2str(fdeltaMin)))
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Check 1: check that for max delta constraints are less then epsilon
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % initialize Xsol
        fXsolLSSP = rand(fK,fN);
        % set current delta to fdeltaMax
        fdeltaHlp = fdeltaMax;
        
        % find solution corresponding for fdeltaMax
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % external iterations
        for fiterE = 1:fNiterE
            
            % matrix columns of which are v_l in row decomposition of constrain
            % function, v_l depend on previous value of solution matrix in
            % external loop, i.e. fXsolLSSP
            fV = zeros(fN,fK);
            for fiN = 1:fN
                fV(fiN,:) = (fY(:,fiN)-fA{fiN}*fXsolLSSP(:,fiN))'*fA{fiN};
            end;
            
            % calculate constant term for internal objective function, this
            % term depends on previous value of solution matrix in external
            % loop, i.e. fXsolLSSP
            fC = 0;
            for fiN = 1:fN
                fC = fC + (fY(:,fiN)-fA{fiN}*fXsolLSSP(:,fiN))'*(fY(:,fiN)-fA{fiN}*fXsolLSSP(:,fiN));
            end;
            
            % initialize internal solution matrix fXsolLSSPin
            fXsolLSSPin = rand(fK,fN);
            
            % internal iterations
            for fiterI = 1:fNiterI
                % vector that contains Eucledian norms of the rows of the
                % solution matrix on previous internal iteration
                fprevIntNorm = sqrt(sum(fXsolLSSPin.^2,2));
                
                % update each row separatly
                for fiK = 1:fK
                    % auxiliary matrix
                    fH = fprevIntNorm(fiK)*fM{fiK}+En./(2*fdeltaHlp);
                    % update of the row
                    fXsolLSSPin(fiK,:) = fprevIntNorm(fiK)*((fH)\(fV(:,fiK)+fM{fiK}*fXsolLSSP(fiK,:)'));
                    % make row nonegative
                    fXsolLSSPin(fiK,:) = fXsolLSSPin(fiK,:).*(fXsolLSSPin(fiK,:)>0);
                end;
            end; % internal iterations
            
            % update of external solution
            fXsolLSSP = fXsolLSSPin;
            
            % calculate the constrains for current deltaMax
            fconstrainHlpMax = 0;
            for fiN = 1:fN
                fconstrainHlpMax = fconstrainHlpMax + norm(fA{fiN}*fXsolLSSP(:,fiN)-fY(:,fiN),2)^2;
            end;
        end; % external iterations
        
        % check that constraints for deltaMax are less then epsilon
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if fconstrainHlpMax > fepsilonLSSP
%             disp('constr for max delta > eps');
            fconstrDeltaMaxFlag = 0;
        else
            fconstrDeltaMaxFlag = 1;
%             disp('constr for max delta < eps');
        end;
        %~~~~~~~~~~~~~~~~~~~~~~% MAX %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Check 2: check that for min delta constraints are greater then epsilon
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % initialize Xsol
        fXsolLSSP = rand(fK,fN);
        % set current delta to min
        fdeltaHlp = fdeltaMin;
        
        % find solution corresponding for fdeltaMin
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % external iterations
        for fiterE = 1:fNiterE
            
            % matrix columns of which are v_l in row decomposition of constrain
            % function, v_l depend on previous value of solution matrix in
            % external loop, i.e. fXsolLSSP
            fV = zeros(fN,fK);
            for fiN = 1:fN
                fV(fiN,:) = (fY(:,fiN)-fA{fiN}*fXsolLSSP(:,fiN))'*fA{fiN};
            end;
            
            % calculate constant term for internal objective function, this
            % term depends on previous value of solution matrix in external
            % loop, i.e. fXsolLSSP
            fC = 0;
            for fiN = 1:fN
                fC = fC + (fY(:,fiN)-fA{fiN}*fXsolLSSP(:,fiN))'*(fY(:,fiN)-fA{fiN}*fXsolLSSP(:,fiN));
            end;
            
            % initialize internal solution matrix fXsolLSSPin
            fXsolLSSPin = rand(fK,fN);
            
            % internal iterations
            for fiterI = 1:fNiterI
                % vector that contains Eucledian norms of the rows of the
                % solution matrix on previous internal iteration
                fprevIntNorm = sqrt(sum(fXsolLSSPin.^2,2));
                
                % update each row separatly
                for fiK = 1:fK
                    % auxiliary matrix
                    fH = fprevIntNorm(fiK)*fM{fiK}+En./(2*fdeltaHlp);
                    % update of the row
                    fXsolLSSPin(fiK,:) = fprevIntNorm(fiK)*((fH)\(fV(:,fiK)+fM{fiK}*fXsolLSSP(fiK,:)'));
                    % make row nonegative
                    fXsolLSSPin(fiK,:) = fXsolLSSPin(fiK,:).*(fXsolLSSPin(fiK,:)>0);
                end;
            end; % internal iterations
            
            % update of external solution
            fXsolLSSP = fXsolLSSPin;
            
            % calculate the constrains for current deltaMax
            fconstrainHlpMin = 0;
            for fiN = 1:fN
                fconstrainHlpMin = fconstrainHlpMin + norm(fA{fiN}*fXsolLSSP(:,fiN)-fY(:,fiN),2)^2;
            end;
        end; % external iterations
        
        % check that constraints for deltaMin are greater then epsilon
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if fconstrainHlpMin  > fepsilonLSSP
%             disp('constr for min delta > eps');
            fconstrDeltaMinFlag = 1;
        else
%             disp('constr for min delta < eps');
            fconstrDeltaMinFlag = 0;
        end;
        %~~~~~~~~~~~~~~~~~~~~~~% MIN %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        % check that constraints for deltaMax are less then epsilon and
        % constraints for deltaMin are greater then epsilon
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % if not - then change the constraints range and check again
        if fconstrDeltaMaxFlag && fconstrDeltaMinFlag
            % the correct range for delta is found - exit while loop
            fdeltaRangeFoundFlag = 1;
        else
            % if deltaMax is too small - increase deltaMax
            if fconstrDeltaMaxFlag == 0
                % if constraints for deltaMin are greater then epsilon and
                % simultaneously constraints for deltaMax are greater then
                % epsilon then move deltaMin to deltaMax position
                if fconstrDeltaMinFlag
                    fdeltaMin = fdeltaMax;
                end;
                fdeltaMax = fdeltaMax*1e1;
            end;
            % if deltaMin is too big - decrease deltaMin
            if fconstrDeltaMinFlag == 0
                % if constraints for deltaMin is less then epsilon and
                % simultaneously constraints for deltaMax is less then
                % epsilon then move deltaMax to deltaMin position
                if fconstrDeltaMaxFlag
                    fdeltaMax = fdeltaMin;
                end;
                fdeltaMin = fdeltaMin*1e-1;
            end;
        end;
    end; % while
    
    % display DeltaMin and DeltaMax that were found
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     disp('**********************************');
%     disp(strcat('fDeltaMax = ',num2str(fdeltaMax)))
%     disp(strcat('fDeltaMin = ',num2str(fdeltaMin)))
%     if fdeltaRangeFoundFlag
%         disp('**********************************');
%         disp('      Range for delta is found ');
%         disp('**********************************');
%     else
%         disp('**********************************');
%         disp('Exit delta range search loop (max number of steps)');
%         disp('**********************************');
%     end;
    %~~~~~~~~~~~~~~~~~~~~~~~ Find range for Delta ~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %                     Find solution matrix
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % initialize Xsol
    fXsolLSSP = fXsolInitial;
    
    % delta iterations
    for fiDelta = 1:fNiterDelta
        % set Lagrangian parameter delta
        fdelta(fiDelta)=(fdeltaMax+fdeltaMin)/2;
        
        % NB!!! initialize Xsol randomly for every external iteration or
        % comment this line to start new external iteration from the previous
        % value of Xsol
        % fXsolLSSP = rand(fK,fN);
        
        % array where number of internal iterations is saved
        fInternalIterN = zeros(1,fNiterE);
        
        % external iterations
        for fiterE = 1:fNiterE
            
            % matrix columns of which are v_l in row decomposition of constrain
            % function, v_l depend on previous value of solution matrix in
            % external loop, i.e. fXsolLSSP
            fV = zeros(fN,fK);
            for fiN = 1:fN
                fV(fiN,:) = (fY(:,fiN)-fA{fiN}*fXsolLSSP(:,fiN))'*fA{fiN};
            end;
            
            % calculate constant term for internal objective function, this
            % term depends on previous value of solution matrix in external
            % loop, i.e. fXsolLSSP
            fC = 0;
            for fiN = 1:fN
                fC = fC + (fY(:,fiN)-fA{fiN}*fXsolLSSP(:,fiN))'*(fY(:,fiN)-fA{fiN}*fXsolLSSP(:,fiN));
            end;
            
            % initialize internal solution matrix fXsolLSSPin
            %fXsolLSSPin = rand(fK,fN);
            fXsolLSSPin = fXsolLSSP+rand(fK,fN);
            
            % initialize internal objective function for each each internal iteration
            % fNiterIMax = 30;
            fobjInternal = zeros(1,fNiterIMax);
            
            % flag for while loop for internal iterations
            fInternalLoopExit = 0;
            % counter of internal iterations
            fInternalIterCntr = 0;
            
            while fInternalLoopExit == 0 && fInternalIterCntr<fNiterIMax
                % update internal iteration counter
                fInternalIterCntr = fInternalIterCntr + 1;
                % vector that contains Eucledian norms of the rows of the
                % solution matrix on previous internal iteration
                fprevIntNorm = sqrt(sum(fXsolLSSPin.^2,2));
                % update each row separatly
                for fiK = 1:fK
                    % auxiliary matrix
                    fH = fprevIntNorm(fiK)*fM{fiK}+En./(2*fdelta(fiDelta));
                    % update of the row
                    fXsolLSSPin(fiK,:) = fprevIntNorm(fiK)*((fH)\(fV(:,fiK)+fM{fiK}*fXsolLSSP(fiK,:)'));
                    % make row nonegative
                    fXsolLSSPin(fiK,:) = fXsolLSSPin(fiK,:).*(fXsolLSSPin(fiK,:)>0);
                    % save values of internal objective function
                    fobjInternal(fInternalIterCntr) = fobjInternal(fInternalIterCntr) + fdelta(fiDelta)*((fXsolLSSPin(fiK,:)-fXsolLSSP(fiK,:))*fM{fiK}*(fXsolLSSPin(fiK,:)-fXsolLSSP(fiK,:))'-2*(fXsolLSSPin(fiK,:)-fXsolLSSP(fiK,:))*fV(:,fiK)) + norm(fXsolLSSPin(fiK,:),2);
                end;
                fobjInternal(fInternalIterCntr) = fobjInternal(fInternalIterCntr) + fdelta(fiDelta)*fC;
                % calculate current value of the external objective function for current internal solution
                fconstrainLSSPhlp = 0;
                for fiN = 1:fN
                    % calculate the constrains
                    fconstrainLSSPhlp = fconstrainLSSPhlp + norm(fA{fiN}*fXsolLSSPin(:,fiN)-fY(:,fiN),2)^2;
                end;
                fobjFuncLSSPhlp = sum(sqrt(sum(fXsolLSSPin.*fXsolLSSPin,2)))+ fdelta(fiDelta)*(fconstrainLSSPhlp - fepsilonLSSP);
                % compare fobjFuncLSSPhlp with the objective on the last
                % external iteration. NB!!! objective function should
                % decrease!!!
                if fiterE > 1
                    if fobjFuncLSSP(fiDelta,fiterE-1) >= fobjFuncLSSPhlp
                        % means that obj function decreased
                        fInternalLoopExit = 1;
                    else
                        % objective function didn't decrese
                        fInternalLoopExit = 0;
                    end;
                end;
            end; % while loop for internal iterations
            
            % display internal objective function for current
            % delta nd current external iteration
%             figure(999)
%             plot(fobjInternal);
%             grid on
%             xlabel('internal iterations')
%             ylabel('internal objective function')
%             title('internal objective function');
           
                             
            
            % save number of internal iterations needed to decrease the
            % external objective function
            fInternalIterN(fiterE) = fInternalIterCntr;
            
            %             % internal iterations
            %             for fiterI = 1:fNiterI
            %                 % vector that contains Eucledian norms of the rows of the
            %                 % solution matrix on previous internal iteration
            %                 fprevIntNorm = sqrt(sum(fXsolLSSPin.^2,2));
            %
            %                 % update each row separatly
            %                 for fiK = 1:fK
            %                     % auxiliary matrix
            %
            %                     fH = fprevIntNorm(fiK)*fM{fiK}+En./(2*fdelta(fiDelta));
            %
            %                     % update of the row
            %                     fXsolLSSPin(fiK,:) = fprevIntNorm(fiK)*((fH)\(fV(:,fiK)+fM{fiK}*fXsolLSSP(fiK,:)'));
            %                     % make row nonegative
            %                     fXsolLSSPin(fiK,:) = fXsolLSSPin(fiK,:).*(fXsolLSSPin(fiK,:)>0);
            %
            %                     % save values of internal objective function
            %                     fobjInternal(fiK,fiterI) = fdelta(fiDelta)*((fXsolLSSPin(fiK,:)-fXsolLSSP(fiK,:))*fM{fiK}*(fXsolLSSPin(fiK,:)-fXsolLSSP(fiK,:))'-2*(fXsolLSSPin(fiK,:)-fXsolLSSP(fiK,:))*fV(:,fiK) + fC) + norm(fXsolLSSPin(fiK,:),2);
            %                 end;
            %             end; % internal iterations
            
            % update of external solution
            fXsolLSSP = fXsolLSSPin;
            
            fconstrainLSSPhlp = 0;
            for fiN = 1:fN
                % calculate the constrains
                fconstrainLSSPhlp = fconstrainLSSPhlp + norm(fA{fiN}*fXsolLSSP(:,fiN)-fY(:,fiN),2)^2;
            end;
            
            % external objective function for current iter and current iDelta, i.e. Lagrangian
            fobjFuncLSSP(fiDelta,fiterE) = sum(sqrt(sum(fXsolLSSP.*fXsolLSSP,2)))+ fdelta(fiDelta)*(fconstrainLSSPhlp - fepsilonLSSP);
        end; % external iterations
        
        % display external objective function for current delta
%         fig = figure('Name','external objective function');
%         plot(fobjFuncLSSP(fiDelta,:));
%         grid on
%         xlabel('external iteration')
%         ylabel('objective function')
%         title('external objective function');
%         hgsave(strcat('ext obj funct for delta iter_',num2str(fiDelta)));
%         saveas(fig,strcat('ext obj funct for delta iter_',num2str(fiDelta)),'jpg');
%         close(fig);
        
        
        % calculate the constraint for current iteration delta (i.e. on the last external iteration)
        fconstrainLSSP(fiDelta) = fconstrainLSSPhlp;
        
        % calculate new range for delta
        if fconstrainLSSP(fiDelta) > fepsilonLSSP
            % then increase delta
            fdeltaMin = fdelta(fiDelta);
            if mod(fiDelta,5) == 0
                disp(['constr > epsilon @ iter = ',num2str(fiDelta)]);
                disp(strcat('deltaMINnow =  ',num2str(fdeltaMin)));
                disp(strcat('deltaMAXwnow =  ',num2str(fdeltaMax)));
                disp('   ');
            end
        else
            % then decrease delta
            fdeltaMax = fdelta(fiDelta);
            if mod(fiDelta,5) == 0
                disp(['constr < epsilon @ iter = ',num2str(fiDelta)]);
                disp(strcat('deltaMINnow =  ',num2str(fdeltaMin)));
                disp(strcat('deltaMAXwnow =  ',num2str(fdeltaMax)));
                disp('   ');
            end
        end;
        
        % save value of 12 norm of the solution for current value of Delta
        %fXsolLSSP12norm(fiDelta,:) = sqrt(sum(fXsolLSSP.*fXsolLSSP,2));
        
        % display number of internal iterations needed to decrease the
        % objective
        %disp('   ');
        %disp('Number of internal iterations needed to decrease the objective for current delta');
        %disp('*****************************************************');
        %disp(fInternalIterN');
        
    end; % for iDelta
    
%     % threshold the 12 norm of the solution matrix XsolLSSP
%     findThr = fXsolLSSP12norm(fiDelta,:)<=fepsilonThreshold;
%     fXsolLSSP12norm(fiDelta,findThr) = 0;
%     
%     % null those rows of the solution matrix XsolLSSP which were thresholded by
%     % the epsilonThreshold in XsolLSSP12norm (we need this for MSE analysis)
%     fXsolLSSP(findThr==1,:) = 0;
    
end; %{if zero point is in epsilon ball}