function [fXsolMLNOSP,fXsolMLNOSP12norm,fobjFuncMLNOSP] = fMLNOSP(fA,fArowsum,fY,fXsolInitial,fNiterMLNOSP,fepsilonThreshold)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function finds ML (no sparsity enforces) solution and thresholds 12 norm
% with epsilonThreshold, then nulls those rows of X that were thresholded 
% in 12norm
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% INPUT: A - cell(1,N) that contains mixing matrices A's
%        Arowsum - matrix, ith row of which is sum of rows of A{i}
%        Y - matrix of Poisson counts
%        XsolInitial - initial matrix X solution: one initial X solution
%        for all algorithms
%        NiterMLNOSP - number of iterations for this MLNOSP algorithm 
%        epsilonThreshold - epsilon for thresholding 12 norm of the 
%        solution mtrix XsolMLNOSP
% OUTPUT: XsolMLNOSP - solution matrix 
%         XsolMLNOSP12norm - 12norm of the solution matrix XsolMLNOSP
%         objFuncMLNOSP - objective function that is minimized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fM,fN] = size(fY);
% initialize objective function
fobjFuncMLNOSP = zeros(fNiterMLNOSP,1);
% initialize Xsol
fXsolMLNOSP = fXsolInitial;

% find Xsol iteratively
for fiter = 1:fNiterMLNOSP
    
    for fiN = 1:fN % update each column of Xsol separately
        
        % indices where y_i(j)~=0 
        fiHlp5 = find(fY(:,fiN)~=0);
        % A_i*x_i term
        fHlp6 = fA{fiN}*fXsolMLNOSP(:,fiN);
        fHlp1 = zeros(fM,1);
        fHlp1(fiHlp5) = fY(fiHlp5,fiN)./(fHlp6(fiHlp5));
        fHlp2 = fHlp1'*fA{fiN};
        
        % update column fiN of the solution matrix             
        fXsolMLNOSP(:,fiN) = (fHlp2.*(1./fArowsum(fiN,:))).*(fXsolMLNOSP(:,fiN)');
        % projection of column on positive semi-axis, since Problem Domain is R^{K*N}_{+}
        fXsolMLNOSfiNPind = (fXsolMLNOSP(:,fiN)<0);
        fXsolMLNOSP(fXsolMLNOSfiNPind,fiN) = 0;
        
        % objective function update for current iteration
        % A_i*x_i term
        fHlp6 = fA{fiN}*fXsolMLNOSP(:,fiN);
        % \sum_j y_i(j)*log((A_i*x_i)_j) for non-zero y_i(j) term
        fHlp7 = fY(fiHlp5,fiN)'*log(fHlp6(fiHlp5));
        % \sum_j y_i(j)*log(y_i(j)) for non-zero y_i(j) term
        fHlp8 = fY(fiHlp5,fiN)'*log(fY(fiHlp5,fiN));
        % \sum_j y_i(j) term
        fHlp9 = sum(fY(:,fiN));
        % save value of objective function for current iteration
        fobjFuncMLNOSP(fiter) = fobjFuncMLNOSP(fiter) + fArowsum(fiN,:)*fXsolMLNOSP(:,fiN) - fHlp7 + fHlp8 - fHlp9;
    end;
end;

% calculate 12norm of the solution matrix XsolMLNOSP
fXsolMLNOSP12norm = sqrt(sum(fXsolMLNOSP.*fXsolMLNOSP,2));

% threshold 12 norm of the solution matrix XsolMLNOSP
findThr = fXsolMLNOSP12norm<=fepsilonThreshold;
fXsolMLNOSP12norm(findThr) = 0;

% null those rows of the solution matrix XsolMLNOSP which were thresholded by
% the epsilonThreshold in XsolMLNOSP12norm (we need this for MSE analysis)
fXsolMLNOSP(findThr==1,:) = 0;

