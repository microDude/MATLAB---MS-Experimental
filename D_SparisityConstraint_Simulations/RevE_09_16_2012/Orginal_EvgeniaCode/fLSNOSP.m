function [fXsolLSNOSP,fXsolLSNOSP12norm,fobjFuncLSNOSP] = fLSNOSP(fA,fY,fXsolInitial,fNiterLSNOSP,fepsilonThreshold)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function finds LS (no sparsity enforces) solution
% using separation of variables trick and thresholds its 12 norm of the
% rows using epsilonThreshold, then nulls those rows of the solution matrix
% that were thresholded in the 12norm
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% INPUT: A - cell(1,N) that contains mixing matrices A's
%        Y - matrix of Poisson counts
%        XsolInitial - initial matrix X solution: one initial X solution
%        for all algorithms
%        NiterLSNOSP - number of iteration for findinf solution matrix X
%        epsilonThreshold - epsilon for thresholding solution matrix
%        XsolLSNOSP
% OUTPUT: XsolLSNOSP - solution matrix 
%         XsolLSNOSP12norm - 12norm of the solution matrix XsolMLNOSP
%         objFuncLSNOSP - objective function that is minimized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fK,~] = size(fXsolInitial);
[~,fN] = size(fY);
% initialize XsolLSNOSP
fXsolLSNOSP = fXsolInitial;
% initialize objective function
fobjFuncLSNOSP = zeros(fNiterLSNOSP,1);
% vector that contains max eig of A{i}'*A{i}
fmaxeigA = zeros(fN,1);
for fiN = 1:fN
    fmaxeigA(fiN) = max(eig(fA{fiN}'*fA{fiN}));
end;

% find solution iteratively
for fiter = 1:fNiterLSNOSP
    % auxiliary variable for objective function calculation
    fobjFuncLSNOSPhlp = 0;
    
    for fiN = 1:fN % update each element of Xsol separately
        % auxiliary vector
        
        % $f_h = \mathbf{A}_i^T \mathbf{A}_i \mathbf{x}_i - \mathbf{A}_i^T \mathbf{y}$
        fh = fA{fiN}'*fA{fiN}*fXsolLSNOSP(:,fiN) - fA{fiN}'*fY(:,fiN);
        for fiK = 1:fK
            % update each element separately, ()_+
            fXsolLSNOSP(fiK,fiN) = max(0,(fXsolLSNOSP(fiK,fiN) - (1/fmaxeigA(fiN))*fh(fiK)));
        end;
        fobjFuncLSNOSPhlp = fobjFuncLSNOSPhlp + norm(fA{fiN}*fXsolLSNOSP(:,fiN)-fY(:,fiN),2)^2;
    end;
    
    % value of objective function on current iteration
    fobjFuncLSNOSP(fiter) = fobjFuncLSNOSPhlp;
end;

% calculate 12norm of the solution matrix XsolLLNOSP
fXsolLSNOSP12norm = sqrt(sum(fXsolLSNOSP.*fXsolLSNOSP,2));

% threshold the 12 norm of the solution matrix XsolLSNOSP
findThr = fXsolLSNOSP12norm<=fepsilonThreshold;
fXsolLSNOSP12norm(findThr) = 0;

% null those rows of the solution matrix XsolLSNOSP which were thresholded by
% the epsilonThreshold in XsolLSNOSP12norm (we need this for MSE analysis)
fXsolLSNOSP(findThr==1,:) = 0;

