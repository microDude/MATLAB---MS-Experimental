function [fXsolLSNOSPsupFull,fXsolLSNOSPsupFull12norm,fobjFuncLSNOSPsup] = fLSNOSPonSupport(fA,fY,fNiterLSNOSPsup,fepsilonThreshold,fSupport)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function finds LS (no sparsity enforces) solution on support
% using separation of variables trick and thresholds its 12 norm of the rows
% using epsilonThreshold, then nulls those rows of X that were thresholded
% in 12norm
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% INPUT: A - cell(1,N) that contains mixing matrices A's
%        Y - matrix of Poisson counts
%        NiterLSNOSPsup - number of iteration to find solution matrix X
%        epsilonThreshold - epsilon for thresholding solution matrix
%        XsolLSNOSPsup
%        Support - support of X, i.e. indices of nonzero positions of 2
%        norm of the rows of X
% OUTPUT: XsolLSNOSPsupFull - full solution matrix
%         fXsolLSNOSPsupFull12norm - full 12norm of the solution matrix
%         objFuncLSNOSPsup - objective function that is minimized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get values of the dimensions N, K from the data
[~,fN] = size(fY);
[~,fK] = size(fA{1});

% if support is empty then set solution matrix to zero matrix, norm to
% zeros and objective function to zeros
if sum(fSupport) == 0
    disp('   ');
    disp('Support is empty');
    disp('~~~~~~~~~~~~~~~~');
    % solution matrix
    fXsolLSNOSPsupFull = zeros(fK,fN);
    % 12 norm of the solution matrix
    fXsolLSNOSPsupFull12norm = zeros(fK,1);
    % objective function
    fobjFuncLSNOSPsup = zeros(fNiterLSNOSPsup,1);
    
else % otherwise run the algorithm
    
    % indices of support
    fSupportIndex = find(fSupport == 1);
    % length of the support
    fS = length(fSupportIndex);
    
    % initialize XsolLSNOSPsup on support
    fXsolLSNOSPsup = rand(fS,fN);
    % initialize full solution matrix
    fXsolLSNOSPsupFull = zeros(fK,fN);
    % initialize objective function
    fobjFuncLSNOSPsup = zeros(fNiterLSNOSPsup,1);
    % vector that contains max eig of A{i}'*A{i} on support
    fmaxeigA = zeros(fN,1);
    for fiN = 1:fN
        fA{fiN} = fA{fiN}(:,fSupportIndex);
        fmaxeigA(fiN) = max(eig(fA{fiN}'*fA{fiN}));
    end;
    
    % find solution on support iteratively
    for fiter = 1:fNiterLSNOSPsup
        % auxiliary variable for objective function calculation
        fobjFuncLSNOSPhlp = 0;
        
        for fiN = 1:fN % update each element of Xsol separately
            % auxiliary vector
            fh = fA{fiN}'*fA{fiN}*fXsolLSNOSPsup(:,fiN) - fA{fiN}'*fY(:,fiN);
            for fiS = 1:fS
                % update each element separately, ()_+
                fXsolLSNOSPsup(fiS,fiN) = max(0,(fXsolLSNOSPsup(fiS,fiN) - (1/fmaxeigA(fiN))*fh(fiS)));
            end;
            fobjFuncLSNOSPhlp = fobjFuncLSNOSPhlp + norm(fA{fiN}*fXsolLSNOSPsup(:,fiN)-fY(:,fiN),2)^2;
        end;
        % value of objective function on current iteration
        fobjFuncLSNOSPsup(fiter) = fobjFuncLSNOSPhlp;
    end;
    
    % calculate 12norm of the solution matrix XsolLLNOSPsup
    fXsolLSNOSP12norm = sqrt(sum(fXsolLSNOSPsup.*fXsolLSNOSPsup,2));
    
    % threshold the 12 norm of the solution matrix XsolLSNOSPsup
    findThr = fXsolLSNOSP12norm<=fepsilonThreshold;
    fXsolLSNOSP12norm(findThr) = 0;
    
    % null those rows of the solution matrix XsolLSNOSPsup which were
    % thresholded by the epsilonThreshold in XsolLSNOSP12norm (we need this
    % for MSE analysis)
    fXsolLSNOSPsup(findThr==1,:) = 0;
    
    % update full soluion matrix XsolLSNOSPsupFull
    fXsolLSNOSPsupFull(fSupportIndex,:) = fXsolLSNOSPsup;
    
    % update full 12 norm of the solution
    fXsolLSNOSPsupFull12norm = zeros(fK,1);
    fXsolLSNOSPsupFull12norm(fSupportIndex) = fXsolLSNOSP12norm;
    
end; %{if support is empty}