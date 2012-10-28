function [fXsolMLNOSPsupFull,fXsolMLNOSPsupFull12norm,fobjFuncMLNOSPsup] = fMLNOSPonSupport(fA,fArowsum,fY,fNiterMLNOSPsup,fepsilonThreshold,fSupport)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function finds ML (no sparsity enforces) solution on support
% and thresholds its 12 norm of the rows using epsilonThreshold, then nulls
% those rows of X that were thresholded in 12norm
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% INPUT: A - cell(1,N) that contains mixing matrices A's
%        Arowsum - matrix, ith row of which is sum of rows of A{i}
%        Y - matrix of Poisson counts
%        NiterMLNOSPsup - number of iterations for this MLNOSP on support
%        algorithm
%        epsilonThreshold - epsilon for thresholding 12 norm of the
%        solution mtrix XsolMLNOSP
%        Support - support of X, i.e. indices of nonzero positions of 2
%        norm of the rows of X
% OUTPUT: XsolMLNOSPsupFull - full solution matrix
%         XsolMLNOSPsupFull12norm - 12norm of the full solution matrix
%         objFuncMLNOSPsup - objective function that is minimized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fM,fN] = size(fY);
[~,fK] = size(fA{1});

% if support is empty then set solution matrix to zero matrix, norm to
% zeros and objective function to zeros
if sum(fSupport) == 0
    disp('   ');
    disp('Support is empty');
    disp('~~~~~~~~~~~~~~~~');
    % solution matrix
    fXsolMLNOSPsupFull = zeros(fK,fN);
    % 12 norm of the solution matrix
    fXsolMLNOSPsupFull12norm = zeros(fK,1);
    % objective function
    fobjFuncMLNOSPsup = zeros(fNiterMLNOSPsup,1);
    
else % otherwise run the algorithm
    
    % dimension of the support,ie. row-sparsity of X found by some other method
    fsupind = fSupport~=0;
    fs = sum(fsupind);
    % initialize Xsol
    fXsolMLNOSPsupport = rand(fs,fN);
    % initialize objective function
    fobjFuncMLNOSPsup = zeros(fNiterMLNOSPsup,1);
    
    % find Xsol iteratively
    for fiter = 1:fNiterMLNOSPsup
        
        % update each column of Xsol separately
        for fiN = 1:fN
            
            % indices where y_i(j)~=0
            fiHlp5 = find(fY(:,fiN)~=0);
            % A_i*x_i term
            fHlp6 = fA{fiN}(:,fsupind)*fXsolMLNOSPsupport(:,fiN);
            fHlp1 = zeros(fM,1);
            fHlp1(fiHlp5) = fY(fiHlp5,fiN)./(fHlp6(fiHlp5));
            fHlp2 = fHlp1'*fA{fiN}(:,fsupind);
            
            % update column fiN of the solution matrix fXsolMLNOSPsupport
            fXsolMLNOSPsupport(:,fiN) = (fHlp2.*(1./fArowsum(fiN,fsupind))).*(fXsolMLNOSPsupport(:,fiN)');
            % projection of column on positive semi-axis, since Problem Domain is R^{K*N}_{+}
            fXsolMLNOSPsupportfiNPind = (fXsolMLNOSPsupport(:,fiN)<0);
            fXsolMLNOSPsupport(fXsolMLNOSPsupportfiNPind,fiN) = 0;
            
            % objective function update for current iteration
            % A_i*x_i term
            fHlp6 = fA{fiN}(:,fsupind)*fXsolMLNOSPsupport(:,fiN);
            % \sum_j y_i(j)*log((A_i*x_i)_j) for non-zero y_i(j) term
            fHlp7 = fY(fiHlp5,fiN)'*log(fHlp6(fiHlp5));
            % \sum_j y_i(j)*log(y_i(j)) for non-zero y_i(j) term
            fHlp8 = fY(fiHlp5,fiN)'*log(fY(fiHlp5,fiN));
            % \sum_j y_i(j) term
            fHlp9 = sum(fY(:,fiN));
            % save value of objective function for current iteration
            fobjFuncMLNOSPsup(fiter) = fobjFuncMLNOSPsup(fiter) + fArowsum(fiN,fsupind)*fXsolMLNOSPsupport(:,fiN) - fHlp7 + fHlp8 - fHlp9;
        end;
    end;
    
    % calculate 12norm of the solution matrix XsolMLNOSP
    fXsolMLNOSPsup12norm = sqrt(sum(fXsolMLNOSPsupport.*fXsolMLNOSPsupport,2));
    % threshold the 12 norm of the solution matrix XsolMLNOSP
    findThr = fXsolMLNOSPsup12norm<=fepsilonThreshold;
    fXsolMLNOSPsup12norm(findThr) = 0;
    
    % full solution matrix with non-zero rows on support and zero rows
    % elsewhere, with null rows that were thresholded in 12 norm
    fXsolMLNOSPsupFull = zeros(fK,fN);
    fXsolMLNOSPsupFull(fsupind,:) = fXsolMLNOSPsupport;
    fXsolMLNOSPsupFull(findThr == 1,:) = 0;
    
    % 12 norm of the full solution
    fXsolMLNOSPsupFull12norm = zeros(fK,1);
    fXsolMLNOSPsupFull12norm(fsupind,1) = fXsolMLNOSPsup12norm;
    
end; %{if support is empty}
