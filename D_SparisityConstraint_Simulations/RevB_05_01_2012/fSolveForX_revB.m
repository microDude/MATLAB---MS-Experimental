function [ Xsol_est ] = fSolveForX_revB(A,y,lambda,iterMax,Xsol)
% fSolveForX_revA(f_k,f_n,f_A,f_Y,lambda,iterMax,f_epsThresholding)
% Returns a iteration run for X-estimated for a given lambda value of the
% Langragian multiplier.

% Find dimensions of y
[m,n] = size(y);   % Length of Spectra x Number of spectra per experiment

% random initialization of Xsol
%Xsol = 0+(1-0).*rand(m,n);

for iter = 1:iterMax
   
    f_D = diag(sqrt(sum(Xsol.^2,2)));
       
    f_Dhat = sqrt(f_D);

    % Loop through each captured spectra per experiment
    for f_i = 1:n
        y_i = y(:,f_i);
        
        Xsol(:,f_i) = f_Dhat*( (lambda*f_Dhat*(A'*A)*f_Dhat+eye(m))\(lambda*f_Dhat*(A'*y_i)) );
    end
   
%     if (mod(iter,40) == 0)
%         plot(sum(Xsol,2));title(['iter = ',num2str(iter)]);
%         pause(0.1);
%     end
    %iter
end

% tresholding the solution, get rid of all negatives
f_ind = (sum(Xsol,2) <= 2*eps);
Xsol(f_ind,:) = 0;

% return the result
Xsol_est = Xsol;

end
