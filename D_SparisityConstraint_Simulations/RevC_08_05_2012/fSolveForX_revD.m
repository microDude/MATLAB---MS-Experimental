function [ Xsol_est ] = fSolveForX_revC(A,y,lambda,iterMax,Xsol,debugFlag)
% fSolveForX_revA(f_k,f_n,f_A,f_Y,lambda,iterMax,f_epsThresholding)
% Returns a iteration run for X-estimated for a given lambda value of the
% Langragian multiplier.

% Find dimensions of y
[m,n] = size(y);   % Length of Spectra x Number of spectra per experiment

% random initialization of Xsol
%Xsol = 0+(1-0).*rand(m,n);

for iter = 1:iterMax
   
    f_D = sqrt(sum(Xsol.^2,2));
       
    f_Dhat = diag(sqrt(f_D));

    % Loop through each captured spectra per experiment
    for f_i = 1:n
        y_i = y(:,f_i);
        
        Xsol(:,f_i) = f_Dhat*( (2*lambda*f_Dhat*(A{1,f_i}'*A{1,f_i})*f_Dhat+eye(m))\(2*lambda*f_Dhat*(A{1,f_i}'*y_i)) );
    end
   
    if (debugFlag && (mod(iter,40) == 0))
        plot(sum(Xsol,2));title(['iter = ',num2str(iter)]);
        pause(0.1);
    end
end

% tresholding the solution, get rid of all negatives

%f_ind = (sum(Xsol,2) <= 2*eps);
%Xsol(f_ind,:) = 0;

% return the result
Xsol_est = Xsol;

end
