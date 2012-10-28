% **************************
% plotting Script
%
% **************************

close all;

% Plotting parameters
[numN,numDELTA,numSNR,numSTAT] = size(XSolGrid); % How many grid points are there

% Calculate the error cell matrix if needed
% Error Parameters
xERRORAVG = zeros(numSNR,numDELTA);

% Compute Solution Error
for iN=1:numN
    for iDELTA=1:numDELTA
        for iSNR=1:numSNR
            % Generate a binary vector for xTRUE
            xTRUE = XTrueGrid{iN,iDELTA,iSNR} > 0;
            cumulative = zeros(numSTAT,1);
            for iSTAT=1:numSTAT
                % Save the Card(x) difference distance
                xEST = XSolGrid{iN,iDELTA,iSNR,iSTAT} > 0; % Convert to a binary vector
                %xEST = XSolGrid{iN,iDELTA,iSNR,iSTAT} > (max(XSolGrid{iN,iDELTA,iSNR,iSTAT})*2e-2); % Convert to a binary vector
                xCARDDIFF = nnz(abs(xEST - xTRUE));
                cumulative(iSTAT) = xCARDDIFF > 0;
                % Debug
%                 figure(99);
%                 stem(xTRUE,'b');hold on;stem(2*xEST,'r');title(['\Delta = ',num2str(iDELTA),', SNR = ',...
%                     num2str(iSNR),', iterN = ',num2str(iSTAT),', xCARDDIFF = ',num2str(xCARDDIFF),', CumulativeError = ',num2str(sum(cumulative))]);hold off;
%                 pause(0.1);
            end
            xERRORAVG(iSNR,iDELTA) = sum(cumulative)/numSTAT;
            display(['\Delta = ',num2str(iDELTA),', SNR = ',...
                    num2str(iSNR)]);
            display(['  xERRORAVR = ',num2str(xERRORAVG(iSNR,iDELTA))]);
            %pause(0.1);
        end
    end
end

% Plot the results
contourf(xERRORAVG);title({'CDF for Prob(error)';'SNR vs. \delta Spacing Distance'});
%imagesc(xERRORAVG);title({'CDF for Prob(error)';'SNR vs. \delta Spacing Distance'});
set(gca,'XTick',1:1:numDELTA);
set(gca,'XTickLabel',{delta});
set(gca,'YTick',1:1:numSNR);
set(gca,'YTickLabel',{SNR});
xlabel('\delta Spacing Distance');
ylabel('SNR(dB)');
colormap(jet);  % Sets the number of colors to display
colorbar; % Shows the color bar

