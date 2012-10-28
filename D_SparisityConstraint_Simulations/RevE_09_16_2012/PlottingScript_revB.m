% **************************
% plotting Script
%
% **************************

close all;

% Plotting parameters
CalculateError = 1;                      % Flag to calculate error, yes or no
[numN,numDELTA,numSNR] = size(XSolGrid); % How many grid points are there
fixed_n.flag = true;                     % Do you want to use a fixed number of columns?
fixed_n.value = 4;                       % If yes, how many columns would you like?
% Choose either:
%   n 
%   delta
%   SNR
axis.row = 'SNR';
axis.column = 'delta';

% Calculate the error cell matrix if needed
if (CalculateError == 1)
    % Error Parameters
    xERROR = cell(numN,numDELTA,numSNR);
    
    % Compute Solution Error
    for iN=1:numN
        for iDELTA=1:numDELTA
            for iSNR=1:numSNR
                xTRUE = XTrueGrid{iN,iDELTA,iSNR};
                xTRUE = xTRUE > 0; % Convert to a binary vector
                % Save the Card(x) difference distance
                xEST = XSolGrid{iN,iDELTA,iSNR} > 0; % Convert to a binary vector
                xCARDDIFF = nnz(abs(xEST - xTRUE));
                xERROR{iN,iDELTA,iSNR} = xCARDDIFF;
            end
        end
    end
end

% Create a R^3 plot of the desired parameter space
switch axis.row
    case 'n'
        switch axis.column
            case 'delta'
                plotERROR = zeros(numN,numDELTA);
                for iP=1:numN
                    for iQ=1:numDELTA
                        plotERROR(iP,iQ) = mean([xERROR{iP,iQ,:}]);
                    end
                end
                numTick = ceil(max(max(plotERROR)) + 1);% account for zero error, 0,1,2,...
                contourf(plotERROR,numTick);title('Error: n vs. \Delta spacing');
                set(gca,'YTick',1:1:numN);
                set(gca,'YTickLabel',{n});
                set(gca,'XTick',1:1:numDELTA);
                set(gca,'XTickLabel',{delta});
                ylabel('number columns of X & y');
                xlabel('\Delta spacing');
            case 'SNR'
                plotERROR = zeros(numN,numSNR);
                for iP=1:numN
                    for iQ=1:numSNR
                        plotERROR(iP,iQ) = mean([xERROR{iP,:,iQ}]);
                    end
                end
                numTick = ceil(max(max(plotERROR)) + 1);% account for zero error, 0,1,2,...
                contourf(plotERROR,numTick);title('Error: n vs. SNR');
                set(gca,'YTick',1:1:numN);
                set(gca,'YTickLabel',{n});
                set(gca,'XTick',1:1:numSNR);
                set(gca,'XTickLabel',{SNR});
                ylabel('number columns of X & y');
                xlabel('SNR(dB)');
            otherwise
                warning('Unexpected plot type');
        end
    case 'delta'
        switch axis.column
            case 'n'
                plotERROR = zeros(numDELTA,numN);
                for iP=1:numDELTA
                    for iQ=1:numN
                        plotERROR(iP,iQ) = mean([xERROR{iQ,iP,:}]);
                    end
                end
                numTick = ceil(max(max(plotERROR)) + 1);% account for zero error, 0,1,2,...
                contourf(plotERROR,numTick);title('Error: \Delta spacing vs. n');
                set(gca,'YTick',1:1:numDELTA);
                set(gca,'YTickLabel',{delta});
                set(gca,'XTick',1:1:numN);
                set(gca,'XTickLabel',{n});
                ylabel('\Delta spacing');
                xlabel('number columns of X & y');
            case 'SNR'
                plotERROR = zeros(numDELTA,numSNR);
                
                for iP=1:numDELTA
                    for iQ=1:numSNR
                        if fixed_n.flag
                            plotERROR(iP,iQ) = xERROR{fixed_n.value,iP,iQ};
                        else
                            plotERROR(iP,iQ) = mean([xERROR{:,iP,iQ}]);
                        end
                    end
                end
                numTick = ceil(max(max(plotERROR)) + 1);% account for zero error, 0,1,2,...
                contourf(plotERROR,numTick);title('Error: \Delta spacing vs. SNR');
                set(gca,'YTick',1:1:numDELTA);
                set(gca,'YTickLabel',{delta});
                set(gca,'XTick',1:1:SNR);
                set(gca,'XTickLabel',{SNR});
                ylabel('\Delta spacing');
                xlabel('SNR(dB)');
            otherwise
                warning('Unexpected plot type');
        end
    case 'SNR'
        switch axis.column
            case 'n'
                plotERROR = zeros(numSNR,numN);
                for iP=1:numSNR
                    for iQ=1:numN
                        plotERROR(iP,iQ) = mean([xERROR{iQ,:,iP}]);
                    end
                end
                numTick = ceil(max(max(plotERROR)) + 1);% account for zero error, 0,1,2,...
                contourf(plotERROR,numTick);title('Error: SNR vs. n');
                set(gca,'YTick',1:1:numSNR);
                set(gca,'YTickLabel',{SNR});
                set(gca,'XTick',1:1:numN);
                set(gca,'XTickLabel',{n});
                ylabel('SNR (dB)');
                xlabel('number columns of X & y');
            case 'delta'
                plotERROR = zeros(numSNR,numDELTA);
                for iP=1:numSNR
                    for iQ=1:numDELTA
                        if fixed_n.flag
                            plotERROR(iP,iQ) = xERROR{fixed_n.value,iQ,iP};
                        else
                            plotERROR(iP,iQ) = mean([xERROR{:,iQ,iP}]);
                        end
                    end
                end
                numTick = ceil(max(max(plotERROR)) + 1);% account for zero error, 0,1,2,...
                contourf(plotERROR,numTick);title('Error: SNR vs. \Delta spacing');
                set(gca,'YTick',1:1:numSNR);
                set(gca,'YTickLabel',{SNR});
                set(gca,'XTick',1:1:numDELTA);
                set(gca,'XTickLabel',{delta});
                ylabel('SNR (dB)');
                xlabel('\Delta spacing');
            otherwise
                warning('Unexpected plot type');
        end
end

% Change the colorbar
colormap(jet(numTick));  % Sets the number of colors to display
hcb = colorbar; % Shows the color bar
YColorBar = get(hcb,'YLim'); % Gets the current min, max color values
set(hcb,'YTickLabelMode','manual','YTick',YColorBar(1):YColorBar(2)/numTick:YColorBar(2),'YTickLabel',{0:1:numTick-1,''}); % 
%set(hcb,'YTick',[0,1,2,3,4],'YTickLabel',{'0','1','2','3'});


% Display values of plot error
plotERROR
    
