% **************************
% plotting Script
%
% **************************

% Plotting parameters
CalculateError = 0;                                     % Flag to calculate error, yes or no
[numN,numDELTA,numSNR,numAlphaSample] = size(XSolGrid); % How many grid points are there
eps3 = 1e-8;
% Choose either:
%   n 
%   delta
%   SNR
%   AlphaSample
axis.row = 'n';
axis.column = 'SNR';
setCOLORMAP = 'jet';

% Calculate the error cell matrix if needed
if (CalculateError == 1)
    % Error Parameters
    xERROR = cell(numN,numDELTA,numSNR,numAlphaSample);
    
    % Compute Solution Error
    for iN=1:numN
        for iDELTA=1:numDELTA
            for iSNR=1:numSNR
                xTRUE = XTrueGrid{iN,iDELTA,iSNR};
                xTRUE = xTRUE > eps3; % Convert to a binary vector
                for iVARIETY=1:numAlphaSample
                    % Save the Card(x) difference distance
                    xEST = XSolGrid{iN,iDELTA,iSNR,iVARIETY} > eps3; % Convert to a binary vector
                    xCARDDIFF = nnz(xEST - xTRUE);
                    xERROR{iN,iDELTA,iSNR,iVARIETY} = xCARDDIFF;
                end
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
                        plotERROR(iP,iQ) = mean([xERROR{iP,iQ,:,:}]);
                    end
                end
                contourf(plotERROR);title('Error: n vs. \Delta spacing');
                set(gca,'YTick',1:1:numN);
                set(gca,'YTickLabel',{n});
                set(gca,'XTick',1:1:numDELTA);
                set(gca,'XTickLabel',{delta});
                ylabel('number columns of X & y');
                xlabel('\Delta spacing');
                colormap(setCOLORMAP);
                colorbar;
            case 'SNR'
                plotERROR = zeros(numN,numSNR);
                for iP=1:numN
                    for iQ=1:numSNR
                        plotERROR(iP,iQ) = mean([xERROR{iP,:,iQ,:}]);
                    end
                end
                contourf(plotERROR);title('Error: n vs. SNR');
                set(gca,'YTick',1:1:numN);
                set(gca,'YTickLabel',{n});
                set(gca,'XTick',1:1:numSNR);
                set(gca,'XTickLabel',{SNR});
                ylabel('number columns of X & y');
                xlabel('SNR(dB)');
                colormap(setCOLORMAP);
                colorbar;
            case 'AlphaSample'
                plotERROR = zeros(numN,numAlphaSample);
                for iP=1:numN
                    for iQ=1:numAlphaSample
                        plotERROR(iP,iQ) = mean([xERROR{iP,:,:,iQ}]);
                    end
                end
                contourf(plotERROR);title('Error: n vs. Variety of A');
                set(gca,'YTick',1:1:numN);
                set(gca,'YTickLabel',{n});
                set(gca,'XTick',1:1:numAlphaSample);
                set(gca,'XTickLabel',AlphaSampleVariety);
                ylabel('number columns of X & y');
                xlabel('Variety of A');
                colormap(setCOLORMAP);
                colorbar;
            otherwise
                warning('Unexpected plot type');
        end
    case 'delta'
        switch axis.column
            case 'n'
                plotERROR = zeros(numDELTA,numN);
                for iP=1:numDELTA
                    for iQ=1:numN
                        plotERROR(iP,iQ) = mean([xERROR{iQ,iP,:,:}]);
                    end
                end
                contourf(plotERROR);title('Error: \Delta spacing vs. n');
                set(gca,'YTick',1:1:numDELTA);
                set(gca,'YTickLabel',{delta});
                set(gca,'XTick',1:1:numN);
                set(gca,'XTickLabel',{n});
                ylabel('\Delta spacing');
                xlabel('number columns of X & y');
                colormap(setCOLORMAP);
                colorbar;
            case 'SNR'
                plotERROR = zeros(numDELTA,numSNR);
                for iP=1:numDELTA
                    for iQ=1:numSNR
                        plotERROR(iP,iQ) = mean([xERROR{:,iP,iQ,:}]);
                    end
                end
                contourf(plotERROR);title('Error: \Delta spacing vs. SNR');
                set(gca,'YTick',1:1:numDELTA);
                set(gca,'YTickLabel',{delta});
                set(gca,'XTick',1:1:SNR);
                set(gca,'XTickLabel',{SNR});
                ylabel('\Delta spacing');
                xlabel('SNR(dB)');
                colormap(setCOLORMAP);
                colorbar;
            case 'AlphaSample'
                plotERROR = zeros(numDELTA,numAlphaSample);
                for iP=1:numDELTA
                    for iQ=1:numAlphaSample
                        plotERROR(iP,iQ) = mean([xERROR{:,iP,:,iQ}]);
                    end
                end
                contourf(plotERROR);title('Error: \Delta spacing vs. Variety of A');
                set(gca,'YTick',1:1:numDELTA);
                set(gca,'YTickLabel',{delta});
                set(gca,'XTick',1:1:numAlphaSample);
                set(gca,'XTickLabel',AlphaSampleVariety);
                ylabel('\Delta spacing');
                xlabel('Variety of A');
                colormap(setCOLORMAP);
                colorbar;
            otherwise
                warning('Unexpected plot type');
        end
    case 'SNR'
        switch axis.column
            case 'n'
                plotERROR = zeros(numSNR,numN);
                for iP=1:numSNR
                    for iQ=1:numN
                        plotERROR(iP,iQ) = mean([xERROR{iQ,:,iP,:}]);
                    end
                end
                contourf(plotERROR);title('Error: SNR vs. n');
                set(gca,'YTick',1:1:numSNR);
                set(gca,'YTickLabel',{SNR});
                set(gca,'XTick',1:1:numN);
                set(gca,'XTickLabel',{n});
                ylabel('SNR (dB)');
                xlabel('number columns of X & y');
                colormap(setCOLORMAP);
                colorbar;
            case 'delta'
                plotERROR = zeros(numSNR,numDELTA);
                for iP=1:numSNR
                    for iQ=1:numDELTA
                        plotERROR(iP,iQ) = mean([xERROR{:,iQ,iP,:}]);
                    end
                end
                contourf(plotERROR);title('Error: SNR vs. \Delta spacing');
                set(gca,'YTick',1:1:numSNR);
                set(gca,'YTickLabel',{SNR});
                set(gca,'XTick',1:1:numDELTA);
                set(gca,'XTickLabel',{delta});
                ylabel('SNR (dB)');
                xlabel('\Delta spacing');
                colormap(setCOLORMAP);
                colorbar;
            case 'AlphaSample'
                plotERROR = zeros(numSNR,numAlphaSample);
                for iP=1:numSNR
                    for iQ=1:numAlphaSample
                        plotERROR(iP,iQ) = mean([xERROR{:,:,iP,iQ}]);
                    end
                end
                contourf(plotERROR);title('Error: SNR vs. Variety of A');
                set(gca,'YTick',1:1:numSNR);
                set(gca,'YTickLabel',{SNR});
                set(gca,'XTick',1:1:numAlphaSample);
                set(gca,'XTickLabel',AlphaSampleVariety);
                ylabel('SNR (dB)');
                xlabel('Variety of A');
                colormap(setCOLORMAP);
                colorbar;
            otherwise
                warning('Unexpected plot type');
        end
    case 'AlphaSample'
        numiP = numAlphaSample;
        switch axis.column
            case 'n'
                plotERROR = zeros(numAlphaSample,numN);
                for iP=1:numAlphaSample
                    for iQ=1:numN
                        plotERROR(iP,iQ) = mean([xERROR{iQ,:,:,iP}]);
                    end
                end
                contourf(plotERROR);title('Error: Variety of A vs. n');
                set(gca,'YTick',1:1:numAlphaSample);
                set(gca,'YTickLabel',AlphaSampleVariety);
                set(gca,'XTick',1:1:numN);
                set(gca,'XTickLabel',{n});
                ylabel('Variety of A');
                xlabel('number columns of X & y');
                colormap(setCOLORMAP);
                colorbar;
            case 'delta'
                plotERROR = zeros(numAlphaSample,numDELTA);
                for iP=1:numAlphaSample
                    for iQ=1:numDELTA
                        plotERROR(iP,iQ) = mean([xERROR{:,iQ,:,iP}]);
                    end
                end
                contourf(plotERROR);title('Error: Variety of A vs. \Delta spacing');
                set(gca,'YTick',1:1:numAlphaSample);
                set(gca,'YTickLabel',AlphaSampleVariety);
                set(gca,'XTick',1:1:numDELTA);
                set(gca,'XTickLabel',{delta});
                ylabel('Variety of A');
                xlabel('\Delta spacing');
                colormap(setCOLORMAP);
                colorbar;
            case 'SNR'
                plotERROR = zeros(numAlphaSample,numSNR);
                for iP=1:numAlphaSample
                    for iQ=1:numSNR
                        plotERROR(iP,iQ) = mean([xERROR{:,:,iQ,iP}]);
                    end
                end
                contourf(plotERROR);title('Error: Variety of A vs. SNR');
                set(gca,'YTick',1:1:numAlphaSample);
                set(gca,'YTickLabel',AlphaSampleVariety);
                set(gca,'XTick',1:1:numSNR);
                set(gca,'XTickLabel',{SNR});
                ylabel('Variety of A');
                xlabel('SNR (dB)');
                colormap(setCOLORMAP);
                colorbar;
            otherwise
                warning('Unexpected plot type');
        end
    otherwise
        warning('Unexpected plot type');
end
