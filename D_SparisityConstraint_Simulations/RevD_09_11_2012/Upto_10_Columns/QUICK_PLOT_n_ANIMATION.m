%* Quick Plot

[numN,numDELTA,numSNR,numAlphaSample] = size(XSolGrid); % How many grid points are there
numN = 9;
setCOLORMAP = 'bone';

% Create a quick mean() plot
for iN=1:numN
    
    plotERROR = zeros(numSNR,numDELTA);
    for iP=1:numSNR
        for iQ=1:numDELTA
            plotERROR(iP,iQ) = [xERROR{iN,iQ,iP,2}];
        end
    end
    contourf(plotERROR);title(['Error: SNR vs. \Delta spacing for, n = ',num2str(n(iN)),'columns']);
    set(gca,'YTick',1:1:numSNR);
    set(gca,'YTickLabel',{SNR});
    set(gca,'XTick',1:1:numDELTA);
    set(gca,'XTickLabel',{delta});
    ylabel('SNR(dB)');
    xlabel('\Delta spacing');
    colormap(setCOLORMAP);
    colorbar;
    pause;
end

