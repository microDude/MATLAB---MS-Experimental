%* Quick Plot

[numN,numDELTA,numSNR,numAlphaSample] = size(XSolGrid); % How many grid points are there
numN = 9;
setCOLORMAP = 'bone';

% Create a quick mean() plot
for iSNR=1:numSNR
    
    plotERROR = zeros(numN,numDELTA);
    for iP=1:numN
        for iQ=1:numDELTA
            plotERROR(iP,iQ) = [xERROR{iP,iQ,iSNR,2}];
        end
    end
    contourf(plotERROR);title(['Error: n vs. \Delta spacing for, SNR = ',num2str(SNR(iSNR)),'dB']);
    set(gca,'YTick',1:1:numN);
    set(gca,'YTickLabel',{2:1:10});
    set(gca,'XTick',1:1:numDELTA);
    set(gca,'XTickLabel',{delta});
    ylabel('number columns of X & y');
    xlabel('\Delta spacing');
    colormap(setCOLORMAP);
    colorbar;
    pause;
end

