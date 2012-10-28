%* Quick Plot

[numN,numDELTA,numSNR,numAlphaSample] = size(XSolGrid); % How many grid points are there
numN = 9;
setCOLORMAP = 'bone';

% Create a quick mean() plot
if 1
    plotERROR = zeros(numN,numDELTA);
    for iP=1:numN
        for iQ=1:numDELTA
            plotERROR(iP,iQ) = mean([xERROR{iP,iQ,:,:}]);
        end
    end
    contourf(plotERROR);title('MEAN() Error: n vs. \Delta spacing');
    set(gca,'YTick',1:1:numN);
    set(gca,'YTickLabel',{2:1:10});
    set(gca,'XTick',1:1:numDELTA);
    set(gca,'XTickLabel',{delta});
    ylabel('number columns of X & y');
    xlabel('\Delta spacing');
    colormap(setCOLORMAP);
    colorbar;
else
    plotERROR = zeros(numN,numDELTA);
    for iP=1:numN
        for iQ=1:numDELTA
            plotERROR(iP,iQ) = [xERROR{iP,iQ,1,2}];
        end
    end
    contourf(plotERROR);title('Error: n vs. \Delta spacing');
    set(gca,'YTick',1:1:numN);
    set(gca,'YTickLabel',{2:1:10});
    set(gca,'XTick',1:1:numDELTA);
    set(gca,'XTickLabel',{delta});
    ylabel('number columns of X & y');
    xlabel('\Delta spacing');
    colormap(setCOLORMAP);
    colorbar;
end
