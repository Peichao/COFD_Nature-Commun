
%% User input, file location
% clear
clc
disk = 'X:';
aniName = '';
unit = ''; % 
dataRoot = strcat(disk, '\',aniName,'\2P_analysis\U',unit, '\_Summary\');
dataFolder = strcat(dataRoot, 'Multiplane\2. Correlation of Cone weights based on ISI\');
ResultfileName = ''; % 
hueaucThres = 0.8;
%%
cd(dataFolder)
resultFolder = strcat(dataFolder, ResultfileName, '\');
fileName = strcat(aniName, '_U', unit,'_', ResultfileName, '_');

if ~isfolder(resultFolder)
    mkdir(resultFolder);
end 
%% 
dklmap = load('');
matFile = dir(fullfile(dataFolder, '*.mat'));
pl0 = load(matFile(1).name);
pl1 = load(matFile(2).name);

pltable = [struct2table(pl0.result);struct2table(pl1.result)];

hueData0 = readtable(strcat(dataRoot,'DataExport\', aniName, '_', unit, '*.csv'));
hue0 = hueData0.maxhue;
cmap0 = dklmap.colors((hue0+1),1:3);

visResp0 = hueData0.visResp;
ax0 = hueData0.hueaxauc;
di0 = hueData0.huediauc;
hueSelt0 = ((ax0 > hueaucThres) | (di0 > hueaucThres)) & (visResp0 > 0);

hueData1 = readtable(strcat(dataRoot,'DataExport\', aniName, '_', unit, '*.csv'));
hue1 = hueData1.maxhue;
cmap1 = dklmap.colors((hue1+1),1:3);

visResp1 = hueData1.visResp;
ax1 = hueData1.hueaxauc;
di1 = hueData1.huediauc;
hueSelt1 = ((ax1 > hueaucThres) | (di1 > hueaucThres)) & (visResp1 > 0);

pltable.maxhue = [hue0;hue1];
pltable.cmap = [cmap0;cmap1];
pltable.hueselect =[hueSelt0;hueSelt1];

%%
nbins = 50;
pos1 = [0.2 0.2 0.7 0.7];
pos2 = [0.2 0.1 0.7 0.1];
pos3 = [0.10 0.2 0.1 0.7];

% X Y label
labelFront = 30;
xlabelPosition = [0 -1.3];
ylabelPosition = [-1.3 0];

dotSize = 40;
% dotTransparency = result.COFD;


% X Y axis
aRange = [-1 1];
bRange = [0 200];  % Range
axThickness = 6.0;
lnThickness = 6.0;
% atickLabel = {'','4','8'};
% btickLabel = {'0    ','4','8'};
atickLabel = {};
btickLabel = {};
% atickLabel = {'','100','200'};
% btickLabel = {'0    ','100','200'};

% X Y ticks
tickFront = 25;
x_ticks = -1:0.5:1;
y_ticks = -1:0.5:1;
x1_ticks = 0:bRange(2)/2:bRange(2);
y1_ticks = 0:bRange(2)/2:bRange(2);

% axis & line color
axColor = [0 0 0];
lnColor = [0 0 0];

selected = find((pltable.COFD == 1) & (pltable.hueselect == 1));
cellNum =  length(selected);
for rp = 1:3
    
    switch rp
        case 1
            a = pltable.Lscale;
            b = pltable.Mscale;
            xName = 'L Cone';
            yName = 'M Cone';
            imgName = strcat(fileName, ['LM cone_cellNum', num2str(cellNum)]);
        case 2
            a = pltable.Lscale;
            b = pltable.Sscale;
            xName = 'L Cone';
            yName = 'S Cone';
            imgName = strcat(fileName, ['LS cone_cellNum', num2str(cellNum)]);
        case 3  
            a = pltable.Mscale;
            b = pltable.Sscale;
            xName = 'M Cone';
            yName = 'S Cone';
            imgName = strcat(fileName, ['MS cone_cellNum', num2str(cellNum)]);
    end
    
    % Initilize figure
    f = figure;
    f.InnerPosition = [100 10 1300 1300];  % Define drawable region
    colormap jet

    % First subplot in the figure
    ax1 = subplot('Position',pos1);
    scatter(a(selected), b(selected), dotSize, pltable.cmap(selected,1:3),'filled')  % Plot as dots

%     histogram2(a(selected), b(selected), nbins, 'DisplayStyle','tile','ShowEmptyBins','off');
%     normN = histcounts2(a(selected),b(selected),nbins)./totalPixel;   
%     if rp == 1
%         result.pLM = normN;
%     elseif rp == 2
%         result.pLS = normN;
%     elseif rp == 3
%         result.pMS = normN;
%     end
    
%     xlabel(xName, 'FontSize', labelFront,'Position', xlabelPosition)
%     ylabel(yName, 'FontSize', labelFront,'Position', ylabelPosition)
    % hXLbl=xlabel('XLabel','Position',[Xlb Ylb],'VerticalAlignment','top','HorizontalAlignment','center'); 

    xlim(aRange)
    ylim(aRange)
    xticks(x_ticks)
    yticks(y_ticks)
    xticklabels({})
    yticklabels({})
    set(ax1,'linewidth', axThickness)
    set(ax1,'Layer','top')
    % ax.XAxisLocation = 'origin';
    % ax.YAxisLocation = 'origin';
    ax1.XColor = axColor;
    ax1.YColor = axColor;
    axis square
    grid off
    box off
    
    % Plot Colorbar
%     cb = colorbar;
%     cb.Position = cbPosition;
%     cb.FontSize = cbFontSize;
%     cb.Color = cbColor;
    
    % Plot indicator line
    line([0 0], [-1 1], 'LineWidth', lnThickness, 'LineStyle', '--', 'Color', lnColor);
    line([-1 1], [0 0], 'LineWidth', lnThickness, 'LineStyle', '--', 'Color', lnColor);

    % add second subplot
    hold on;
    ax2 = subplot('Position',pos2);
    hL = histogram(a(selected), nbins);

    xlim(aRange)
    ylim(bRange)
    ax2.XAxisLocation = 'top';
    set(ax2,'linewidth',axThickness)
    ax2.XColor = axColor;
    ax2.YColor = axColor;
    ax2.FontSize = tickFront;
    xticks(x_ticks);
%     xticklabels({'      -1','-0.5','0   ','0.5','1  '})
    xticklabels({})

    yticks(y1_ticks);
    yticklabels(atickLabel)
    grid off
    box off
    set(ax2,'view',[0 -90])

    % add third subplot
    hold on;
    ax3 = subplot('Position',pos3);
    hM = histogram(b(selected), nbins, 'Orientation','horizontal');
    
    xlim(bRange)
    ylim(aRange)
    ax3.YAxisLocation = 'right';
    set(ax3,'linewidth',axThickness)
    ax3.XColor = axColor;
    ax3.YColor = axColor;
    ax3.FontSize = tickFront;
    yticks(y_ticks)
%     yticklabels({'',' -0.5',' 0',' 0.5',' 1'})
    yticklabels({})
    xticks(x1_ticks);
    xticklabels(btickLabel)
    grid off
    box off
    set(gca,'view',[180 -90])

    % Save figure
    
%     title([aniName, ' ', hemiSph, ' ', area, ' ', ResultfileName, ' LM'])
%     savefig([resultFolder, imgName, '.fig'])
%     saveas(gcf, [resultFolder, imgName, '.png'])
%     saveas(gcf, [resultFolder, imgName, '.svg'])
end

% save([resultFolder, fileName, 'result.mat'], 'result', '-v7.3');
% close all;
