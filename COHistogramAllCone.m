
% Plot the pixel value of CO image in the regions of ISI response as histogram.


%% Inputs
clear
aniName = 'AF5';
hemiSph = 'Left';
area = 'V1'; % 'V1' or 'V2'
dataFolder = 'D:\AF5\ISI_analysis\3. CO histogram';
ResultfileName = 'ConeNoncone'; % CO, Liso; Miso, Siso, RG, BY, RB, or other name

labelFront = 15;

COImage = imread(fullfile(dataFolder,'CO5_no holes.tif'));  % 8bit CO image that there is no holes 
% COImage = uint8(mean(imread(fullfile(dataFolder,'xx.png')),3));

%AF5
areaMask = logical(mean(imread(fullfile(dataFolder, 'V1.png')),3));  % Mask of V1
areaCOBleach = logical(ones(size(areaMask)));
areaMaskbv = logical(mean(imread(fullfile(dataFolder, 'AF5_V1_BV.png')),3));  % Mask of V1
areaMaskL = logical(mean(imread(fullfile(dataFolder, 'AF5_Left_Liso_V1_FinalMask.tif')),3));  % Mask of L cone response region
areaMaskM = logical(mean(imread(fullfile(dataFolder, 'AF5_Left_Miso_V1_FinalMask.tif')),3)); % Mask of M cone response region
areaMaskS = logical(mean(imread(fullfile(dataFolder, 'AF5_Left_Siso_V1_FinalMask.tif')),3)); % Mask of S cone response region

%AF4
% areaMask = logical(mean(imread(fullfile(dataFolder, 'V1.png')),3));  % Mask of V1
% areaCOBleach = logical(ones(size(areaMask)));
% areaMaskbv = logical(mean(imread(fullfile(dataFolder, 'AF4_Left_V1_Mask_15_0.tif')),3));  % Mask of V1
% areaMaskL = logical(mean(imread(fullfile(dataFolder, 'AF4_Left_Liso_V1_FinalMask.tif')),3));  % Mask of L cone response region
% areaMaskM = logical(mean(imread(fullfile(dataFolder, 'AF4_Left_Miso_V1_FinalMask.tif')),3)); % Mask of M cone response region
% areaMaskS = logical(mean(imread(fullfile(dataFolder, 'AF4_Left_Siso_V1_FinalMask.tif')),3)); % Mask of S cone response region



% %AE9
% areaCOBleach = logical(mean(imread(fullfile(dataFolder, 'CO_bleach mask.png')),3));  % Mask of V1
% areaMask = logical(mean(imread(fullfile(dataFolder, 'V1 mask_crop.tif')),3));  % Mask of V1
% areaMaskbv = logical(mean(imread(fullfile(dataFolder, 'AE9_Left_V1_Mask_30_crop.tif')),3));  % Mask of V1
% areaMaskL = logical(mean(imread(fullfile(dataFolder, 'AE9_Left_Liso_V1_FinalMask.tif')),3));  % Mask of L cone response region
% areaMaskM = logical(mean(imread(fullfile(dataFolder, 'AE9_Left_Miso_V1_FinalMask.tif')),3)); % Mask of M cone response region
% areaMaskS = logical(mean(imread(fullfile(dataFolder, 'AE9_Left_Siso_V1_FinalMask.tif')),3)); % Mask of S cone response region

% %AF3
% areaMask = logical(mean(imread(fullfile(dataFolder, 'V1 mask.tif')),3));  % Mask of V1
% areaCOBleach = logical(ones(size(areaMask)));
% areaMaskbv = logical(mean(imread(fullfile(dataFolder, 'AF3_Left_V1_Mask_30.tif')),3));  % Mask of V1
% areaMaskL = logical(mean(imread(fullfile(dataFolder, 'AF3_Left_Liso_V1_FinalMask.tif')),3));  % Mask of L cone response region
% areaMaskM = logical(mean(imread(fullfile(dataFolder, 'AF3_Left_Miso_V1_FinalMask.tif')),3)); % Mask of M cone response region
% areaMaskS = logical(mean(imread(fullfile(dataFolder, 'AF3_Left_Siso_V1_FinalMask.tif')),3)); % Mask of S cone response region

% %AE8 
% areaMask = logical(mean(imread(fullfile(dataFolder, 'V1 mask.png')),3));  % Mask of V1
% areaCOBleach = logical(mean(imread(fullfile(dataFolder, 'COBleach_mask.png')),3));  % Mask of CO bleach
% areaMaskbv = logical(mean(imread(fullfile(dataFolder, 'AE8_Left_V1_Mask_18.tif')),3));  % Mask of V1
% areaMaskL = logical(mean(imread(fullfile(dataFolder, 'AE8_Left_Liso_V1_FinalMask.tif')),3));  % Mask of L cone response region
% areaMaskM = logical(mean(imread(fullfile(dataFolder, 'AE8_Left_Miso_V1_FinalMask.tif')),3)); % Mask of M cone response region
% areaMaskS = logical(mean(imread(fullfile(dataFolder, 'AE8_Left_Siso_V1_FinalMask.tif')),3)); % Mask of S cone response region


areaMaskCone = logical(areaMaskL + areaMaskM + areaMaskS);
areaMaskNoncone = logical((1-areaMaskCone).*areaMaskbv);

areaMask = areaMask .* areaCOBleach;
areaMaskCone = areaMaskCone .* areaCOBleach;
areaMaskNoncone = areaMaskNoncone .* areaCOBleach;

resultFolder = strcat(dataFolder, '\', ResultfileName, '_', area, '\');
fileName = strcat(aniName, '_', hemiSph,'_', area, '_', ResultfileName, '_');

if ~isfolder(resultFolder)
    mkdir(resultFolder);    
end 
cd(dataFolder)

%% Filtering
COImage = imgaussfilt(COImage,5); % Smoothing image, filling small holes

%% Image normalization [0 1]
COImage = mat2gray(COImage);
COImage = abs(COImage-1);   % CO intensity, higher number means closer to CO blob
%% Store pixel value selected from CO image based on the ON OFF masks of cones
result.CO = zeros(1, nnz(areaMask));
result.cone = zeros(1, nnz(areaMaskCone));
result.noncone = zeros(1, nnz(areaMaskNoncone));

n = 1;
for i = 1:size(COImage, 1)
    for j = 1:size(COImage, 2)
        if areaMask(i, j) == 1
            result.CO(n) = COImage(i, j);
            n = n+1;
        end
    end
end

n = 1;
for i = 1:size(COImage, 1)
    for j = 1:size(COImage, 2)
        if areaMaskCone(i, j) == 1
            result.cone(n) = COImage(i, j);
            n = n+1;
        end
    end
end

n = 1;
for i = 1:size(COImage, 1)
    for j = 1:size(COImage, 2)
        if areaMaskNoncone(i, j) == 1
            result.noncone(n) = COImage(i, j);
            n = n+1;
        end
    end
end

%%
% Get the median
midCO = median(result.CO);
midCone = median(result.cone);
midNoncone = median(result.noncone);

figure
h0 = histogram(result.CO);
xline(midCO, '--', 'LineWidth', 1, 'Color', [0 0 0], 'Alpha', 1)
hold on
h1 = histogram(result.cone);
xline(midCone, '--', 'LineWidth', 1, 'Color',[1 0 0], 'Alpha', 0.5)
hold on
h2 = histogram(result.noncone);
xline(midNoncone, '--', 'LineWidth', 1, 'Color',[0.5 0.5 0.5], 'Alpha', 0.5)

h0.Normalization = 'probability';
h0.DisplayStyle = 'stairs';
h0.EdgeColor = [0 0 0];
h0.LineWidth = 0.5;
h0.LineStyle = '--';
% h0.FaceColor = [0.2 0.2 0.2];
% h0.FaceAlpha = 1;
h0. BinLimits = [0 1];
h0.BinWidth = 0.05;
result.COVal = h0.Values;

h1.Normalization = 'probability';
h1.EdgeColor = 'none';
h1.FaceColor = [1 0 0];
h1.FaceAlpha = 0.5;
h1. BinLimits = [0 1];
h1.BinWidth = 0.05;
result.ConeVal = h1.Values;

h2.Normalization = 'probability';
h2.EdgeColor = 'none';
h2.FaceColor = [0.5 0.5 0.5];
h2.FaceAlpha = 0.5;
h2. BinLimits = [0 1];
h2.BinWidth = 0.05;
result.nonConeVal = h2.Values;

xlabel('CO Intensity', 'FontSize', labelFront)
ylabel('Proportion of Pixel Number', 'FontSize', labelFront)

legend([h0 h1 h2],'V1 region', 'Cone domain region', 'Non-cone domain region','Location', 'northeast' )
% title([aniName, ' ', hemiSph, ' ', area, ' ', ResultfileName, ' CO'])

savefig([resultFolder, fileName, 'CO.fig'])
saveas(gcf, [resultFolder, fileName, 'CO.png'])
saveas(gcf, [resultFolder, fileName, 'CO.svg'])
% fig2svg(gcf,[resultFolder, fileName, 'CO.svg'])

%% Statistical test of two distributions
[h,p,ci,stats]=ttest2(result.ConeVal, result.nonConeVal, 'Vartype','unequal'); % "Welch's unequal variances t-test" (unequal variances t-test)
result.wt_h = h;
result.wt_p = p;
[p,h,stats] = signrank(result.ConeVal, result.nonConeVal);  % Wilcoxon signed rank testcollapse
result.ws_h = h;
result.ws_p = p;

save([resultFolder, fileName, 'result.mat'], 'result', '-v7.3');

%% Plot the ratio, nomalize by the distribution of CO intensity
figure;
b0 = bar(result.ConeVal ./result.COVal);
hold on;
b1 = bar(result.nonConeVal ./result.COVal);
hold on; 
hline = refline([0 1]);
hline.Color = 'k';
hline.LineStyle = '--';
hline.LineWidth = 1;

b0.EdgeColor = 'none';
b0.FaceColor = [1 0 0];
b0.FaceAlpha = 0.5;
% b0.BarLayout = 'histc';

b1.EdgeColor = 'none';
b1.FaceColor = [0.5 0.5 0.5];
b1.FaceAlpha = 0.5;
% b1.BarLayout = 'histc';

tickFront = 12;
axisThickness = 1;
ax = gca;
set(ax,'linewidth',axisThickness)
ax.FontSize = tickFront;

xlabel('CO Intensity', 'FontSize',labelFront)
ylabel('Ratio of distribution probability', 'FontSize',labelFront)

xticks([0:4:20])
xticklabels({'0','0.2','0.4','0.6','0.8', '1'})
legend([b0, b1], 'Cone domain region', 'Non-cone domain region','Location', 'northwest' )
savefig([resultFolder, fileName, 'CO2.fig'])
saveas(gcf, [resultFolder, fileName, 'CO2.png'])
saveas(gcf, [resultFolder, fileName, 'CO2.svg'])
hold off;
%%
binNum = 20;
cone = result.cone;
nonCone = result.noncone;
total = sort([cone,nonCone]);
totalNum = length(total);
% binEdge = flip(totalNum:-totalNum/binNum:1);
% binX = total(binEdge);
binX = 0:1/binNum:1;
% binX = binX(2:end);

binCount = zeros(2,binNum);
for bi = 1:binNum
    cone(cone<=binX(bi+1)) = nan;
    nonCone(nonCone<=binX(bi+1)) = nan;
    binCount(1,bi) = sum(isnan(cone))/(sum(isnan(cone))+sum(isnan(nonCone)));
    binCount(2,bi) = sum(isnan(nonCone))/(sum(isnan(cone))+sum(isnan(nonCone)));
    cone = cone(~isnan(cone));
    nonCone = nonCone(~isnan(nonCone));
end
figure
b=bar(binCount','stacked');
b(1).EdgeColor = 'none';
b(1).FaceColor = [1 0 0];
b(1).FaceAlpha = 0.5;
b(2).EdgeColor = 'none';
b(2).FaceColor = [0.5 0.5 0.5];
b(2).FaceAlpha = 0.5;

tickFront = 12;
axisThickness = 1;
ax = gca;
set(ax,'linewidth',axisThickness)
ax.FontSize = tickFront;


xlabel('CO Intensity', 'FontSize',labelFront)
ylabel('Proportion in Each Bin', 'FontSize',labelFront)

ylim([0 1.2])
xticks([0:4:20])
xticklabels({'0','0.2','0.4','0.6','0.8', '1'})
legend([b(1), b(2)], 'Cone domain region', 'Non-cone domain region','Location', 'northwest' )
savefig([resultFolder, fileName, 'CO3.fig'])
saveas(gcf, [resultFolder, fileName, 'CO3.png'])
saveas(gcf, [resultFolder, fileName, 'CO3.svg'])
%% Plot
% 
% dotSize = 20;
% dotTransparency = 0.2;
% labelFront = 10;
% tickFront = 12;
% axisThickness = 2;
% 
% figure
% scatter(result.L, result.M, dotSize,'filled', 'MarkerFaceAlpha', dotTransparency)
% 
% xlabel('L cone', 'FontSize',labelFront)
% ylabel('M cone', 'FontSize',labelFront)
% 
% xlim([-1.1 1.1])
% ylim([-1.1 1.1])
% xticks([-1:0.5:1])
% yticks([-1:0.5:1])
% xticklabels({'-1','-0.5','0','0.5','1'})
% yticklabels({'-1','-0.5','','0.5','1'})
% ax = gca;
% set(ax,'linewidth',axisThickness)
% ax.FontSize = tickFront;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% axis square
% savefig([resultFolder, fileName, 'LMcone.fig'])
% 
% figure
% scatter(result.L, result.S, dotSize,'filled', 'MarkerFaceAlpha', dotTransparency)
% 
% xlabel('L cone', 'FontSize',labelFront)
% ylabel('S cone', 'FontSize',labelFront)
% 
% xlim([-1.1 1.1])
% ylim([-1.1 1.1])
% xticks([-1:0.5:1])
% yticks([-1:0.5:1])
% xticklabels({'-1','-0.5','0','0.5','1'})
% yticklabels({'-1','-0.5','','0.5','1'})
% ax = gca;
% set(ax,'linewidth',axisThickness)
% ax.FontSize = tickFront;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% axis square
% savefig([resultFolder, fileName, 'LScone.fig'])
% 
% figure
% scatter(result.M, result.S, dotSize,'filled', 'MarkerFaceAlpha', dotTransparency)
% 
% xlabel('M cone', 'FontSize',labelFront)
% ylabel('S cone', 'FontSize',labelFront)
% 
% xlim([-1.1 1.1])
% ylim([-1.1 1.1])
% xticks([-1:0.5:1])
% yticks([-1:0.5:1])
% xticklabels({'-1','-0.5','0','0.5','1'})
% yticklabels({'-1','-0.5','','0.5','1'})
% ax = gca;
% set(ax,'linewidth',axisThickness)
% ax.FontSize = tickFront;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% axis square
% savefig([resultFolder, fileName, 'MScone.fig'])
%% Save current script
Version='4';
fileNameAndLocation='C:\Users\lyr19\Dropbox\Github\ISI_analysis_NHP\ISI_analysis_Peichao\Matlab_data_process\12_Spatial_pattern';
% fileNameAndLocation = mfileName('fullpath');   % mfileName('fullpath') does not work in current version!!
% fileNameAndLocation = which(mfileName);
newbackup= strcat(resultFolder,strcat('COHistogramAllCone_',Version, '.m'));
currentfile=strcat(fileNameAndLocation, '\COHistogramAllCone.m');
copyfile(currentfile,newbackup);

