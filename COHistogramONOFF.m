
% Plot the pixel value of CO image in the regions of ISI response as histogram.


%% Inputs
clear;
aniName = '';
hemiSph = '';
area = ''; % 'V1' or 'V2'
dataFolder = '';
coneType = ''; % L, M, S, WB

labelFront = ;

switch coneType 
    case 'L'
        ResultfileName = 'Lonoff';
        maskOnName = 'Liso_on.png';
        maskOffName = 'Liso_off.png';
        legend1 = 'L-cone ON region';
        legend2 = 'L-cone OFF region';
    case 'M'
        ResultfileName = 'Monoff';
        maskOnName = 'Miso_on.png';
        maskOffName = 'Miso_off.png';
        legend1 = 'M-cone ON region';
        legend2 = 'M-cone OFF region';
    case 'S'
        ResultfileName = 'Sonoff';
        maskOnName = 'Siso_on.png';
        maskOffName = 'Siso_off.png';
        legend1 = 'S-cone ON region';
        legend2 = 'S-cone OFF region';
    case 'WB'
        ResultfileName = 'WBonoff';
        maskOnName = 'WB_on.png';
        maskOffName = 'WB_off.png';
        legend1 = 'Wihte Black ON region';
        legend2 = 'Wihte Black OFF region';    
end


COImage = imread(fullfile(dataFolder,'CO4_no holes.tif'));  % 8bit CO image that there is no holes 

areaMask = logical(mean(imread(fullfile(dataFolder,'V1 mask.tif')),3));  % Mask of cortical region
areaMaskOn = logical(mean(imread(fullfile(dataFolder, maskOnName)),3));  % Mask of ON response region
areaMaskOff = logical(mean(imread(fullfile(dataFolder, maskOffName)),3)); % Mask of OFF response region
areaMaskBloodvessel = logical(mean(imread(fullfile(dataFolder,'AF3_Left_V1_Mask_30.tif')),3)); % Mask of OFF response region

resultFolder = strcat(dataFolder, '\', ResultfileName, '_', area, '\');
fileName = strcat(aniName, '_', hemiSph,'_', area, '_', ResultfileName, '_');

if ~isfolder(resultFolder)
    mkdir(resultFolder);    
end 
cd(dataFolder)

%% Filtering
COImage = imgaussfilt(COImage,5); % Smoothing image, filling small holes
% figure 
% imshow(COImage)
%% Image normalization [0 1]
COImage = mat2gray(COImage);
COImage = abs(COImage-1);   % CO intensity, higher number means closer to CO blob

%% Store pixel value selected from CO image based on the ON OFF masks of cones
areaMaskOn = areaMaskOn.*areaMaskBloodvessel;
areaMaskOff = areaMaskOff.*areaMaskBloodvessel;

result.CO = zeros(1, nnz(areaMask));
result.On = zeros(1, nnz(areaMaskOn));
result.Off = zeros(1, nnz(areaMaskOff));

% Region in areaMask
n = 1;
for i = 1:size(COImage, 1)
    for j = 1:size(COImage, 2)
        if areaMask(i, j) == 1
            result.CO(n) = COImage(i, j);
            n = n+1;
        end
    end
end

% ON regions
n = 1;
for i = 1:size(COImage, 1)
    for j = 1:size(COImage, 2)
        if areaMaskOn(i, j) == 1
            result.On(n) = COImage(i, j);
            n = n+1;
        end
    end
end

% OFF regions
n = 1;
for i = 1:size(COImage, 1)
    for j = 1:size(COImage, 2)
        if areaMaskOff(i, j) == 1
            result.Off(n) = COImage(i, j);
            n = n+1;
        end
    end
end

save([resultFolder, fileName, 'result.mat'], 'result', '-v7.3');

%%

midCO = median(result.CO);
midOn = median(result.On);
midOff = median(result.Off);


BinWidth = 0.05;

figure
h0 = histogram(result.CO);
xline(midCO, '--', 'LineWidth', 1, 'Color', [0 0 0], 'Alpha', 1);
hold on
h1 = histogram(result.On);
xline(midOn, '--', 'LineWidth', 1, 'Color', [0.5 0 0], 'Alpha', 0.5);
hold on
h2 = histogram(result.Off);
xline(midOff, '--', 'LineWidth', 1, 'Color', [0 0 0.5], 'Alpha', 0.5);

h0.Normalization = 'probability';
h0.DisplayStyle = 'stairs';
h0.EdgeColor = [0 0 0];
h0.LineWidth = 1;
h0.LineStyle = '--';
% h0.FaceColor = [0.2 0.2 0.2];
% h0.FaceAlpha = 1;
h0. BinLimits = [0 1];
h0.BinWidth = BinWidth;

h1.Normalization = 'probability';
h1.EdgeColor = 'none';
h1.FaceColor = [0.5 0 0];
h1.FaceAlpha = 0.5;
h1. BinLimits = [0 1];
h1.BinWidth = BinWidth;

h2.Normalization = 'probability';
h2.EdgeColor = 'none';
h2.FaceColor = [0 0 0.5];
h2.FaceAlpha = 0.5;
h2. BinLimits = [0 1];
h2.BinWidth = BinWidth;

ylim([0 0.15])
xlabel('CO Intensity', 'FontSize', labelFront)
ylabel('Proportion of Pixel Number', 'FontSize', labelFront)

[lgd,legendIcons]=legend([h0 h1 h2],'V1 region', legend1, legend2, 'Location', 'northeast');
% title([aniName, ' ', hemiSph, ' ', area, ' ', ResultfileName, ' CO'])


savefig([resultFolder, fileName, 'CO.fig'])
saveas(gcf, [resultFolder, fileName, 'CO.svg'])
saveas(gcf, [resultFolder, fileName, 'CO.png'])
