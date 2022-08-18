

% Peichao Li

% Plot the histogram of cone ON/OFF map based on the normalized pixel value.
% Plot the bivariate histogram of ON/OFF map of two types of cones, e.g. L vs M.

%% Inputs
% clear
% aniName = 'AF4';
% hemiSph = 'Left';
% area = 'V1'; % 'V1' or 'V2'
% % isimapFolder = 'O:\AF4\ISI_analysis\0. original ISI maps';
% dataFolder = 'O:\AF4\ISI_analysis\2. Spatial correlation of COFD';
% ResultfileName = 'Achromatic OnOff Relation'; % 
% 
% % isiData=dir(fullfile(isimapFolder, '*_meta.mat'));
% % image = load(fullfile(isiData.folder, isiData.name));
% % isiMap = {image.ISImeta.map};
% % 
% % Image1 = isiMap(ismember({image.ISImeta.mapName},'Liso'));
% % Image2 = isiMap(ismember({image.ISImeta.mapName},'Miso'));
% % Image3 = isiMap(ismember({image.ISImeta.mapName},'Siso'));
% 
% Image1 = imread(fullfile(dataFolder,'AF4_Left_u000_009_Liso_freq1_1st_0_Clip1.5_10_Phase_Filted.tif'));  % Should be 2D array in uint8; I use ImageJ to convert image from rgb to 8bits.
% Image2 = imread(fullfile(dataFolder,'AF4_Left_u000_010_Miso_freq1_2nd_0_Clip1.5_9_Phase_Filted.tif'));
% Image3 = imread(fullfile(dataFolder,'AF4_Left_u000_011_Siso_freq1_2nd_0_Clip1.5_9_Phase_Filted.tif'));
% Image4 = imread(fullfile(dataFolder,'AF4_Left_u000_012_WB_freq1_2nd_0_Clip1.5_9_Phase_Filted.tif'));
% 
% areaMask = logical(imread(fullfile(dataFolder,'AF4_Left_LMSA_V1_DomainMask.tif')));
% bvMask = logical(mean(imread(fullfile(dataFolder,'AF4_Left_V1_Mask_15.tif')),3));
% areaMask = areaMask .* bvMask;
%% Inputs
clear
aniName = 'AF5';
hemiSph = 'Left';
area = 'V1'; % 'V1' or 'V2'
isimapFolder = '/media/peichao/PL_NHP_AE9/AF5/ISI_analysis/0. original ISI maps';
dataFolder = '/media/peichao/PL_NHP_AE9/AF5/ISI_analysis/2. Spatial relation';
ResultfileName = 'Achromatic OnOff Relation'; % 

% isiData=dir(fullfile(isimapFolder, '*_meta.mat'));
% image = load(fullfile(isiData.folder, isiData.name));
% isiMap = {image.ISImeta.map};
% 
% Image1 = isiMap(ismember({image.ISImeta.mapName},'Liso'));
% Image2 = isiMap(ismember({image.ISImeta.mapName},'Miso'));
% Image3 = isiMap(ismember({image.ISImeta.mapName},'Siso'));

Image1 = imread(fullfile(dataFolder,'Liso_aligned.tif'));  % Should be 2D array in uint8; I use ImageJ to convert image from rgb to 8bits.
Image2 = imread(fullfile(dataFolder,'Miso.tif'));
Image3 = imread(fullfile(dataFolder,'Siso_aligned.tif'));
Image4 = imread(fullfile(dataFolder,'AF5_Left_u000_045_WB_freq1_1st_0_Clip1.5_9_Phase_Filted.tif'));

areaMask = logical(imread(fullfile(dataFolder,'AF5_Left_LMSA_V1_DomainMask.tif')));
bvMask = logical(mean(imread(fullfile(dataFolder,'V1_BV.png')),3));
areaMask = areaMask .* bvMask;

%%
% clear
% aniName = 'AF3';
% hemiSph = 'Left';
% area = 'V1'; % 'V1' or 'V2'
% dataFolder = 'J:\AF3\ISI_analysis\2. Spatial relation';
% ResultfileName = 'OnOff Relation'; % 
% 
% Image1 = imread(fullfile(dataFolder,'AF3_Left_u000_025_Liso_freq1_2nd_31_Clip1.5_24_Phase_Filted.tif'));  % Should be 2D array in uint8; I use ImageJ to convert image from rgb to 8bits.
% Image2 = imread(fullfile(dataFolder,'AF3_Left_u000_026_Miso_freq1_1st_0_Clip1.5_32_Phase_Filted.tif'));
% Image3 = imread(fullfile(dataFolder,'AF3_Left_u000_027_Siso_freq1_2nd_0_Clip1.5_21_Phase_Filted.tif'));
% 
% areaMask = logical(imread(fullfile(dataFolder,'AF3_Left_LMS_V1_FinalMask.tif')));

%%
bRange = [0 10000];  % Range
%%
resultFolder = strcat(dataFolder, '/', ResultfileName, '_', area, '/');
fileName = strcat(aniName, '_', hemiSph,'_', area, '_', ResultfileName, '_');

if ~isfolder(resultFolder)
    mkdir(resultFolder);    
end 
cd(dataFolder)


%%

Image1 = norm_to_uint8(OIClip(double(Image1), 1, 2)); % Smoothing image
Image2 = norm_to_uint8(OIClip(double(Image2), 1, 2)); % Smoothing image
Image3 = norm_to_uint8(OIClip(double(Image3), 1, 2)); % Smoothing image
Image4 = norm_to_uint8(OIClip(double(Image4), 1, 2)); % Smoothing image
% tileNum = [8, 8];   % The total number of image tiles is equal to M*N.
% ClipLimit = 0.04;   % 0-1. higher number results in more contrast.
% Image1 = adapthisteq(Image1,'NumTiles',tileNum,'Distribution','Exponential', 'ClipLimit', ClipLimit);
% Image2 = adapthisteq(Image2,'NumTiles',tileNum,'Distribution','Exponential', 'ClipLimit', ClipLimit);
% Image3 = adapthisteq(Image3,'NumTiles',tileNum,'Distribution','Exponential', 'ClipLimit', ClipLimit);
%% Filtering
Image1 = imgaussfilt(Image1,5); % Smoothing image
Image2 = imgaussfilt(Image2,5); % Smoothing image
Image3 = imgaussfilt(Image3,5); % Smoothing image
Image4 = imgaussfilt(Image4,5); % Smoothing image
% imshow(Image1)
% figure
% imshow(Image2)
% figure
% imshow(Image3)

%% Generate matrix to store pixel value of the same position from two images

result.L = zeros(1, nnz(areaMask));
result.M = result.L;
result.S = result.L;
result.A = result.L;
n = 1;
for i = 1:size(Image1, 1)
    for j = 1:size(Image1, 2)
        if areaMask(i, j) == 1
            result.L(n) = Image1(i, j);
            result.M(n) = Image2(i, j);
            result.S(n) = Image3(i, j);
            result.A(n) = Image4(i, j);
            n = n+1;
        end
    end
end

%% Normalization to [-1 1]
result.Lscale = 2 * mat2gray(result.L) - 1;
result.Mscale = 2 * mat2gray(result.M) - 1;
result.Sscale = 2 * mat2gray(result.S) - 1;
result.Ascale = 2 * mat2gray(result.A) - 1;

%% Histogram plotting parameters
% dotSize = 20;
% dotTransparency = 1;

% Histogram 
nbins = 50;

%
pos1 = [0.2 0.2 0.7 0.7];
pos2 = [0.2 0.1 0.7 0.1];
pos3 = [0.10 0.2 0.1 0.7];

% X Y label
labelFront = 30;
xlabelPosition = [0 -1.3];
ylabelPosition = [-1.3 0];

% X Y axis
aRange = [-1 1];
% bRange = [0 6000];  % Range
axThickness = 6.0;
lnThickness = 6.0;
% atickLabel = {'','4','8'};
% btickLabel = {'0    ','4','8'};
atickLabel = {};
btickLabel = {};
% atickLabel = {'','5+E3','1+E4'};
% btickLabel = {'0    ','5+E3','1+E4'};

% X Y ticks
tickFront = 25;
x_ticks = -1:0.5:1;
y_ticks = -1:0.5:1;
x1_ticks = 0:bRange(2)/2:bRange(2);
y1_ticks = 0:bRange(2)/2:bRange(2);

% axis & line color
axColor = [0 0 0];
lnColor = [0 0 0];

% Colorbar
cbPosition = [0.91 0.2 0.02 0.7];
cbFontSize = 20;
cbColor = [0.5 0.5 0.5];

% Calculate the normalized histocounts/proportion
totalPixel = nnz(areaMask);
result.pL = histcounts(result.Lscale,nbins)./totalPixel;
result.pM = histcounts(result.Mscale,nbins)./totalPixel;
result.pS = histcounts(result.Sscale,nbins)./totalPixel;
result.pA = histcounts(result.Ascale,nbins)./totalPixel;
result.totalPixel = totalPixel;
%% Plot histogram2 
for rp = 1:3
    
    switch rp

        case 1 
            a = result.Lscale;
            b = result.Ascale;
            xName = 'L Cone';
            yName = 'Achromatic';
            imgName = strcat(fileName, 'LAchromatic');
        case 2
            a = result.Mscale;
            b = result.Ascale;
            xName = 'M Cone';
            yName = 'Achromatic';
            imgName = strcat(fileName, 'MAchromatic');
        case 3
            a = result.Sscale;
            b = result.Ascale;
            xName = 'S Cone';
            yName = 'Achromatic';
            imgName = strcat(fileName, 'SAchromatic');
    end
    
    % Initilize figure
    f = figure;
    f.InnerPosition = [100 10 1300 1300];  % Define drawable region
    colormap jet

    % First subplot in the figure
    ax1 = subplot('Position',pos1);
    % scatter(result.Lscale, result.Mscale, dotSize, dotColor1,'filled', 'MarkerFaceAlpha', dotTransparency)  % Plot as dots
    histogram2(a, b, nbins, 'DisplayStyle','tile','ShowEmptyBins','off');
    normN = rot90(histcounts2(a,b,nbins)./totalPixel);
    orignN = rot90(histcounts2(a,b,nbins));

    samNperQua = nbins^2/4;
    samMatrix = zeros(samNperQua,4);
    samMatrix(:,1) = reshape(orignN(1:25, 26:50), [],1);
    samMatrix(:,2) = reshape(orignN(1:25, 1:25), [],1);
    samMatrix(:,3) = reshape(orignN(26:50, 1:25), [],1);
    samMatrix(:,4) = reshape(orignN(26:50, 26:50), [],1);
    
    samMatrixp(:,1) = reshape(normN(1:25, 26:50), [],1);
    samMatrixp(:,2) = reshape(normN(1:25, 1:25), [],1);
    samMatrixp(:,3) = reshape(normN(26:50, 1:25), [],1);
    samMatrixp(:,4) = reshape(normN(26:50, 26:50), [],1);
    % 1st & 3rd quadrant
    q13 = sum(cat(1,samMatrix(:,1), samMatrix(:,3)));
    q13p = sum(cat(1,samMatrixp(:,1), samMatrixp(:,3)));
    q13pd = cat(1,samMatrixp(:,1), samMatrixp(:,3));
    % 2nd & 4th quadrant
    q24 = sum(cat(1,samMatrix(:,2), samMatrix(:,4)));
    q24p = sum(cat(1,samMatrixp(:,2), samMatrixp(:,4)));
    q24pd = cat(1,samMatrixp(:,2), samMatrixp(:,4));
%  Using the Fisher's exact test (fishertest in Matlab)
    xtable = [q13,q24;q24,q13];
    [h,p_ft] = fishertest(xtable);
    
 % Using the K-S test (kstest2 in Matlab)
    qall = cat(2,q13pd,q24pd);
    aa = hist(qall, 20);
    [h,p_ks] = kstest2(aa(:,1),aa(:,2));


    if rp == 1
        result.pLA = normN;
        result.pLA_qua = qall;
        result.pLA_pft = p_ft;
        result.pLA_pks = p_ks;
        result.pLAsamesign = q13p;
        result.pLAdiffsign = q24p;
    elseif rp == 2
        result.pMA = normN;
        result.pMA_qua = qall;
        result.pMA_pft = p_ft;
        result.pMA_pks = p_ks;
        result.pMAsamesign = q13p;
        result.pMAdiffsign = q24p;
    elseif rp == 3
        result.pSA = normN;
        result.pSA_qua = qall;
        result.pSA_pft = p_ft;
        result.pSA_pks = p_ks;
        result.pSAsamesign = q13p;
        result.pSAdiffsign = q24p;
    end
    
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
    hL = histogram(a, nbins);

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
    hM = histogram(b, nbins, 'Orientation','horizontal');
    
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
    savefig([resultFolder, imgName, '.fig'])
    saveas(gcf, [resultFolder, imgName, '.png'])
    saveas(gcf, [resultFolder, imgName, '.svg'])
end

save([resultFolder, fileName, 'result.mat'], 'result', '-v7.3');
close all;
%% Save current script
Version='2';
fileNameAndLocation='/home/peichao/Dropbox/Github/ISI_analysis_NHP/ISI_analysis_Peichao/Matlab_data_process/12_Spatial_pattern';
% fileNameAndLocation = mfileName('fullpath');   % mfileName('fullpath') does not work in current version!!
% fileNameAndLocation = which(mfileName);
newbackup= strcat(resultFolder,strcat('SpatialRelation_ConeAchromatic_',Version, '.m'));
currentfile=strcat(fileNameAndLocation, '/SpatialRelation_ConeAchromatic.m');
copyfile(currentfile,newbackup);

