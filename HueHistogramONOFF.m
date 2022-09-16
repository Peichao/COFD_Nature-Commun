
% Plot the pixel value of CO image in the regions of ISI response as histogram.


%% Inputs
clear;
aniName = '';
hemiSph = '';
area = ''; % 
dataFolder = '';
coneType = 'L'; % L, M, S, cone

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
end

resultFolder = strcat(dataFolder, '\', ResultfileName, '_', area, '\');

if ~isfolder(resultFolder)
    mkdir(resultFolder);    
end 
cd(dataFolder)

% areaMask = logical(mean(imread(fullfile(dataFolder,'V1 mask_crop.tif')),3));  % Mask of ON response region
areaMaskOn = logical(mean(imread(fullfile(dataFolder, maskOnName)),3));  % Mask of ON response region
areaMaskOff = logical(mean(imread(fullfile(dataFolder, maskOffName)),3)); % Mask of OFF response region
areaMaskBloodvessel = logical(imread(fullfile(dataFolder,'*.bmp'))); % Mask of OFF response region

hueFiles = dir(strcat(dataFolder, '\', '*u000*.tif'));
numFiles = size(hueFiles, 1);

 for k = 1:numFiles
        disp(strcat('Analyzing Hue ', num2str(k)));
        HueImage = imread(fullfile(hueFiles(k).folder, hueFiles(k).name));  % 8bit CO image that there is no holes 
        hueName = hueFiles(k).name(19:20);
        fileName = strcat(aniName, '_', hemiSph,'_', area, '_', ResultfileName, '_', hueName, '_');
        %% Filtering
        HueImage = imgaussfilt(HueImage,5); % Smoothing image, filling small holes
        % figure 
        % imshow(HueImage)
        %% Image normalization [0 1]
        HueImage = mat2gray(HueImage)*2-1;
        % HueImage = abs(HueImage-1);   % CO intensity, higher number means closer to CO blob

        %% Store pixel value selected from CO image based on the ON OFF masks of cones
        areaMaskOn = areaMaskOn.*areaMaskBloodvessel;
        areaMaskOff = areaMaskOff.*areaMaskBloodvessel;

        % result.Hue = zeros(1, nnz(areaMask));
        result.On = zeros(1, nnz(areaMaskOn));
        result.Off = zeros(1, nnz(areaMaskOff));

        % Region in areaMask
        % n = 1;
        % for i = 1:size(HueImage, 1)
        %     for j = 1:size(HueImage, 2)
        %         if areaMask(i, j) == 1
        %             result.Hue(n) = HueImage(i, j);
        %             n = n+1;
        %         end
        %     end
        % end

        % ON regions
        n = 1;
        for i = 1:size(HueImage, 1)
            for j = 1:size(HueImage, 2)
                if areaMaskOn(i, j) == 1
                    result.On(n) = HueImage(i, j);
                    n = n+1;
                end
            end
        end

        % OFF regions
        n = 1;
        for i = 1:size(HueImage, 1)
            for j = 1:size(HueImage, 2)
                if areaMaskOff(i, j) == 1
                    result.Off(n) = HueImage(i, j);
                    n = n+1;
                end
            end
        end

        save([resultFolder, fileName, 'result.mat'], 'result', '-v7.3');

        %% Plot
        BinWidth = 0.1;

        switch hueName
            case 'RB'
                color1 = [0.5 0 0];
                color2 = [0 0 0.5];
            case 'RG'
                color1 = [0.5 0 0];
                color2 = [0 0.5 0];
            case 'YB'
                color1 = [0.5 0.5 0];
                color2 = [0 0 0.5];
        end
      
        
        figure
        % h0 = histogram(result.Hue);
        % hold on
        h1 = histogram(result.On);  % ON phase 
        hold on
        h2 = histogram(result.Off);  % OFF phase

        % h0.Normalization = 'probability';
        % h0.DisplayStyle = 'stairs';
        % h0.EdgeColor = [0 0 0];
        % h0.LineWidth = 1;
        % h0.LineStyle = '--';
        % % h0.FaceColor = [0.2 0.2 0.2];
        % % h0.FaceAlpha = 1;
        % h0. BinLimits = [0 1];
        % h0.BinWidth = BinWidth;

        h1.Normalization = 'probability';
        h1.EdgeColor = 'none';
        h1.FaceColor = color1;  % This color is the color of hue
        h1.FaceAlpha = 0.5;
        h1. BinLimits = [-1 1];
        h1.BinWidth = BinWidth;

        h2.Normalization = 'probability';
        h2.EdgeColor = 'none';
        h2.FaceColor = color2;  % Gray
        h2.FaceAlpha = 0.5;
        h2. BinLimits = [-1 1];
        h2.BinWidth = BinWidth;

        % xlabel('CO Intensity', 'FontSize', labelFront)
        % ylabel('Proportion of Pixel Number', 'FontSize', labelFront)

        % legend([h1 h2], legend1, legend2, 'Location', 'northwest')
        % title([aniName, ' ', hemiSph, ' ', area, ' ', ResultfileName, ' CO'])


        savefig([resultFolder, fileName, 'fig.fig'])
        saveas(gcf, [resultFolder, fileName, 'fig.png'])
 end


