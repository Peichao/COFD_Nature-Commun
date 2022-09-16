
% Plot the pixel value of CO image in the regions of ISI response as histogram.


%% Inputs
clear
aniName = '';
hemiSph = '';
area = ''; % 'V1' or 'V2'
dataFolder = '';
ResultfileName = ''; % CO, Liso; Miso, Siso, RG, BY, RB, or other name

labelFront = 25;

COImage = imread(fullfile(dataFolder,'*.tif'));  % 8bit CO image that there is no holes 

% 
areaMask = logical(mean(imread(fullfile(dataFolder, '*.png')),3));  % Mask of V1
areaMaskbv = logical(mean(imread(fullfile(dataFolder, '*.png')),3));  % Mask of V1

%%
resultFolder = strcat(dataFolder, '\', ResultfileName, '_', area, '\');
fileName = strcat(aniName, '_', hemiSph,'_', area, '_', ResultfileName, '_');

if ~isfolder(resultFolder)
    mkdir(resultFolder);    
end 
cd(dataFolder)

%% Filtering
COImage = imgaussfilt(COImage,5); % Smoothing image, filling small holes

%% Image normalization [0 1]
COImage = mat2gray(COImage);  % CO intensity, small number means closer to CO blob

%% 

maskFiles = dir(strcat(dataFolder, '\coneMask\', '*.png'));
numFiles = size(maskFiles, 1);

areaMaskCone = zeros(size(areaMask));

for k = 1:numFiles

    maskFiles(k).name

    mask = imread(fullfile(maskFiles(k).folder, maskFiles(k).name));  % 8bit CO image that there is no holes 
    mask = logical(mean(mask,3));

    imgName = split(maskFiles(k).name, '.');       
    cone = imgName{1,1};

    if strcmp(cone, 'V1')
        countMask = mask;
    elseif strcmp(cone, 'WB')
        countMask = areaMaskbv .* mask;
    else
        countMask = areaMaskbv .* mask;
        areaMaskCone = areaMaskCone + mask;
    end

    idx = find(countMask);
    % Store pixel value selected from CO image based on the ON OFF masks of cones
    result.(cone) = COImage(idx);
end

areaMaskCone = logical(areaMaskCone) .* areaMaskbv;
idx = find(areaMaskCone);
result.Allcone = COImage(idx);
areaMaskNoncone = logical((1-areaMaskCone).*areaMaskbv);
idx = find(areaMaskNoncone);
result.Noncone = COImage(idx);

save([resultFolder, fileName, 'result.mat'], 'result', '-v7.3');

%%
edges = [0:0.25:1];
fns = fieldnames(result);
barmtr = zeros(length(fns), length(edges)-1);
barName = {'V1', 'Non Cone', 'Cone', 'Lon','Loff','Mon','Moff', 'Son', 'Soff', 'Aon','Aoff'};

for ii = 1:length(fns)

    data = result.(fns{ii});
    [N,ed]=histcounts(data,edges);

    v = matlab.lang.makeValidName(fns{ii});
    eval([v ' = N./length(data)']);
    
    if strcmp(fns{ii}, 'V1')
        jj=1;
    elseif strcmp(fns{ii}, 'Noncone')
        jj=2;
    elseif strcmp(fns{ii}, 'Allcone')
        jj=3;
    elseif strcmp(fns{ii}, 'Liso_on')
        jj=4;
    elseif strcmp(fns{ii}, 'Liso_off')
        jj=5;
    elseif strcmp(fns{ii}, 'Miso_on')
        jj=6;
    elseif strcmp(fns{ii}, 'Miso_off')
        jj=7;
    elseif strcmp(fns{ii}, 'Siso_on')
        jj=8;
    elseif strcmp(fns{ii}, 'Siso_off')
        jj=9;
    elseif strcmp(fns{ii}, 'WB_on')
        jj=10;
    elseif strcmp(fns{ii}, 'WB_off')
        jj=11; 
    end
    barmtr(jj,:) = eval(v);
end

%% Plot the ratio, nomalize by the distribution of CO intensity

f=figure; 
f.InnerPosition = [10 10 1000 600];  % Define drawable region

ba=bar(barmtr, 'stacked','FaceColor','flat', 'Barwidth',0.5);

ba(1).CData = [1 1 1]*0.1;
ba(2).CData = [1 1 1]*0.3;
ba(3).CData = [1 1 1]*0.7;
ba(4).CData = [1 1 1]*0.9;


hold on
hline = refline([0 0.5]);
hline.Color = 'k';
hline.LineStyle = '--';
hline.LineWidth = 1;


tickFront = 20;
axisThickness = 1;
ax = gca;
set(ax,'linewidth',axisThickness)
ax.FontSize = tickFront;

% xlabel({}, 'FontSize',labelFront)
ylabel('Distribution Proportion', 'FontSize', labelFront)
ylim([0 1])
yticks([0:0.25:1])
xticklabels(barName)
xtickangle(45)
% legend([ba(1), ba(2), ba(3), ba(4)], '25%', '50%','75%','100%','Location', 'northwest' )


savefig([resultFolder, fileName, 'CO.fig'])
saveas(gcf, [resultFolder, fileName, 'CO.png'])
saveas(gcf, [resultFolder, fileName, 'CO.svg'])
hold off;



