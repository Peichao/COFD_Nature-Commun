% Peichao
% pull CO intensities based on 2P cell mask (.segment)


%% User input, file location
clear
clc
disk = 'X:';
aniName = '';
unit = ''; % 'V1' or 'V2'
dataRoot = strcat(disk, '\',aniName,'',unit, '');
dataFolder = strcat(dataRoot, 'Multiplane\5. Hue preference & CO intensity\');

colorSpace = 'DKL';
hueaucThres = 0.8;
labelFront = 25;

ResultfileName = ''; % 
resultName = strcat(aniName,'_U', unit,'_hueaucThres', num2str(hueaucThres), '_', ResultfileName);


%% Load data
cd(dataFolder)
COImage = uint8(mean(imread(fullfile(dataFolder,'*.png')),3));  % 8bit CO image that there is no holes
segFiles = dir(fullfile(dataFolder,'**\','*.segment'));
segNum = length(segFiles);

%% Filtering
COImage = imgaussfilt(COImage,40); % Smoothing image
% Image normalization [0 1]
COImage = mat2gray(COImage); % CO intensity, small number means closer to CO blob


%% STA proportion & Plot correlation figures

for ii = 1:segNum

%     ii=1
    segName = segFiles(ii).name
    planeId = ii-1
    
    Hue = readtable(fullfile(dataRoot, 'DataExport\', strcat(aniName, '_', unit,'_00', num2str(planeId), '_sum.csv')));
    seg=load(fullfile(segFiles(ii).folder,segName), '-mat'); % load segmentation
    
    cellNum = length(seg.vert);
    cellMask = seg.mask;
 
    ishue = Hue.visResp & ((Hue.hueaxauc>hueaucThres) | (Hue.huediauc>hueaucThres));
    
    result.(strcat('plane',num2str(planeId))).('CO') = zeros(1, sum(ishue));
    result.(strcat('plane',num2str(planeId))).('Hue') = zeros(1,sum(ishue));
    
    count = 1;

    for cell = 1:cellNum
        
        if ishue(cell)
            result.(strcat('plane',num2str(planeId))).('CO')(count) = mean(COImage(find(cellMask == cell)));
            result.(strcat('plane',num2str(planeId))).('Hue')(count) = Hue.maxhue(cell);
            count = count + 1;
        end
        
    end
end

%% Summary of two planes

result.Hue = [result.plane0.Hue,result.plane1.Hue];
result.CO = [result.plane0.CO,result.plane1.CO];

hueNum = length(unique(result.Hue));
for ii = 1:hueNum
    
    result.distr.(strcat(colorSpace,num2str(ii))) = result.CO(result.Hue == 360/hueNum*(ii-1));
end
save(fullfile(dataFolder,  [resultName, '.mat']), 'result', '-v7.3');

%% Organize data for bar plotting
edges = [0:0.25:1];
fns = fieldnames(result.distr);
barmtr = zeros(length(fns), length(edges)-1);
barName = {'0', '30', '60', '90','120','150','180', '210', '240', '270','300','330'};

for ii = 1:length(fns)

    data = result.distr.(fns{ii});
    [N,ed]=histcounts(data,edges);

    v = matlab.lang.makeValidName(fns{ii});
    eval([v ' = N./length(data)']);
    
    if strcmp(fns{ii}, 'DKL1')
        jj=1;
    elseif strcmp(fns{ii}, 'DKL2')
        jj=2;
    elseif strcmp(fns{ii}, 'DKL3')
        jj=3;
    elseif strcmp(fns{ii}, 'DKL4')
        jj=4;
    elseif strcmp(fns{ii}, 'DKL5')
        jj=5;
    elseif strcmp(fns{ii}, 'DKL6')
        jj=6;
    elseif strcmp(fns{ii}, 'DKL7')
        jj=7;
    elseif strcmp(fns{ii}, 'DKL8')
        jj=8;
    elseif strcmp(fns{ii}, 'DKL9')
        jj=9;
    elseif strcmp(fns{ii}, 'DKL10')
        jj=10;
    elseif strcmp(fns{ii}, 'DKL11')
        jj=11;
    elseif strcmp(fns{ii}, 'DKL12')
        jj=12; 
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


savefig([dataFolder, resultName, 'CO.fig'])
saveas(gcf, [dataFolder, resultName, 'CO.png'])
saveas(gcf, [dataFolder, resultName, 'CO.svg'])
hold off;





