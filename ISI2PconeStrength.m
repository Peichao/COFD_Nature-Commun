% Peichao
% pull ISI signals based on mask (.segment)


%% User input, file location
clear
clc
disk = 'X:';
aniName = '';
unit = ''; % 'V1' or 'V2'
dataRoot = strcat(disk, '\',aniName,'\2P_analysis\U',unit, '\_Summary\');
dataFolder = strcat(dataRoot, 'Multiplane\2. Correlation of Cone weights based on ISI\');
ResultfileName = ''; % 

%% Load data
cd(dataFolder)
imgFiles = dir(fullfile(dataFolder,'*iso*.png'));
segFiles = dir(fullfile(dataFolder,'**\','*.segment'));
maskISI = logical(mean(imread(fullfile(dataFolder, strcat(aniName,'_U',unit,'_LMS_mask.png'))),3));  % Mask of LMS COFDs

numImg = length(imgFiles);
numSeg = length(segFiles);

for ii = 1: numImg
    v=genvarname(imgFiles(ii).name(1:4));
    eval([v ' = mean(imread(fullfile(imgFiles(ii).folder,imgFiles(ii).name)),3);'])
end

%% Filtering
L = imgaussfilt(Liso,20); % Smoothing image
M = imgaussfilt(Miso,20); % Smoothing image
S = imgaussfilt(Siso,20); % Smoothing image

for ii = 1:numSeg
    segName = segFiles(ii).name
    
    load(fullfile(segFiles(ii).folder,segName), '-mat'); % load segmentation
    
    resultName = strcat(segName(1:end-8),'_', ResultfileName);
    
    ncell = max(mask(:));
    
    result.L = zeros(ncell,1);
    result.M = result.L;
    result.S = result.L;
    result.COFD = result.L;
    
    for(jj=1:ncell)
        result.L(jj) = mean(L(find(mask==jj)));  
        result.M(jj) = mean(M(find(mask==jj)));
        result.S(jj) = mean(S(find(mask==jj)));
        if sum(maskISI(find(mask==jj)))>0
            result.COFD(jj) = 1;
        elseif sum(maskISI(find(mask==jj))) == 0
            result.COFD(jj) = 0;
        end
    end

    % Normalization to [-1 1]
    result.Lscale = 2 * mat2gray(result.L) - 1;
    result.Mscale = 2 * mat2gray(result.M) - 1;
    result.Sscale = 2 * mat2gray(result.S) - 1;

    save(fullfile(dataFolder,  [resultName, '.mat']), 'result', '-v7.3');
end

