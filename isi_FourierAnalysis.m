%%
% Expt info.
clear
clc
animal = '';
hemiSph = '';
unit = '';
expt ='';  % 25 24 23 22 21 20 19 18 17 16
ResultFileName = ''; % eg.: Liso, Miso, Siso, WB, RG, YB, RB, GB
dataDrive = '';  %

% Stimulus parameters for Fourier analysis
freq = 1;  % 1: First harmonic response; 2: Second harmonic response
stimCycle = 10;   % cycles in stim 48cyc (12.5sec); 10cyc(60sec)
stimTimeperCycle = 60;  % 10sec(0.1Hz), 60sec(1/60hz) modulation in temporal domain
startPhase = '1st';   % '1st' or '2nd'
frameShift = 0; % 0-9
delay = 0;
%
if freq == 1
    stimCycle  = stimCycle-1; 
elseif freq == 2
    stimCycle  = stimCycle*2-1;
    stimTimeperCycle = stimTimeperCycle./2;  
end


% Normlization 
normFlag = 0;  % 0: no normalization; 1: normalization

% Spatial filtering after Fourier analysis not before, and only Highpass
filterFlag = 1;   % 0: no filtering; 1: filtering
% % HP = 400;  % high pass
tileNum = [8, 8];   % The total number of image tiles is equal to M*N.
ClipLimit = 0.04;   % 0-1. higher number results in more contrast.

% % sizedum = 2.5*HP;
% % H = -fspecial('gaussian', sizedum, HP);
% % H(round(sizedum/2),round(sizedum/2)) = 1+H(round(sizedum/2),round(sizedum/2));

%Clipping
clip_method=1;    %  0: no clipping;   1: clipping to +-SD (value) (on each side of the median);    2: clipping to +-SD (value) with mask (immask)   
clip_value=1.5;     %  sd value for clipping, usually 2 or 1.5

%% Find & Make folders
dataFolder = strcat(dataDrive, animal, '/ISI_data/', 'u', unit, '_e', expt, '/');
analyzerFolder = strcat(dataDrive, animal, '/AnalyzerFiles/', animal);

resultFolder = strcat(dataFolder, animal, '_', hemiSph, '_', ResultFileName, '_freq', num2str(freq), '_', startPhase, '_', num2str(frameShift), '/');
fileName = strcat(animal, '_', hemiSph, '_u', unit, '_', expt, '_freq', num2str(freq), '_',startPhase, '_');
mapName = strcat(animal, '_', hemiSph, '_u', unit, '_', expt, '_', ResultFileName, '_freq', num2str(freq), '_', startPhase, '_', num2str(frameShift), '_', 'Clip',num2str(clip_value));
if ~isfolder(resultFolder)
    mkdir(resultFolder);    
end 

%% Read data, Calculate frame size & number
movFiles = dir(strcat(dataFolder, '*.mj2'));
numFiles = size(movFiles, 1);

analyzerFile = dir(fullfile(analyzerFolder, strcat(animal, '_', unit, '_', expt, '.analyzer')));
[trials, numFiles, stimTime] = looper2(analyzerFile, numFiles);    % compare trial numbers in analyzer and folder to see if all trials are completed, them organize trial no. and parameter.

v = VideoReader(fullfile(movFiles.folder, movFiles.name));
vidWidth = v.Width;
vidHeight = v.Height;
% vidFrames = v.Duration * v.FrameRate;

cameraFrameTime = 1/v.FrameRate;
preStimTime = stimTime(1) -delay;   % sec
durStimTime = stimTime(2);   % sec

cameraFrameperCycle = stimTimeperCycle * v.FrameRate;
% numSample = stimCycle * cameraFrameperCycle;   % Number of samples/frames during stimulus presentation
cameraFramepreStim = v.FrameRate * preStimTime; % Number of samples/frames during pre-stimulus presentation
% cameraFramedurStim = v.FrameRate * durStimTime; % Number of samples/frames during stimulus presentation

%% Fourier analysis
cameraFrameSeq = cameraFrameTime:cameraFrameTime:stimTimeperCycle;
cameraFrameAng = cameraFrameSeq/stimTimeperCycle*2*pi;

imgPre = zeros(vidHeight, vidWidth);
for i=1:cameraFramepreStim
    img = double(read(v,i));
    imgPre = imgPre+img;
end
imgPre=imgPre/cameraFramepreStim;

acc = zeros(vidHeight, vidWidth);
imgAll = zeros(vidHeight, vidWidth);
if strcmp(startPhase,'2nd')
    f = round(cameraFramepreStim+cameraFrameperCycle*0.25);
    imaginary = 1j;
elseif strcmp(startPhase,'1st')
    f = round(cameraFramepreStim+cameraFrameperCycle*0.75);
    imaginary = -1j;   % Fourier Transform. Using negitive frequency, bright spot is ON response. PL
end

bw =  ones(size(imgPre));
h = waitbar(0,'Processing...');
frameNum = 1;

for cycNum = 1:stimCycle
    for k = 1:cameraFrameperCycle
        img = 4096*double(read(v, f+k+frameShift))./imgPre;  
        acc = acc + exp(imaginary*cameraFrameAng(k)).*img;   
        imgAll = imgAll+img;
        frameNum = frameNum +1;
    end
    f = f+cameraFrameperCycle;
    waitbar(cycNum/stimCycle,h);

    F0 = imgAll./ (f-cameraFramepreStim);
    acc1 = acc - F0*cycNum*sum(exp(-1j*cameraFrameAng)); %Subtract f0 leakage
    acc2 = 2*acc1 ./ (f-cameraFramepreStim);
    acc2(isnan(acc2))= 0;
    
    if normFlag == 1
        re = OIClip(real(acc2), 1, 1);
        im = OIClip(imag(acc2), 1, 2.5);
        imInput = cat(3,re, im);
        imOut =  ImNorm(imInput, bw);  % Image normalize by diviation
        re2 = imOut(:,:,1);
        im2 = imOut(:,:,2);
        acc2 = re2 + 1i*im2;
    end
    
    Phase = angle(acc2);  %range: -Pi to Pi
    Mag = abs(acc2); % >0
    result.phase{cycNum} = Phase;
    result.mag{cycNum} = Mag;

    imwrite(norm_to_uint8(OIClip(Phase, clip_method, clip_value)), [resultFolder, mapName, '_', num2str(cycNum), '_Phase.tif'])
    if filterFlag == 1
        imwrite(adapthisteq(norm_to_uint8(OIClip(Phase, clip_method, clip_value)), 'NumTiles',tileNum,'Distribution','Exponential', 'ClipLimit', ClipLimit), [resultFolder, mapName, '_', num2str(cycNum), '_Phase_Filted.tif'])
    end
end

result.F0 = F0;
result.allImg = imgAll;
result.preImg = imgPre;

save([resultFolder, fileName, 'result.mat'], 'result', '-v7.3');
imwrite(norm_to_uint8(OIClip(Mag, clip_method, clip_value)), [resultFolder, mapName, '_Mag.tif'])

figure
imshow(norm_to_uint8(Phase))
figure
imshow(adapthisteq(norm_to_uint8(Phase),'NumTiles',[8 8],'Distribution','Exponential'));
close(h);


%% Save results


