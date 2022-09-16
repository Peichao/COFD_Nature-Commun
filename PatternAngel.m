

% Peichao Li

% Plot the histogram of cone ON/OFF map based on the normalized pixel value.
% Plot the bivariate histogram of ON/OFF map of two types of cones, e.g. L vs M.

%% Inputs
clear
aniName = 'AF4';
hemiSph = 'Left';
area = 'V1'; % 'V1' or 'V2'
dataFolder = 'O:\AF4\ISI_analysis\7. Spatial correlation of COFD Hue\OnOff Relation_V1_angel12';
ResultfileName = 'Hue_Cone_Angle'; % 
fileName = strcat(aniName, '_', hemiSph,'_', area, '_', ResultfileName);
%% Make folders
resultFolder = strcat(dataFolder, '\', ResultfileName, '_', area, '\');
fileName = strcat(aniName, '_', hemiSph,'_', area, '_', ResultfileName);

if ~isfolder(resultFolder)
    mkdir(resultFolder);
end 
cd(dataFolder);

%% Read in images

imgFiles = dir(strcat(dataFolder, '\', '*.png'));
numFiles = size(imgFiles, 1);
hueAngle = linspace(0,330, numFiles/3);
% hueAngle = [0 30 60 75 80 83 90 95 101 105 110 120 150 180 210 240 255 260 263 270 275 281 285 290 300 330];
result.Lhue = zeros(1, numFiles/3);
result.Mhue = zeros(1, numFiles/3);
result.Shue = zeros(1, numFiles/3);

for k = 1:numFiles

        imgFiles(k).name
  
        img = imread(fullfile(imgFiles(k).folder, imgFiles(k).name));  % 8bit CO image that there is no holes 

        img = img(2:end-1,2:end-1);
        img = imgaussfilt(img,30); % Smoothing image, filling small holes
        sz=(size(img,1)-1)/2;
        img1 = img(1:sz,:);
        img2 = img(sz+1:end,:);
        img1 = mat2gray(img1);
        img2 = mat2gray(img2);
        img = cat(1,img1,img2);
        
        
%         img = adapthisteq(img,'clipLimit',0.01, 'NumTiles',[16 16],'Distribution','Exponential');
%         [xx,yy] = meshgrid(-sz:sz,-sz:sz);   % PL: sholud be the same size as hartley space (-kxmax:kxmax, -kymax:kymax)
%         zz = xx+1i*yy;    
%         za = abs(zz).*exp(1i*angle(zz)*2);
%         za=flip(round(rad2deg(angle(zz)))+180,2);
%         za=mod(za,180)+1;
%         zzm = mean(za(id));
%         zzm = abs(zzm)*exp(1i*angle(zzm)/2);
%         ori = rad2deg(angle(zzm))+90;
        
        
        % estimate ori
        bw = img>(max(img(:))*.90);
        [row,col] = find(bw);
        
        id1 = find(row < sz);
        row1=mean(row(id1));
        col1=mean(col(id1));
        
        id2 = find(row > sz);
        row2=mean(row(id2));
        col2=mean(col(id2));
        ori = atand((row2-row1)/(col1-col2));
        ori = mod(ori+180,180);
        
        
        imgName = split(imgFiles(k).name, '_');       
        ang = imgName{4};
        cone = imgName{6}(1);
        idx = find(hueAngle == str2double(ang));        
        
        if strcmp(cone, 'L')
            result.Lhue(idx) = ori;
        elseif strcmp(cone, 'M')
            result.Mhue(idx) = ori;
        elseif strcmp(cone, 'S')
            result.Shue(idx) = ori;
        end  
end

%%
save([resultFolder, fileName, '_result.mat'], 'result', '-v7.3');

Lhue = result.Lhue;
Mhue = result.Mhue;
Shue = result.Shue;
fL = fit(hueAngle', Lhue', 'smoothingspline');
fM = fit(hueAngle', Mhue', 'smoothingspline');
fS = fit(hueAngle', Shue', 'smoothingspline');



Lwidth = 1;
% X Y axis
aRange = [-10 330];
bRange = [40 140];  % Range
axThickness = 1;
lnThickness = 2;
atickLabel = {'0','30','60', '90','120','150', '180','210','240','270','300','330'};
btickLabel = {'0    ','3000','6000'};

% X Y ticks
tickFront = 30;
x_ticks = 0:30:330;
y_ticks = 40:20:140;
dotSize = 5;
% axis & line color
axColor = [0 0 0];
lnColor = [0 0 0];

f=figure;
f.InnerPosition = [10 10 3000 300];  % Define drawable region

L1=plot(fL);
hold on;
S1=scatter(hueAngle', Lhue', dotSize, 'MarkerEdgeColor', [1 0 0],...
    'MarkerFaceColor', [1 0 0], 'LineWidth',0.5);
 
hold on;
L2=plot(fM);
hold on;
S2=scatter(hueAngle', Mhue', dotSize, 'MarkerEdgeColor', [0 1 0],...
    'MarkerFaceColor', [0 1 0], 'LineWidth',0.5);

hold on;
L3=plot(fS);
hold on;
S3=scatter(hueAngle', Shue', dotSize, 'MarkerEdgeColor', [0 0 1],...
    'MarkerFaceColor', [0 0 1], 'LineWidth',0.5);

xlim(aRange)
ylim(bRange)
xticks(x_ticks)
yticks(y_ticks)
xticklabels(atickLabel)
% xticklabels({})
% yticklabels({})
set(gca,'FontSize',tickFront)
set(gca, 'linewidth',axThickness)

set(L1,'color','r','LineWidth', Lwidth)
set(L2,'color','g','LineWidth', Lwidth)
set(L3,'color','b','LineWidth', Lwidth)

grid off
box off


% savefig([resultFolder, fileName, '.fig'])
export_fig (gcf, [resultFolder, fileName, '.png'])
fig2svg([resultFolder, fileName, '.svg'])



%% Save current script
Version='1';
fileNameAndLocation='C:\Users\lyr19\Dropbox\Github\ISI_analysis_NHP\ISI_analysis_Peichao\Matlab_data_process\12_Spatial_pattern';
% fileNameAndLocation = mfileName('fullpath');   % mfileName('fullpath') does not work in current version!!
% fileNameAndLocation = which(mfileName);
newbackup= strcat(resultFolder,strcat('PatternAngel_',Version, '.m'));
currentfile=strcat(fileNameAndLocation, '\PatternAngel.m');
copyfile(currentfile,newbackup);

