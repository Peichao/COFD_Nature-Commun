clear
dataFolder = '/media/peichao/PL_NHP2_AE7AE8AF4/AF4/COFD_manuscript/';
cd(dataFolder);

matFile = dir(fullfile(dataFolder, '*V1_ConeNoncone_result.mat'));
numFile = length(matFile);
COFD=zeros(numFile,20);
nonCOFD=zeros(numFile,20);
pvalue=zeros(numFile,1);
COFDmedInt = zeros(numFile,1);
nonCOFDmedInt = zeros(numFile,1);
for i = 1:numFile
%     i=1
    load(fullfile(matFile(i).folder, matFile(i).name))
%     COFD(i,:) = result.ConeVal;
%     nonCOFD(i,:) = result.nonConeVal;
    COFDmedInt(i) = median(result.cone);
    nonCOFDmedInt(i) = median(result.noncone);
%     [H,p,CI,STATS] = ttest(COFD(i,:)', nonCOFD(i,:)', 0.05, 'right', 1);
%     [p,h,stats] = signrank(COFD(i,:)', nonCOFD(i,:)', 'tail', 'right');
    [p,h,stats] = ranksum(result.cone', result.noncone', 'tail', 'right');
%     [h,p,ks2stat] = kstest2(COFD(i,:)', nonCOFD(i,:)');
%     [h,p,stat] = signtest(COFD(i,:)', nonCOFD(i,:)','Tail','right');
    pvalue(i) = p;
    
end

% aa=cat(1,COFD,nonCOFD)';
% bb={'COFD','COFD','COFD','COFD','COFD','nonCOFD','nonCOFD','nonCOFD','nonCOFD','nonCOFD'};
% [p,tbl,stats]=anova-2(aa,bb);


figure
plot(mean(COFD,1))
hold on
plot(mean(nonCOFD,1))

% [h,p,ci,stats]=ttest2(COFD',nonCOFD');



% [h,p,stat] = signtest(COFDmedInt, nonCOFDmedInt,'Tail','right')

[h,p,stat] = signrank(COFDmedInt, nonCOFDmedInt,'Tail','right')

% [p,h,stats] = kstest2(result.cone', result.noncone', 'tail', 'larger')
