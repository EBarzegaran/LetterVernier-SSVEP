
clear; clc;

PATH = '/Users/kohler/code/git/mrC';
addpath(genpath(PATH));
path2 = '/Users/kohler/code/git/gardner/mrTools/mrUtilities/mrFlatMesh';
addpath(genpath(path2));
%% Read the inverses
Path = '/Volumes/svndl/mrC_Projects/VernierLetters/source';
InvName = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangkgsROIsCorr.inv';
[Inverses,SubIDInv] = mrC.Simulate.ReadInverses(Path,InvName);
% Read  morph maps
for s = 1:numel(SubIDInv)
    mapMtx{s} = makeDefaultCortexMorphMap(SubIDInv{s},SubIDInv{2});
    InvMapped{s} = mapMtx{s}*Inverses{s}.';
end
% mapMtx*Inverses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SubIDInv{strcmp(SubIDInv,'skeri0003')} = 'nl-0033';
SubIDInv{strcmp(SubIDInv,'skeri0004')} = 'nl-0034';

%% Load ROIs
%Path = '/Users/babylab/Documents/Elham/simulation/mrCSimulate/Examples/Dataset/FwdProject';
[ROIs,subIDROI] = mrC.Simulate.GetRoiClass(Path);
% correct the labels
ROIs{strcmp(subIDROI,'skeri0003')}.subID = 'nl-0033';
subIDROI{strcmp(subIDROI,'skeri0003')} = 'nl-0033';
ROIs{strcmp(subIDROI,'skeri0004')}.subID = 'nl-0034';
subIDROI{strcmp(subIDROI,'skeri0004')} = 'nl-0034';
Ords = cell2mat(arrayfun(@(x) find(strcmp(subIDROI,SubIDInv{x})),1:numel(SubIDInv),'uni',false));
ROIs = ROIs(Ords);

% Prepare Wang and KGS atlases
WangROIs = cellfun(@(x) x.getAtlasROIs('wang'),ROIs,'uni',false); 
KgsROIs = cellfun(@(x) x.getAtlasROIs('kgs'),ROIs,'uni',false); 
WKROIs = arrayfun(@(x) WangROIs{x}.mergROIs(KgsROIs{x}),1:numel(ROIs),'uni',false); 
WKROIsM = cellfun(@(x) x.ROI2mat(20484),WKROIs,'uni',false);
ROILabel = WKROIs{1}.getFullNames('noatlas');
InvROIs = arrayfun(@(x) Inverses{x}*WKROIsM{x}./(sum(WKROIsM{x})+eps),1:numel(Inverses),'uni',false);
%% Load in the odd and even waves
Task = {'Letter','Vernier'};
Harms = {'Odd','Even','All'};
load(fullfile('ResultData','EvenOddWaves.mat'));
SubSelect = cell2mat(arrayfun(@(x) find(strcmp(SubIDs,SubIDInv{x})),1:numel(SubIDInv),'uni',false));

for ts = 1:numel(Task)
    for h = 1:numel(Harms)
        for cond = 1:5
            RWave_temp = squeeze(RWave_all.(Task{ts}).(Harms{h})(:,:,cond,SubSelect));
            
            RWaveSource_temp = arrayfun(@(x) squeeze(RWave_temp(:,:,x))*InvMapped{x}',1:numel(SubSelect),'uni',false);
            RWaveSource.(Task{ts}).(Harms{h})(:,:,:,cond) = cat(3,RWaveSource_temp{:});  
            
            RWaveROI_temp = arrayfun(@(x) squeeze(RWave_temp(:,:,x))*InvROIs{x},1:numel(SubSelect),'uni',false);
            RWaveROI.(Task{ts}).(Harms{h})(:,:,:,cond) = cat(3,RWaveROI_temp{:});  
        end
    end
end

%%
cond = 3;
ts = 1;
h = 1;
FIG = figure;
subplot(2,1,1),imagesc(mean(RWave_all.(Task{ts}).(Harms{h})(:,:,cond,SubSelect),4))
set(gca,'ytick',1:10:140,'yticklabel',2.81:28.1:(2.81*140));

subplot(2,1,2),imagesc(mean(RWaveROI.(Task{ts}).(Harms{h})(:,:,:,cond),3))
set(gca,'xtick',1:2:64,'xticklabel',ROILabel(1:2:64));
xtickangle(90)
set(gca,'ytick',1:10:140,'yticklabel',2.81:28.1:(2.81*140));
set(FIG,'unit','inch','position',[5 5 10 8]);

%
% FIG = figure;
% Matx =mean(RWaveROI.(Task{ts}).(Harms{h})(:,:,:,cond),3);
% FMatx = fft(Matx);
% subplot(2,1,1),bar(abs(FMatx(h+1,1:2:end))');
% set(gca,'xtick',1:32,'xticklabel',ROILabel(1:2:64));
% xtickangle(90);xlim([1 32])
% subplot(2,1,2),bar(unwrap(angle((FMatx(h+1,1:2:end))))')
% set(gca,'xtick',1:32,'xticklabel',ROILabel(1:2:64));
% xtickangle(90);xlim([1 32])
% set(FIG,'unit','inch','position',[15 5 10 4]);

%%
SMdata = smoothdata(mean(RWaveROI.(Task{ts}).(Harms{h})(:,:,:,cond),3),1);
[~,MI] = max(SMdata);
Times = 2.381:2.381:(140*2.381);
figure,
subplot(2,1,1),imagesc(SMdata);
subplot(2,1,2),plot(Times(MI(1:2:end)));xlim([1 31])
set(gca,'xtick',1:32,'xticklabels',ROILabel(1:2:end)); xtickangle(90)
%%
figure,mrC.Simulate.VisualizeSourceRoi2('nl-0014',[],'Wang',45:50)
figure,mrC.Simulate.VisualizeSourceRoi2('nl-0014',[],'kgs',[],[],'L')

%%
Temp = cat(3,InvROIs{:});
FIG = figure;
for i =1:30
   subplot(5,6,i),mrC.plotOnEgi(mean(mean(Temp(:,(i*2)-1:i*2,:),3),2)); axis tight
   title(ROILabel{(i*2)-1})
end
set(FIG,'unit','inch','position',[5 5 20 20]);
%% convert to source space
for s = 1:numel(SubIDInv)
    SourceWaveT(:,:,s) = squeeze(Test_Wave(:,:,s))*InvMapped{s}.';
end

% videos 
cond = 3;
M = max(max(abs(mean(SourceWaveT(:,:,:),3))));
vidfile = VideoWriter(['Figures/' Task{2} '_Cond' num2str(cond) '_Source.mp4'],'MPEG-4');
vidfile.FrameRate = 10;
open(vidfile);
for i = 1:140
    figure; mrC.Simulate.VisualizeSourceData('nl-0014',mean(SourceWaveT(i,:,:),3),[],'coolhotcortex');
    caxis([-M/2 M/2]);
    pause(.1);
    F(i)= getframe(gcf);
    writeVideo(vidfile,F(i));
    close;
end
close(vidfile);
close all