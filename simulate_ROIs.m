% This script generates an example simulated EEG with SSVEP signals
% for this script to run correctly you need three paths:
    % mrCpath: The latest mrC package
    % ProjectPath: Pointing to the mrC project folder of an SPECIFIC 
                % subject (this version of simulation package is not for group level data)
    % AnatomyPath: pointing to the folder where anatomy data is ()including
                % freesurfer files are and meshes and ROIs are (check the Example folder)
                
% Elham Barzegaran 3/14/2018

%% Add latest mrC
clear;clc
mrCFolder = '/Users/kohler/code/git/mrC';
addpath(genpath(mrCFolder));


addpath(genpath('/Users/kohler/code/git/surfing'));% this tool can be found in github

ProjectPath = '/Volumes/svndl/mrC_Projects/VernierLetters/source';


%% simulated signal generation

[outSignal, FundFreq, SF]= mrC.Simulate.ModelSeedSignal('signalType','SSVEP','signalFreq',[3 3],'signalHarmonic',{[1],[1]},'signalPhase',{[pi/2], [pi]});

%% load ROIs
Rois = mrC.Simulate.GetRoiClass(ProjectPath);
kgsRois = cellfun(@(x) x.getAtlasROIs('kgs'),Rois,'UniformOutput',false);% KGS ROI
wangRois = cellfun(@(x) x.getAtlasROIs('wang'),Rois,'UniformOutput',false);% % wang ROI
bensonRoisCent = cellfun(@(x) x.getAtlasROIs('benson',[0 4]),Rois,'UniformOutput',false);% % wang ROI
bensonRoisSurr = cellfun(@(x) x.getAtlasROIs('benson',[4 12]),Rois,'UniformOutput',false);% % wang ROI

savepath = '/Users/babylab/Documents/Elham/text-scramble/Simulations';
FigFolders = '/Users/babylab/Documents/Elham/text-scramble/Figs/Simulate';
noise.mu=3;
RoisIOG = cellfun(@(x) x.searchROIs('VWFA1','kgs','L'),Rois,'UniformOutput',false);% % wang ROI
[~,EEGAxx,~,masterList,subIDs] = mrC.Simulate.SimulateProject(ProjectPath,'signalArray',outSignal(:,1),'noiseParams',noise,'rois',RoisIOG,'SavePath',savepath,'cndNum',1);%'roiType','kgs','roiList',RoiList([1])'); % IOG_L
%[~,EEGAxx,~,masterList,subIDs] = mrC.Simulate.SimulateProject(ProjectPath,'signalArray',outSignal(:,1),'noiseParams',noise,'SavePath',savepath,'cndNum',1);%'roiType','kgs','roiList',RoiList([1])'); % IOG_L


%% One ROI
% savepath = '/Users/babylab/Documents/Elham/text-scramble/Simulations';
% FigFolders = '/Users/babylab/Documents/Elham/text-scramble/Figs/Simulate';
% noise.mu=3;
% [~,EEGAxx,~,masterList,subIDs] = mrC.Simulate.SimulateProject(ProjectPath,'signalArray',outSignal(:,1),'noiseParams',noise,'rois',Roi1(1),'SavePath',savepath,'cndNum',1);%'roiType','kgs','roiList',RoiList([1])'); % IOG_L
% SavePath = fullfile(FigFolders,'IOG');
% mrC.Simulate.PlotSSVEPonEGI(EEGAxx, 'Amplitude' , SavePath,'IOG',FundFreq,masterList,subIDs,'Individuals');%'Average')
% mrC.Simulate.PlotSSVEPonEGI(EEGAxx, 'Phase' , SavePath,'IOG',FundFreq,masterList,subIDs,'Individuals');%'Average')
% 
% [~,EEGAxx,~,masterList,subIDs] = mrC.Simulate.SimulateProject(ProjectPath,'signalArray',outSignal(:,1),'noiseParams',noise,'roiType','kgs','rois',Roi1(7)','SavePath',savepath,'cndNum',2); % VWFA1_L
% SavePath = fullfile(FigFolders,'VWFA1');
% mrC.Simulate.PlotSSVEPonEGI(EEGAxx, 'Amplitude' , SavePath,'VWFA1',FundFreq,masterList,subIDs,'Individuals');%'Average')
% mrC.Simulate.PlotSSVEPonEGI(EEGAxx, 'Phase' , SavePath,'VWFA1',FundFreq,masterList,subIDs,'Individuals');%'Average')
% 
% noise.mu=3;
% [~,EEGAxx,~,masterList,subIDs] = mrC.Simulate.SimulateProject(ProjectPath,'signalArray',outSignal(:,1),'noiseParams',noise,'roiType','wang','rois',Roi2(15)','SavePath',savepath,'cndNum',3); % LO1_L
% SavePath = fullfile(FigFolders,'LO1');
% mrC.Simulate.PlotSSVEPonEGI(EEGAxx, 'Amplitude' , SavePath,'LO1',FundFreq,masterList,subIDs,'Individuals');%'Average')
% mrC.Simulate.PlotSSVEPonEGI(EEGAxx, 'Phase' , SavePath,'LO1',FundFreq,masterList,subIDs,'Individuals');%'Average')
% 
% [~,EEGAxx,~,masterList,subIDs] = mrC.Simulate.SimulateProject(ProjectPath,'signalArray',outSignal(:,1),'noiseParams',noise,'roiType','wang','rois',Roi2(17)','SavePath',savepath,'cndNum',4); % LO2_L
% SavePath = fullfile(FigFolders,'LO2');
% mrC.Simulate.PlotSSVEPonEGI(EEGAxx, 'Amplitude' , SavePath,'LO2',FundFreq,masterList,subIDs,'Individuals');%'Average')
% mrC.Simulate.PlotSSVEPonEGI(EEGAxx, 'Phase' , SavePath,'LO2',FundFreq,masterList,subIDs,'Individuals');%'Average')

%% Two ROIs (one V1)
% 
% [~,EEGAxx,~,masterList,subIDs] = mrC.Simulate.SimulateProject(ProjectPath,'signalArray',outSignal,'noiseParams',noise,'roiType','wang','rois',Roi2([15 31])','roiSize',30,'SavePath',savepath,'cndNum',5); % LO1_L + V1v_L
% SavePath = fullfile(FigFolders,'LO1-V1v');
% mrC.Simulate.PlotSSVEPonEGI(EEGAxx, 'Amplitude' , SavePath,'LO1-V1v',FundFreq,masterList,subIDs,'Individuals');%'Average')
% mrC.Simulate.PlotSSVEPonEGI(EEGAxx, 'Phase' , SavePath,'LO1-V1v',FundFreq,masterList,subIDs,'Individuals');%'Average')
% 
% [~,EEGAxx,~,masterList,subIDs] = mrC.Simulate.SimulateProject(ProjectPath,'signalArray',outSignal,'noiseParams',noise,'roiType','wang','rois',Roi2([17 31])','roiSize',30,'SavePath',savepath,'cndNum',6); % LO2_L + V1v_L
% SavePath = fullfile(FigFolders,'LO2-V1v');
% mrC.Simulate.PlotSSVEPonEGI(EEGAxx, 'Amplitude' , SavePath,'LO2-V1v',FundFreq,masterList,subIDs,'Individuals');%'Average')
% mrC.Simulate.PlotSSVEPonEGI(EEGAxx, 'Phase' , SavePath,'LO2-V1v',FundFreq,masterList,subIDs,'Individuals');%'Average')
% 
% [~,EEGAxx,~,masterList,subIDs] = mrC.Simulate.SimulateProject(ProjectPath,'signalArray',outSignal,'noiseParams',noise,'roiType','wang','rois',[Roi1(1),Roi2(31)]','roiSize',30,'SavePath',savepath,'cndNum',7); % IOG_L + V1v_L
% SavePath = fullfile(FigFolders,'IOG-V1v');
% mrC.Simulate.PlotSSVEPonEGI(EEGAxx, 'Amplitude' , SavePath,'IOG-V1v',FundFreq,masterList,subIDs,'Individuals');%'Average')
% mrC.Simulate.PlotSSVEPonEGI(EEGAxx, 'Phase' , SavePath,'IOG-V1v',FundFreq,masterList,subIDs,'Individuals');%'Average')
% 
% [~,EEGAxx,~,masterList,subIDs] = mrC.Simulate.SimulateProject(ProjectPath,'signalArray',outSignal,'noiseParams',noise,'roiType','wang','rois',[Roi1(7),Roi2(31)]','roiSize',30,'SavePath',savepath,'cndNum',8); % VWFA1_L + V1v_L
% SavePath = fullfile(FigFolders,'VWFA1-V1v');
% mrC.Simulate.PlotSSVEPonEGI(EEGAxx, 'Amplitude' , SavePath,'VWFA1-V1v',FundFreq,masterList,subIDs,'Individuals');%'Average')
% mrC.Simulate.PlotSSVEPonEGI(EEGAxx, 'Phase' , SavePath,'VWFA1-V1v',FundFreq,masterList,subIDs,'Individuals');%'Average')

%% Two ROIs and different phases

% for p = 0:19
%     
%     [outSignal, FundFreq, SF]= mrC.Simulate.ModelSeedSignal('signalType','SSVEP','signalFreq',[3 3],'signalHarmonic',{[1],[1]},'signalPhase',{[pi/2], [(pi/2)+(p*pi/10)]});
%     
%     [EEGData,EEGAxx,sourceDataOrigin,masterList,subIDs] = mrC.Simulate.SimulateProject(ProjectPath,'signalArray',outSignal,'noiseParams',noise,'roiType','wang','rois',[Roi1(1),Roi2(31)]','roiSize',30); % IOG_L + V1v_L
%     close all;
%     %k = input('Press Enter');
%     SavePath = '/Users/babylab/Documents/Elham/text-scramble/Figs/Simulate/IOG-V1v-phase';
%     %%% CHANGE THESE PLOTS SO THAT THEY WILL AUTOMATICALLY SAVE THE RESULTS
%     mrC.Simulate.PlotSSVEPonEGI(EEGAxx, 'Amplitude' , SavePath,[num2str(rad2deg((p*pi/10))) 'deg'],FundFreq,masterList,subIDs,'Average','Simple');close all;
%     mrC.Simulate.PlotSSVEPonEGI(EEGAxx, 'Phase' , SavePath,[num2str(rad2deg((p*pi/10))) 'deg'],FundFreq,masterList,subIDs,'Average','Simple');close all;
% 
% end

%% Create resolution matrices
Cecc = 2.5;
Secc = 7.5;
%[CT1 ] = mrC.Simulate.ResolutionMatrices(ProjectPath,'inverse','mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangkgsROIsCorr.inv','rois',Rois);%mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangROIsCorr.inv
[CT3,~,ROISrc,List,SubIDs]= mrC.Simulate.ResolutionMatrices(ProjectPath,'inverse','mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangROIsCorr.inv','rois',Rois,'roiType','benson','eccRange',[Cecc Secc]);%
%[CT3,~,ROISrc,List,SubIDs]= mrC.Simulate.ResolutionMatrices(ProjectPath,'inverse','mneInv_bem_gcv_regu_F1_1_2_3_4_5_6.inv','rois',Rois,'roiType','benson','eccRange',[Cecc Secc]);%


MCT1 = mean(cat(3,CT1{:}),3);%MCT1([2 14 27],:)=[];MCT1(:,[2 14 27])=[];
MCT2 = mean(cat(3,CT2{:}),3);%MCT2([2 14 27],:)=[];MCT2(:,[2 14 27])=[];
MCT3 = mean(cat(3,CT3{:}),3);%MCT3([2 14 27],:)=[];MCT3(:,[2 14 27])=[];

MCT1 = MCT1./repmat(max(MCT1,[],2),[1 length(MCT1)]);
MCT2 = MCT2./repmat(max(MCT2,[],2),[1 length(MCT2)]);
MCT3 = MCT3./repmat(max(MCT3,[],2),[1 length(MCT3)]);

%List = [cellfun(@(x) [x.Name '_' x.Hemi],Roi2,'UniformOutput',false) cellfun(@(x) [x.Name '_' x.Hemi],Roi1,'UniformOutput',false)]; 
figure,
subplot(1,3,1),imagesc((MCT1));colorbar;caxis([-1 1]);colormap(jmaColors('coolhot'))
set(gca,'ytick',1:numel(List),'yticklabel',List)
subplot(1,3,2),imagesc((MCT2));colorbar;caxis([-1 1]);colormap(jmaColors('coolhot'))
set(gca,'ytick',1:numel(List),'yticklabel',List)
subplot(1,3,3),imagesc((MCT3));colorbar;caxis([-1 1]);
set(gca,'ytick',1:numel(List),'yticklabel',List)

%%

for s = 1: numel(SubIDs)
    mapMtx{s} = makeDefaultCortexMorphMap(SubIDs{s},SubIDs{2});
    morphSrc{s} = mapMtx{s}*ROISrc{s}.';
end

%%

MSrc = mean(cat(3,morphSrc{:}),3);

Roistemp = Rois{2};
RoisC = Roistemp.getAtlasROIs('benson',[0 Cecc]);
RoisS = Roistemp.getAtlasROIs('benson',[Cecc Secc]);
RoisCS = RoisC.mergROIs(RoisS);

roiChunk = RoisCS.ROI2mat(size(MSrc,1));

plotCentSurr(MSrc, roiChunk, SubIDs, List,'ventral')
plotCentSurr(MSrc, roiChunk, SubIDs, List,'anterior')
%


%%
for cent = 2:5
    for sur = 8:12
        display (['CENTER ' num2str(cent) ' SURROUND:' num2str(sur)]);
        [CT,CTnorm,~,Lists]= mrC.Simulate.ResolutionMatrices(ProjectPath,'inverse','mneInv_bem_gcv_regu_F1_1_2_3_4_5_6.inv','rois',Rois,'roiType','benson','eccRange',[cent sur]);%
        CTM = nanmean(cat(3,CT{:}),3);
        CTnormM = nanmean(cat(3,CTnorm{:}),3);
        ctm(:,:,cent,sur) = CTM./repmat(max(CTM,[],2),[1 length(CTM)]);
        ctnormm(:,:,cent,sur) = CTnormM./repmat(max(CTnormM,[],2),[1 length(CTnormM)]);
    end
end

%% calculate the crosstalk error
ctm2 = ctm(:,:,2:5,8:12);
for roi = 1:size(ctm2,1)
    cterr(roi,:,:) = squeeze(sum(abs(ctm2(roi,[1:roi-1 roi+1:end],:,:)),2)./ctm2(roi,roi,:,:));
    cterr2(roi,:,:) = squeeze(sum(abs(ctm2(roi,[1:roi-1 roi+1:end],:,:)),2));

end

%%
figure,
subplot(1,2,1), plot(2:5,squeeze(mean(cterr(13:end,:,3))),'b','linewidth',2)
hold on; plot(2:5,squeeze(mean(cterr(1:12,:,3))),'r','linewidth',2)
%legend('Surround ROI','Center ROI');
ylabel('Cross Talk Error','fontsize',14,'fontweight','bold')
xlabel('Center ecc','fontsize',14,'fontweight','bold')
title('Varying center size');

subplot(1,2,2), plot(8:12,squeeze(mean(cterr(13:end,2,:))),'b','linewidth',2)
hold on; plot(8:12,squeeze(mean(cterr(1:12,2,:))),'r','linewidth',2)
h = legend('Surround ROI','Center ROI');set(h,'fontsize',14);
xlabel('Surround ecc','fontsize',14,'fontweight','bold')
title('Varying surround size');

ylim([3.5 8.5]);

%%
figure,
subplot(1,2,1), imagesc(squeeze(mean(cterr(1:12,:,:),1)));
set(gca,'ytick',1:4,'yticklabel',2:5,'xtick',1:5,'xticklabel',8:12);
title('Center Error','fontsize',14);xlabel('Surround Size','fontsize',14);ylabel('Center Size','fontsize',14)
caxis([3.5 10])

subplot(1,2,2),imagesc(squeeze(mean(cterr(13:end,:,:),1)));
set(gca,'ytick',1:4,'yticklabel',2:5,'xtick',1:5,'xticklabel',8:12);
title('Surround Error','fontsize',14);xlabel('Surround Size','fontsize',14);
caxis([3.5 10])
colormap('hot');

%%



