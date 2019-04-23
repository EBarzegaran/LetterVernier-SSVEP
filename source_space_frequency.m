% This script is for preliminary analysis and visualization of vernier-text
% scramble task in source space
% Elham Barzegaran, 5.22.2018
% modified: EB, 4.19.2019

clear;clc;

PATH = '/Users/kohler/code/git/mrC';
addpath(genpath(PATH));
path2 = '/Users/kohler/code/git/gardner/mrTools/mrUtilities/mrFlatMesh';
addpath(genpath(path2));

addpath(genpath('/Users/kohler/code/git/sweepAnalysis/functions/helper'));

%% Convert spectrum data into source space
Path = '/Volumes/svndl/mrC_Projects/VernierLetters/source';
Inverse = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangROIsCorr.inv';
Inverse = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangkgsROIsCorr.inv';
[outData,~,FreqFeatures,subIDs] = mrC.SourceBrain(Path,Inverse,'domain','frequency','doSmooth' , true);%,'template','nl-0014');
% load morph maps
for s = 1:numel(subIDs)
    mapMtx{s} = makeDefaultCortexMorphMap(subIDs{s},subIDs{2});
end

%% Organize the data
Task = {'Letter','Vernier'};
Harms = {'h1F1','h2F1'};
for ts = 1:2
    for cond = (1:5)+((ts-1)*5)
        FreqData.(Task{ts})(:,:,cond-((ts-1)*5),:) = cat(3,outData{cond,:});
    end
end
clear outData;

%% plot amplitude and phase group level

anatDir = '/volumes/svndl/anatomy';
direction='ventral';
Finds  = [7 13 19 25];
for h = 1:2
    FIG = figure;
    set(FIG,'unit','inch','position',[3 5 20 8]);
    for ts = 1:2       
        M = max(max(abs(squeeze(mean(FreqData.(Task{ts})(:,Finds(h),:,:),4)))));
        for cond = 1:5
            FData_temp = squeeze(FreqData.(Task{ts})(:,Finds(h),cond,:));
            FData = arrayfun(@(x) mapMtx{x}*FData_temp(:,x),1:numel(subIDs),'uni',false);
            subplot(2,5,cond+((ts-1)*5)),mrC.Simulate.VisualizeSourceData('nl-0014',abs(mean(cat(2,FData),2)),anatDir);
            caxis([-0 M]);
            axis tight;
        end
    end
end

close all;
%% Roi analysis: Estimate amplitude and phase for each Roi

% Load ROIs
[ROIs,subIDROI] = mrC.Simulate.GetRoiClass(Path);
Ords = cell2mat(arrayfun(@(x) find(strcmp(subIDROI,subIDs{x})),1:numel(subIDs),'uni',false));
ROIs = ROIs(Ords);

% Prepare Wang and KGS atlases
WangROIs = cellfun(@(x) x.getAtlasROIs('wang'),ROIs,'uni',false); 
KgsROIs = cellfun(@(x) x.getAtlasROIs('kgs'),ROIs,'uni',false); 
WKROIs = arrayfun(@(x) WangROIs{x}.mergROIs(KgsROIs{x}),1:numel(ROIs),'uni',false); 
WKROIsM = cellfun(@(x) x.ROI2mat(20484),WKROIs,'uni',false);
ROILabel = WKROIs{1}.getFullNames('noatlas');
InvROIs = arrayfun(@(x) Inverses{x}*WKROIsM{x}./(sum(WKROIsM{x})+eps),1:numel(Inverses),'uni',false);


%%
ROIselect = [3:10 15:18 29:60];
ROILabelSel = ROILabel(ROIselect);
for h = 1:2
    for ts = 1:2
       FIG = figure;
       set(FIG,'unit','inch','position',[3 5 20 8]);
       
       for cond = 1:5
           FData_temp = squeeze(FreqData.(Task{ts})(:,Finds(h),cond,:));
           %FDataRoi = arrayfun(@(x) mapMtx{x}*FData_temp(:,x),1:numel(subIDs),'uni',false);
           FDataROI_temp = arrayfun(@(x) FData_temp(:,x).'*WKROIsM{x}./sum(WKROIsM{x}),1:numel(subIDs),'uni',false);
           FDataROI.(Task{ts}).(Harms{h})(:,:,cond) = cat(1,FDataROI_temp{:});
           
       end
       
       PData = squeeze(mean(FDataROI.(Task{ts}).(Harms{h})(:,ROIselect,:),1));
       subplot(2,1,1),bar(abs(PData(1:2:end,:)));
       set(gca,'xtick',1:22,'xticklabel',ROILabelSel(1:2:end));
       title([Task{ts} '_' Harms{h}])
       ylim([0 7])
       
       subplot(2,1,2),bar(abs(PData(2:2:end,:)));
       set(gca,'xtick',1:22,'xticklabel',ROILabelSel(2:2:end));
       ylim([0 7])
    end
end

