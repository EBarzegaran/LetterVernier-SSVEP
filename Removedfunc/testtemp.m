path2 = '/Users/kohler/code/git/gardner/mrTools/mrUtilities/mrFlatMesh';
addpath(genpath(path2));

addpath(genpath('/Users/kohler/code/git/sweepAnalysis/functions/helper'));

%% Convert spectrum data into source space
Path = '/Volumes/svndl/mrC_Projects/VernierLetters/source';
% Inverse = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangkgsROIsCorr.inv';
Inverse = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangkgsROIsCorr_DepthWeight.inv';
[outData,~,FreqFeatures,subIDs] = mrC.SourceBrain(Path,Inverse,'domain','frequency','doSmooth' , true);%,'template','nl-0014');
%% load morph maps
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
%InvROIs = arrayfun(@(x) Inverses{x}*WKROIsM{x}./(sum(WKROIsM{x})+eps),1:numel(Inverses),'uni',false);

Inverse1 = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangkgsROIsCorr_DepthWeight.inv';
[Inverses,subIDs_Inverse] = mrC.Simulate.ReadInverses(Path,Inverse1);

%%
% 
% for i = 1:64
%     for subj = 1:numel(subIDs)
%         TData = squeeze(FreqData.Letter(WKROIsM{subj}(:,i)>0,7,3,subj));
%         if ~isempty(TData) && numel(TData)>2
%             %[Z_est, confidence_radii, p(subj,i), t2circ(subj,i)] =t2circ_1tag(TData.');
%             [~,~,zSNR(subj,i),errorEllipse] = fitErrorEllipse([real(TData(:)) imag(TData(:))],[],0);
%             ROImat = Inverses{subj}(:,WKROIsM{subj}(:,i)>0);
%             sourcenum(subj,i) = sum(WKROIsM{subj}(:,i)>0);
%             Temp = triu(corr(ROImat));
%             Temp(1:length(Temp)+1:end)=0;
%             ROIMCorr(subj,i) = sum(Temp(:))/sum(Temp(:)~=0);
%             abROIMCorr(subj,i) = sum(abs(Temp(:)))/sum(Temp(:)~=0);
%         else
%             ROIMCorr(subj,i) = 0;
%             zSNR(subj,i) = nan;
% %             p(subj,i)=1;
% %             t2circ(subj,i)=0;
%         end
%     end
% end
% ROIMCorr(isnan(ROIMCorr))=0;

%%

for subj = 1:numel(subIDs)
   Cmat = corr(Inverses{subj});
   for j = 1:64
        for i = 1:64
            ROImat = Cmat(WKROIsM{subj}(:,j)>0,WKROIsM{subj}(:,i)>0);
            if ~isempty(ROImat)
                if i==j
                    Temp = triu(ROImat);
                    Temp(1:length(Temp)+1:end)=0;
                    ROIMCorr(subj,i,j) = sum(Temp(:))/sum(Temp(:)~=0);
                    %abROIMCorr(subj,i,j) = sum(abs(Temp(:)))/sum(Temp(:)~=0);
                else
                    ROIMCorr(subj,i,j) = sum(ROImat(:))/sum(ROImat(:)~=0);
                end
            else
                ROIMCorr(subj,i,j) = nan;
            end
        end
    end
end

%%
MM = squeeze(mean(ROIMCorr(:,:,:),1));
%MM = MM./(repmat(max(MM),[64 1]));
I = imagesc(MM);
set(gca,'xtick',1:numel(ROILabel),'xticklabel',ROILabel,'position',get(gca,'position')+[0 .07 0 0]);
set(gca,'ytick',1:numel(ROILabel),'yticklabel',ROILabel,'position',get(gca,'position')+[0 .07 0 0]);
xtickangle(90);
colormap(jmaColors('coolhot'));caxis([-1 1])
set(gca,'position',get(gca,'position')-[0 .1 0 0])
%%
