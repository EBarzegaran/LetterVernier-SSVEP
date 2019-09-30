
clear; clc;

PATH = '/Users/kohler/code/git/mrC';
addpath(genpath(PATH));
path2 = '/Users/kohler/code/git/gardner/mrTools/mrUtilities/mrFlatMesh';
addpath(genpath(path2));
%% Read the inverses
Path = '/Volumes/svndl/mrC_Projects/VernierLetters/source';
InvName = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangkgsROIsCorr_DepthWeight.inv';
InvName = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangkgsROIsCorr.inv';
[Inverses,SubIDInv] = mrC.Simulate.ReadInverses(Path,InvName);
% Read  morph maps
for s = 1:numel(SubIDInv)
    mapMtx{s} = makeDefaultCortexMorphMap(SubIDInv{s},SubIDInv{2});
    InvMapped{s} = mapMtx{s}*Inverses{s}.';
    
    %smoothing
    fromCtx = readDefaultCortex(SubIDInv{s});
    fromCtx.uniqueVertices = fromCtx.vertices;
    fromCtx.uniqueFaceIndexList = fromCtx.faces;
    [fromCtx.connectionMatrix] = findConnectionMatrix(fromCtx);
    fromCtx.connectionMatrix = fromCtx.connectionMatrix + speye(length(fromCtx.connectionMatrix));
    sumNeighbours=sum(fromCtx.connectionMatrix,2); % Although it should be symmetric, we specify row-summation
    smMtx=bsxfun(@rdivide,fromCtx.connectionMatrix,sumNeighbours);
    smoothMtx{s} = smMtx*smMtx;
    
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
%% load pre-calculate RCA results
Task = {'Letter','Vernier'};
Harms = {'H1F1','H2F1','H3F1','H4F1'};
Finds  = [7 13 19 25];
analHarms = [1 2];

RedoRCA = false;
if RedoRCA || ~exist(fullfile('ResultData','GroupRCA.mat'),'file')
   
else
    load(fullfile('ResultData','GroupRCA.mat'));
end

%% All subject IDs and select the subjects with inverse
PDiva_Path = '/Volumes/Denali_DATA1/Elham/EEG_Textscamble/';
Subjfolders = subfolders(PDiva_Path,0);
Subjfolders =  Subjfolders(cellfun(@(x) ~isempty(x),strfind(Subjfolders,'nl-')));
SubIDs = cellfun(@(x) x(1:7),Subjfolders,'UniformOutput',false);
SubSelect = cell2mat(arrayfun(@(x) find(strcmp(SubIDs,SubIDInv{x})),1:numel(SubIDInv),'uni',false));

%% Reverse of RCA
CompName = {'RC1','RC2'};
for ts = 1:numel(Task)
    for h = 1:2
        clear OutAxx;
        W = W_all.(Task{ts}).(Harms{h});
        dims = 128;
        WI = pinv(W);
        InAxx = decompAxx_all.(Task{ts}).(Harms{h});
        OutAxx = InAxx ;
        for comp =1:2
            temp = WI(comp,:)'*reshape(permute(InAxx.Cos(:,comp,:),[2,1,3]),size(InAxx.Cos(:,comp,:),2),[]);
            OutAxx.Cos = permute(reshape(temp,dims,size(InAxx.Cos(:,comp,:),1),size(InAxx.Cos(:,comp,:),3)),[2,1,3]);
            temp = WI(comp,:)'*reshape(permute(InAxx.Sin(:,comp,:),[2,1,3]),size(InAxx.Sin(:,comp,:),2),[]);
            OutAxx.Sin = permute(reshape(temp,dims,size(InAxx.Sin(:,comp,:),1),size(InAxx.Sin(:,comp,:),3)),[2,1,3]);
            OutAxx.Amp = abs(OutAxx.Cos +1i *OutAxx.Sin);
            temp = WI(comp,:)'*reshape(permute(InAxx.Wave(:,comp,:),[2,1,3]),size(InAxx.Wave(:,comp,:),2),[]);
            OutAxx.Wave = permute(reshape(temp,dims,size(InAxx.Wave(:,comp,:),1),size(InAxx.Wave(:,comp,:),3)),[2,1,3]);
            R_axx.(Task{ts}).(Harms{h}).(CompName{comp}) = OutAxx;
        end
    end
end

%% reshape the RCs to prepare for source estimation
Condnum = 5; % number of conditions per task
Subnum = numel(SubIDs);
 for ts = 1:numel(Task)
     for f = 1:2
        if f==1
            Flips = [1 -1]; 
        else
            Flips = [1 1];
        end
        for comp = 1:2
            numcomp = 128;
            TCos = squeeze(mean(reshape(squeeze(R_axx.(Task{ts}).(Harms{analHarms(f)}).(CompName{comp}).Cos(Finds(analHarms(f)),:,1:end-70)),[numcomp,16,Condnum,Subnum-1]),2));% last subject has 14 trails
            TSin = squeeze(mean(reshape(squeeze(R_axx.(Task{ts}).(Harms{analHarms(f)}).(CompName{comp}).Sin(Finds(analHarms(f)),:,1:end-70)),[numcomp,16,Condnum,Subnum-1]),2));

            TCos(:,:,Subnum) = squeeze(mean(reshape(squeeze(R_axx.(Task{ts}).(Harms{analHarms(f)}).(CompName{comp}).Cos(Finds(analHarms(f)),:,end-69:end)),[numcomp,14,Condnum,1]),2));
            TSin(:,:,Subnum) = squeeze(mean(reshape(squeeze(R_axx.(Task{ts}).(Harms{analHarms(f)}).(CompName{comp}).Sin(Finds(analHarms(f)),:,end-69:end)),[numcomp,14,Condnum,1]),2));
            TCmplx.(Task{ts}).(Harms{analHarms(f)}).(CompName{comp}) = Flips(ts)*TCos+(TSin*Flips(ts)*1i);
        end
     end
 end

%% calculate sources of Rcs for subjects with inverse solution
FigPath = 'Figures';
HemiName = {'Left','Right'};
load 'resultData/LogMar_Val.mat';
LM = num2cell(round(logMAR_letter,2));
Legend.(Task{1}) = cellfun(@(x) ['LogMAR = ' num2str(x)],LM,'uni',false);
LM = num2cell(round(logMAR_ver,2));
Legend.(Task{2}) =  cellfun(@(x) ['LogMAR = ' num2str(x)],LM,'uni',false);
colors = winter(5);

plotType = 1;% 1 for ROI and 2 for powermaps
Space = 2; % 1 for amplitude and 2 for phase
for h = 1:1
    for ts = 1:numel(Task)
         for f = 1:2
             FIG = figure;
             set(FIG,'unit','inch','position',[10 10 15 8])
             for comp = 1:2
                 Temp = TCmplx.(Task{ts}).(Harms{analHarms(f)}).(CompName{comp})(:,:,SubSelect);
                 for sub = 1:size(Temp,3)
                     SourceTemp(:,:,sub) = squeeze(Temp(:,:,sub).'*InvMapped{sub}.'*smoothMtx{sub});
                     ROITemp(:,:,sub) = squeeze(Temp(:,:,sub).'*Inverses{sub}*smoothMtx{sub}'*WKROIsM{sub});
                 end
                 if plotType == 1
                     subplot(2,1,comp) 
                     if Space ==1
                         B = bar(abs(mean(ROITemp(:,h:2:end,:),3))','FaceColor','flat');
                     else
                         B = bar(wrapTo2Pi(angle(mean(ROITemp(:,h:2:end,:),3)))','FaceColor','flat');
                     end
                     title([Task{ts} ' _ ' Harms{f}(2:end) ' _ ' CompName{comp}]);
                     set(gca,'xtick',1:32,'xticklabel',ROILabel(h:2:end));
                     xtickangle(90);
                     if Space==1
                        ylim([0 .0003])
                     else
                        ylim([0 2*pi])
                     end
                    for k = 1:5
                        B(k).CData = colors(k,:);
                    end
                 else
%                      for cond =1:5
%                          S = subplot(2,5,cond+(comp-1)*5);
%                          mrC.Simulate.VisualizeSourceData('nl-0014',abs(mean(SourceTemp(cond,:,:),3)),[],[],'ventral');
%                          caxis([0 max(abs(mean(SourceTemp(cond,:,:),3)))/1.5])
%                          set(S,'position',get(S,'position')+[-.05 0 0.08 0])
%                          colormap('jet')
%                     end
                    S = subplot(2,2,1+(comp-1)*2);
                    %mrC.Simulate.VisualizeSourceData('nl-0014',abs(mean(SourceTemp(5,:,:),3)),[],[],'ventral');
                    mrC.plotOnEgi(wrapTo2Pi(angle(mean(Temp(:,5,:),3))));
                    set(S,'position',get(S,'position')+[-.05 0 0.1 0.08])
                    %caxis([0 max(abs(mean(SourceTemp(cond,:,:),3)))/1.2])

                    S = subplot(2,2,2+(comp-1)*2);
                    %mrC.Simulate.VisualizeSourceData('nl-0014',abs(mean(SourceTemp(5,:,:),3)),[],[]);
                    set(S,'position',get(S,'position')+[-.05 0 0.08 0])
                    %caxis([0 max(abs(mean(SourceTemp(cond,:,:),3)))/1.2])
                    colormap('jet')
                 end
             end
             if plotType == 1
                 legend(Legend.(Task{ts}));
                 export_fig(FIG,fullfile(FigPath,['SourceEstimation_RCA_' Task{ts} ' _ ' Harms{f}(2:end) '_' HemiName{h}]),'-pdf')
                 close;
             else
                 
             end
             
         end
    end
end


