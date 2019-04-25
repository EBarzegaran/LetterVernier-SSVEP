
clear; clc;

PATH = '/Users/kohler/code/git/mrC';
addpath(genpath(PATH));
path2 = '/Users/kohler/code/git/gardner/mrTools/mrUtilities/mrFlatMesh';
addpath(genpath(path2));
%% Read the inverses
Path = '/Volumes/svndl/mrC_Projects/VernierLetters/source';
InvName = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangkgsROIsCorr_DepthWeight.inv';
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
for h = 1:2
    for ts = 1:numel(Task)
         for f = 1:2
             FIG = figure;
             set(FIG,'unit','inch','position',[10 10 15 8])
             for comp = 1:2
                 Temp = TCmplx.(Task{ts}).(Harms{analHarms(f)}).(CompName{comp})(:,:,SubSelect);
                 for sub = 1:size(Temp,3)
                     ROITemp(:,:,sub) = squeeze(Temp(:,:,sub).'*Inverses{sub}*WKROIsM{sub});
                 end
                 subplot(2,1,comp)
                 bar(abs(mean(ROITemp(:,h:2:end,:),3))')
                 title([Task{ts} ' _ ' Harms{f}(2:end) ' _ ' CompName{comp}]);
                 set(gca,'xtick',1:32,'xticklabel',ROILabel(h:2:end));
                 xtickangle(90);
                 %ylim([0 .0003])
                 ylim([-pi pi])
             end
             %export_fig(FIG,fullfile(FigPath,['SourceEstimation_RCA_' Task{ts} ' _ ' Harms{f}(2:end) '_' HemiName{h}]),'-pdf')
             %close;
         end
    end
end


