
clear; clc;

PATH = '/Users/elhamb/Documents/Codes/Git/mrC';
addpath(genpath('/Users/elhamb/Documents/Codes/Git/EEGSourceSim'));
addpath(genpath(PATH));

%% 
 Path = '/Users/elhamb/Documents/Data/TextScramble_exp1/VernierLetters/source';
Inverse1 = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangkgsROIsCorr.inv';
Inverse2 = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangkgsROIsCorr_DepthWeight.inv';

ProjectPath = fullfile(Path);
%% Select subjects with inverses
Lateral = {'V','D','VD'};
l = 3;

[Inverse,subIDs_Inverse] = mrC.Simulate.ReadInverses(ProjectPath,Inverse1);
subIDs_Inverse = subIDs_Inverse(cellfun(@(x) ~isempty(x),Inverse));
clear Inverse;

% Pre-select ROIs
[RoiList,subIDs] = ESSim.Simulate.GetRoiClass(ProjectPath,[],subIDs_Inverse);% 13 subjects with Wang atlab 
Wang_RoiList = cellfun(@(x) {x.getAtlasROIs('wang')},RoiList);
kgs_RoiList = cellfun(@(x) {x.getAtlasROIs('kgs')},RoiList);
Wangkgs_RoiList = arrayfun(@(x) Wang_RoiList{x}.mergROIs(kgs_RoiList{x}),1:numel(Wang_RoiList),'uni',false);

switch l
    case 1 % Ventral V1-V3 and IOG and VWFA
        ROIind = [31:32 35:36 43:44 51:52 57:58];          
        Wangkgs_RoiList = cellfun(@(x) {x.selectROIs(ROIind)},Wangkgs_RoiList);
        ROI_Final = Wangkgs_RoiList;
    case 2 % Dorsal V1-V3 and IOG and VWFA
        ROIind = [[31:32 35:36 43:44]-2 [51:52 57:58]];          
        Wangkgs_RoiList = cellfun(@(x) {x.selectROIs(ROIind)},Wangkgs_RoiList);
        ROI_Final = Wangkgs_RoiList;
    case 3 % Both dorsa and Ventral
        ROIind = [29:36 41:44 51:52 57:58];                
        Wangkgs_RoiList = cellfun(@(x) {x.selectROIs(ROIind)},Wangkgs_RoiList);
        Merges = {[1 3],[2 4], [5 7], [6 8], [9 11], [10 12], 13, 14, 15, 16}; % Merge ROIs
        Temp_Names = {'V1','V1','V2','V2','V3','V3','IOG','IOG','VWFA','VWFA'};
        ROI_Final = cellfun(@(x) x.mergeIndROIs(Merges,Temp_Names),Wangkgs_RoiList,'uni',false);

end
 
 ROILabel = ROI_Final{1}.getFullNames('noatlas');
%% Generate Resolution matrices
ResultPath = 'ResultData';
FilePath = fullfile(ResultPath,['LocalizationExampleData_Paper' Lateral{l} '.mat']);
do_new_data_generation = true;

if ~exist(FilePath,'file') || do_new_data_generation
    [CrossTalk1,Error1,ROISource1,ScalpData,LIST,subIDs] = ESSim.Simulate.ResolutionMatrices(ProjectPath,'subSelect',subIDs,...
        'rois',ROI_Final,'doAUC',true,'inverse',Inverse1,'roiType','all');
    
%     [CrossTalk2,Error2,ROISource2,ScalpData,LIST,subIDs] = ESSim.Simulate.ResolutionMatrices(ProjectPath,'subSelect',subIDs,...
%         'rois',ROI_Final,'doAUC',true,'inverse',Inverse2,'roiType','all');
    save(FilePath,'CrossTalk1','Error1','ROISource1','ScalpData','LIST','subIDs','ROILabel');
else
    load(FilePath);
end

%% 
FS = 12;
FigPath = fullfile('Figures','SourceSpace');
HemiName = {'left','right'};

% Plot Cross Talk Matrices
CrossTalk1 = cellfun(@(x) x./repmat(max(x),[size(x,1) 1]),CrossTalk1,'uni',false);% normalize
% CrossTalk2 = cellfun(@(x) x./repmat(max(x),[size(x,1) 1]),CrossTalk2,'uni',false);% normalize
CT1 = (cat(3,CrossTalk1{:}));
% CT2 = (cat(3,CrossTalk2{:}));


FIG = figure;
for i = 1
    eval(['CT_plot = ((mean(CT' num2str(i) ',3)));']);
    %CT_plot(1:length(CT_plot)+1:end)=0;
    imagesc(abs(CT_plot));

   colormap(jmaColors('coolhot'));
   colormap('gray')
    %caxis([-max(abs(CTMM(:))) max(abs(CTMM(:)))]);
    caxis([0 1]/1.5);

    set(gca,'ytick',1:numel(ROILabel),'yticklabel',ROILabel,'xtick',1:numel(ROILabel),'xticklabel',ROILabel,'fontsize',FS);
    if exist('xtickangle'), xtickangle(90); end
    xlabel('receiving ROI','fontsize',12);
    if i==1
        %title('MN+FACE');
        title('Cross Talk Matrix')
    end
    colorbar
end
set(FIG,'PaperPosition',[1 1 8.5 8]);
set(gcf, 'Color', 'w');
set(FIG,'Units','Inch')
set(FIG,'Position',[1 1 8.5 8]);
mkdir(fullfile(FigPath,'simulation'))
print(fullfile(FigPath,'simulation','CrossTalk_All'),'-dtiff','-r300');
close;

%%
for i = 1:2
    eval(['CTM(:,:,i) = mean(cat(3,CrossTalk' num2str(i) '{:}),3);']);
end

for r = 1:64
    FIG= figure;
    indR = mod(r,2);
    if indR==0,indR=2;end
    
    bar(squeeze(CTM(r,indR:2:end,:)));
    ylim([-1 1]);
    
    set(gca,'xtick',1:numel(ROILabel(indR:2:end)),'xticklabel',ROILabel(indR:2:end),'fontsize',10);
    if exist('xtickangle'), xtickangle(90); end
    title(ROILabel(r))
    legend('MN+FACE','WMN+FACE')
    
    set(FIG,'PaperPosition',[1 1 12 4]);
    %print(fullfile(FigPath,['SourceEstimation_Simulation_Error_' HemiName{Hemi}]),'-r300','-dtiff')
    set(gcf, 'Color', 'w');
    set(FIG,'Units','Inch')
    set(FIG,'Position',[1 1 12 4]);
    print(fullfile(FigPath,'simulation',['CrossTalk_' ROILabel{r}]),'-dtiff','-r300');
    close;
end

%% AUC relative Enegry, Focalization Error
% prepare data


for Hemi = 1:2% 1 for left  and 2 for right
    ROINum = 64;
    for i = 1:2
        eval(['Error = Error' num2str(i) ';']);
        for sub = 1:numel(subIDs)
            TP = Error{sub}.TP;
            FP = Error{sub}.FP;
            FN = Error{sub}.FN;
            Recall = TP./(TP+FN);
            Precision = TP./(TP+FP);
            for r = 1:size(TP,1)
                AUC(r,sub) = trapz(Recall(r,:),Precision(r,:));
            end
            Relative(:,sub) = Error{sub}.Relative;
            Focal(:,sub) = Error{sub}.Focalization;
        end

        AUC2(:,:,i) = (AUC(Hemi:2:end,:));%+AUC(2:2:end,:))/2; % Average over left and right
        Focal2(:,:,i) = (Focal(Hemi:2:end,:));%+Focal(2:2:end,:))/2; % Average over left and right
        Relative2(:,:,i) = (Relative(Hemi:2:end,:));%+Relative(2:2:end,:))/2; % Average over left and right
    end

    FS = 11;
    FIG3 = figure;
    C = [.6 .6 .6];%Colors;%

    % AREA UNDER CURVE
    AUCM = squeeze(mean(AUC2,2));
    AUCS = squeeze(std(AUC2,[],2))./sqrt(size(AUC2,2));
    S1 = subplot(3,1,1); 
    B = bar(AUCM);
    set(B(1),'FaceColor',[0.3,.3,.3])
    set(B(2),'FaceColor',[.7,.7,.7])
    %set(gca,'xtick',1:numel(LIST2)*2,'xticklabel',XTLabel,'fontsize',FS-1,'ytick',0.2:.2:1,'yticklabel',arrayfun(@(x) num2str(x),0.2:.2:1,'uni',false));
    set(gca,'xtick',1:numel(ROILabel(Hemi:2:end)),'xticklabel',ROILabel(Hemi:2:end),'fontsize',10);
    xtickangle(90);
    ylim([0.0 0.85]);
    hold on;
    XT = [(1:size(AUCM,1))-.15; (1:size(AUCM,1))+.15];
    errorbar(reshape(XT,1,numel(XT)),reshape(AUCM',1,numel(AUCM)),reshape(AUCS',1,numel(AUCS)),'.k');
    YL = ylabel('AUCPR','fontsize',FS,'fontweight','bold'); YLP = get(YL,'position');
    set(YL,'position',get(YL,'position')+[-0.2 0 0]); 
    set(S1,'position',get(S1,'position')+[0 0 0 .04]);

    legend({'MN + FACE','WMN + FACE'},'location','northeast')

    % RELATIVE ENERGY
    RelativeM = squeeze(mean(Relative2,2));
    RelativeS = squeeze(std(Relative2,[],2))./sqrt(size(Relative2,2));
    S2 = subplot(3,1,2); 
    B = bar(RelativeM);
    set(B(1),'FaceColor',[0.3,.3,.3])
    set(B(2),'FaceColor',[.7,.7,.7])
    %set(gca,'xtick',1:numel(LIST2)*2,'xticklabel',XTLabel,'fontsize',FS-1,'ytick',0:.1:.3,'yticklabel',arrayfun(@(x) num2str(x),0:.1:.3,'uni',false));
    set(gca,'xtick',1:numel(ROILabel(Hemi:2:end)),'xticklabel',ROILabel(Hemi:2:end),'fontsize',10);
    xtickangle(90);
    ylim([0 .45]);
    hold on;
    XT = [(1:size(RelativeM,1))-.15; (1:size(RelativeM,1))+.15];
    errorbar(reshape(XT,1,numel(XT)),reshape(RelativeM',1,numel(RelativeM)),reshape(RelativeS',1,numel(RelativeS)),'.k');
    YL = ylabel('Relative Energy','fontsize',FS,'fontweight','bold');
    YLPT = get(YL,'position');
    set(YL,'position',[YLP(1) YLPT(2:end)]); 
    set(YL,'position',get(YL,'position')+[-0.15 0 0]);
    set(S2,'position',get(S2,'position')+[0 -.02 0 .04]);



    % FOCALIZATION ERROR
    FocalM = squeeze(mean(Focal2,2));
    FocalS = squeeze(std(Focal2,[],2))./sqrt(size(Focal2,2));
    S3 = subplot(3,1,3); 
    B = bar(FocalM);
    set(B(1),'FaceColor',[0.3,.3,.3])
    set(B(2),'FaceColor',[.7,.7,.7])
    %set(gca,'xtick',1:numel(LIST2)*2,'xticklabel',XTLabel,'fontsize',FS-1,'ytick',0.2:.2:1,'yticklabel',arrayfun(@(x) num2str(x),0.2:.2:1,'uni',false));
    set(gca,'xtick',1:numel(ROILabel(Hemi:2:end)),'xticklabel',ROILabel(Hemi:2:end),'fontsize',10);
    xtickangle(90);
    ylim([0.2 .9]);
    hold on;
    XT = [(1:size(FocalM,1))-.15; (1:size(FocalM,1))+.15];
    errorbar(reshape(XT,1,numel(XT)),reshape(FocalM',1,numel(FocalM)),reshape(FocalS',1,numel(FocalS)),'.k');
    YL = ylabel('Focalization Error','fontsize',FS,'fontweight','bold');
    YLPT = get(YL,'position');
    set(YL,'position',[YLP(1) YLPT(2:end)]); 
    set(S3,'position',get(S3,'position')+[0 -.04 0 .04]);


    set(FIG3,'PaperPosition',[1 1 12 8]);
    print(fullfile(FigPath,['SourceEstimation_Simulation_Error_' HemiName{Hemi}]),'-r300','-dtiff')
    set(gcf, 'Color', 'w');
    set(FIG3,'Units','Inch')
    set(FIG3,'Position',[1 1 12 8]);
    export_fig(FIG3,fullfile(FigPath,['SourceEstimation_Simulation_Error_' HemiName{Hemi}]),'-pdf')
close;
end

%% make scalp plots

Temp = cat(3,ScalpData{:});
FIG = figure;
for i =1:30
   subplot(5,6,i),mrC.plotOnEgi(mean(mean(Temp((i*2)-1:(i*2),:,:),1),3)); axis tight
   title(ROILabel{(i*2)}(1:end-1))
end
set(FIG,'unit','inch','position',[5 5 20 20]);
