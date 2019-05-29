clear; clc;

%% 
Path = '/Volumes/svndl/mrC_Projects/VernierLetters/source';
Inverse1 = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangkgsROIsCorr.inv';
Inverse2 = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangkgsROIsCorr_DepthWeight.inv';

ProjectPath = fullfile(Path);
%% Select subjects with inverses
[Inverse,subIDs_Inverse] = mrC.Simulate.ReadInverses(ProjectPath,Inverse1);
subIDs_Inverse = subIDs_Inverse(cellfun(@(x) ~isempty(x),Inverse));
clear Inverse;

% Pre-select ROIs
[RoiList,subIDs] = mrC.Simulate.GetRoiClass(ProjectPath,[],subIDs_Inverse);% 13 subjects with Wang atlab 
Wang_RoiList = cellfun(@(x) {x.getAtlasROIs('wang')},RoiList);
kgs_RoiList = cellfun(@(x) {x.getAtlasROIs('kgs')},RoiList);
Wangkgs_RoiList = arrayfun(@(x) Wang_RoiList{x}.mergROIs(kgs_RoiList{x}),1:numel(Wang_RoiList),'uni',false);
ROIind = [15:18 25:52 57:59];% select a few ROIs
Wangkgs_RoiList = cellfun(@(x) {x.selectROIs(ROIind)},Wangkgs_RoiList);
ROILabel = Wangkgs_RoiList{1}.getFullNames('noatlas');
%% Generate Resolution matrices
ResultPath = 'ResultData';
FilePath = fullfile(ResultPath,'LocalizationExampleData_Paper.mat');
do_new_data_generation = false;

if ~exist(FilePath,'file') || do_new_data_generation
    [CrossTalk1,Error1,ROISource1,~,~,~] = mrC.Simulate.ResolutionMatrices(ProjectPath,'subSelect',subIDs,...
        'rois',Wangkgs_RoiList,'doAUC',true,'inverse',Inverse1,'roiType','all');
    
    [CrossTalk2,Error2,ROISource2,ScalpData,LIST,subIDs] = mrC.Simulate.ResolutionMatrices(ProjectPath,'subSelect',subIDs,...
        'rois',Wangkgs_RoiList,'doAUC',true,'inverse',Inverse2,'roiType','all');
    save(FilePath,'CrossTalk1','Error1','ROISource1','CrossTalk2','Error2','ROISource2','ScalpData','LIST','subIDs');
else
    load(FilePath);
end

%% 
FigPath = fullfile('Figures','SourceSpace');
HemiName = {'left','right'};

% Plot Cross Talk Matrices
CrossTalk1 = cellfun(@(x) x./repmat(max(x),[size(x,1) 1]),CrossTalk1,'uni',false);% normalize
CrossTalk2 = cellfun(@(x) x./repmat(max(x),[size(x,1) 1]),CrossTalk2,'uni',false);% normalize
CT1 = (cat(3,CrossTalk1{:}));%CT1 = (CT1(1:2:end,1:2:end,:)+CT1(2:2:end,2:2:end,:))./2;
CT2 = (cat(3,CrossTalk2{:}));%CT2 = (CT2(1:2:end,1:2:end,:)+CT2(2:2:end,2:2:end,:))./2;
CT1 = CT1(ROIind,ROIind);
CT2 = CT2(ROIind,ROIind);

FIG = figure;
for i = 1:2
    %eval(['CTM1 = mean(cat(3,CrossTalk' num2str(i) '{:}),3);']);
    S = subplot(1,2,i);
    eval(['imagesc((mean(CT' num2str(i) ',3)));']);
    colormap(jmaColors('coolhot'));
    %colormap('gray')
    %caxis([-max(abs(CTMM(:))) max(abs(CTMM(:)))]);
    caxis([-1 1]);

    set(gca,'ytick',1:numel(ROILabel),'yticklabel',ROILabel,'xtick',1:numel(ROILabel),'xticklabel',ROILabel,'fontsize',10);
    if exist('xtickangle'), xtickangle(90); end
    xlabel('receiving ROI','fontsize',12);
    if i==1
        title('MN+FACE');
    else
        title('WMN+FACE');
    end
end
set(FIG,'PaperPosition',[1 1 20 8]);
set(gcf, 'Color', 'w');
set(FIG,'Units','Inch')
set(FIG,'Position',[1 1 20 9]);
print(fullfile(FigPath,'simulation',['CrossTalk_All']),'-dtiff','-r300');
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
