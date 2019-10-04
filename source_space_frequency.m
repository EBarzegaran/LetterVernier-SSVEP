% This script is for preliminary analysis and visualization of vernier-text
% scramble task in source space
% Elham Barzegaran, 5.22.2018
% modified: EB, 4.19.2019

clear;clc;

PATH = '/Users/elhamb/Documents/Codes/Git/mrC';
addpath(genpath('/Users/elhamb/Documents/Codes/Git/EEGSourceSim'));
addpath(genpath(PATH));
path2 = '/Users/elhamb/Documents/Codes/NonGit/mrTools/mrUtilities/mrFlatMesh';
addpath(genpath(path2));

addpath(genpath('/Users/kohler/code/git/sweepAnalysis/functions/helper'));

%% Convert spectrum data into source space
Path = '/Users/elhamb/Documents/Data/TextScramble_exp1/VernierLetters/source';
Inverse = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangkgsROIsCorr.inv';
%Inverse = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangkgsROIsCorr_DepthWeight.inv';
[outData,~,FreqFeatures,subIDs] = mrC.SourceBrain(Path,Inverse,'domain','frequency','doSmooth' , true);%,'template','nl-0014');

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
Lateral = {'V','D','VD'};
l = 3;

% Load ROIs
[ROIs,subIDROI] = ESSim.Simulate.GetRoiClass(Path);
Ords = cell2mat(arrayfun(@(x) find(strcmp(subIDROI,subIDs{x})),1:numel(subIDs),'uni',false));
ROIs = ROIs(Ords);

% Load Wang and KGS atlases
WangROIs = cellfun(@(x) x.getAtlasROIs('wang'),ROIs,'uni',false); 
KgsROIs = cellfun(@(x) x.getAtlasROIs('kgs'),ROIs,'uni',false); 
WKROIs = arrayfun(@(x) WangROIs{x}.mergROIs(KgsROIs{x}),1:numel(ROIs),'uni',false); 



switch l
    case 1 % Ventral V1-V3 and IOG and VWFA
        ROIind = [31:32 35:36 43:44 51:52 57:58];          
        Wangkgs_RoiList = cellfun(@(x) {x.selectROIs(ROIind)},WKROIs);
        ROI_Final = Wangkgs_RoiList;
    case 2 % Dorsal V1-V3 and IOG and VWFA
        ROIind = [[31:32 35:36 43:44]-2 [51:52 57:58]];          
        Wangkgs_RoiList = cellfun(@(x) {x.selectROIs(ROIind)},WKROIs);
        ROI_Final = Wangkgs_RoiList;
    case 3 % Both dorsa and Ventral
        ROIind = [29:36 41:44 51:52 57:58];                
        Wangkgs_RoiList = cellfun(@(x) {x.selectROIs(ROIind)},WKROIs);
        Merges = {[1 3],[2 4], [5 7], [6 8], [9 11], [10 12], 13, 14, 15, 16}; % Merge ROIs
        Temp_Names = {'V1','V1','V2','V2','V3','V3','IOG','IOG','VWFA','VWFA'};
        ROI_Final = cellfun(@(x) x.mergeIndROIs(Merges,Temp_Names),Wangkgs_RoiList,'uni',false);

end
 

WKROIsM = cellfun(@(x) x.ROI2mat(20484),ROI_Final,'uni',false);
ROILabelSel = ROI_Final{1}.getFullNames('noatlas');
%InvROIs = arrayfun(@(x) Inverses{x}*WKROIsM{x}./(sum(WKROIsM{x})+eps),1:numel(Inverses),'uni',false);

%% plot amps
Finds  = [7];

for h = 1:1
    for ts = 1:2
       FIG = figure;
       set(FIG,'unit','inch','position',[3 5 20 8]);
       
       for cond = 1:5
           FData_temp = squeeze(FreqData.(Task{ts})(:,Finds(h),cond,:));
           FDataROI_temp = arrayfun(@(x) (FData_temp(:,x)).'*WKROIsM{x}./sum(WKROIsM{x}),1:numel(subIDs),'uni',false);
           FDataROI.(Task{ts}).(Harms{h})(:,:,cond) = cat(1,FDataROI_temp{:});
       end
       
       PData = squeeze(mean(FDataROI.(Task{ts}).(Harms{h})(:,:,:),1));
       %Results.('Amp').(Task{ts}).(Harms{h}) = abs(PData);
       subplot(2,1,1),bar(abs(PData(1:2:end,:)));
       set(gca,'xtick',1:22,'xticklabel',ROILabelSel(1:2:end));
       title([Task{ts} '_' Harms{h}])
       ylim([0 8])
       
       subplot(2,1,2),bar(abs(PData(2:2:end,:)));
       set(gca,'xtick',1:22,'xticklabel',ROILabelSel(2:2:end));
       ylim([0 8])
    end
end


%% plot phase

for h = 1:1
    for ts = 1:2
       FIG = figure;
       set(FIG,'unit','inch','position',[3 5 20 8]);
       
       for cond = 1:5
           FData_temp = squeeze(FreqData.(Task{ts})(:,Finds(h),cond,:));
           %FDataRoi = arrayfun(@(x) mapMtx{x}*FData_temp(:,x),1:numel(subIDs),'uni',false);
           FDataROI_temp = arrayfun(@(x) FData_temp(:,x).'*WKROIsM{x}./sum(WKROIsM{x}),1:numel(subIDs),'uni',false);
           FDataROI.(Task{ts}).(Harms{h})(:,:,cond) = cat(1,FDataROI_temp{:});
           
       end
       
       PData = squeeze(nanmean(FDataROI.(Task{ts}).(Harms{h})(:,:,:),1));
       PData2 = squeeze(nanstd(FDataROI.(Task{ts}).(Harms{h})(:,:,:),[],1))/sqrt(10);
       
       Results.(Task{ts}).(Harms{h}).M = (PData);
       Results.(Task{ts}).(Harms{h}).S = (PData2);
       %Results.('Phase').(Task{ts}).(Harms{h})(Results.('Phase').(Task{ts}).(Harms{h})>5.5) = 0;
       subplot(2,1,1),bar(wrapTo2Pi(angle(Results.(Task{ts}).(Harms{h}).M(1:2:end,:))));
       set(gca,'xtick',1:22,'xticklabel',ROILabelSel(1:2:end));
       title([Task{ts} '_' Harms{h}])
       ylim([0 2*pi])
       
       subplot(2,1,2),bar(abs(PData(2:2:end,:)));
       set(gca,'xtick',1:22,'xticklabel',ROILabelSel(2:2:end));
       ylim([0 7])
    end
end

%close all;
save(fullfile('ResultData',['ROISourceResults_freq_' Lateral{l}]),'Results','ROILabelSel');


%% sort ROIs according to phase
InvName = 'MN';
Cs = distinguishable_colors(18);
Colors(1:2:36,:) = Cs;
Colors(2:2:36,:) = Cs;
Colors_amp = winter(5);
Colors_phase = autumn(5);

ROILabelSel2 = ROILabelSel(1:end);
HemName = {'left','right'};
for ts = 1:2
    for h = 1:2
        FIG = figure;
        set(FIG,'PaperPosition',[5 5 20 12]);
        set(gcf, 'Color', 'w');
        set(FIG,'Units','Inch')
        set(FIG,'Position',[5 5 20 12]);
        R = Results.(Task{ts}).(Harms{h}).M;
        RS = Results.(Task{ts}).(Harms{h}).S;
        RM = mean(R,2);
        for hem = 1:2
            subplot(2,3,1+((hem-1)*3)), % plot the complex values

            ind = 1;
            clear hnd
            for roi= hem:2:numel(ROILabelSel2)
                if abs(RM(roi))>=nanmean(abs(RM(:)))
                    hnd(ind)=line([0 real(RM(roi))],[0 imag(RM(roi))],'Color',Colors(roi,:),'linewidth',2);
                    RoiInd(ind) = roi;
                    ind = ind+1;
                end
            end
            Leg = legend(hnd,ROILabelSel2(RoiInd));
            set(Leg,'position',get(Leg,'position')+[.03 0 0 0])
            xlim([-max(abs(RM)) max(abs(RM))]);
            ylim([-max(abs(RM)) max(abs(RM))]);
            axis equal
            title([HemName{hem} ' hemisphere'],'fontsize',12)

            S1 = subplot(4,2,2+(hem-1)*4); % plot the amplitudes
            B = bar(abs(squeeze(R(hem:2:end,:))),'FaceColor','flat');
            hold on;
            %errorbar()
            for k = 1:5
                 B(k).CData = Colors_amp(k,:);
            end
            ylabel('Amplitude (\mu V)')
            ylim([0 max(max((abs(R))))]*1.1);
            if hem==1
                T = title([Task{ts} ' - ' Harms{h}(2:end) ' - ' InvName],'fontsiz',16);
                set(T,'position',get(T,'position')+[-10 1 0])
            end
            
            S2 = subplot(4,2,4+(hem-1)*4); % plot the amplitudes
            B = bar(-1*wrapTo2Pi(angle(squeeze(R(hem:2:end,:)))),'FaceColor','flat');
            for k = 1:5
                 B(k).CData = Colors_phase(k,:);
            end
            
            set(S1,'position',get(S1,'position')+[-.1 -.02 .1 .02],'xticklabel',[]);
            set(S2,'position',get(S2,'position')+[-.1 .02 .1 .02],'xtick',1:numel(ROILabelSel2(hem:2:numel(ROILabelSel2))),...
                'xticklabel',ROILabelSel2(hem:2:numel(ROILabelSel2)),'ytick',-2*pi:pi/2:0,'yticklabel',{'2\pi','','\pi','','0'});
            xtickangle(90);
            ylim([-2*pi+eps 0]);
            ylabel('Phase (Radian)')

        end 
        print(fullfile('Figures','SourceSpace',[Task{ts} '_' Harms{h}(2:end) '_sources_' InvName '.tif']),'-dtiff','-r300');
        close;
    end
end
% Then take into account the amplitude

