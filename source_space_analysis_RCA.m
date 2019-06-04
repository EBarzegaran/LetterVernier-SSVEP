% The goal of this script is to
% (1) Create a project with axx trial files
% (2) Read in axx trial files from subjects
% (3) Apply decomposition approaches on the subjects

clear;
clc;
addpath(genpath('/Users/kohler/code/git/mrC'));

%% load inverses and ROI files
Path = '/Volumes/svndl/mrC_Projects/VernierLetters/source';
% Load ROIs
[ROIs,subIDROI] = mrC.Simulate.GetRoiClass(Path);

% Prepare Wang and KGS atlases
WangROIs = cellfun(@(x) x.getAtlasROIs('wang'),ROIs,'uni',false); 
KgsROIs = cellfun(@(x) x.getAtlasROIs('kgs'),ROIs,'uni',false); 
WKROIs = arrayfun(@(x) WangROIs{x}.mergROIs(KgsROIs{x}),1:numel(ROIs),'uni',false); 
ROIind = [15:18 25:52 57:60];% select a few ROIs
WKROIs = cellfun(@(x) {x.selectROIs(ROIind)},WKROIs);
WKROIsM = cellfun(@(x) x.ROI2mat(20484),WKROIs,'uni',false);

ROILabel = WKROIs{1}.getFullNames('noatlas');
%InvROIs = arrayfun(@(x) Inverses{x}*WKROIsM{x}./(sum(WKROIsM{x})+eps),1:numel(Inverses),'uni',false);

% InvName = 'WMN';
% Inverse1 = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangkgsROIsCorr_DepthWeight.inv';

InvName = 'MN';
Inverse1 = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangkgsROIsCorr.inv';

[Inverses,subIDs_Inverse] = mrC.Simulate.ReadInverses(Path,Inverse1);

%% correct the IDs
ROIs{strcmp(subIDROI,'skeri0003')}.subID = 'nl-0033';
subIDROI{strcmp(subIDROI,'skeri0003')} = 'nl-0033';
ROIs{strcmp(subIDROI,'skeri0004')}.subID = 'nl-0034';
subIDROI{strcmp(subIDROI,'skeri0004')} = 'nl-0034';

%% Load Axx trial files

PDiva_Path = '/Volumes/Denali_DATA1/Elham/EEG_Textscamble/';
%mrC_Path = '/Volumes/svndl/mrC_Projects/VernierLetters/';

% Read files from thge original folder
Subjfolders = subfolders(PDiva_Path,0);
Subjfolders =  Subjfolders(cellfun(@(x) ~isempty(x),strfind(Subjfolders,'nl-')));
SubIDs = cellfun(@(x) x(1:7),Subjfolders,'UniformOutput',false);

Ords = cell2mat(arrayfun(@(x) find(strcmp(SubIDs,subIDROI{x})),1:numel(subIDROI),'uni',false));
% ROIs = ROIs(Ords);
SubIDs = SubIDs(Ords);

for S = 1:numel(Ords)
    sub = Ords(S);
axx_trialFiles = subfiles(fullfile(PDiva_Path,Subjfolders{sub},'Exp_MATL_HCN_128_Avg','Axx*_trials.mat'),1);
    for Cond=1:length(axx_trialFiles)
        axxStrct = matfile(axx_trialFiles{Cond});
        outData{S,Cond} = mrC.axx.loadobj(axxStrct);
    end
end
clear axx_trialFiles axxStrct Cond Sub;
%save(fullfile('ResultData','FFT_Trial_Data'),'outData','SubIDs');
%% calculate source signal and apply RCA to summerize the ROI data
Harms = {'nF1','nF2','all'};
Task = {'Letter','Vernier'};

Freqs = 0:outData{1,1}.dFHz:outData{1,1}.dFHz*(outData{1,1}.nFr-1);
FreqIdx(:,1) = [7 19];
FreqIdx(:,2) = [13 25];
for sub = 1:numel(SubIDs)
    disp(['Subject:' SubIDs{sub}])
    for ts = 1:2
        disp(['Task #' Task{ts}])
        for roi = 1:numel(ROILabel)
            %disp(['ROI #' num2str(roi)])
            Sens_data = MergeAxx(outData(sub,(1:5)+((ts-1)*5)));
            InvM = Inverses{sub}(:,WKROIsM{sub}(:,roi)>0);
            Source_data = AXXtosource(Sens_data,InvM);
            for h = 1:2
                [data,W,A_all,d.(Task{ts}).(Harms{h}).('Complx'){sub,roi}] = mrC.SpatialFilters.RCA(Source_data,'freq_range',Freqs(FreqIdx(:,h)));
                WI = pinv(W);
                decompAxx.(Task{ts}).(Harms{h}).('Complx'){sub,roi} = (data.Cos+(1i*data.Sin)).*repmat(mean(WI,2)',[size(data.Cos,1) 1 size(data.Cos,3)]);
                decompAxx.(Task{ts}).(Harms{h}).('Wave'){sub,roi} = data.Wave.*repmat(mean(WI,2)',[size(data.Wave,1) 1 size(data.Wave,3)]);
                
                sourceAxx.(Task{ts}).(Harms{h}).('Complx'){sub,roi} = (Source_data.Cos+(1i*Source_data.Sin));
                sourceAxx.(Task{ts}).(Harms{h}).('Wave'){sub,roi} = Source_data.Wave;
                
            end
        end
    end
end
%%
Cs = distinguishable_colors(18);
Colors(1:2:36,:) = Cs;
Colors(2:2:36,:) = Cs;
Fs = [7 13];
HemName = {'Left','Right'};
Colors_amp = winter(5);
Colors_phase = autumn(5);

for ts = 1:2
    for h = 1:2
        for comp = 1:2
            clear temp; 
            FIG = figure;
            set(FIG,'PaperPosition',[5 5 20 12]);
            set(gcf, 'Color', 'w');
            set(FIG,'Units','Inch')
            set(FIG,'Position',[5 5 20 12]);
            for sub = 1:numel(SubIDs)
                for roi = 1:numel(ROILabel)
                 t = mean(decompAxx.(Task{ts}).(Harms{h}).('Wave'){sub,roi},3);
                 t_block = (mean(reshape(decompAxx.(Task{ts}).(Harms{h}).('Wave'){sub,roi},size(t,1),size(t,2),16,5),3));

                 tc = mean(decompAxx.(Task{ts}).(Harms{h}).('Complx'){sub,roi},3);
                 tc_block = (mean(reshape(decompAxx.(Task{ts}).(Harms{h}).('Complx'){sub,roi},size(tc,1),size(tc,2),16,5),3));

                if ~isempty(t)
    %                  temp(sub,roi,:) = mean(t(:,1:min(2,size(t,2))),2);
    %                  tempc(sub,roi,:) = mean(tc(:,1:min(2,size(t,2))),2);
                     temp(sub,roi,:) = mean(t(:,min(comp,size(t,2))),2);
                     tempc(sub,roi,:) = mean(tc(:,min(comp,size(t,2))),2);               
                     temp_bl(sub,roi,:,:) = mean(t_block(:,min(comp,size(t_block,2)),:),2);
                     tempc_bl(sub,roi,:,:) = mean(tc_block(:,min(comp,size(tc_block,2)),:),2);
                else
                    temp(sub,roi,:) = 0;
                    tempc(sub,roi,:) = 0;
                    temp_bl(sub,roi,:,:) = 0;
                    tempc_bl(sub,roi,:,:) = 0;
                end
                end
            end
            Wave_Average = squeeze(mean(temp,1));
            Cmplx_Average = squeeze(mean(tempc,1));

            for hem = 1:2
                subplot(2,3,1+((hem-1)*3)), % plot the complex values

                ind = 1;
                clear hnd
                for roi= hem:2:numel(ROILabel)
                    if abs(Cmplx_Average(roi,Fs(h)))>=nanmean(abs(Cmplx_Average(:,Fs(h))))
                        hnd(ind)=line([0 real(Cmplx_Average(roi,Fs(h)))],[0 imag(Cmplx_Average(roi,Fs(h)))],'Color',Colors(roi,:),'linewidth',2);
                        RoiInd(ind) = roi;
                        ind = ind+1;
                    end
                end
                Leg = legend(hnd,ROILabel(RoiInd));
                set(Leg,'position',get(Leg,'position')+[.03 0 0 0])
                xlim([-max(abs(Cmplx_Average(:,Fs(h)))) max(abs(Cmplx_Average(:,Fs(h))))]);
                ylim([-max(abs(Cmplx_Average(:,Fs(h)))) max(abs(Cmplx_Average(:,Fs(h))))]);
                axis equal
                title([HemName{hem} ' hemisphere'],'fontsize',12)

                S1 = subplot(4,2,2+(hem-1)*4); % plot the amplitudes
                B = bar(abs(squeeze(mean(tempc_bl(:,hem:2:end,Fs(h),:),1))),'FaceColor','flat');
                for k = 1:5
                     B(k).CData = Colors_amp(k,:);
                end
                ylabel('Amplitude (\mu V)')
                ylim([0 max(max(mean(abs(tempc_bl(:,:,Fs(h),:)))))]);
                if hem==1
                    T = title([Task{ts} ' - ' Harms{h}(2:end) ' - ' InvName],'fontsiz',16);
                    %set(T,'position',get(T,'position')+[-10 1 0])
                end
                
                S2 = subplot(4,2,4+(hem-1)*4); % plot the amplitudes
                B = bar(-1*wrapTo2Pi(angle(squeeze(mean(tempc_bl(:,hem:2:end,Fs(h),:),1)))),'FaceColor','flat');
                for k = 1:5
                     B(k).CData = Colors_phase(k,:);
                end
                set(S1,'position',get(S1,'position')+[-.1 -.02 .1 .02],'xticklabel',[]);
                set(S2,'position',get(S2,'position')+[-.1 .02 .1 .02],'xtick',1:numel(ROILabel(hem:2:numel(ROILabel))),...
                    'xticklabel',ROILabel(hem:2:numel(ROILabel)),'ytick',-2*pi:pi/2:0,'yticklabel',{'2\pi','','\pi','','0'});
                xtickangle(90);
                ylim([-2*pi+eps 0]);
                ylabel('Phase (Radian)')

            end
            print(fullfile('Figures','SourceSpace','Source1_ROIRCA',[Task{ts} '_' Harms{h} '_sources_RC' num2str(comp) '_' InvName '.tif']),'-dtiff','-r300');
            close;
        end
    end
end

