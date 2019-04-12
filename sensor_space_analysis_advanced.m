% The goal of this script is to
% (1) Create a project with axx trial files
% (2) Read in axx trial files from subjects
% (3) Apply decomposition approaches on the subjects

clear;
clc;
addpath(genpath('/Users/kohler/code/git/mrC'));
%% Load Axx trial files

PDiva_Path = '/Volumes/Denali_DATA1/Elham/EEG_Textscamble/';
%mrC_Path = '/Volumes/svndl/mrC_Projects/VernierLetters/';

% Read files from thge original folder
Subjfolders = subfolders(PDiva_Path,0);
Subjfolders =  Subjfolders(cellfun(@(x) ~isempty(x),strfind(Subjfolders,'nl-')));
SubIDs = cellfun(@(x) x(1:7),Subjfolders,'UniformOutput',false);

for Sub = 1:numel(Subjfolders)
axx_trialFiles = subfiles(fullfile(PDiva_Path,Subjfolders{Sub},'Exp_MATL_HCN_128_Avg','Axx*_trials.mat'),1);
    for Cond=1:length(axx_trialFiles)
        axxStrct = matfile(axx_trialFiles{Cond});
        outData{Sub,Cond} = mrC.axx.loadobj(axxStrct);
    end
end
clear axx_trialFiles axxStrct Cond Sub;
%save(fullfile('ResultData','FFT_Trial_Data'),'outData','SubIDs');

%% RCA on individual subjects and prepare data for group level RCA
%load(fullfile('ResultData','FFT_Trial_Data'));
Freqs = 0:outData{1,1}.dFHz:outData{1,1}.dFHz*(outData{1,1}.nFr-1);

for Sub = 1:numel(SubIDs)
    % merge letter and vernier conditions for 
    axx_Letter{Sub} = MergeAxx(outData(Sub,1:5));
    [decompAxx_Letter{Sub},~,A_Letter{Sub},~] = mrC.SpatialFilters.PCA(axx_Letter{Sub},'freq_range',Freqs([7 19]));
    axx_Letter{Sub} = MergeAxx(outData(Sub,1:5));
     
    axx_Vernier{Sub} = MergeAxx(outData(Sub,6:10));
    [decompAxx_Vernier{Sub},~,A_Vernier{Sub},~] = mrC.SpatialFilters.PCA(axx_Vernier{Sub},'freq_range',Freqs([7 19]));
    axx_Vernier{Sub} = MergeAxx(outData(Sub,6:10));
end


%% plot individual and average ASD results

Mode = 'Amp';%'Phase'
Ords = [1:10];
Finds  = [7 13 19 25];
Harms = {'1F1','2F1','3F1','4F1'};
%Lims = [1.5 1.5 .4 .45];
Sublist = num2cell(1:16);
Sublist{end+1} = 1:16;
load ResultData/LogMar_Val.mat
logMAR = [logMAR_letter logMAR_ver];
% figure params
FS = 12;
SizeInc = .05;

for Sub = 17: numel(Sublist)
    for i = 1:numel(Finds)%4
        FIG = figure;
        set(FIG,'unit','inch','position',[17 10 16 6])
        set(FIG,'unit','inch','paperposition',[17 10 16 6])
        clear Datamean Lim;
        for Cond = 10:-1:1
            axx_allsub = MergeAxx(outData(Sublist{Sub},Ords(Cond)));
            FFT_allsub = axx_allsub.Cos+axx_allsub.Sin*1j;
            Datamean{Cond} = abs(mean(FFT_allsub,3));
            Lim(Cond) = max(mean(Datamean{Cond}([Finds(i) Finds(i)],:)));
        end
        LimC = [max(Lim(1:5)) max(Lim(6:10))];
        for Cond = 10:-1:1
            S = subplot(2,5,Cond); mrC.Simulate.plotOnEgi(mean(Datamean{Cond}([Finds(i) Finds(i)],:))); axis tight equal;
            set(S,'position',get(S,'position')+[-SizeInc/2 -SizeInc/2-.01 SizeInc SizeInc]);
            title(['logMAR = ' num2str(round(logMAR(Cond),2))],'fontsize',FS,'fontweight','normal');
            
            if (Cond == 10)||(Cond==5)               
                SP = get(S,'position');
                h = colorbar;
                set(h,'ylim',[0 LimC(ceil(Cond/5))]);
                set(S,'position',SP);
                set(h,'position',get(h,'position')+[ .02 0 0 .0],'fontsize',FS);
                set(get(h,'label'),'String','1F1 ASD \muV','position',get(get(h,'label'),'position')+[-3.7 0 0])
            end
            
            caxis([0 LimC(ceil(Cond/5))]);
            if Cond==6
                text(-1.8,-.2,'Vernier','fontsize',FS+4,'fontweight','bold','rotation',90);
            elseif Cond ==1
                text(-1.8,-.2,'Letters','fontsize',FS+4,'fontweight','bold','rotation',90);
            end
         end
        
         colormap(jmaColors('hotcortex'));%'jet')

         if Sub<=numel(SubIDs)
             Subname = SubIDs{Sub};
         else
             Subname = 'AverageAll';
         end
         print(FIG,['../Presentation/TopoMap_individuals/TopoMap_' num2str(i) 'f1_' Subname '_' Mode '_hot'],'-r300','-dtiff');
         close all;
    end
end

%% GROUP LEVEL ANALYSIS
%SELECT the harmonic to do analysis 
analFreq = 1;
%%%%%%
%% Group Level RCA

fRCA = Finds(analFreq);
[decompAxx_Letter_all,~,A_Letter_all,D_Letter_all] = mrC.SpatialFilters.RCA(MergeAxx(axx_Letter),'freq_range',Freqs(fRCA));%,
[decompAxx_Vernier_all,~,A_Vernier_all,D_Vernier_all] = mrC.SpatialFilters.RCA(MergeAxx(axx_Vernier),'freq_range',Freqs(fRCA)); 

%% Plot RC components calculated on all subjects

% figure params
Cols = brewermap(2,'Set2');

Vinv = 1;
Linv = 1;
caxisRV = [-1*max(max(A_Vernier_all(:,1:2))) max(max(A_Vernier_all(:,1:2)))];
caxisRL = [-1*max(max(A_Letter_all(:,1:2))) max(max(A_Letter_all(:,1:2)))];
NCOMP  = 2;
FIG=figure;
set(FIG,'unit','inch','position',[17 10 8 14])
set(FIG,'unit','inch','paperposition',[17 10 8 14])
for comp = 1:NCOMP 
    subplot(4,2,comp*2),mrC.Simulate.plotOnEgi(A_Vernier_all(:,comp)*Vinv);axis tight equal;caxis(caxisRV);
    title(['RC' num2str(comp)],'fontsize',FS,'Color',Cols(comp,:));%,'fontweight','normal');%colorbar;
    if comp==1
        text(-.4, 2,'Vernier','fontsize',FS+4,'fontweight','bold')
    end
    
    subplot(4,2,comp*2-1),mrC.Simulate.plotOnEgi(A_Letter_all(:,comp)*Linv);axis tight equal;caxis(caxisRL);
    title(['RC' num2str(comp)],'fontsize',FS,'Color',Cols(comp,:));%,'fontweight','normal');%colorbar;
    if comp==1
        text(-.4, 2,'Letter','fontsize',FS+4,'fontweight','bold')
    end
end
%
Subnum = numel(SubIDs);
Condnum = 5;
nF1 = fRCA;
Ver_cos = squeeze(mean(reshape(squeeze(decompAxx_Vernier_all.Cos(nF1,:,1:end-70)),[numel(D_Vernier_all),16,Condnum,Subnum-1]),2));% last subject has 14 trails
Ver_sin = squeeze(mean(reshape(squeeze(decompAxx_Vernier_all.Sin(nF1,:,1:end-70)),[numel(D_Vernier_all),16,Condnum,Subnum-1]),2));

Ver_cos(:,:,Subnum) = squeeze(mean(reshape(squeeze(decompAxx_Vernier_all.Cos(nF1,:,end-69:end)),[numel(D_Vernier_all),14,Condnum,1]),2));
Ver_sin(:,:,Subnum) = squeeze(mean(reshape(squeeze(decompAxx_Vernier_all.Sin(nF1,:,end-69:end)),[numel(D_Vernier_all),14,Condnum,1]),2));
Ver_cmplx = Vinv*Ver_cos+(Ver_sin*1i*Vinv);


Let_cos = squeeze(mean(reshape(squeeze(decompAxx_Letter_all.Cos(nF1,:,1:end-70)),[numel(D_Letter_all),16,Condnum,Subnum-1]),2));% last subject has 14 trails
Let_sin = squeeze(mean(reshape(squeeze(decompAxx_Letter_all.Sin(nF1,:,1:end-70)),[numel(D_Letter_all),16,Condnum,Subnum-1]),2));

Let_cos(:,:,Subnum) = squeeze(mean(reshape(squeeze(decompAxx_Letter_all.Cos(nF1,:,end-69:end)),[numel(D_Letter_all),14,Condnum,1]),2));
Let_sin(:,:,Subnum) = squeeze(mean(reshape(squeeze(decompAxx_Letter_all.Sin(nF1,:,end-69:end)),[numel(D_Letter_all),14,Condnum,1]),2));
Let_cmplx = Linv*Let_cos+(Let_sin*1i*Linv);

for comp = 1:NCOMP 
    subplot(4,2,6), %plot(logMAR_ver,mean(abs(Ver_cmplx(comp,:,:)),3),'-o','Color',Cols(comp,:),'linewidth',1.5);hold on;
    errorbar(logMAR_ver,mean(abs(Ver_cmplx(comp,:,:)),3),std(squeeze(abs(Ver_cmplx(comp,:,:))),[],2)/sqrt(16),'Color',Cols(comp,:),'MarkerFaceColor',Cols(comp,:),'linewidth',1.5);hold on;
    xlim([0 .9])
    ylim([.2 .9])
    ylabel('1F1 Amplitude \muV')
    
    subplot(4,2,5), 
    errorbar(logMAR_letter,mean(abs(Let_cmplx(comp,:,:)),3),std(squeeze(abs(Let_cmplx(comp,:,:))),[],2)/sqrt(16),'Color',Cols(comp,:),'MarkerFaceColor',Cols(comp,:),'linewidth',1.5);hold on;
    xlim([0.1 1.4])
    ylim([.2 .9])
    ylabel('1F1 Amplitude \muV')
    
    subplot(4,2,8), 
    Ver_ang = wrapTo360(rad2deg(angle(Ver_cmplx(comp,:,:))));
    errorbar(logMAR_ver,mean(Ver_ang,3),std(squeeze(Ver_ang),[],2)/sqrt(16),'Color',Cols(comp,:),'MarkerFaceColor',Cols(comp,:),'linewidth',1.5);hold on;
    xlim([0 .9])
    ylim([0 360])
    xlabel('LogMAR');
    ylabel('1F1 Phase degrees')
    
    S = subplot(4,2,7);
    Let_ang = wrapTo360(rad2deg(angle(Let_cmplx(comp,:,:))));
    errorbar(logMAR_letter,mean(Let_ang,3),std(squeeze(Let_ang),[],2)/sqrt(16),'Color',Cols(comp,:),'MarkerFaceColor',Cols(comp,:),'linewidth',1.5);hold on;
    %plot(logMAR_letter,mean(wrapTo360(rad2deg(angle(Let_cmplx(comp,:,:)))),3),'-o','Color',Cols(comp,:),'linewidth',1.5);hold on;
    xlim([0.1 1.4])
    ylim([0 360])
    xlabel('LogMAR')
    ylabel('1F1 Phase degrees')
end

print(FIG,['../Presentation/TopoMap_individuals/RCA_' Harms{Finds} '_AverageAll'],'-r300','-dtiff');
close all;

%% Plot individual-level RCs 
% sublist = [1:16];
% figure;
% for s = 1:numel(sublist)
%     sub = sublist(s);
%     subplot(2,numel(sublist),s), mrC.Simulate.plotOnEgi(A_Letter{sub}(:,1)); axis tight;
%     M = max(max(abs(A_Letter{sub}(:,1:2))));
%     caxis([-M M]);
%     title(SubIDs{sub})
%     subplot(2,numel(sublist),s+numel(sublist)), mrC.Simulate.plotOnEgi(A_Letter{sub}(:,2)); axis tight
%     caxis([-M M]);
% end
% set(gcf,'unit','inch','position',[5 10 34 4.5])
% 
% figure;
% for s = 1:numel(sublist)
%     sub = sublist(s);
%     subplot(2,numel(sublist),s), mrC.Simulate.plotOnEgi(A_Vernier{sub}(:,1)); axis tight 
%     M = max(max(abs(A_Vernier{sub}(:,1:2))));
%     caxis([-M M]);
%     title(SubIDs{sub})
%     subplot(2,numel(sublist),s+numel(sublist)), mrC.Simulate.plotOnEgi(A_Vernier{sub}(:,2)); axis tight
%     caxis([-M M]);
% end
% set(gcf,'unit','inch','position',[5 5 34 4.5])


%% functions
function out_axx = MergeAxx(axxlist)
out_axx = axxlist{1};
    for C = 2:numel(axxlist)
        out_axx = out_axx.MergeTrials(axxlist{C});
    end
end

