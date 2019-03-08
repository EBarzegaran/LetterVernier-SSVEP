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
    axx_Letter{Sub} = MergeAxx(outData(Sub,4:5));
    [decompAxx_Letter{Sub},~,A_Letter{Sub},~] = mrC.SpatialFilters.PCA(axx_Letter{Sub},'freq_range',Freqs([7 19]));
    axx_Letter{Sub} = MergeAxx(outData(Sub,1:5));
     
    axx_Vernier{Sub} = MergeAxx(outData(Sub,6:10));
    [decompAxx_Vernier{Sub},~,A_Vernier{Sub},~] = mrC.SpatialFilters.PCA(axx_Vernier{Sub},'freq_range',Freqs([7 19]));
    axx_Vernier{Sub} = MergeAxx(outData(Sub,6:10));
end


%% Group Level RCA
fRCA = 19;
[decompAxx_Letter_all,~,A_Letter_all,D_Letter_all] = mrC.SpatialFilters.RCA(MergeAxx(axx_Letter),'freq_range',Freqs(fRCA));%,
[decompAxx_Vernier_all,~,A_Vernier_all,D_Vernier_all] = mrC.SpatialFilters.RCA(MergeAxx(axx_Vernier),'freq_range',Freqs(fRCA));

%% plot average ASD results
Ords = [1 6 2 7 3 8 4 9 5 10];
% for sub = 1:11
Finds  = [7 13 19 25];
Lims = [1 1.5 .4 .45];
for i = 3:3
    FIG = figure;
    for Cond = 1:10
        axx_allsub = MergeAxx(outData(:,Ords(Cond)));
        FFT_allsub = axx_allsub.Cos+axx_allsub.Sin*1j;
        AMP = abs(mean(FFT_allsub,3));
        %AMP = mean(abs(FFT_allsub),3);
        %subplot(5,2,Cond), mrC.Simulate.plotOnEgi(mean(AMP([7 19],:))); axis tight;caxis([-1.5 1.5])
        subplot(5,2,Cond), mrC.Simulate.plotOnEgi(mean(AMP([Finds(i) Finds(i)],:))); axis tight;caxis([-Lims(i)/2 Lims(i)])
        if Cond ==1
            title('Letter');
        elseif Cond==2
            title('Vernier');
        end
    end
    h = colorbar;
    set(h,'ylim',[0 Lims(i)]);
    set(gcf,'unit','inch','position',[17 10 8 15])
     set(gcf,'unit','inch','paperposition',[17 10 8 15])
    print(FIG,['../Presentation/TopoMap_' num2str(i) 'f1'],'-r300','-dtiff');
    
end

 
%% Plot RC components calculated on all subjects
caxisR2 = [-1.5 1.5];
caxisR1 = caxisR2;
NCOMP  = 2;
FIG=figure;
for comp = 1:NCOMP 
    subplot(NCOMP,2,comp*2),mrC.Simulate.plotOnEgi(A_Vernier_all(:,comp));axis tight;%caxis(caxisR2);
    title(['Vernier-Comp' num2str(comp)]);colorbar;
    
    subplot(NCOMP,2,comp*2-1),mrC.Simulate.plotOnEgi(A_Letter_all(:,comp));axis tight;%caxis(caxisR2);
    title(['Letter-Comp' num2str(comp)]);colorbar;
end
set(gcf,'unit','inch','position',[1 10 8 4*NCOMP]);
set(gcf,'unit','inch','paperposition',[1 10 8 4*NCOMP]);

print(FIG,['../Presentation/RCA_TopoMap_' num2str((fRCA-1)/6) 'f1'],'-r300','-dtiff');

%% Plot group level RCs and their 2D phase plots
Subnum = numel(SubIDs);
Condnum = 5;
nF1 = fRCA;
Ver_cos = squeeze(mean(reshape(squeeze(decompAxx_Vernier_all.Cos(nF1,:,1:end-70)),[128,16,Condnum,Subnum-1]),2));% last subject has 14 trails
Ver_sin = squeeze(mean(reshape(squeeze(decompAxx_Vernier_all.Sin(nF1,:,1:end-70)),[128,16,Condnum,Subnum-1]),2));

Ver_cos(:,:,Subnum) = squeeze(mean(reshape(squeeze(decompAxx_Vernier_all.Cos(nF1,:,end-69:end)),[128,14,Condnum,1]),2));
Ver_sin(:,:,Subnum) = squeeze(mean(reshape(squeeze(decompAxx_Vernier_all.Sin(nF1,:,end-69:end)),[128,14,Condnum,1]),2));

Colors = {'k','g','b'};
FIG=figure;
for cond = 1:Condnum
    subplot(5,2,2*cond);
    for comp = 1:2
        if comp==1
            l(comp) = line([0 mean(Ver_cos(comp,cond,:),3)]*-1,[0 mean(Ver_sin(comp,cond,:),3)]*-1,'color',Colors{comp},'LineWidth',1.5);hold on;
            [ampDiff,phaseDiff,zSNR,errorEllipse] = fitErrorEllipse([squeeze(Ver_cos(comp,cond,:)) squeeze(Ver_sin(comp,cond,:))]*-1,'SEM');
        else
            l(comp) = line([0 mean(Ver_cos(comp,cond,:),3)],[0 mean(Ver_sin(comp,cond,:),3)],'color',Colors{comp},'LineWidth',1.5);hold on;
            [ampDiff,phaseDiff,zSNR,errorEllipse] = fitErrorEllipse([squeeze(Ver_cos(comp,cond,:)) squeeze(Ver_sin(comp,cond,:))],'95CI');
        end
        
        plot(errorEllipse(:,1), errorEllipse(:,2),'k-','LineWidth',1.5,'color',Colors{comp}) ;
    end
    if cond ==1
        title('Vernier');
        legend(l,{'comp1','comp2'});
    end
    xlim([-1 1]);
    ylim([-1 1]);
end

Let_cos = squeeze(mean(reshape(squeeze(decompAxx_Letter_all.Cos(nF1,:,1:end-70)),[128,16,Condnum,Subnum-1]),2));% last subject has 14 trails
Let_sin = squeeze(mean(reshape(squeeze(decompAxx_Letter_all.Sin(nF1,:,1:end-70)),[128,16,Condnum,Subnum-1]),2));

Let_cos(:,:,Subnum) = squeeze(mean(reshape(squeeze(decompAxx_Letter_all.Cos(nF1,:,end-69:end)),[128,14,Condnum,1]),2));
Let_sin(:,:,Subnum) = squeeze(mean(reshape(squeeze(decompAxx_Letter_all.Sin(nF1,:,end-69:end)),[128,14,Condnum,1]),2));

for cond = 1:5
    subplot(5,2,2*cond-1);
    for comp = 1:2
        l(comp) = line([0 mean(Let_cos(comp,cond,:),3)],[0 mean(Let_sin(comp,cond,:),3)],'color',Colors{comp},'LineWidth',1.5);hold on;
        [ampDiff,phaseDiff,zSNR,errorEllipse] = fitErrorEllipse([squeeze(Let_cos(comp,cond,:)) squeeze(Let_sin(comp,cond,:))],'95CI');
        plot(errorEllipse(:,1), errorEllipse(:,2),'k-','LineWidth',1.5,'color',Colors{comp}) ;
    end   
    if cond ==1
        title('Letter')
        legend(l,{'comp1','comp2'});
    end
    xlim([-1 1]);
    ylim([-1 1]);
end

set(gcf,'unit','inch','position',[9 10 8 15]);
print(FIG,['../Presentation/RCA_phase_' num2str((fRCA-1)/6) 'f1'],'-r300','-dtiff');
%% Plot individual-level RCs 
sublist = [1:16];
figure;
for s = 1:numel(sublist)
    sub = sublist(s);
    subplot(2,numel(sublist),s), mrC.Simulate.plotOnEgi(A_Letter{sub}(:,1)); axis tight;
    M = max(max(abs(A_Letter{sub}(:,1:2))));
    caxis([-M M]);
    title(SubIDs{sub})
    subplot(2,numel(sublist),s+numel(sublist)), mrC.Simulate.plotOnEgi(A_Letter{sub}(:,2)); axis tight
    caxis([-M M]);
end
set(gcf,'unit','inch','position',[5 10 34 4.5])

figure;
for s = 1:numel(sublist)
    sub = sublist(s);
    subplot(2,numel(sublist),s), mrC.Simulate.plotOnEgi(A_Vernier{sub}(:,1)); axis tight 
    M = max(max(abs(A_Vernier{sub}(:,1:2))));
    caxis([-M M]);
    title(SubIDs{sub})
    subplot(2,numel(sublist),s+numel(sublist)), mrC.Simulate.plotOnEgi(A_Vernier{sub}(:,2)); axis tight
    caxis([-M M]);
end
set(gcf,'unit','inch','position',[5 5 34 4.5])


%% functions
function out_axx = MergeAxx(axxlist)
out_axx = axxlist{1};
    for C = 2:numel(axxlist)
        out_axx = out_axx.MergeTrials(axxlist{C});
    end
end

