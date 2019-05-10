% The goal of this script is to
% (1) Create a project with axx trial files
% (2) Read in axx trial files from subjects
% (3) Apply decomposition approaches on the subjects

clear;
clc;
addpath(genpath('/Users/kohler/code/git/mrC'));
%% Load Axx trial files

PDiva_Path = '/Volumes/Denali_DATA1/2018_TextScramble/experiment2';
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
%%

%axx_allsub.Cos+axx_allsub.Sin*1j;
Sublist = num2cell(1:numel(SubIDs));
Sublist{end+1} = 1:numel(SubIDs);

f1Idx = outData{1,1}.i1F1; % index of 1f1
f2Idx = outData{1,1}.i1F2;% index of 1f2

nF1C = [1 2 4 5 ]; % nF1 clean harmonics 2 4 5 7 8 
nF2 = [1 2 ]; % nF2 harmonics 2 3 

nFC_Idx{1} = nF1C * f1Idx +1;
nFC_Idx{2} = nF2 * f2Idx +1;



%%
Condnum = size(outData,2);
SubNames = SubIDs;SubNames{end+1} = 'Average_all';

for Sub = 1: numel(Sublist)
    FIG = figure;
    for Cond = 1: Condnum
        for f = 1:2
            if numel(Sublist{Sub})==1
                if isempty(outData{Sublist{Sub},Cond})
                    continue
                end
            end
            
            axx_allsub = MergeAxx(outData(Sublist{Sub},Cond));
            FFT_allsub = axx_allsub.Cos+axx_allsub.Sin*1j;
            Datamean{Cond} = abs(mean(FFT_allsub,3));
            Lim(Cond,f) = max(mean(Datamean{Cond}(nFC_Idx{f},:)));
            S = subplot(2,size(outData,2),Cond+(f-1)*Condnum); mrC.Simulate.plotOnEgi(mean(Datamean{Cond}(nFC_Idx{f},:))); axis tight equal;
            caxis([0 Lim(Cond,f)]);
            colorbar
            if Cond ==1
                if f==1
                    TX1 =text(-1.8,-.6,'nF1 clean','fontsize',14,'rotation',90);
                    title('Letter Flip Task')    
                else
                    TX2 =text(-1.8,-.3,'nF2','fontsize',14,'rotation',90);
                end
            elseif Cond==2
                if f==1
                    title('Color Task')
                end
            elseif Cond==3
                if f==1
                    title('Fixation - No Taks')
                end
                for c = 1:Condnum
                    subplot(2,size(outData,2),c+(f-1)*Condnum), caxis([0 max(Lim(:,f))])
                end
            end

        end
    end
    colormap('jet')
    set(FIG,'unit','inch','Paperposition',[1 1 12 6],'position',[1 1 12 6])
    print(FIG,['Figures/EXPERIMENT2/TopoMap_' SubNames{Sub}],'-r300','-dtiff');
    close all;
end
   
%%
Condnum = size(outData,2);
Task = {'FlipLetter','ColorTask','Fixation'};
Subnum = numel(SubIDs);
elecnum = 128;
Harms = {'nF1C','nF2'};



for ts = 1:3
    Maxx{ts} = MergeAxx(outData(:,ts));
    for h = 1:2
        TCmplx = zeros(840,128,Subnum); % fft time window = 2*FS
        TCos = squeeze(mean(reshape(squeeze(Maxx{ts}.Cos(nFC_Idx{h},:,:)),[numel(nFC_Idx{h}),elecnum,40,Subnum]),3));% last subject has 14 trails
        TSin = squeeze(mean(reshape(squeeze(Maxx{ts}.Sin(nFC_Idx{h},:,:)),[numel(nFC_Idx{h}),elecnum,40,Subnum]),3));

        TCmplx(nFC_Idx{h},:,:) = TCos+(TSin*1i);
        TCmplx(end:-1:422,:,:) = conj(TCmplx(2:420,:,:));
        RWave = ifft(mean(conj(TCmplx),3),840,1);
        
        RWaveT =  RWave(1:420,:)+RWave(421:end,:) ;

        RWave_all.(Task{ts}).(Harms{h}) = RWaveT/2*840/2;
        clear RWaveT RWave TCos TSin TCmplx;
    end
end

%% make the video
Time = (1:420)*2.38095;
FIG2 = figure;
set(FIG2,'unit','inch','Paperposition',[2 2 12 6],'position',[5 4 12 6])
Wave(:,:,1) = RWave_all.(Task{1}).(Harms{1});
Wave(:,:,2) = RWave_all.(Task{2}).(Harms{1});
Wave(:,:,3) = RWave_all.(Task{3}).(Harms{1});
Wave(:,:,4) = RWave_all.(Task{1}).(Harms{2});
Wave(:,:,5) = RWave_all.(Task{2}).(Harms{2});
Wave(:,:,6) = RWave_all.(Task{3}).(Harms{2});

vidfile = VideoWriter(['Figures/EXPERIMENT2/Average_timecourse.mp4'],'MPEG-4');
vidfile.FrameRate = 10;
open(vidfile);

for i = 1:420
    for sub = 1:6
        subplot(2,3,sub)
        mrC.plotOnEgi(Wave(i,:,sub));
        colormap('jet');
        caxis([-max(max(Wave(:,:,sub))) max(max(Wave(:,:,sub)))]/1.2);
        if sub ==1
            title('Letter Flip Task')
            if exist('TX1')
                delete(TX1)
            end
            TX1 =text(-1.8,-.6,'nF1 clean','fontsize',14,'rotation',90);
        elseif sub==2
            title('Color Task')
        elseif sub==3
            title('Fixation-No Task')
        elseif sub==5
            if exist('TX')
                delete(TX)
            end
            TX =text(-0.8,-1.8,['time = ' num2str(round(Time(i),2))],'fontsize',14);
        elseif sub==4
            if exist('TX2')
                delete(TX2)
            end
            TX2 =text(-1.8,-.3,'nF2','fontsize',14,'rotation',90);
        end
    end
    
    pause(.01)
    F(i)= getframe(gcf);
    writeVideo(vidfile,F(i));
end
close(vidfile);
close all

%% RCA analysis on individuals
SUBJ = 1;
Freqs = 0:outData{1,1}.dFHz:outData{SUBJ,1}.dFHz*(outData{1,1}.nFr-1);
Cols = brewermap(5,'Dark2');
Cols = Cols(2:5,:);
numcom = 3;


for ts = 1:3
    Maxx{ts} = MergeAxx(outData(:,ts));
    for h = 1:2
        [decompAxx_ind.(Task{ts}).(Harms{h}),~,A_ind.(Task{ts}).(Harms{h}),D] = mrC.SpatialFilters.RCA(Maxx{ts},'freq_range',Freqs(nFC_Idx{h}));
    end
    FIG = figure;
    for comp = 1:numcom
        subplot(2,numcom+1,comp),mrC.plotOnEgi(A_ind.(Task{ts}).(Harms{1})(:,comp));axis tight equal
        if comp==1
            TX2 =text(-1.8,-.3,Harms{1},'fontsize',14,'rotation',90);
        end
        title(['RC' num2str(comp)],'color',Cols(comp,:));
        % Plot the RC 2D phases
        C = mean(decompAxx_ind.(Task{ts}).(Harms{1}).Cos(nFC_Idx{1},1:numcom,:),3);
        S = mean(decompAxx_ind.(Task{ts}).(Harms{1}).Sin(nFC_Idx{1},1:numcom,:),3);
        subplot(2,numcom+1,numcom+1),
        line([0 C(1,comp)],[0 S(1,comp)],'Color',Cols(comp,:),'linewidth',2);hold on;
        M1(comp) = max(max(max(abs(C(1:2,:)))),max(max(abs(S(1:2,:)))))*1.1;

        axis equal
        if comp ==numcom
            M = max(M1);
            xlim([-M M]);ylim([-M M]);
        end
        %%%%%%%

        subplot(2,numcom+1,comp+numcom+1),mrC.plotOnEgi(A_ind.(Task{ts}).(Harms{2})(:,comp));axis tight
        if comp==1
            TX2 =text(-1.8,-.3,Harms{2},'fontsize',14,'rotation',90);
        end

        % plot the RC 2D phases
        C = mean(decompAxx_ind.(Task{ts}).(Harms{2}).Cos(nFC_Idx{2},1:numcom,:),3);
        S = mean(decompAxx_ind.(Task{ts}).(Harms{2}).Sin(nFC_Idx{2},1:numcom,:),3);
        subplot(2,numcom+1,(numcom+1)*2),
        line([0 C(1,comp)],[0 S(1,comp)],'Color',Cols(comp,:),'linewidth',2);hold on;
        M2(comp) = max(max(max(abs(C(1:2,:)))),max(max(abs(S(1:2,:)))))*1.1;
        %xlim([-M M]);ylim([-M M]);
        axis equal
        if comp ==numcom
            M = max(M2);
            xlim([-M M]);ylim([-M M]);
        end
        %%%%%%
    end
    colormap('jet')
    set(FIG,'unit','inch','Paperposition',[1 1 12 6],'position',[1 1 12 6])
    print(FIG,['Figures/EXPERIMENT2/RCA_AverageAll_' Task{ts}],'-r300','-dtiff');
    close all;
end





