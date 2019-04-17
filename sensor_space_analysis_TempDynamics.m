
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
Task = {'Letter','Vernier'};

for Sub = 1:numel(SubIDs)
    % merge letter and vernier conditions for 
    axx.(Task{1}){Sub} = MergeAxx(outData(Sub,1:5));
    %[decompAxx_ind.(Task{1}){Sub},~,A_ind.(Task{1}){Sub},~] = mrC.SpatialFilters.PCA(axx.(Task{1}){Sub},'freq_range',Freqs([7 19]));
     
    axx.(Task{2}){Sub} = MergeAxx(outData(Sub,6:10));
    %[decompAxx_ind.(Task{2}){Sub},~,A_ind.(Task{2}){Sub},~] = mrC.SpatialFilters.PCA(axx.(Task{2}){Sub},'freq_range',Freqs([7 19]));
end

clear Subjfolders

%% time course of even and odd harmonics by reconstructing the fourier coefs
elecnum= 128;
Condnum = 5;
Subnum = numel(SubIDs);


FreqIdxall = 6*(1:16)+1;
Harms = {'Odd','Even','All'};
OEFreqIdx = [FreqIdxall(1:2:end);FreqIdxall(2:2:end)];
for ts = 1:numel(Task)
    Maxx.(Task{ts}) = MergeAxx(axx.(Task{ts}));
    for h = 1:2
        TCmplx = zeros(840,128,5,16); % fft time window = 2*FS
        TCos = squeeze(mean(reshape(squeeze(Maxx.(Task{ts}).Cos(OEFreqIdx(h,:),:,1:end-70)),[size(OEFreqIdx,2),elecnum,16,Condnum,Subnum-1]),3));% last subject has 14 trails
        TSin = squeeze(mean(reshape(squeeze(Maxx.(Task{ts}).Sin(OEFreqIdx(h,:),:,1:end-70)),[size(OEFreqIdx,2),elecnum,16,Condnum,Subnum-1]),3));

        TCos(:,:,:,Subnum) = squeeze(mean(reshape(squeeze(Maxx.(Task{ts}).Cos(OEFreqIdx(h,:),:,end-69:end)),[size(OEFreqIdx,2),elecnum,14,Condnum,1]),3));
        TSin(:,:,:,Subnum) = squeeze(mean(reshape(squeeze(Maxx.(Task{ts}).Sin(OEFreqIdx(h,:),:,end-69:end)),[size(OEFreqIdx,2),elecnum,14,Condnum,1]),3));
        TCmplx(OEFreqIdx(h,:),:,:,:) = TCos+(TSin*1i);
        TCmplx(end:-1:422,:,:,:) = conj(TCmplx(2:420,:,:,:));
        RWave = ifft(conj(TCmplx),840,1);
        Idxs = 0:140:840;
        RWaveT =  RWave(Idxs(1)+1:Idxs(1+1),:,:,:) ;
        
        for i = 1:numel(Idxs)-1
            RWaveT = RWaveT+RWave(Idxs(i)+1:Idxs(i+1),:,:,:);
        end
        RWave_all.(Task{ts}).(Harms{h}) = RWaveT/6*840/2;
        clear RWaveT RWave TCos TSin TCmplx;
    end
    MWave = squeeze(mean(reshape(squeeze(Maxx.(Task{ts}).Wave(:,:,1:end-70)),[140,elecnum,16,Condnum,Subnum-1]),3));% last subject has 14 trails
    MWave(:,:,:,Subnum) = squeeze(mean(reshape(squeeze(Maxx.(Task{ts}).Wave(:,:,end-69:end)),[140,elecnum,14,Condnum,1]),3));
    RWave_all.(Task{ts}).(Harms{3}) = MWave;
end
%% load stim images
load('../Stimuli/SloanTextScrambles_devel.mat','img');
load('../Stimuli/Vernier_Example.mat');

Img.(Task{1}) = double(squeeze(img(:,:,1,5,1,:)));
Img.(Task{2}) = cat(3,Im1,Im2);
clear Im1 Im2 img;
%% Prepare videos
ImNames.(Task{1})={'Letters','Scrambled Letters'};
ImNames.(Task{2})={'Uniform','Segmented'};
Generatevids = true;
StimPlot = 0;
if Generatevids
    for ts = 1:numel(Task)
        for cond = 1:Condnum
            if StimPlot
                vidfile = VideoWriter(['Figures/' Task{ts} '_Cond' num2str(cond) '_Stim.mp4'],'MPEG-4');
            else
                vidfile = VideoWriter(['Figures/' Task{ts} '_Cond' num2str(cond) '.mp4'],'MPEG-4');
            end
            vidfile.FrameRate = 10;
            open(vidfile);
            FIG = figure;
            if StimPlot
                set(FIG,'unit','inch','position',[5 5 10 10]);
            else
                set(FIG,'unit','inch','position',[5 5 10 5])
            end
                
            for h = 1:numel(Harms)
                M(h) = max(max(abs(squeeze(mean(RWave_all.(Task{ts}).(Harms{h})(:,:,cond,:),4)))));
            end
            for i = 1:140
                for h = 1:2%numel(Harms)
                    subplot(StimPlot+1,2,h);
                    mrC.plotOnEgi(squeeze(mean(RWave_all.(Task{ts}).(Harms{h})(i,:,cond,:),4))*-1); colormap('jet'); caxis([-M(h) M(h)]);
                    T = title([Harms{h} ' harmonics'],'fontsize',12);
                    if h==2, H = text(-2.5,-1.8,['Time = ' num2str(round(2.381*i,1)) ' ms'],'fontsize',14); end
                end
                % plot stimulus
                if StimPlot
                    S = subplot(2,1,2);
                    imshow(Img.(Task{ts})(:,:,ceil(i/70))); 
                    axis equal;
                    set(S,'position',get(S,'position')+[0.02 0 -.08 -.08])
                    title('Stimulus','fontsize',14)
                    
                    % plot time points
                    axes('Position',[0.1 0.42 0.85 0.1],'next','add'); axis off               
                    line([0 .45],[.5 .5],'linewidth',10,'color',[.2 .4 .4])               
                    line([.45 .9],[.5 .5],'linewidth',10,'color',[.5 .5 .8])
                    A = annotation('rectangle',[(i/140)*.85+.1 .44 .01 .05],'FaceColor','k','FaceAlpha',1);
                    axes('Position',[0.1 0.40 0.85 0.1],'next','add'); axis off
                    text(.2,.3,ImNames.(Task{ts}){1},'fontsize',14,'color',[.2 .4 .4]);
                    text(.6,.3,ImNames.(Task{ts}){2},'fontsize',14,'color',[.5 .5 .8]);
                else
                    % plot time points
                    axes('Position',[0.1 0.1 0.85 0.1],'next','add'); axis off               
                    line([0 .45],[.5 .5],'linewidth',10,'color',[.2 .4 .4])               
                    line([.45 .9],[.5 .5],'linewidth',10,'color',[.5 .5 .8])
                    A = annotation('rectangle',[(i/140)*.85+.1 .115 .01 .05],'FaceColor','k','FaceAlpha',1);
                    axes('Position',[0.1 0.07 0.85 0.1],'next','add'); axis off
                    text(.2,.3,ImNames.(Task{ts}){1},'fontsize',14,'color',[.2 .4 .4]);
                    text(.65,.3,ImNames.(Task{ts}){2},'fontsize',14,'color',[.5 .5 .8]);
                end
                % save the frame
                F(i)= getframe(gcf);
                writeVideo(vidfile,F(i));
                pause(.01)
                delete(H)
                delete(A)
            end
            close(vidfile);
            close all
        end
    end
end