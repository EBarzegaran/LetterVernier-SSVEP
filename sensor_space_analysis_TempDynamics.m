% The goal of this script is to
% (1) Create a project with axx trial files
% (2) Read in axx trial files from subjects
% (3) Visualize and analyize temporal dynamics

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
%% time course
% reconstruct the fourier
clear Ver_cmplx Ver_recon_even Ver_recon_odd;
FreqIdxall = 6*(1:16)+1;
OddFreqIdx = FreqIdxall(1:2:end);
EvenFreqIdx = FreqIdxall(2:2:end);
Maxx_Letter = MergeAxx(axx_Letter);
Letter_cmplx_odd = zeros(840,128,5,16);
Letter_cmplx_even = zeros(840,128,5,16);
elecnum= 128;
Condnum = 5;
Subnum = numel(SubIDs);
Vinv = 1;
Linv = 1;

for FI = 1: numel(OddFreqIdx)
    % odd harmanics
    clear Letter_cos Letter_sin;
    Letter_cos = squeeze(mean(reshape(squeeze(Maxx_Letter.Cos(OddFreqIdx(FI),:,1:end-70)),[elecnum,16,Condnum,Subnum-1]),2));% last subject has 14 trails
    Letter_sin = squeeze(mean(reshape(squeeze(Maxx_Letter.Sin(OddFreqIdx(FI),:,1:end-70)),[elecnum,16,Condnum,Subnum-1]),2));

    Letter_cos(:,:,Subnum) = squeeze(mean(reshape(squeeze(Maxx_Letter.Cos(OddFreqIdx(FI),:,end-69:end)),[elecnum,14,Condnum,1]),2));
    Letter_sin(:,:,Subnum) = squeeze(mean(reshape(squeeze(Maxx_Letter.Sin(OddFreqIdx(FI),:,end-69:end)),[elecnum,14,Condnum,1]),2));
    Letter_cmplx_odd(OddFreqIdx(FI),:,:,:) = Vinv*Letter_cos+(Letter_sin*1i*Vinv);
    
    % even harmaonics
    clear Letter_cos Letter_sin;
    Letter_cos = squeeze(mean(reshape(squeeze(Maxx_Letter.Cos(EvenFreqIdx(FI),:,1:end-70)),[elecnum,16,Condnum,Subnum-1]),2));% last subject has 14 trails
    Letter_sin = squeeze(mean(reshape(squeeze(Maxx_Letter.Sin(EvenFreqIdx(FI),:,1:end-70)),[elecnum,16,Condnum,Subnum-1]),2));

    Letter_cos(:,:,Subnum) = squeeze(mean(reshape(squeeze(Maxx_Letter.Cos(EvenFreqIdx(FI),:,end-69:end)),[elecnum,14,Condnum,1]),2));
    Letter_sin(:,:,Subnum) = squeeze(mean(reshape(squeeze(Maxx_Letter.Sin(EvenFreqIdx(FI),:,end-69:end)),[elecnum,14,Condnum,1]),2));
    Letter_cmplx_even(EvenFreqIdx(FI),:,:,:) = Vinv*Letter_cos+(Letter_sin*1i*Vinv);
    
end
Letter_cmplx_odd(end:-1:422,:,:,:) = conj(Letter_cmplx_odd(2:420,:,:,:));
Letter_recon_odd = ifft(conj(Letter_cmplx_odd),840,1);

Letter_cmplx_even(end:-1:422,:,:,:) = conj(Letter_cmplx_even(2:420,:,:,:));
Letter_recon_even = ifft(conj(Letter_cmplx_even),840,1);
for i = 1:(840/140)
    Letter_recon_oddM(:,:,:,:,i) = Letter_recon_odd((i-1)*140+1:i*140,:,:,:);
    Letter_recon_evenM(:,:,:,:,i) = Letter_recon_even((i-1)*140+1:i*140,:,:,:);
end
Letter_recon_oddM = mean(Letter_recon_oddM,5);
Letter_recon_evenM = mean(Letter_recon_evenM,5);
% all wave
MWave = squeeze(mean(reshape(squeeze(Maxx_Letter.Wave(:,:,1:end-70)),[140,elecnum,16,Condnum,Subnum-1]),3));% last subject has 14 trails
MWave(:,:,:,Subnum) = squeeze(mean(reshape(squeeze(Maxx_Letter.Wave(:,:,end-69:end)),[140,elecnum,14,Condnum,1]),3));

clear Letter_cmplx_even Letter_cmplx_odd Letter_cos Letter_sin Letter_recon_even Letter_recon_odd
%%
Condnum = 5;
vidfile = VideoWriter(['Figures/Letter_Cond' num2str(Condnum) '.mp4'],'MPEG-4');
vidfile.FrameRate = 10;
open(vidfile);
FIG = figure;
set(FIG,'unit','inch','position',[5 5 15 5])
M1 = max(max(abs(squeeze(mean(Letter_recon_evenM(:,:,Condnum,1:16),4)))));
M2 = max(max(abs(squeeze(mean(Letter_recon_oddM(:,:,Condnum,1:16),4)))));
M3 = max(max(abs(squeeze(mean(MWave(:,:,Condnum,:),4)))));
for i = 1:140
    subplot(1,3,1);
    mrC.plotOnEgi(squeeze(mean(Letter_recon_evenM(i,:,Condnum,1:16),4))*-1); colormap('jet'); caxis([-M1 M1]);
    title('Even harmonics','fontsize',12)
    subplot(1,3,2);
     mrC.plotOnEgi(squeeze(mean(Letter_recon_oddM(i,:,Condnum,1:16),4))*-1); colormap('jet'); caxis([-M2 M2]);
    title('Odd harmonics','fontsize',12)
    h = text(-.5,-1.5,['Time = ' num2str(round(2.381*i,1)) ' ms'],'fontsize',12);
    subplot(1,3,3);
    mrC.plotOnEgi(squeeze(mean(MWave(i,:,Condnum,:),4))*-1); colormap('jet'); caxis([-M3 M3]);
    title('all harmonics','fontsize',12)
    F(i)= getframe(gcf);
    writeVideo(vidfile,F(i));
    pause(.25)
    delete(h)
end

close(vidfile);
close all
%%

%% functions
function out_axx = MergeAxx(axxlist)
out_axx = axxlist{1};
    for C = 2:numel(axxlist)
        out_axx = out_axx.MergeTrials(axxlist{C});
    end
end