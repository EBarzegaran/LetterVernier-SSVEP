% The goal of this script is to
% (1) Create a project with axx trial files
% (2) Read in axx trial files from subjects
% (3) Apply decomposition approaches on the subjects - and do statistics
% using bootstrapping

clear;
clc;
addpath(genpath('/Users/elhamb/Documents/Codes/Git/mrC'));
addpath(genpath('/Users/elhamb/Documents/Codes/Git/EEGSourceSim'));
addpath(genpath(fileparts(mfilename('fullpath'))))
%% Load Axx trial files

PDiva_Path  =   '/Users/elhamb/Documents/Data/TextScramble_exp1';
%mrC_Path = '/Volumes/svndl/mrC_Projects/VernierLetters/';

% Read files from thge original folder
Subjfolders =   subfolders(PDiva_Path,0);
Subjfolders =   Subjfolders(cellfun(@(x) ~isempty(x),strfind(Subjfolders,'nl-')));
SubIDs      =   cellfun(@(x) x(1:7),Subjfolders,'UniformOutput',false);

for Sub     =   1:numel(Subjfolders)
axx_trialFiles  =   subfiles(fullfile(PDiva_Path,Subjfolders{Sub},'Exp_MATL_HCN_128_Avg','Axx*_trials.mat'),1);
    for Cond    =   1:length(axx_trialFiles)
        axxStrct            =   matfile(axx_trialFiles{Cond});
        outData{Sub,Cond}   =   mrC.axx.loadobj(axxStrct);
    end
end
clear axx_trialFiles axxStrct Cond Sub;
%save(fullfile('ResultData','FFT_Trial_Data'),'outData','SubIDs');

%% PCA on individual subjects and prepare data for group level RCA
%load(fullfile('ResultData','FFT_Trial_Data'));
Freqs   =   0:outData{1,1}.dFHz:outData{1,1}.dFHz*(outData{1,1}.nFr-1);
Task    =   {'Letter','Vernier'};

for Sub = 1:numel(SubIDs)
    % merge letter and vernier conditions for 
    axx.(Task{1}){Sub}  =   MergeAxx(outData(Sub,1:5));
    axxM.(Task{1}){Sub} =   MergeAxx(cellfun(@(x) x.AverageTrials(),outData(Sub,1:5),'uni',false));
    %[decompAxx_ind.(Task{1}){Sub},~,A_ind.(Task{1}){Sub},~] = mrC.SpatialFilters.RCA(axx.(Task{1}){Sub},'freq_range',Freqs([7 19]));
     
    axx.(Task{2}){Sub}  =   MergeAxx(outData(Sub,6:10));
    axxM.(Task{2}){Sub} =   MergeAxx(cellfun(@(x) x.AverageTrials(),outData(Sub,6:10),'uni',false));
    %[decompAxx_ind.(Task{2}){Sub},~,A_ind.(Task{2}){Sub},~] = mrC.SpatialFilters.PCA(axx.(Task{2}){Sub},'freq_range',Freqs([7 19]));
end

clear Subjfolders;% outData;


%% GROUP LEVEL ANALYSIS -  Bootstrapping
%SELECT the harmonic to do analysis 

% (1) fist make the boostrap samples
nboot = 500;
clear Boots;
%for b = 1:nboot, Boots(b,:) = randsample(18,15);end
Boots = randi(18,nboot,10);

analHarms = [1];
Harms = {'H1F1','H2F1'};
for b = 1:nboot
    for f = 1:numel(analHarms)
        tic
        fRCA = Finds(analHarms(f));
        for ts = 1:numel(Task)

            [~,~,A_temp] = mrC.SpatialFilters.RCA(MergeAxx(axx.(Task{ts})(Boots(b,:))),'freq_range',Freqs(fRCA));
            A_all.(Task{ts}).(Harms{analHarms(f)}){b} = A_temp(:,1:5);
        end
        toc
    end
end

save('ResultData\GroupRCA_bootstrap_500','A_all')
%%
A_all_c = cat(3,A_all.Vernier.H1F1{:});
A1 = squeeze(A_all_c(:,1,:));
A1 = [squeeze(A_all_c(:,1,:)) squeeze(A_all_c(:,2,:))];
A_all.Letter.H1F1 = cellfun(@(x) x(:,1:5),A_all.Letter.H1F1,'uni',false);
A_all_c2 = cat(3,A_all.Letter.H1F1{:});
A2 = [squeeze(A_all_c2(:,1,:)) squeeze(A_all_c2(:,2,:))];
A12 = [A1 A2];
A12 = A12 ./ mean(abs(A12),1); % different kind of normalizations can be applied

imagesc(corr(A12))
[R,P] = (corr(A12));
H = imagesc(abs(R));
set(H,'alphadata',(P)<.01/800);
colormap('jet');
caxis([0 1]);

% or clusterize the scalp topographies
%%

A_all_c = cat(3,A_all.Vernier.H1F1{:});
A1 = squeeze(A_all_c(:,1,:));
A1 = [squeeze(A_all_c(:,1,:)) squeeze(A_all_c(:,2,:))];
A_all.Letter.H1F1 = cellfun(@(x) x(:,1:5),A_all.Letter.H1F1,'uni',false);
A_all_c2 = cat(3,A_all.Letter.H1F1{:});
A2 = [squeeze(A_all_c2(:,1,:)) squeeze(A_all_c2(:,2,:))];
A12 = [A1 A2];
A12 = A12 ./ mean(abs(A12),1);

figure;
for i = 1:4
    ind = (i-1:i)*nboot;
    CC = corr(A12(:,ind(1)+1:ind(2)));
    A12(:,ind(1)+1:ind(2)) =  A12(:,ind(1)+1:ind(2)) .*sign(CC(1,:));
    for j  = i+1:4
        ind2 = (j-1:j)*nboot;
        for e = 1:128
            [~,P(i,j,e),~,stat] = ttest2(A12(e,ind(1)+1:ind(2)),A12(e,ind2(1)+1:ind2(2)));
            Tval(i,j,e) = stat.tstat;
        end
    end
    
    subplot(1,4,i),ESSim.Simulate.plotOnEgi(mean(A12(:,ind(1)+1:ind(2)),2))
%     A_temp(:,i) = mean(A12(:,ind(1)+1:ind(2)),2);
end


figure,
subplot(1,2,1),imagesc((squareform(pdist(A12'))));
subplot(1,2,2),imagesc(abs(corr(A12)));
caxis([0 1]);


%%
Ac12 = A12 ./ sqrt(mean((A12).^2,1));
Dist_temp = (squareform(pdist(Ac12'.^2)));
inds = arrayfun(@(i) ((i-1)*nboot+1):(i*nboot),1:4,'uni',false);
inds = cat(1,inds{:});

for i = 1:4
    for j = 1:4
        DD = Dist_temp(inds(i,:),inds(j,:));
        subplot(4,4,(i-1)*4+j);
        x = (log10(DD(:)));
        x(isinf(x))=[];
        [N,edges] = histcounts(x);
        %hist(log10(DD(:)),100);
        plot(edges(1:end),[N 0],'linewidth',1.5);

        vline(mean(x),'k--');
        xlim([-.4 1.5])
        
        DD2 = Dist_temp(inds(i,:),inds(i,:));
        DD3 = Dist_temp(inds(j,:),inds(j,:));
        x = log10(mean(DD));
        y = log10(mean(DD2));
        z = log10(mean(DD3));
        x(isinf(x))=[];
        y(isinf(y))=[];
        [h,p(i,j),~,stats] = ttest2(x,[y z]);
        TVAL(i,j) = stats.tstat;

    end
end

figure,imagesc(TVAL);

%%
figure;
Ac12 = A12;
%Ac12(Ac12>0)=0;
PD = pdist(Ac12'.^2,'spearman');
imagesc(abs(corr(Ac12.^2)))
Z = linkage(PD,'average');
c = cophenet(Z,PD);
T = cluster(Z,'maxclust',10);
subplot(1,2,1), plot(T);
 
subplot(1,2,2), [H,T2,outperm] = dendrogram(Z,'ColorThreshold',.42);

figure;
for c = 1:10
    subplot(2,5,c)
    ESSim.Simulate.plotOnEgi(mean(Ac12(:,T==c),2))
end

%%
figure,
for i = 1:70
   subplot(7,10,i)
   ESSim.Simulate.plotOnEgi(mean(Ac12(:,i+210).^2,2))
end
