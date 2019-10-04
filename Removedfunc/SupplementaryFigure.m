% Supplementary Figure

clear; clc;
PATH = '/Users/elhamb/Documents/Codes/Git/mrC';
addpath(genpath('/Users/elhamb/Documents/Codes/Git/EEGSourceSim'));
addpath(genpath(PATH));
load(fullfile('ResultData','GroupRCA.mat'));

Lateral = {'V','D','VD'};
l = 3;

%% set parameters and load RCA results
FS = 12;
FigPath = fullfile('Figures','SourceSpace');

FIG = figure;
set(FIG,'Units','Inch','PaperPosition',[1 1 25 13],'Position',[1 1 25 13],'Color', 'w');

% subplot A: RCA results
PDiva_Path  =   '/Users/elhamb/Documents/Data/TextScramble_exp1';
%mrC_Path = '/Volumes/svndl/mrC_Projects/VernierLetters/';

% Read files from thge original folder
Subjfolders =   subfolders(PDiva_Path,0);
Subjfolders =   Subjfolders(cellfun(@(x) ~isempty(x),strfind(Subjfolders,'nl-')));
SubIDs      =   cellfun(@(x) x(1:7),Subjfolders,'UniformOutput',false);
%
load ResultData/LogMar_Val.mat
Task            =   {'Letter','Vernier'};
logMARs.(Task{1}) = logMAR_letter;
logMARs.(Task{2}) = logMAR_ver;
Subnum          =   numel(SubIDs);% number of subjects
Condnum         =   5; % number of conditions per task
NCOMP           =   2;
Cols            =   brewermap(4,'Set1');
Cols            =   Cols(3:4,:);
analHarms       =   [1 2];
Harms           =   {'H1F1','H2F1','H3F1','H4F1'};
Finds           =   [7 13 19 25];
f               =   1; % harmonic to be analyzed
ts              =   1; % Task to be analyzed
Flips           =   [1 -1]; 

caxisS          =   [(max(abs(A_all.(Task{1}).(Harms{analHarms(f)})(:,1:2)))) ];


for comp = NCOMP:-1:1
    S(comp) = subplot(4,6,comp);mrC.Simulate.plotOnEgi(A_all.(Task{ts}).(Harms{analHarms(f)})(:,comp)*Flips(ts));axis tight equal;caxis([-caxisS(ts,comp) caxisS(ts,comp)]);
    axis tight
    [~,IndElec] = sort(A_all.(Task{ts}).(Harms{analHarms(f)})(:,comp)*Flips(ts),'descend');
    Elec_max(f,ts,comp,:) = IndElec; 
    title(['RC' num2str(comp)],'fontsize',FS,'Color',Cols(comp,:));%,'fontweight','normal');%colorbar;
    set(S(comp),'position',get(S(comp),'position')+[0.0+((comp-1)*-.01) 0 0 0]);
    %SPP = get(S(comp+(ts-1)*2),'position');
    %h=colorbar;
    %set(h,'location','westoutside');
    %set(S(comp+(ts-1)*2),'position',SPP)
    %set(h,'position',get(h,'position')+[.02 0 0 0],'fontsize',FS);
    colormap(jmaColors('coolhotcortex'))
    freezeColors
end

% estimate the phases and amplitudes
fRCA = Finds(analHarms(f));
numcomp = numel(D_all.(Task{ts}).(Harms{analHarms(f)}));
TCos = squeeze(mean(reshape(squeeze(decompAxx_all.(Task{ts}).(Harms{analHarms(f)}).Cos(fRCA-1:fRCA+1,:,1:end-70)),[3,numcomp,16,Condnum,Subnum-1]),3));% last subject has 14 trails
TSin = squeeze(mean(reshape(squeeze(decompAxx_all.(Task{ts}).(Harms{analHarms(f)}).Sin(fRCA-1:fRCA+1,:,1:end-70)),[3,numcomp,16,Condnum,Subnum-1]),3));

TCos(:,:,:,Subnum) = squeeze(mean(reshape(squeeze(decompAxx_all.(Task{ts}).(Harms{analHarms(f)}).Cos(fRCA-1:fRCA+1,:,end-69:end)),[3,numcomp,14,Condnum,1]),3));
TSin(:,:,:,Subnum) = squeeze(mean(reshape(squeeze(decompAxx_all.(Task{ts}).(Harms{analHarms(f)}).Sin(fRCA-1:fRCA+1,:,end-69:end)),[3,numcomp,14,Condnum,1]),3));
TCmplx.(Task{ts}).(Harms{analHarms(f)}) = Flips(ts)*TCos+(TSin*1i*Flips(ts));


for comp = 1:NCOMP 

    SP(comp,1) = subplot(8,6,comp+12); %plot(logMAR_ver,mean(abs(Ver_cmplx(comp,:,:)),3),'-o','Color',Cols(comp,:),'linewidth',1.5);hold on;
    errorbar(logMARs.(Task{ts}),mean(squeeze(abs(TCmplx.(Task{ts}).(Harms{analHarms(f)})(2,comp,:,:))),2),std(squeeze(abs(TCmplx.(Task{ts}).(Harms{analHarms(f)})(2,comp,:,:))),[],2)/sqrt(18),'Color','k','MarkerFaceColor','k','linewidth',1.5);hold on;box off
    xlim([0 max(logMARs.(Task{ts}))+.15])
    ylim([.0 .95])
    if comp==1
        ylabel('Amplitude','fontsize',FS-6);
    end
    Noise = mean(squeeze(mean(abs(TCmplx.(Task{ts}).(Harms{analHarms(f)})([1 3],comp,:,:)),1)),2);
    fill([logMARs.(Task{ts}) logMARs.(Task{ts})(end:-1:1)],[Noise' zeros(size(logMARs.(Task{ts})))+.00],[.7 .7 .7],'facealpha',.7,'linestyle','none');

    SP(comp,2) = subplot(8,6,comp+18);
    TAng = wrapTo360(rad2deg(squeeze(angle(TCmplx.(Task{ts}).(Harms{analHarms(f)})(2,comp,:,:)))));
    TAng_all(f,ts,comp,:,:) = TAng;
    errorbar(logMARs.(Task{ts}),mean(TAng,2),std(squeeze(TAng),[],2)/sqrt(18),'Color',[.7 0 0],'MarkerFaceColor',[.7 0 0],'linewidth',1.5);hold on;box off
    xlim([0 max(logMARs.(Task{ts}))+.15])
    ylim([0 450])
    xlabel('LogMAR','fontsize',FS-4);
    if comp ==1
        ylabel('Phase','fontsize',FS-6);
    end
    xtickangle(90)

end
for comp =1:2
set(SP(comp,1),'position',get(SP(comp,1),'position')+[0.01+((comp-1)*-.01) .01 -0.02 0],'xtick',round(logMAR_letter,2),'xticklabel',[],'fontsize',FS,'linewidth',1.2,'layer','top')
set(SP(comp,2),'position',get(SP(comp,2),'position')+[0.01+((comp-1)*-.01) .035 -0.02 0],'fontsize',FS,'ytick',0:90:360,'yticklabel',{'0','','\pi','','2\pi'},'xtick',round(logMAR_letter,2),'xticklabel',round(logMAR_letter,2),'linewidth',1.2,'layer','top');

end

axes('Position',[0.2 0.97 0.1 0.1],'Box','off'); axis off
text (0,0,'Reliable Components','fontsize',FS+4,'fontweight','bold')

%% subplot B: forward projections
ResultPath = 'ResultData';
FilePath = fullfile(ResultPath,['LocalizationExampleData_Paper' Lateral{l} '.mat']);

if exist(FilePath)
    load(FilePath);
    load(FilePath,'ROILabel');
else
    error('Run CompareInverses.m first')
end

ScalpData = cat(3,ScalpData{:});% the first 10 is ROI and the second 10 is subject

% VWFA
VWFA = mean(mean(ScalpData(9:10,:,:),1),3);
S(1) = subplot(4,6,9);
ESSim.Simulate.plotOnEgi(VWFA);
title('VWFA','fontsize',FS);
axis tight
freezeColors

% IOG
IOG = mean(mean(ScalpData(7:8,:,:),1),3);
S(2) = subplot(4,6,3);
ESSim.Simulate.plotOnEgi(IOG);
title('IOG','fontsize',FS);
axis tight
freezeColors

% V1v
V1 = mean(mean(ScalpData(1:2,:,:),1),3);
S(3)=subplot(4,6,4);
ESSim.Simulate.plotOnEgi(V1);
title('V1d','fontsize',FS);
axis tight
freezeColors

% V3v
V2 = mean(mean(ScalpData(5:6,:,:),1),3);
S(4) = subplot(4,6,10);
ESSim.Simulate.plotOnEgi(V2);
title('V3v','fontsize',FS);
axis tight
freezeColors

set(S(1),'position',get(S(1),'position')+[.01 0 0 0])
set(S(2),'position',get(S(2),'position')+[.01 0 0 0])
set(S(3),'position',get(S(3),'position')+[-.01 0 0 0])
set(S(4),'position',get(S(4),'position')+[-.01 0 0 0])

axes('Position',[0.472 0.97 0.1 0.1],'Box','off'); axis off
text (0,0,'Forward Projections','fontsize',FS+4,'fontweight','bold')
% subplot C: Cross talks

S = subplot(2,3,3);
set(S,'position',get(S,'position')+[.01 0 0 0]);
% HemiName = {'left','right'};

%% Plot Cross Talk Matrices
CrossTalk1 = cellfun(@(x) x./repmat(max(x),[size(x,1) 1]),CrossTalk1,'uni',false);% normalize
CT1 = (cat(3,CrossTalk1{:}));

for i = 1
    eval(['CT_plot = ((mean(CT' num2str(i) ',3)));']);
    imagesc((CT_plot));

    colormap(jmaColors('coolhot'));
    %colormap('gray')
    caxis([-1 1]);

    set(gca,'ytick',1:numel(ROILabel),'yticklabel',ROILabel,'xtick',1:numel(ROILabel),'xticklabel',ROILabel,'fontsize',FS);
    if exist('xtickangle'), xtickangle(90); end
    xlabel('Receiving ROI','fontsize',FS+2);
    ylabel('Seed ROI','fontsize',FS+2);
    colorbar
end
freezeColors
%print(fullfile(FigPath,'simulation',['CrossTalk_All']),'-dtiff','-r300');

axes('Position',[0.75 0.97 0.1 0.1],'Box','off'); axis off
text (0,0,'Cross-Talk Matrix','fontsize',FS+4,'fontweight','bold')

%% subplot D: Source spectrum data

load(fullfile('ResultData',['ROISourceResults_freq_' Lateral{l}]));
Harms           =   {'h1F1','h2F1','h3F1','h4F1'};
for roi = 1:numel(ROILabelSel)

    SP(roi,1) = subplot(4,10,roi+20); %plot(logMAR_ver,mean(abs(Ver_cmplx(comp,:,:)),3),'-o','Color',Cols(comp,:),'linewidth',1.5);hold on;
    errorbar(logMARs.(Task{ts}),abs(Results.(Task{ts}).(Harms{analHarms(f)}).M(roi,:)),Results.(Task{ts}).(Harms{analHarms(f)}).S(roi,:),'Color','k','MarkerFaceColor','k','linewidth',1.5);hold on;box off
    xlim([0 max(logMARs.(Task{ts}))+.15])
    ylim([-1 8])
    if roi==1
        ylabel('Amplitude','fontsize',FS);
    end
    title(ROILabelSel{roi},'fontsize',FS)
    
    SP(roi,2) = subplot(4,10,roi+30); %plot(logMAR_ver,mean(abs(Ver_cmplx(comp,:,:)),3),'-o','Color',Cols(comp,:),'linewidth',1.5);hold on;
    errorbar(logMARs.(Task{ts}),wrapTo360(rad2deg(angle(Results.(Task{ts}).(Harms{analHarms(f)}).M(roi,:)))),Results.(Task{ts}).(Harms{analHarms(f)}).S(roi,:),'Color',[.7 0 0],'MarkerFaceColor','k','linewidth',1.5);hold on;box off
    xlim([0 max(logMARs.(Task{ts}))+.15])
    xtickangle(90);
    ylim([0 360])
    if roi==1
        ylabel('Phase','fontsize',FS);
        xlabel('logMAR','fontsize',FS)
    end
    
    
    
%     Noise = mean(squeeze(mean(abs(TCmplx.(Task{ts}).(Harms{analHarms(f)})([1 3],comp,:,:)),1)),2);
%     fill([logMARs.(Task{ts}) logMARs.(Task{ts})(end:-1:1)],[Noise' zeros(size(logMARs.(Task{ts})))+.00],[.7 .7 .7],'facealpha',.7,'linestyle','none');
% 
%     SP(comp,2) = subplot(8,6,comp+18);
%     TAng = wrapTo360(rad2deg(squeeze(angle(TCmplx.(Task{ts}).(Harms{analHarms(f)})(2,comp,:,:)))));
%     TAng_all(f,ts,comp,:,:) = TAng;
%     errorbar(logMARs.(Task{ts}),mean(TAng,2),std(squeeze(TAng),[],2)/sqrt(18),'Color',[.7 0 0],'MarkerFaceColor',[.7 0 0],'linewidth',1.5);hold on;box off
%     xlim([0 max(logMARs.(Task{ts}))+.15])
%     ylim([0 450])
%     xlabel('LogMAR','fontsize',FS-4);
%     if comp ==1
%         ylabel('Phase','fontsize',FS-6);
%     end
%     xtickangle(90)

end


for comp =1:roi
    if comp<=1
        set(SP(comp,1),'position',get(SP(comp,1),'position')+[0.01+((comp-1)*-.005) -.08 0.01 -.01],'xtick',round(logMAR_letter,2),'xticklabel',[],'fontsize',FS,'linewidth',1.2,'layer','top')
        set(SP(comp,2),'position',get(SP(comp,2),'position')+[0.01+((comp-1)*-.005) -.02 0.01 -.01],'fontsize',FS,'ytick',0:90:360,'yticklabel',{'0','','\pi','','2\pi'},'xtick',round(logMAR_letter,2),'xticklabel',round(logMAR_letter,2),'linewidth',1.2,'layer','top');
    else
        set(SP(comp,1),'position',get(SP(comp,1),'position')+[0.01+((comp-1)*-.005) -.08 0.01 -.01],'xtick',round(logMAR_letter,2),'xticklabel',[],'yticklabel',[],'fontsize',FS,'linewidth',1.2,'layer','top')
        set(SP(comp,2),'position',get(SP(comp,2),'position')+[0.01+((comp-1)*-.005) -.02 0.01 -.01],'fontsize',FS,'ytick',0:90:360,'yticklabel',[],'xtick',round(logMAR_letter,2),'xticklabel',[],'linewidth',1.2,'layer','top');

    end
end
% 
axes('Position',[0.475 0.47 0.1 0.1],'Box','off'); axis off
text (0,0,'Source Estimations','fontsize',FS+4,'fontweight','bold')

export_fig(FIG,fullfile('Figures','SourceSpace','Supplementary_fig_v1'),'-pdf')
