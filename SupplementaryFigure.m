% Supplementary Figure

clear; clc;
PATH = '/Users/elhamb/Documents/Codes/Git/mrC';
addpath(genpath('/Users/elhamb/Documents/Codes/Git/EEGSourceSim'));
addpath(genpath(PATH));
load(fullfile('ResultData','GroupRCA.mat'));

Lateral = {'V','D','VD'};
l = 2;

%% set parameters and load RCA results
FS = 14;
FigPath = fullfile('Figures','SourceSpace');

FIG = figure;
set(FIG,'Units','Inch','PaperPosition',[1 1 15 8],'Position',[1 1 15 8],'Color', 'w');

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

%
for comp = NCOMP:-1:1
    S(comp) = subplot(4,6,(comp-1)*6+1);mrC.Simulate.plotOnEgi(A_all.(Task{ts}).(Harms{analHarms(f)})(:,comp)*Flips(ts));axis tight equal;caxis([-caxisS(ts,comp) caxisS(ts,comp)]);
    axis tight
    [~,IndElec] = sort(A_all.(Task{ts}).(Harms{analHarms(f)})(:,comp)*Flips(ts),'descend');
    Elec_max(f,ts,comp,:) = IndElec; 
    title(['RC' num2str(comp)],'fontsize',FS,'Color',Cols(comp,:));%,'fontweight','normal');%colorbar;
    set(S(comp),'position',get(S(comp),'position')+[0.0+(-.04) (comp-1)*.017 0.01 0.01]);
    
    colormap(jmaColors('coolhotcortex'))
    freezeColors
end

axes('Position',[0.02 0.97 0.1 0.1],'Box','off'); axis off
text (0,0,'A','fontsize',FS+6,'fontweight','bold');

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

for roi = 1:5
    SP(roi,1) = subplot(4,6,roi+1);
    ESSim.Simulate.plotOnEgi(mean(ScalpData((roi-1)*2+1,:,:),3));
    title(ROILabel{roi*2}(1:end-2),'fontsize',FS);
    axis tight;
    freezeColors;
    
    SP(roi,2) = subplot(4,6,roi+7);
    ESSim.Simulate.plotOnEgi(mean(ScalpData(roi*2,:,:),3));
    axis tight;
    freezeColors;
end

for roi = 1:5
    set(SP(roi,1),'position',get(SP(roi,1),'position')+[.027 0 0.01 0.01]);
    set(SP(roi,2),'position',get(SP(roi,2),'position')+[.027 .017 0.01 0.01]);
end

axes('Position',[0.25 0.97 0.1 0.1],'Box','off'); axis off
text (0,0,'B','fontsize',FS+4,'fontweight','bold')

axes('Position',[0.27 0.83 0.1 0.1],'Box','off'); axis off
text (0,0,'Left','fontsize',FS,'fontweight','bold','rotation',90)

axes('Position',[0.27 0.63 0.1 0.1],'Box','off'); axis off
text (0,0,'Right','fontsize',FS,'fontweight','bold','rotation',90)

%% Plot Cross Talk Matrices
SP2 = subplot(4,6,13);

CrossTalk1 = cellfun(@(x) x./repmat(max(x),[size(x,1) 1]),CrossTalk1,'uni',false);% normalize
CT1 = (cat(3,CrossTalk1{:}));

CT_plot = (mean(CT1,3));
imagesc((CT_plot));

colormap(jmaColors('coolhot'));
%colormap('jet')
caxis([-1 1]);

set(gca,'ytick',1:numel(ROILabel),'yticklabel',ROILabel,'xtick',1:numel(ROILabel),'xticklabel',ROILabel,'fontsize',FS);
if exist('xtickangle'), xtickangle(90); end
xlabel('Receiving ROI','fontsize',FS+2);
ylabel('Seed ROI','fontsize',FS+2);
CB = colorbar;
CB.Ticks = [-1 1];
freezeColors

axes('Position',[0.02 0.47 0.1 0.1],'Box','off'); axis off
text (0,0,'C','fontsize',FS+4,'fontweight','bold')

set(SP2,'position',get(SP2,'position')+[-.06 -.15 .02 .08]);

%% subplot D: Source spectrum data

load(fullfile('ResultData',['ROISourceResults_freq_' Lateral{l}]));
Harms           =   {'h1F1','h2F1','h3F1','h4F1'};
Colors_amp = winter(5);
for roi = 1:numel(ROILabelSel)/2

    SP(roi,1) = subplot(4,6,roi+13); 
    B = bar(logMARs.(Task{ts}),abs(Results.(Task{ts}).(Harms{analHarms(f)}).M((roi-1)*2+1,:)),'FaceColor','flat');
    hold on;
    B.CData  = Colors_amp;
    
    errorbar(logMARs.(Task{ts}),abs(Results.(Task{ts}).(Harms{analHarms(f)}).M((roi-1)*2+1,:)),Results.(Task{ts}).(Harms{analHarms(f)}).S((roi-1)*2+1,:),'.','Color','k','MarkerFaceColor','k','linewidth',1.5);hold on;box off
    xlim([0 max(logMARs.(Task{ts}))+.15])
    ylim([-1 7])
    title(ROILabelSel{(roi-1)*2+1}(1:end-2),'fontsize',FS)
   
    SP(roi,2) = subplot(4,6,roi+19);
    B = bar(logMARs.(Task{ts}),abs(Results.(Task{ts}).(Harms{analHarms(f)}).M((roi)*2,:)),'FaceColor','flat');
    hold on;
    B.CData  = Colors_amp;
    
    errorbar(logMARs.(Task{ts}),abs(Results.(Task{ts}).(Harms{analHarms(f)}).M((roi)*2,:)),Results.(Task{ts}).(Harms{analHarms(f)}).S((roi)*2,:),'.','Color','k','MarkerFaceColor','k','linewidth',1.5);hold on;box off
    xlim([0 max(logMARs.(Task{ts}))+.15])
    ylim([-1 7])
    if roi==1
        ylabel('Amplitude','fontsize',FS);
        xlabel('LogMAR');
    end
    xtickangle(90);

end


for comp =1:roi
    if comp<=1
        set(SP(comp,1),'position',get(SP(comp,1),'position')+[.027 -.04  0.01 -.01],'xtick',round(logMAR_letter,2),'xticklabel',[],'fontsize',FS,'linewidth',1.2,'layer','top')
        set(SP(comp,2),'position',get(SP(comp,2),'position')+[0.027 0.03 0.01 -.02],'xtick',round(logMAR_letter,2),'fontsize',FS,'linewidth',1.2,'layer','top')
    else
        set(SP(comp,1),'position',get(SP(comp,1),'position')+[0.027 -.04  0.01 -.01],'xtick',round(logMAR_letter,2),'xticklabel',[],'yticklabel',[],'fontsize',FS,'linewidth',1.2,'layer','top')
        set(SP(comp,2),'position',get(SP(comp,2),'position')+[0.027 0.03 0.01 -.01],'xtick',round(logMAR_letter,2),'yticklabel',[],'fontsize',FS,'linewidth',1.2,'layer','top')
    end
end
% 
axes('Position',[0.25 0.47 0.1 0.1],'Box','off'); axis off
text (0,0,'D','fontsize',FS+4,'fontweight','bold')

axes('Position',[0.25 0.32 0.1 0.1],'Box','off'); axis off
text (0,0,'Left','fontsize',FS,'fontweight','bold','rotation',90)

axes('Position',[0.25 0.17 0.1 0.1],'Box','off'); axis off
text (0,0,'Right','fontsize',FS,'fontweight','bold','rotation',90)

export_fig(FIG,fullfile('Figures','SourceSpace',['Supplementary_fig_' Lateral{l}]),'-pdf')
