% The goal of this script is to
% (1) Create a project with axx trial files
% (2) Read in axx trial files from subjects
% (3) Apply decomposition approaches on the subjects

clear;
clc;
addpath(genpath('/Users/elhamb/Documents/Codes/Git/mrC'));
addpath(genpath('/Users/elhamb/Documents/Codes/Git/EEGSourceSim'));
addpath(genpath(fileparts(mfilename('fullpath'))))
%% Load Axx trial files

PDiva_Path = '/Users/elhamb/Documents/Data/TextScramble_exp1';
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

%% PCA on individual subjects and prepare data for group level RCA
%load(fullfile('ResultData','FFT_Trial_Data'));
Freqs = 0:outData{1,1}.dFHz:outData{1,1}.dFHz*(outData{1,1}.nFr-1);
Task = {'Letter','Vernier'};

for Sub = 1:numel(SubIDs)
    % merge letter and vernier conditions for 
    axx.(Task{1}){Sub} = MergeAxx(outData(Sub,1:5));
    axxM.(Task{1}){Sub} =MergeAxx(cellfun(@(x) x.AverageTrials(),outData(Sub,1:5),'uni',false));
    %[decompAxx_ind.(Task{1}){Sub},~,A_ind.(Task{1}){Sub},~] = mrC.SpatialFilters.RCA(axx.(Task{1}){Sub},'freq_range',Freqs([7 19]));
     
    axx.(Task{2}){Sub} = MergeAxx(outData(Sub,6:10));
    axxM.(Task{2}){Sub} =MergeAxx(cellfun(@(x) x.AverageTrials(),outData(Sub,6:10),'uni',false));
    %[decompAxx_ind.(Task{2}){Sub},~,A_ind.(Task{2}){Sub},~] = mrC.SpatialFilters.PCA(axx.(Task{2}){Sub},'freq_range',Freqs([7 19]));
end

clear Subjfolders;% outData;
%% plot individual and average ASD results

Ords = 1:10;
Finds  = [7 13 19 25];

Sublist = num2cell(1:numel(SubIDs));
Sublist{end+1} = 1:numel(SubIDs);
load ResultData/LogMar_Val.mat
logMAR = [logMAR_letter logMAR_ver];
% figure params
FS = 9;
SizeInc = .05;
if false
    for Sub =  numel(Sublist): numel(Sublist)
        for i = 1:numel(Finds)%4
            FIG = figure;
            set(FIG,'unit','inch','position',[17 10 7 2.36],'color','w')
            set(FIG,'unit','inch','paperposition',[17 10 7 2.36])
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
                    set(get(h,'label'),'String',[num2str(i) 'F ASD \muV'],'position',get(get(h,'label'),'position')+[-5.4 0 0])
                    if Cond ==5
                        set(get(h,'label'),'position',get(get(h,'label'),'position')+[-0.25 0 0])
                    end
                end

                caxis([0 LimC(ceil(Cond/5))]);
                if Cond==6
                    text(-1.8,-.2,'Vernier','fontsize',FS,'fontweight','bold','rotation',90);
                elseif Cond ==1
                    text(-1.8,-.2,'Letter','fontsize',FS,'fontweight','bold','rotation',90);
                end
             end

             %colormap('jet')
             colormap(jmaColors('hotcortex'))

             if Sub<=numel(SubIDs)
                 Subname = SubIDs{Sub};
             else
                 Subname = 'AverageAll';
             end
             %print(FIG,['Figures/TopoMap_individuals/TopoMap_' num2str(i) 'f1_' Subname '_Amp'],'-r300','-dtiff');
             export_fig(FIG,['Figures/TopoMap_individuals/TopoMap_' num2str(i) 'f1_' Subname '_Amp'],'-pdf');
   
             close all;
        end
    end
end

% clear decompAxx_Letter_all decompAxx_Vernier_all A_Vernier_all A_Letter_all D_Letter_all D_Vernier_all;
%% GROUP LEVEL ANALYSIS
%SELECT the harmonic to do analysis 
analHarms = [1 2];
Harms = {'H1F1','H2F1','H3F1','H4F1'};
RedoRCA = false;
% Group Level RCA
if RedoRCA || ~exist(fullfile('ResultData','GroupRCA.mat'),'file')
    for f = 1:numel(analHarms)
        fRCA = Finds(analHarms(f));
        for ts = 1:numel(Task)
            [decompAxx_all.(Task{ts}).(Harms{analHarms(f)}),W_all.(Task{ts}).(Harms{analHarms(f)}),A_all.(Task{ts}).(Harms{analHarms(f)}),D_all.(Task{ts}).(Harms{analHarms(f)})] = mrC.SpatialFilters.RCA(MergeAxx(axx.(Task{ts})),'freq_range',Freqs(fRCA));
        end
    end
    save(fullfile('ResultData','GroupRCA'),'decompAxx_all','W_all','A_all','D_all','-v7.3');
else
    load(fullfile('ResultData','GroupRCA.mat'));
end

%% Plot RC components calculated on all subjects

% figure params
Subnum = numel(SubIDs);% number of subjects
Condnum = 5; % number of conditions per task
NCOMP  = 2;
Cols = brewermap(4,'Set1');
Cols = Cols(3:4,:);
sublabels = {'A','D','B','E','C','F'};

logMARs.(Task{1}) = logMAR_letter;
logMARs.(Task{2}) = logMAR_ver;
FS = 9;
if true
    for f = 1:numel(analHarms)
        if f==1
            Flips = [1 -1]; 
        else
            Flips = [1 1];
        end
        caxisS = [(max(abs(A_all.(Task{1}).(Harms{analHarms(f)})(:,1:2)))); (max(abs(A_all.(Task{2}).(Harms{analHarms(f)})(:,1:2)))) ];

        FIG=figure;
        set(FIG,'unit','inch','position',[10 5 7 6],'color','w');
        set(FIG,'unit','inch','paperposition',[17 10 7 6])
        for comp = NCOMP:-1:1
            for ts = 1:numel(Task)
                S(comp+(ts-1)*2) = subplot(3,NCOMP*2,comp+(ts-1)*NCOMP);mrC.Simulate.plotOnEgi(A_all.(Task{ts}).(Harms{analHarms(f)})(:,comp)*Flips(ts));axis tight equal;caxis([-caxisS(ts,comp) caxisS(ts,comp)]);
                [~,IndElec] = sort(A_all.(Task{ts}).(Harms{analHarms(f)})(:,comp)*Flips(ts),'descend');
                Elec_max(f,ts,comp,:) = IndElec; 
                title(['RC' num2str(comp)],'fontsize',FS,'Color',Cols(comp,:));%,'fontweight','normal');%colorbar;
                set(S(comp+(ts-1)*2),'position',get(S(comp+(ts-1)*2),'position')+[-0.02+ts*.032 -.07 -.02 -.02]);
                SPP = get(S(comp+(ts-1)*2),'position');
                h=colorbar;
                set(h,'location','westoutside');
                set(S(comp+(ts-1)*2),'position',SPP)
                set(h,'position',get(h,'position')+[.02 0 0 0],'fontsize',FS);
                if comp==1
                    text(1, 2.4,Task{ts},'fontsize',FS+2,'fontweight','bold')
                    SP(ts) = S(comp+(ts-1)*2);
                end
            end

            colormap(jmaColors('coolhotcortex'))
        end

        %
        fRCA = Finds(analHarms(f));
        for ts = 1:numel(Task)
            numcomp = numel(D_all.(Task{ts}).(Harms{analHarms(f)}));
            TCos = squeeze(mean(reshape(squeeze(decompAxx_all.(Task{ts}).(Harms{analHarms(f)}).Cos(fRCA-1:fRCA+1,:,1:end-70)),[3,numcomp,16,Condnum,Subnum-1]),3));% last subject has 14 trails
            TSin = squeeze(mean(reshape(squeeze(decompAxx_all.(Task{ts}).(Harms{analHarms(f)}).Sin(fRCA-1:fRCA+1,:,1:end-70)),[3,numcomp,16,Condnum,Subnum-1]),3));

            TCos(:,:,:,Subnum) = squeeze(mean(reshape(squeeze(decompAxx_all.(Task{ts}).(Harms{analHarms(f)}).Cos(fRCA-1:fRCA+1,:,end-69:end)),[3,numcomp,14,Condnum,1]),3));
            TSin(:,:,:,Subnum) = squeeze(mean(reshape(squeeze(decompAxx_all.(Task{ts}).(Harms{analHarms(f)}).Sin(fRCA-1:fRCA+1,:,end-69:end)),[3,numcomp,14,Condnum,1]),3));
            TCmplx.(Task{ts}).(Harms{analHarms(f)}) = Flips(ts)*TCos+(TSin*1i*Flips(ts));
        end

        for comp = 1:NCOMP 
            for ts = 1:numel(Task)
                SP(ts+2) = subplot(3,2,ts+2); %plot(logMAR_ver,mean(abs(Ver_cmplx(comp,:,:)),3),'-o','Color',Cols(comp,:),'linewidth',1.5);hold on;
                errorbar(logMARs.(Task{ts}),mean(squeeze(abs(TCmplx.(Task{ts}).(Harms{analHarms(f)})(2,comp,:,:))),2),std(squeeze(abs(TCmplx.(Task{ts}).(Harms{analHarms(f)})(2,comp,:,:))),[],2)/sqrt(18),'Color',Cols(comp,:),'MarkerFaceColor',Cols(comp,:),'linewidth',1.5);hold on;box off
                xlim([0 max(logMARs.(Task{ts}))+.15])
                ylim([.0 .95])
                ylabel('Amplitude [\muV]','fontsize',FS)
                Noise = mean(squeeze(mean(abs(TCmplx.(Task{ts}).(Harms{analHarms(f)})([1 3],comp,:,:)),1)),2);
                fill([logMARs.(Task{ts}) logMARs.(Task{ts})(end:-1:1)],[Noise' zeros(size(logMARs.(Task{ts})))+.00],[.7 .7 .7],'facealpha',.7,'linestyle','none');

                SP(ts+4) = subplot(3,2,ts+4);
                TAng = wrapTo360(rad2deg(squeeze(angle(TCmplx.(Task{ts}).(Harms{analHarms(f)})(2,comp,:,:)))));
                TAng_all(f,ts,comp,:,:) = TAng;
                errorbar(logMARs.(Task{ts}),mean(TAng,2),std(squeeze(TAng),[],2)/sqrt(18),'Color',Cols(comp,:),'MarkerFaceColor',Cols(comp,:),'linewidth',1.5);hold on;box off
                xlim([0 max(logMARs.(Task{ts}))+.15])
                ylim([0 450])
                xlabel('LogMAR','fontsize',FS);
                ylabel('Phase [r]','fontsize',FS);
            end

        end

        set(SP(4),'position',get(SP(4),'position')+[.025 -0.05 -.05 0],'xtick',round(logMAR_ver,2),'xticklabel',round(logMAR_ver,2),'fontsize',FS,'linewidth',1.2,'layer','top')
        set(SP(3),'position',get(SP(3),'position')+[.025 -0.05 -.05 0],'xtick',round(logMAR_letter,2),'xticklabel',round(logMAR_letter,2),'fontsize',FS,'linewidth',1.2,'layer','top')
        set(SP(6),'position',get(SP(6),'position')+[.025 -0.04 -.05 0],'fontsize',FS,'ytick',0:90:360,'yticklabel',{'0','','\pi','','2\pi'},'xtick',round(logMAR_ver,2),'xticklabel',round(logMAR_ver,2),'linewidth',1.2,'layer','top');
        set(SP(5),'position',get(SP(5),'position')+[.025 -0.04 -.05 .0],'fontsize',FS,'ytick',0:90:360,'yticklabel',{'0','','\pi','','2\pi'},'xtick',round(logMAR_letter,2),'xticklabel',round(logMAR_letter,2),'linewidth',1.2,'layer','top'); 

        % ABC labels
        for sub = 6:-1:1
            set(gcf,'currentaxes',SP(sub),'unit','normalized')
            switch sub
                case {1,2}
                    text(-2,2.5,sublabels{sub},'fontsize',FS+5,'fontweight','bold','color','k')
                case {3}
                    text(-0.2,1,sublabels{sub},'fontsize',FS+5,'fontweight','bold','color','k')
                case {4}
                    text(-0.15,1,sublabels{sub},'fontsize',FS+5,'fontweight','bold','color','k')
                case{5}
                    text(-0.2,460,sublabels{sub},'fontsize',FS+5,'fontweight','bold','color','k')
                case{6}
                    text(-0.15,460,sublabels{sub},'fontsize',FS+5,'fontweight','bold','color','k')
            end
        end

%         print(FIG,['Figures/TopoMap_individuals/RCA_' Harms{analHarms(f)} '_AverageAll'],'-r300','-dtiff');
%         export_fig(FIG,['Figures/TopoMap_individuals/RCA_' Harms{analHarms(f)} '_AverageAll'],'-pdf');
        close all;
        clear TSin TCos TCmplx;
    end
    %clear decompAxx_all W_all A_all D_all;
end

%% Estimate temporal parameters of RC1 and RC2
Harms2 = {'1F','2F'};
FS = 11;
tsec = 2.381:2.381:140*2.381;
elecnum = [2 1;1 2];
for f = 2:2
    FIG  = figure;
    set(FIG,'unit','inch','color','w','position',[3 3 8 6])
    for ts = 1:numel(Task)
         
         subplot(2,2,ts),
         [Elloc] = mrC.plotOnEgi(zeros(128,1));
         hold on;
         Ps = mean(Elloc(Elec_max(f,ts,1,1:elecnum(f,ts)),:),1);
         ellipse(.1,.1*elecnum(f,ts),0,Ps(1),Ps(2),Cols(1,:))
         text(Ps(1)-.2,Ps(2)+.3,'RC1','Color',Cols(1,:),'fontsize',FS,'fontweight','bold')
         P1 = Ps(1:2);

         hold on;
         Ps = mean(Elloc(Elec_max(f,ts,2,1:elecnum(f,ts)),:),1);
         ellipse(.1,.1*elecnum(f,ts),pi/2,Ps(1),Ps(2),Cols(2,:))
         text(Ps(1)-.2,Ps(2)+.3,'RC2','Color',Cols(2,:),'fontsize',FS,'fontweight','bold')
         P2 = Ps(1:2);
         axis tight
         text(-0.3, 1.6,Task{ts},'fontsize',FS+2,'fontweight','bold')
         
         % Plot the arrow
         if f ==1
            arrow('start',P1,'stop',P2,'Length',8,'Width',1.5)
         else
             arrow('start',P2,'stop',P1,'Length',8,'Width',1.5)
         end
         colormap(jmaColors('coolhotcortex'));
         caxis([-.99 1])
         %
         PHtemp  = squeeze(TAng_all(f,ts,:,:,:));
         PHdiff = squeeze(diff(PHtemp,[],1))/360*1000/3/f;
         if f==2 % From RC2 to RC1
             PHdiff = PHdiff*-1;
         end
         
         S = subplot(2,2,ts+2);
         errorbar(logMARs.(Task{ts}),mean(PHdiff,2),std(squeeze(PHdiff),[],2)/sqrt(18),'Color','k','MarkerFaceColor',Cols(comp,:),'linewidth',1.5);hold on;box off
         M = max(abs(mean(PHdiff,2))+std(squeeze(PHdiff),[],2)/sqrt(18));
         set(S,'position',get(S,'position')+[0 .1 0 0],'xtick',round(logMARs.(Task{ts}),2),'xticklabel',round(logMARs.(Task{ts}),2),'fontsize',FS,'linewidth',1.2,'layer','top')
         xlim([0 max(logMARs.(Task{ts}))+.15])
         xlabel('LogMAR','fontsize',FS);
        
         
         if f==1
              ylim([20 125]); 
              if ts ==1
                ylabel('Latency RC1-RC2 (mS)')
              end
         else
             ylim([25 55]);
             if ts ==1
                ylabel('Latency RC2-RC1 (mS)')
              end
         end
    end
    
    %print(FIG,['Figures/Temporal_RCA_' Harms2{f} '_AverageAll'],'-r300','-dtiff');axis tight;
    export_fig(FIG,['Figures/Latency_RCA_' Harms2{f} '_AverageAll'],'-pdf');
    close;
end



%% Temporal dynamics of the RCs
Styles = {'-','-'};
FreqIdxall = 6*(1:16)+1;
OddFreqIdx = FreqIdxall(1:2:end);
EvenFreqIdx = FreqIdxall(2:2:end);
OEharms = [OddFreqIdx;EvenFreqIdx];


%GROUP LEVEL ANALYSIS 2
%SELECT the harmonic to do analysis 
analHarms = [1 2];
Harms = {'H1F1','H2F1','H3F1','H4F1'};
RedoRCA = false;
% Group Level RCA
if RedoRCA || ~exist(fullfile('ResultData','GroupRCA2.mat'),'file')
    for f = 1:numel(analHarms)
        for ts = 1:numel(Task)
            ts
            [decompAxx_all.(Task{ts}).(Harms{analHarms(f)}),W_all.(Task{ts}).(Harms{analHarms(f)}),A_all.(Task{ts}).(Harms{analHarms(f)}),D_all.(Task{ts}).(Harms{analHarms(f)})] = mrC.SpatialFilters.RCA(MergeAxx(axxM.(Task{ts})),'freq_range',Freqs(OEharms(f,1:8)));
        end
    end
    save(fullfile('ResultData','GroupRCA2'),'decompAxx_all','W_all','A_all','D_all','-v7.3');
else
    load(fullfile('ResultData','GroupRCA2.mat'));
end
%
Harms2 = {'1F','2F'};
Fac = 1;
FS2 = FS*Fac;


for f = 1:numel(analHarms) % for each harmonics
    if f==1
        Flips = [-1 1]; 
    else
        Flips = [1 -1];
    end
    
    for ts = 1:numel(Task) % for each task
        
         
        Wave_temp = decompAxx_all.(Task{ts}).(Harms{analHarms(f)}).Wave(:,1:NCOMP,:);
        TWave = squeeze(mean(reshape(squeeze(Wave_temp),[140,NCOMP,Condnum,Subnum]),3));% Reshape the waves
        Wave_all.(Task{ts}).(Harms{analHarms(f)}) = TWave*Flips(ts);
        clear Wave_temp TWave;
        
        % separate even and odd
        TCmplx = zeros(840,NCOMP,Condnum,Subnum);
        fRCA=OEharms(f,:);
        TCos = squeeze(reshape(squeeze(decompAxx_all.(Task{ts}).(Harms{analHarms(f)}).Cos(fRCA,1:NCOMP,:)),[numel(fRCA),NCOMP,Condnum,Subnum]));% last subject has 14 trails
        TSin = squeeze(reshape(squeeze(decompAxx_all.(Task{ts}).(Harms{analHarms(f)}).Sin(fRCA,1:NCOMP,:)),[numel(fRCA),NCOMP,Condnum,Subnum]));
       
        TCmplx(fRCA,:,:,:) = Flips(ts)*TCos+(TSin*1i*Flips(ts));
        TCmplx(end:-1:422,:,:,:) = conj(TCmplx(2:420,:,:,:));
        RWave = ifft(conj(TCmplx),840,1);
        Idxs = 0:140:840;
        RWaveT =  RWave(Idxs(1)+1:Idxs(1+1),:,:,:) ;
        
        for i = 1:numel(Idxs)-1
            RWaveT = RWaveT+RWave(Idxs(i)+1:Idxs(i+1),:,:,:);
        end
        RWave_all.(Task{ts}).(Harms{analHarms(f)}) = RWaveT/6*840/2;
     end
     
    FIG = figure;
    set(FIG,'unit','inch','position',[1 10 11 2.8]*Fac,'color','w')
    set(FIG,'unit','inch','paperposition',[1 10 11 2.8]*Fac)
     for ts = 1:numel(Task)
         MWave = mean(RWave_all.(Task{ts}).(Harms{analHarms(f)}),4);
        % Topo plots
         for cmp = 1:2
             S = subplot(2,7,((Condnum+2)*(ts-1))+cmp);
             mrC.plotOnEgi(A_all.(Task{ts}).(Harms{analHarms(f)})(:,cmp)*Flips(ts))
             T = title(['RC' num2str(cmp)],'color',Cols(cmp,:),'fontweight','bold','fontsize',FS-1);
             %T.Position = T.Position-[0 .2 0] ;
             SP  = S.Position +[.025*cmp-.1 0 0 0];
             h=colorbar;
             S.Position = SP;
             %set(h,'location','westoutside');
             set(h,'position',get(h,'position')-[0.00 0 0 0.02])
             
             if cmp==1
                 text(-2,-0.4,Task{ts},'fontsize',FS2+2,'fontweight','bold','rotation',90);
             end
         end
         
         % Temporal dynamics
         for cond = 1:5
             S = subplot(2,7,cond+(ts-1)*(Condnum+2)+2); 
             S.Position = S.Position - [.0 0 0 0];
             line([0 333],[0 0],'color',[.4 .4 .4],'linestyle','--');hold on;
             for c = 1:2
                Cm(c) = plot(2.381:2.381:140*2.381,MWave(:,c,cond),Styles{f},'Color',Cols(c,:),'linewidth',1.5);%hold on;
                set(gca,'fontsize',FS-2,'xtick',0:100:300);
             end
             
             title(['LogMAR = ' num2str(round(logMARs.(Task{ts})(cond),2))],'fontsize',FS2,'fontweight','normal');
             if cond==5 && f==1
%                  L = legend (Cm,'RC1','RC2');
%                  L.FontSize = FS2-4;
             end
             if cond==1 && ts ==2
                 ylabel('Amplitude(\muV)','fontsize',FS2,'fontweight','bold');%,'fontweight','bold');
                 xlabel('Time(mS)','fontsize',FS2,'fontweight','bold');
                 %ylabel(Task{ts},'fontsize',FS2,'fontweight','bold');%,'fontweight','bold');
             end
             xlim([0 333]);
             ylim([-1.4 1.4]);%axis tight;
             
         end
         
         %set(S,'position',get(S,'position')+[-.01 0 .02 0])
     end
     
     print(FIG,['Figures/Temporal_RCA_' Harms2{f} '_AverageAll'],'-r300','-dtiff');axis tight;
     export_fig(FIG,['Figures/Temporal_RCA_' Harms2{f} '_AverageAll'],'-pdf');
     close;
end



