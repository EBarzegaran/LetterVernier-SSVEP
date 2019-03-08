% This script is for preliminary analysis and visualization of vernier-text
% scramble task in sensor space
% Elham Barzegaran, 5.22.2018

clear;clc;

PATH = '/Users/kohler/code/git/mrC';
addpath(genpath(PATH));
addpath('Tcirc');
%% Convert spectrum data into source space
Path = '/Volumes/svndl/mrC_Projects/VernierLetters/source';
CondLabels = {'Text-1','Text-2','Text-3','Text-4','Text-5','Vernier-1','Vernier-2','Vernier-3','Vernier-4','Vernier-5'};

[outData,FreqFeatures,subIDs] = mrC.EEGRead(Path,'domain','frequency');

%% plot the results
plottype = 'amp';%'phase';

for sub = 1:numel(subIDs)
    figure,
    for cond = 1:10
       
        % phase data
        Data1 = angle(outData{cond,sub}(:,(FreqFeatures.i1F1)+1));
        Data1 = wrapTo2Pi(Data1);
        phasdata(:,cond,sub) = Data1;
        
        % amplitude data
        Data2 = abs(outData{cond,sub}(:,(FreqFeatures.i1F1)+1));
        ampdata(:,cond,sub) = Data2;
        
        % plot the results
        h = subplot(2,5,cond);
        switch plottype
            case 'phase'
                mrC.plotOnEgi((Data1));
                colormap('jet');
                caxis([0 2*pi])
            case 'amp'
                mrC.plotOnEgi(Data2);
                colormap(jmaColors('coolhotcortex'));
                %caxis([0 12]);
        end
        set(h,'position',get(h,'position')+[-.08 -.08 .03 .03]);
        axis tight
        title(CondLabels{cond},'fontweight','bold','fontsize',14);
    end
        
end

%% group level stats
for cond = 1:10
    %outDataM{cond,1} = mean(cat(3,outData{cond,:})./repmat(max(abs(cat(3,outData{cond,:}))),[20484 1 1]),3)*10;% normalized over each condition
    %outDataM{cond,1} = mean(cat(3,outData{cond,:})./repmat(max(abs(cat(3,outData{10,:}))),[20484 1 1]),3)*10;% normalized over all condition
    datacond = cat(3,outData{cond,:});

%     for Har = 1:2
%         for s = 1:size(datacond,1)
%             sourcedata = squeeze(datacond(s,(FreqFeatures.i1F1)+1,:));
%             %[ampDiff,phaseDiff,zSNR,errorEllipse] =fitErrorEllipse([real(sourcedata) imag(sourcedata)],[],1);
%             [results] = tSquaredFourierCoefs([real(sourcedata) imag(sourcedata)],[],0.01);
%             p1(s,Har)=results.pVal;
%             tsqrd1(s,Har)= results.tSqrd;
%         end
%     end
%     p{cond,1} = p1;%<(0.05/128);
%     tsqrd{cond,1}= tsqrd1;
    outDataM{cond,1} = mean(datacond,3);% not-normalized
end

% Har = 1;
% for cond = 1:10
%     h = subplot(2,5,cond);
%     mrC.plotOnEgi(double(squeeze(tsqrd{cond,1}(:,Har))));
%     PH = get(h,'position');
%     set(h,'position',PH+[-.08 -.08 .03 .03]);
%     axis tight
%     title(CondLabels{cond},'fontweight','bold','fontsize',14)
%     %caxis([0 max(ampdata(:))/1.2])
% end
% colormap(jmaColors('hotcortex'));

%% group level average

data  = cat(3,outDataM{:});

Har = 1;

figure,
ampdata = squeeze(abs(data(:,(FreqFeatures.i1F1*Har)+1,:)));

for cond = 1:10
    h = subplot(2,5,cond);
    mrC.plotOnEgi(ampdata(:,cond),[],0,[],0);
    PH = get(h,'position');
    set(h,'position',PH+[-.08 -.08 .03 .03]);
    axis tight
    title(CondLabels{cond},'fontweight','bold','fontsize',14)
    caxis([0 max(ampdata(:))/1.2])
end
colormap(jmaColors('hotcortex'));

figure,
phasdata = squeeze(wrapTo2Pi(angle(data(:,(FreqFeatures.i1F1*Har)+1,:))));
for cond = 1:10
    h = subplot(2,5,cond);
    mrC.plotOnEgi(phasdata(:,cond),[],0,[],0);
    PH = get(h,'position');
    set(h,'position',PH+[-.08 -.08 .03 .03]);
    axis tight
    title(CondLabels{cond},'fontweight','bold','fontsize',14)
    caxis([0 2*pi])
end
colormap(jmaColors('phasecolor'));


