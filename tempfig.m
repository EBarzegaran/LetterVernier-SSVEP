
%% subplot B: forward projections
clear; clc
Lateral = {'V','D','VD'};
l = 3;

ResultPath = 'ResultData';
FilePath = fullfile(ResultPath,['LocalizationExampleData_Paper' Lateral{l} '.mat']);

if exist(FilePath)
    load(FilePath);
    load(FilePath,'ROILabel');
else
    error('Run CompareInverses.m first')
end

ScalpData = cat(3,ScalpData{:});% the first 10 is ROI and the second 10 is subject
%%
FIG = figure;
for i =1:5
    subplot(2,5,i),
    ESSim.Simulate.plotOnEgi(mean(ScalpData(2*(i-1)+1,:,:),3));
    TT = strsplit(ROILabel{2*(i-1)+1},'_');
    title(TT{1});
    axis tight;
    if i == 1
        text(-2,.0,'Left','rotation',90,'fontsize',12)
    end
    
    subplot(2,5,i+5),
    ESSim.Simulate.plotOnEgi(mean(ScalpData(2*i,:,:),3));
    axis tight;
    if i == 1
        text(-2,.0,'Right','rotation',90,'fontsize',12);
    end
end

% subplot C: Cross talks
set(FIG,'unit','inch','position',[5 10 15 5],'color','w');
export_fig(FIG,fullfile('Figures','SourceSpace',['Topographies_' Lateral{l}]),'-pdf');

%%

