clear;clc;
Imsize = 1000;
Im1 = zeros(Imsize,Imsize);
Im2 = zeros(Imsize,Imsize);
%load('../Stimuli/SloanTextScrambles_1920x1024_8to80.mat');
load('../Stimuli/SloanTextScrambles_devel.mat');
img = squeeze(img);
%% Bar information
Barnum = 8;
Barsize = Imsize/(Barnum*2);
BarIdx1 = [round(1:Barsize:Imsize) Imsize]; % Index of Y
if mod(numel(BarIdx1),2)~=0
    BarIdx1 = [BarIdx1 Imsize];
end

BarIdx2 = [round(1:Barsize*2:Imsize) Imsize];%Index of X
if mod(numel(BarIdx2),2)~=0
    BarIdx2 = [BarIdx2 Imsize];
end

Voffset = 12;

for i = 1:2:numel(BarIdx1)
    Im1(BarIdx1(i):BarIdx1(i+1),:)=1;
    for j = 1:2:numel(BarIdx2)
        Im2(BarIdx1(i):BarIdx1(i+1),BarIdx2(j):BarIdx2(j+1))=1;
    end
    for j = 2:2:numel(BarIdx2)-1
        Im2((BarIdx1(i)+Voffset):min(BarIdx1(i+1)+Voffset,Imsize),BarIdx2(j):BarIdx2(j+1))=1;
    end
end

%%

FS = 12;
Fig_Ver = figure;
set(Fig_Ver,'unit','inch','Position',[5 5 6 5.7],'color','w','PaperPosition',[5 5 6 5.7]);

subplot(2,2,1),imagesc(Im1); 
colormap('gray');
set(gca,'xtick',[],'ytick',[])
xlabel('Stage 1','fontsize',FS)

T = title('Vernier displacement at 3 Hz','fontsize',FS);
set(T,'position',get(T,'position')+[650 -20 0])

subplot(2,2,2),imagesc(Im2);
set(gca,'xtick',[],'ytick',[])
xlabel('Stage 2','fontsize',FS)


Im3 = img(:,:,7,1,1);
Im4 = img(:,:,7,1,2);

S3 = subplot(2,2,3);imagesc(Im3(:,300:1100));
set(gca,'xtick',[],'ytick',[])
xlabel('Stage 1','fontsize',FS);
%set(S3,'position',get(S3,'position')+[0 .08 0 -.1])

T = title('Letters vs. Scambled letters at 3 Hz','fontsize',FS);
set(T,'position',get(T,'position')+[600 -20 0])

S4 = subplot(2,2,4);imagesc(Im4(:,300:1100));
set(gca,'xtick',[],'ytick',[])
xlabel('Stage 2','fontsize',FS)
%set(S4,'position',get(S4,'position')+[0 .08 0 -.1])
%
%axes('Position',[0.47 .7 0.1 0.1])
annotation('doublearrow',[.48 .55],[.75 .75]);
annotation('doublearrow',[.48 .55],[.28 .28]);
% 
 export_fig(Fig_Ver,'Stim_visualize','-pdf');
 print(Fig_Ver,'Stim_visualize','-r300','-dtiff');

 
