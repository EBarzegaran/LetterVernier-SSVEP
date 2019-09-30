clear;clc;

addpath(genpath('/Users/elhamb/Documents/Codes/Git/EEGSourceSim'));
% POWER DIVA and monitor parameters
ImPixY = 1080;
ImPixX = 1011;
ImSizeY = 30.5; % In cm
SpatFreq = 2;% cpd: cycle per minute
AspectR = 2; % row to columns

pixCM = ImPixY/ImSizeY;
ViewDist = 150; % In cm
% how many pixel in one ARC min
cmARC = atan(deg2rad(1))*150;
pixARC = cmARC*pixCM;
%log10(atan((HCharPix)/pixCM/150/5)*60*180/pi);
HBarnum = round(ImPixY/(pixARC)*SpatFreq*2);% how many pixels for one horizental bar
HBarSize = round((pixARC)/(SpatFreq*2));

VBarnum = round(ImPixX/(pixARC));% How many pixels for one vertical bar (Vernier displacement)
VBarSize = round((pixARC));


% Vernier misalignment condition 4
MisAarcmin = [1.4 2 2.8 4 5.6];
%Im = zeros(ImPixX,ImPixY,numel(MisAarcmin),2);


for VI = 1:numel(MisAarcmin)
    Im1 = zeros(ImPixX,ImPixY);
    Im2 = zeros(ImPixX,ImPixY);
    misApix = round(MisAarcmin(VI)/60*pixARC);

    BarIdxY = [1 round(HBarSize/2:HBarSize:ImPixY) ImPixY];% Indices of Y with a half a bar shift up
    BarIdxX = [round(1:VBarSize:ImPixX) ImPixX];% Indices of Y with a half a bar shift up
    BarIdxYof = [1 round(HBarSize/2-misApix:HBarSize:ImPixY) ImPixY];
    for x = 1:numel(BarIdxX)-1
        if mod(x,2)==1
            for y = 1:2:numel(BarIdxY)-1
                Im2(BarIdxX(x):BarIdxX(x+1),BarIdxY(y):BarIdxY(y+1)) = 1;
                Im1(BarIdxX(x):BarIdxX(x+1),BarIdxY(y):BarIdxY(y+1)) = 1;
            end
        else
            for y = 1:2:numel(BarIdxY)-1
                Im2(BarIdxX(x):BarIdxX(x+1),BarIdxYof(y):BarIdxYof(y+1)) = 1;
                Im1(BarIdxX(x):BarIdxX(x+1),BarIdxY(y):BarIdxY(y+1)) = 1;
            end
        end
    end
    Im2 = Im2';
    Im1 = Im1';
    
    Im(:,:,VI,1) = Im1;
    Im(:,:,VI,2) = Im2;
end
save('Stimuli/Vernier_Stim.mat','Im');

%%
load('Stimuli/SloanTextScrambles_devel.mat');
img = squeeze(img);

VI = 4;
FS = 12;
Fig_Ver = figure;
set(Fig_Ver,'unit','inch','Position',[5 5 6 5.7],'color','w','PaperPosition',[5 5 6 5.7]);

subplot(2,2,1),imagesc(Im(1:768,1:768,VI,1)); 
colormap('gray');
set(gca,'xtick',[],'ytick',[])
xlabel('Stage 1','fontsize',FS)

T = title('Vernier displacement at 3 Hz','fontsize',FS);
set(T,'position',get(T,'position')+[500 -20 0])

subplot(2,2,2),imagesc(Im(1:768,1:768,VI,2));
set(gca,'xtick',[],'ytick',[])
xlabel('Stage 2','fontsize',FS)


Im3 = img(:,:,7,1,1);
Im4 = img(:,:,7,1,2);

S3 = subplot(2,2,3);imagesc(Im4(:,313:1080));%axis equal
set(gca,'xtick',[],'ytick',[])
xlabel('Stage 1','fontsize',FS);
%set(S3,'position',get(S3,'position')+[0 .08 0 -.1])

T = title('Letters vs. Scambled letters at 3 Hz','fontsize',FS);
set(T,'position',get(T,'position')+[500 -20 0])

S4 = subplot(2,2,4);imagesc(Im3(:,313:1080));%axis equal
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

%% Video of the stimulus

vidfile = VideoWriter(['Figures/EXPERIMENT2/Paradigm.mp4'],'MPEG-4');
vidfile.FrameRate = 60;
open(vidfile);

cnd = 7;
ccolors = {'b','r'};
thold = 1;

FIG = figure;

set(FIG,'unit','pixel','position',[10 10 1440 768])
colormap('gray')
Tinds = [100 400];

for t = 20 : 620
    if mod(t,20)==0
        hold off;
        if mod(t,60)~=0
            imagesc(img(:,:,cnd,randi(20),2));
        else
            if rand<.2
                imagesc(img(:,end:-1:1,cnd,randi(20),1));
            else
                imagesc(img(:,:,cnd,randi(20),1));
            end
        end
    end
    
    hold on;    
    if t>(thold+15) % hold on to that color
        holdc=0;
        ind = 1;
    end
    
    if (t>10) && (ismember(t,Tinds+20)) % change the color to red
        ind = 2;
        holdc =1;
        thold = t;
    end
    
    L = 30;
    l(1) = line([1440 1440]/2,[768-L 768+L]/2,'linewidth',2.5,'color',ccolors{ind});
    l(2) = line([1440-L 1440+L]/2,[768 768]/2,'linewidth',2.5,'color',ccolors{ind});
    axis off
    
    
    F(t-19)= getframe(FIG);
    writeVideo(vidfile,F(t-19));
    delete(l);
    
end
 

close(vidfile);
close all;