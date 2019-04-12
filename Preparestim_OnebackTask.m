%% estimate the logMAR for viewing distance of 150 cm instead of 200 cm
load('../Stimuli/SloanTextScrambles_devel.mat');

%%
ObjectVerticalSizePix = hCharPix;%[11 hCharPix(1:9)];

IMG = img;
img = squeeze(IMG(:,:,:,:,:,1));
img = permute(img,[1 2 4 3]);
img =  reshape(img,size(img,1),size(img,2),1,size(img,3),size(img,4));
save('../Stimuli/Sloan_intact_768x1440.mat','img','ObjectVerticalSizePix');
img = squeeze(IMG(:,:,:,:,:,2));
img = permute(img,[1 2 4 3]);
img =  reshape(img,size(img,1),size(img,2),1,size(img,3),size(img,4));
save('../Stimuli/Sloan_scrambled_768x1440.mat','img','ObjectVerticalSizePix');

%%
HCharPix =  hCharPix(1:10);
MAR = 10.^logMAR;
MAR2 = atan((4/3)*tan(MAR/(60*180)*pi))*60*180/pi;
logMAR2 = log10(MAR2); 
% Or
screenHeightCm = 31; % Monitor Upstairs
screenWidthCm = 55;

pixCM = screenHeightPix/screenHeightCm;
logMAR_letter = log10(atan((HCharPix)/pixCM/150/5)*60*180/pi);
logMAR_letter = logMAR_letter(1:2:9);


%%
MAR_ver = [1.4 2 2.8 4 5.6];
logMAR_ver = log10(MAR_ver);
save(fullfile('ResultData','LogMar_Val'),'logMAR_letter','logMAR_ver')

