clear; clc;
addpath(genpath('/Users/babylab/Documents/Elham/git/stimulus_assessment'))

load('Stimuli/SloanTextScrambles_devel.mat','img');
load('Stimuli/Vernier_Stim.mat');
img_Letter = squeeze(img(:,:,:,1:2:end,:,:));
img_Vernier = Im;
clear Im img;
%%

img_Vernier_input = img_Vernier(1:768,1:768,:,:);
img_Letter_input = img_Letter(1:768,313:1080,:,:,:);

exp_fov =  12;%8.5;

%[resp_vals_V, exampleV] = soc_assessment(img_Vernier_input(1:768,1:768,1,1), exp_fov,2);

for cond = 4:5
    for ts = 1:2
        [resp_vals_L{cond,ts}, exampleL] = soc_assessment(squeeze(img_Letter_input(1:400,1:400,cond,1:2,ts)), exp_fov,1);
        [resp_vals_V{cond,ts}, exampleV] = soc_assessment(squeeze(img_Vernier_input(1:400,1:400,1,2)), exp_fov,1);
    end
end

%% Letters
RA = cellfun(@(x) mean(x),resp_vals_L,'uni',false);
RA = arrayfun(@(x) cat(1,RA{:,x}),1:2,'uni',false);

figure,
subplot(1,3,1);bar(RA{1}');
subplot(1,3,2);bar(RA{2}');
subplot(1,3,3);bar(abs(RA{1}'-RA{2}'));ylim([0 .9])

%% Vernier
%RA = cellfun(@(x) mean(x),resp_vals_V,'uni',false);
RA = arrayfun(@(x) cat(1,resp_vals_V{:,x}),1:2,'uni',false);

figure,
subplot(1,3,1);bar(RA{1}');
subplot(1,3,2);bar(RA{2}');
subplot(1,3,3);bar(abs(RA{1}'-RA{2}'));ylim([0 .9])

%%
% for i = 1:size(img,3)
%     i
%     for j = 1:size(img,4)
%         fftimg(:,:,i,j,1)=fft2(img(:,:,i,j,1),100,100);
%         fftimg(:,:,i,j,2)=fft2(img(:,:,i,j,2),100,100);
%     end
% end
% 
% %%
% fftimg(1,1,:,:,:)=0;
% figure,
% for j = 1:5
%     subplot(5,2,j*2-1),imagesc(mean(abs(squeeze(fftimg(1:50,1:50,j*2,1,1))),4));
%     subplot(5,2,j*2),imagesc(mean(abs(squeeze(fftimg(1:50,1:50,j*2,1,2))),4));
% end


%% bandpass filtering

load BPFilter
% prepare the image
res = 256;
 imexample = squeeze(img(:,:,8,1,2));
% imexample = imresize(imexample,[res,res]);
figure, subplot(2,1,1); imagesc(imexample);

% prepare the filter
% temp = zeros(res,res);
% temp (1:size(bpfilter,1),1:size(bpfilter,2))= bpfilter;
% temp = fftshift(abs(fft2(temp)));
% figure; imagesc(temp);

II = conv2(imexample,bpfilter);
subplot(2,1,2); imagesc(II);
colormap('bone')