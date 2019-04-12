clear; clc;
load('../Stimuli/SloanTextScrambles_1920x1024_8to80.mat');

%%
img = squeeze(img);

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