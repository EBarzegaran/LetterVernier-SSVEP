clear; clc;
addpath(genpath('/Users/elhamb/Documents/Codes/NonGit/stimulus_assessment'))

load('Stimuli/SloanTextScrambles_devel.mat','img');
load('Stimuli/Vernier_Stim.mat');
img_Letter = squeeze(img(:,:,:,1:2:end,:,:));
img_Vernier = Im;
clear Im img;
%%

img_Vernier_input = img_Vernier(1:768,1:768,:,:);
img_Letter_input = img_Letter(1:768,313:1080,:,:,:);

exp_fov = 10% atan((768*30.5/1080)/150)*180/pi;%8.5;

%[resp_vals_V, exampleV] = soc_assessment(img_Vernier_input(1:768,1:768,1,1), exp_fov,2);

for cond = 1:5
    for ts = 1:2
        [resp_vals_L{cond,ts}, Resp_img_L{cond,ts}] = soc_assessment(squeeze(img_Letter_input(1:400,1:400,cond,:,ts)), exp_fov,1);
        [resp_vals_V{cond,ts}, Resp_img_V{cond,ts}] = soc_assessment(squeeze(img_Vernier_input(1:400,1:400,cond,ts)), exp_fov,1);
    end
end

%% Letters


for cond = 1:5
    for i = 1:size(Resp_img_L{cond,1}.resp,3)
        R1 = squeeze((Resp_img_L{cond,1}.resp(:,:,i,:)));
        for j = 1:size(Resp_img_L{cond,2}.resp,3)
            R2 = squeeze((Resp_img_L{cond,2}.resp(:,:,j,:)));
            Rdif = R1-R2;
            Rdif = Rdif(5:end-5,5:end-5,:);
            RespDiff(cond,i,j,:) = squeeze(mean(mean(abs(Rdif),1),2));
            RespDiff2(cond,i,j,:) = mean(mean(abs(R1)))-mean(mean(abs(R2)));
        end
    end
    for v = 1:4
        RD = squeeze(RespDiff(cond,:,:,v));
        ResDiffL(cond,v) = mean(RD(:));
        RD2 = squeeze(RespDiff(cond,:,:,v));
        ResDiffL2(cond,v) = mean(RD2(:));
    end
end


RA = cellfun(@(x) mean(x),resp_vals_L,'uni',false);
RA = arrayfun(@(x) cat(1,RA{:,x}),1:2,'uni',false);

figure,
subplot(1,3,1);bar(RA{1}');
subplot(1,3,2);bar(RA{2}');
subplot(1,3,3);bar((RA{1}'-RA{2}'));ylim([0 .9])

clear RespDiff RD;
%% Vernier
RA = cellfun(@(x) mean(x),resp_vals_V,'uni',false);
RA = arrayfun(@(x) cat(1,resp_vals_V{:,x}),1:2,'uni',false);

figure,
subplot(1,3,1);bar(RA{1}');
subplot(1,3,2);bar(RA{2}');
subplot(1,3,3);bar(abs(RA{1}'-RA{2}'));ylim([0 .9])

for cond = 1:5
    for i = 1:size(Resp_img_V{cond,1}.resp,3)
        R1 = squeeze((Resp_img_V{cond,1}.resp(:,:,i,:)));
        for j = 1:size(Resp_img_V{cond,2}.resp,3)
            R2 = squeeze((Resp_img_V{cond,2}.resp(:,:,j,:)));
            Rdif = R1-R2;
            Rdif = Rdif(5:end-5,5:end-5,:);
            RespDiff(cond,i,j,:) = squeeze(mean(mean((Rdif),1),2));
        end
    end
    for v = 1:4
        RD = squeeze(RespDiff(cond,:,:,v));
        ResDiffV(cond,v) = mean(RD(:));
    end
end

%%
figure,

subplot(1,2,1),imagesc(ResDiffL)
subplot(1,2,2),imagesc(ResDiffV)
