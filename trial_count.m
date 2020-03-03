clear;
clc;
addpath(genpath('/Users/elhamb/Documents/Codes/Git/mrC'));
addpath(genpath('/Users/elhamb/Documents/Codes/Git/EEGSourceSim'));
addpath(genpath(fileparts(mfilename('fullpath'))))
%% Load Axx trial files

PDiva_Path  =   '/Users/elhamb/Documents/Data/TextScramble_exp1';
%mrC_Path = '/Volumes/svndl/mrC_Projects/VernierLetters/';

% Read files from thge original folder
Subjfolders =   subfolders(PDiva_Path,0);
Subjfolders =   Subjfolders(cellfun(@(x) ~isempty(x),strfind(Subjfolders,'nl-')));
SubIDs      =   cellfun(@(x) x(1:7),Subjfolders,'UniformOutput',false);

for Sub     =   1:numel(Subjfolders)
rawFiles  =   subfiles(fullfile(PDiva_Path,Subjfolders{Sub},'Exp_MATL_HCN_128_Avg','Raw*.mat'),1);
    for F    =   1:length(rawFiles)
        Data = load(rawFiles{F},'IsEpochOK');
        S(Sub,F,:) = sum(Data.IsEpochOK==0,2);
        
    end
end
%%
S2 = S(:,:,2:11);
sum(S2(:)>20)/sum(S2(:)>=0)

