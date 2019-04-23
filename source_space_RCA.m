
clear; clc;

PATH = '/Users/kohler/code/git/mrC';
addpath(genpath(PATH));
path2 = '/Users/kohler/code/git/gardner/mrTools/mrUtilities/mrFlatMesh';
addpath(genpath(path2));
%% Read the inverses
Path = '/Volumes/svndl/mrC_Projects/VernierLetters/source';
InvName = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangkgsROIsCorr.inv';
[Inverses,SubIDInv] = mrC.Simulate.ReadInverses(Path,InvName);
% Read  morph maps
for s = 1:numel(SubIDInv)
    mapMtx{s} = makeDefaultCortexMorphMap(SubIDInv{s},SubIDInv{2});
    InvMapped{s} = mapMtx{s}*Inverses{s}.';
end
% mapMtx*Inverses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SubIDInv{strcmp(SubIDInv,'skeri0003')} = 'nl-0033';
SubIDInv{strcmp(SubIDInv,'skeri0004')} = 'nl-0034';

%%
Task = {'Letter','Vernier'};
Harms = {'H1F1','H2F1','H3F1','H4F1'};
Finds  = [7 13 19 25];
analHarms = [1 2];

RedoRCA = false;
if RedoRCA || ~exist(fullfile('ResultData','GroupRCA.mat'),'file')
   
else
    load(fullfile('ResultData','GroupRCA.mat'));
end

%%