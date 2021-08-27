close all;
clear;
clc;
addpath(genpath('./'));

dataDir = 'data';
outDir = 'output';

lambR = 1;
lambA = 1;

if ~exist(outDir, 'dir') 
    mkdir(outDir);
end

files = dir(fullfile(dataDir, '*.mat'));
files = {files.name}';
for i = 1 : numel(files)
    fname = fullfile(dataDir, files{i});
    dataStr = strsplit(files{i}, '.');
%     dataStr = strsplit(dataStr{1}, '_'); 
    load(fname);
    data{1} = struct2array(exp);
    data{2} = struct2array(methy);
    data{3} = struct2array(mirna);
%     fprintf('for test only');
    [idx_eg, idx_rc] = MALS_wrapper(data, lambR, lambA);
    sampleNames = fieldnames(exp);
    outFile = sprintf('%s/%s.mat', outDir, dataStr{1});
    save(outFile, 'sampleNames', 'idx_eg', 'idx_rc');
end



