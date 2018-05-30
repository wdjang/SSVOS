%% ===================================================
%  Semi-supervised video object segmentation
%  Code written by Won-Dong Jang (wdjang@mcl.korea.ac.kr)
%  ===================================================

clc
clear all
close all

%% System setting
addpath(genpath('./opticalflow'));
addpath('others');

%% Parameter setting
param_list = set_param;

%% DB setting
data_dir = './dataset';
user_dir = 'annotation';
result_dir = './results';

seq_name = 'girl';

%% Set sequence directory
fprintf('=============================================\n');
fprintf('%s\n',seq_name);
fprintf('=============================================\n');
seq_path = fullfile(data_dir,seq_name);
ann_path = fullfile(data_dir,seq_name,user_dir);
result_path = fullfile(result_dir,seq_name);

%% Perform Semi-Supervised Video Object Segmentation (SSVOS)
overall_tic = tic;
segment_track = SSVOS(seq_path, ann_path, result_path, param_list);
overall_toc = toc(overall_tic);

%% Visualize segmentation result
figure;
for frame_id = 1:length(segment_track)
    imshow(segment_track{frame_id}==1);
    pause(0.1);
end


