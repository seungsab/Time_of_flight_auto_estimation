
%%
clc, clear all, close all


addpath Utils
addpath Data

%%

load(fn_data)

signal = Y{1,N_trial};
Data_type = [fn_data(1:end-4) '_Trial#' num2str(N_trial)];
ch_dist = 0.005; % mm => m
Dist = ch_dist*[1:size(signal,2)];