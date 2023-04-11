
%% Initialze
clc, clear all, close all

color_type = 'rbk';

%% Load and add Utils
addpath Utils
addpath Data


fn_data ='Rawdata_MJY_80.mat';
N_trial = 1;
load(fn_data)

signal = Y{1,N_trial};
Data_type = [fn_data(1:end-4) '_Trial#' num2str(N_trial)];

%% Narray
signal = Y{N_trial};

figure; hold on; grid on

for i=1:16
    plot3(t,i*ones(size(t)),signal(:,i),'linewidth',2);
end

view([-3.7 42.7]); axis tight
axis off

%% Narray
% Ch = 2;
Ch = 8;
% Ch = 16;
% Ch = 42;

% for Ch=1:64
    figure(22); set(gcf,'position',[383 724 920 387])
    clf;
    hold on; grid on
    plot(t,signal(:,Ch),'linewidth',2);
    
    axis tight;
    hold on
    title(num2str(Ch))
    xlabel('Time (Sec)');
    ylabel('Signal');
    set(gca,'fontsize',15,'fontweight','bold')
%     pause
% end