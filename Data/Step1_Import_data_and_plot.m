
%% Initialize
clc, clear all
close all

color_type = [
    0 0.4470 0.7410
    0.8500 0.3250 0.0980
    0.9290 0.6940 0.1250
    0.4940 0.1840 0.5560
    ];

% N_trial // N_array
N_trial = 1; N_array = 2;

%% load time
t = load('time.txt');

%% Load Ground Truth based on visual selection
N_cycle = '80'; load('80cycle_Y_true');
% N_cycle = '260'; load('260cycle_Y_true');
Temp_ToF = Y_true;


%% Load signal
Y = {}; ToF = {};
for i=1:9
    Y_true = Temp_ToF{i};
    
    temp0 = []; temp1 = [];
    for j=1:16
        temp = load([pwd '\60\' num2str(i) '\60' num2str(i) '_' N_cycle '_' num2str(i) '.txt']);
        temp0 = [temp0 temp];
        temp1 = [temp1 Y_true(:,j)'];
    end
    
    Y{1,i} = temp0; ToF{1,i} = temp1;
end



%% Save Result
save(['Rawdata_MJY_' N_cycle],'t','Y','ToF');

%% Plot Some Samples.
% disp(['# Trial: ' num2str(N_trial) '// # Array: ' num2str(N_array)])
% ToF_GT = Y_true{1,N_trial}(N_array,:);
% figure(1); set(gcf,'position',[187.7 117 1034.7 1008])
% cla;
% 
% y = Y{1,N_trial};
% for i=1:size(y,2)
%     
%     subplot(4,1,i); cla
%     a = plot(t,y(:,i),'color',color_type(i,:),'linewidth',1);
%     axis tight; grid on;  hold on
%     
%     [B, IX] = min(abs(t-ToF_GT(i)));
%     plot(t(IX),y(IX,i),'ko','markersize',7,'linewidth',1);
%     set(gca,'fontsize',15,'fontweight','bold')
%     legend(a,['Ch #' num2str(i)],'location','best');
% end

