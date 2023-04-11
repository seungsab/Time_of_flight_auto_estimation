%% Initialize
clc, clearvars -except N_trial fn_data
close all;


cur_fig_save = [pwd '\Result'];
if ~isdir('Result')
    mkdir(cur_fig_save);    
end



%% Load and add Utils
addpath Utils
addpath Data

% load data_JY
% ch_dist = 0.005; % mm => m
% t = time;  signal = Sig1; Data_type = 'No_damage';
% % t = time;  signal = Sig2; Data_type = 'Damaged';
% Dist = ch_dist*[1:size(signal,2)];

fn_data ='Rawdata_MJY_80.mat';
N_trial = 2;

load(fn_data)


signal = Y{1,N_trial};

Data_type = [fn_data(1:end-4) '_Trial#' num2str(N_trial)];
ch_dist = 0.005; % mm => m
Dist = ch_dist*[1:size(signal,2)];


%%
ToF_target = ToF{N_trial};
IND_final = []; Run_time = []; 

figure(1); set(gcf,'position',[321.7 79 984 809.3])
figure(2); set(gcf,'position',[1397.7 212.3 1020.7 390.7])
for Ch=1:size(signal,2)
    n_step = 1;
    tic
    % Min-Max Normalize data between [-1, 1]
    y = signal(:,Ch);
    Y_norm = normalize_sample(y);
    t_norm = normalize_sample(t);
    
    [B_min, IX_min] = min(abs(t-ToF_target(Ch)));
    
    %% Find local peaks!
    % Due to noise components in the signal, local flucuations in the
    % signal generate local peaks in very small intervals
    locs1 = islocalmax(Y_norm,'FlatSelection','center'); locs1 = find(locs1 == 1);
    locs2 = islocalmin(Y_norm,'FlatSelection','center'); locs2 = find(locs2 == 1);
    
    IND_peak = sort(unique([locs1; locs2]));
    
    Run_time(n_step,Ch) = toc; n_step = n_step + 1;
    
    %% Stage #1: K-means clustering
    tic
    X0 =[t_norm' Y_norm];
    
    feature1 = cumsum((Y_norm.^2));
    feature1 = normalize_sample(feature1);
    
    [ipt,residual] = findchangepts(feature1,'Statistic','linear','MaxNumChanges',2);
    figure(33); clf; findchangepts(feature1,'Statistic','linear','MaxNumChanges',2);
    idx(1:ipt(1)) = 1; idx(ipt(1)+1:ipt(2)) = 2; idx(ipt(2)+1:length(Y_norm)) = 3;
    
    IND_Final1 = 1:ipt(1); temp_IND2 = ipt(2);
    Run_time(n_step,Ch) = toc; n_step = n_step + 1;
    
    % Perform K-means with K=3 to divide the RoI into two parts
%     feature(:,1) = t_norm;
%     feature(:,2) = cumsum((X0(:,2).^2));
%     feature = normalize_sample(feature);
%     [idx, ctrs] = kmeans(feature,3,'start',[min(feature); median(feature); max(feature)]);
%     idx1 = find(idx == 1); % Towards Smaller
%     idx2 = find(idx == 2); % Towards Medium
%     idx3 = find(idx == 3); % Towards Larger
    
%     figure(55); clf; gscatter(feature(:,1),feature(:,2),idx); hold on
    
    % Find 1st Segment from the result of K-means
%     [~, sort_max_IND] = sort([max(feature(idx1)), max(feature(idx2)), max(feature(idx3))]);
%     IND_Final1 = find(idx == sort_max_IND(1)); % Selected the first segment (output of Information Source #1)
%     temp_IND2 = find(idx == sort_max_IND(2)); % For plot 
%     Run_time(n_step,Ch) = toc; n_step = n_step + 1;
    
    % Plot results from information source #1
    figure(1);
    subplot(311); cla
    a=gscatter(t,Y_norm,idx); hold on
%     a(sort_max_IND(2)).MarkerSize = 2; a(sort_max_IND(2)).Color = 0.6*[1 1 1];
%     a(sort_max_IND(3)).MarkerSize = 2; a(sort_max_IND(3)).Color = 0.6*[1 1 1];
    plot(t(IND_Final1(end))*[1 1],get(gca,'ylim'),'k:','linewidth',2);
    plot(t(temp_IND2(end))*[1 1],get(gca,'ylim'),'k:','linewidth',2);
    grid on; axis tight;  ylim([-1 1])
    plot(get(gca,'xlim'),median(Y_norm)*[1 1],'k-','linewidth',0.5);
    set(gca,'fontsize',15,'fontweight','bold');
    title(['[Stage #1] Signal Segmentation using K-means clustering (# Ch: ' num2str(Ch) ')']);
    legend off
%     legend(a(sort_max_IND(1)),'1st Segment','location','northwest')
    
    %% Stage #2: Perform Robust COV to extract extracted time segment before excitation
    tic    
    feature1 = cumsum((Y_norm(IND_Final1).^2));
    feature1 = normalize_sample(feature1);
    
    [ipt2,residual] = findchangepts(feature1,'Statistic','linear','MaxNumChanges',2);
    figure(55); clf; findchangepts(feature1,'Statistic','linear','MaxNumChanges',2);
%     [ipt,residual] = findchangepts(feature1,'Statistic','linear','MinThreshold',.01);
%     figure(55); clf; findchangepts(feature1,'Statistic','linear','MinThreshold',.01);
    
    % Generate Label for all local peak to Stage #3 
    IND_Final2 = 1:ipt2(1);
    
    % Find maximum in the extracted time segment
%     [B_max, IX_max] = max(Y_norm(IND_Final2));
%     IND_Final2 = 1:IX_max;
    
    Run_time(n_step,Ch) = toc; n_step = n_step + 1;
    
    figure(1);
    % Plot results from information source #2
    subplot(312); cla;
    plot(t,Y_norm,'color',0.6*[1 1 1]); hold on
    plot(t(IND_Final1(end))*[1 1],get(gca,'ylim'),'k:','linewidth',2);
    plot(t(IND_Final2(end))*[1 1],get(gca,'ylim'),'r:','linewidth',2);
    grid on; axis tight;  ylim([-1 1])
    plot(get(gca,'xlim'),median(Y_norm)*[1 1],'k-','linewidth',0.5);
    set(gca,'fontsize',15,'fontweight','bold');
    title(['[Stage #2]  Segmentation using K-means clustering over 1st Segment']);
    legend off
    
    %% Stage #3: Information fusion and Path Clustering via GOP
    tic
    feature1 = cumsum((Y_norm(IND_Final2).^2));
    feature1 = normalize_sample(feature1);
    
    [ipt3,residual2] = findchangepts(feature1,'Statistic','linear','MaxNumChanges',1);
    figure(77); clf;  findchangepts(feature1,'Statistic','linear','MaxNumChanges',1);
    %     [residual1, residual2, residual1./residual2]
    %     [ipt3,residual] = findchangepts(feature1,'Statistic','linear','MinThreshold',.01);
    %     figure(77); clf; findchangepts(feature1,'Statistic','linear','MinThreshold',.01);
        
    K = length(ipt3); nseg = K+1;
    istart = [1 ipt3]; istop = [ipt3-1 length(feature1)];
    p_coeff = []; S = {}; mu = []; y_hat = []; delta_hat=[];
    for s=1:nseg
        ix = (istart(s):istop(s));
        b1 = 0.5*(s>1); b2 = 0.5*(s<nseg);        
        [p,S{1,s},mu] = polyfit(ix,feature1(ix,1),1);
        [y_hat(:,s),delta_hat(:,s)] = polyval(p,1:length(feature1),S{1,s},mu);
        p_coeff(s,:) = p;
    end
    p_coeff
    figure(99); clf; plot(y_hat(:,1),'b.');hold on; plot(y_hat(:,2),'r:','linewidth',2);
    
    % Generate Label for all local peak to Stage #3
    IND_Final3 = 1:ipt3(1);
    
    % Information Source #1: Perfomr Robust COV over the extracted time segment
    X0 =[t_norm' Y_norm]; X0 = normalize_sample(X0(IND_Final3,:));
    [sig,mu,mah,Extract_INDa] = robustcov(X0,'Method','olivehawkins');
    Extract_INDa(1:ceil(length(Extract_INDa)/2)) = 0;
    
    IND_Outliers = ones(size(Y_norm));
    IND_Outliers(1:IND_Final3(end)) = Extract_INDa;
    IND_Outliers = IND_Outliers(IND_peak);
    
    
    IND_positive = Y_norm(IND_peak) > median(Y_norm);
    
    IND_Final = IND_Outliers.*IND_positive;
    temp_INDEX_final = IND_peak(find(IND_Final==1));
    IND_final(1,Ch) = temp_INDEX_final(1);

    Run_time(n_step,Ch) = toc; n_step = n_step + 1;
    
    
    figure(1);
    subplot(313); cla;
    plot(t,Y_norm,'color',0.6*[1 1 1]); hold on
    gscatter(t(IND_Final3),Y_norm(IND_Final3),Extract_INDa);    
    plot(t(IND_peak),Y_norm(IND_peak),'ko');
    gscatter(t(IND_peak),Y_norm(IND_peak),IND_Final);   
    plot(t(IND_Final1(end))*[1 1],get(gca,'ylim'),'k:','linewidth',2);
    plot(t(IND_Final2(end))*[1 1],get(gca,'ylim'),'r:','linewidth',2);
    plot(t(IND_final(1,Ch))*[1 1],get(gca,'ylim'),'g:','linewidth',2);
    grid on; axis tight;  ylim([-1 1])
    plot(get(gca,'xlim'),median(Y_norm)*[1 1],'k-')
    set(gca,'fontsize',15,'fontweight','bold');
    title('[Stage #3] 2nd Information Fusion (Result of Stage 2 + Sign indices)');
    legend off
    
    if Ch<10
        fn_sv = [Data_type '_0' num2str(Ch)];
    else
        fn_sv = [Data_type '_' num2str(Ch)];
    end
    
%     exportgraphics(figure(1),[cur_fig_save '\' fn_sv '.png'],'Resolution',300)
    
    % Generate gif-Animation
%     drawnow
%     frame = getframe(1);
%     img = frame2im(frame);
%     [imind1, cm1] = rgb2ind(img,256);
    
%     if Ch == 1
%         imwrite(imind1,cm1,[Data_type '.gif'],'Loopcount',5);
%     else
%         imwrite(imind1,cm1,[Data_type '.gif'],'WriteMode','append');
%     end
    
    %% Selected
    figure(2); cla
    plot(t,Y_norm); hold on
    a=plot(t(IX_min)*[1 1],get(gca,'ylim'),'b-','linewidth',2);
    b=plot(t(IND_final(1,Ch))*[1 1],get(gca,'ylim'),'r:','linewidth',2);
    grid on; axis tight;  ylim([-1 1])
    plot(get(gca,'xlim'),mean(Y_norm)*[1 1],'k-')
    title(['Final Result (# Ch:' num2str(Ch) ' in Trial #' num2str(N_trial) ')']);
    set(gca,'fontsize',15,'fontweight','bold');
    legend([a, b],{'Target','Predicted'},'location','best')
    
%     exportgraphics(figure(2),[cur_fig_save '\Final_' fn_sv '.png'],'Resolution',300)
    
    drawnow
    % figure에서의 frame을 가져욤
    frame = getframe(2);
    % 가져온 frame을 image로 변화시킴
    img = frame2im(frame);
    % index화된 이미지로 변화시킴
    [imind2, cm2] = rgb2ind(img,256);
    
    if Ch == 1
        imwrite(imind2,cm2,['Result_' Data_type '.gif'],'Loopcount',5);
    else
        imwrite(imind2,cm2,['Result_' Data_type '.gif'],'WriteMode','append');
    end
end
ToF_hat{1,N_trial} = t(IND_final);

save([Data_type '.mat']);
