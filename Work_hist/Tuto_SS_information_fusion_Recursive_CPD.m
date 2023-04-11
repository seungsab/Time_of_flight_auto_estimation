%% Initialize
clc, clearvars -except N_trial fn_data
close all;


cur_fig_save = [pwd '\Result'];
if ~isdir('Result')
    mkdir(cur_fig_save);
end


color_type = 'rbk';

%% Load and add Utils
addpath Utils
addpath Data

% load data_JY
% ch_dist = 0.005; % mm => m
% t = time;  signal = Sig1; Data_type = 'No_damage';
% % t = time;  signal = Sig2; Data_type = 'Damaged';
% Dist = ch_dist*[1:size(signal,2)];

fn_data ='Rawdata_MJY_80.mat';
N_trial = 1;

load(fn_data)

signal = Y{1,N_trial};

Data_type = [fn_data(1:end-4) '_Trial#' num2str(N_trial)];
ch_dist = 0.005; % mm => m
Dist = ch_dist*[1:size(signal,2)];


%%
ToF_target = ToF{N_trial};
IND_final = []; Run_time = [];

figure(1); set(gcf,'position',[321.7 79 984 809.3])
figure(2); set(gcf,'position',[1397.7 212.3 1020.7 390.7]);
for Ch=1:size(signal,2)
%     figure(11); clf;
    
    Ch = 2;
    gif_on = 1;
    if Ch<10
        fn_sv = [Data_type '_0' num2str(Ch)];
    else
        fn_sv = [Data_type '_' num2str(Ch)];
    end
    
    % Target TOF (Target value => selected by human)
    [B_target, IX_target] = min(abs(t-ToF_target(Ch)));
    
    
    %% Stage #1: Recursive Segmentation using Change Point Detection with # CP = 1 (for whole signal)    
    n_step = 1;
    
    tic    
    % Min-Max Normalize data between [-1, 1]
    y = signal(:,Ch);
    Y_norm = normalize_sample(y);
    t_norm = normalize_sample(t);
    
    IND_Result = 1:length(Y_norm);
    ipt_hist = []; IND_Result_hist = {};
    RMSE_resi = []; adjR2 =[]; BIC_hist =[];
    
    IND_Do = 1;
    while IND_Do
        IND_Result_hist{1,n_step} = IND_Result;
        tic
        % Set Feature vecotrs for whole signal
        feature1 = Y_norm(IND_Result); feature1 = normalize_sample(feature1);
        feature2 = cumsum((Y_norm(IND_Result).^2));  feature2 = normalize_sample(feature2);
        
        % Peform Sinal Segmentation using Change point dection with "Linear" statistics
        [ipt, residue1] = findchangepts(feature1,'Statistic','rms','MaxNumChanges',1);
        
        if ~isempty(ipt)    
            % [ipt,residual(1)] = findchangepts(feature1,'Statistic','linear','MaxNumChanges',0);
            idx(1:ipt) = 0; idx(ipt+1:length(Y_norm)) = 1;
            
            % Construct Linear regression for all data
            mdl0 = fitlm([1:ipt]',feature2(1:ipt));
            y_hat = predict(mdl0,[1:length(feature1)]');
            
        else
            break;            
        end
        
        ipt_hist(n_step) = ipt; RMSE_resi(n_step) = mdl0.RMSE;
        adjR2(n_step) = mdl0.Rsquared.Adjusted;
        BIC_hist(n_step) = mdl0.ModelCriterion.BIC;
        
        IND_Result = 1:ipt;
        
        % Save time for execution
        Run_time(n_step,Ch) = toc;
        
        % Plot Result of the Segmentation
        title_str = ['(' num2str(n_step) ') Signal Segmentation using Change Point Detection (# Ch: ' num2str(Ch) ')'];
        figure(1);
        Plot_Result_CPD
        
        n_step = n_step + 1;
    end
    
    %% Stage #1: Select Optimal Segment based on Maximal Adjust R-squared values
    [B_max, IX_max] = max(adjR2); IND_Result = 1:ipt_hist(IX_max); IX_max
    
    % If the R-squared value is not larger than 0.98, current segment
    % contains some signal component due to excitation (not noise)    
    if B_max<=0.95
        [B_min, IX_min] = min(RMSE_resi(1:IX_max)); IND_Result = 1:ipt_hist(IX_min);
    end
    
    %% Stage #3: Perform Outlier detection based on Robust Covariance estimation (for the Final Segment from Stage #1)
    IND_IND_ROBCOV = ones(size(Y_norm));
    
    % Perfomr Robust COV over the Final Segment
    X0 =[t_norm' Y_norm];
    [sig,mu,mah,IND_outlier] = robustcov(X0(IND_Result,:),'Method','olivehawkins');
    
    robustdist = pdist2(X0, mu,'mahal',sig);
    p = size(X0,2);
    chi2quantile = sqrt(chi2inv(0.975,p)); % Some TOF are chosend with one periodic delay
    IND_IND_ROBCOV = robustdist>chi2quantile;
    IND_IND_ROBCOV(1:ceil(length(IND_outlier)/2)) = 0;
    
    % Plot Result
    figure(1);
    subplot(211); cla;
    plot(t,Y_norm,'color',0.6*[1 1 1]);  hold on; axis tight; grid on
    plot(t(IND_Result(end))*[1 1],get(gca,'ylim'),'k:','linewidth',2);
    gscatter(t(1:length(IND_Result)),Y_norm(1:length(IND_Result)),IND_IND_ROBCOV(1:length(IND_Result)),color_type,'...');
    patch([0 t(length(IND_Result)) t(length(IND_Result)) 0],[-1 -1 1 1],'y','FaceAlpha',0.1)
    set(gca,'fontsize',15,'fontweight','bold');
    xlabel(''); ylabel('Normalized Signal (Ynorm)'); axis tight
    legend off;
    title(['(' num2str(n_step) ') Outlier Detection using Robust MCD Final Segmentation (# Ch: ' num2str(Ch) ')'],'fontsize',17,'fontweight','bold')
    
    subplot(212); cla;  hold on; grid on
    a=gscatter(t(1:length(IND_Result)),Y_norm(1:length(IND_Result)),IND_IND_ROBCOV(1:length(IND_Result)),color_type,'...');
    patch([0 t(length(IND_Result)) t(length(IND_Result)) 0],[-1 -1 1 1],'y','FaceAlpha',0.1)
    xlabel('Time (sec)');ylabel('Normalized Signal (Ynorm)');
    title('Outlier Detection based on Robust Covariance Estimation');
    set(gca,'fontsize',15,'fontweight','bold'); axis tight;
    legend off
    
    Run_time(n_step,Ch) = toc; n_step = n_step + 1;
    
    %% Stage #4: Estimate Time-of-Flgith (TOF) based on results from Stage #1-2    
    % Find local peaks (maximum)
    locs1 = islocalmax(Y_norm,'FlatSelection','center'); locs1 = find(locs1 == 1);
    locs2 = islocalmin(Y_norm,'FlatSelection','center'); locs2 = find(locs2 == 1);
    IND_peak = sort(unique([locs1; locs2]));
    
    % Compute Sign index for local peaks
    IND_positive = Y_norm(IND_peak) > median(Y_norm);
    
    % Find TOF using information fusion
    IND_Final = IND_IND_ROBCOV(IND_peak).*IND_positive;
    
    
    temp_INDEX_final = IND_peak(find(IND_Final==1));
    IND_final(1,Ch) = temp_INDEX_final(1);
    
    Run_time(n_step,Ch) = toc; n_step = n_step + 1;
    
    % Plot Result  
%     figure(1); clf
%     plot(t,Y_norm,'color',0.6*[1 1 1]);  hold on; grid on
%     gscatter(t(IND_peak),Y_norm(IND_peak),IND_Final,'rc');
%     patch([0 t(length(IND_Result)) t(length(IND_Result)) 0],[-1 -1 1 1],'y','FaceAlpha',0.1); axis tight
%     plot(get(gca,'xlim'),median(Y_norm)*[1 1],'k-')
%     a=plot(t(IX_target)*[1 1],get(gca,'ylim'),'b-','linewidth',2);
%     b=plot(t(IND_final(1,Ch))*[1 1],get(gca,'ylim'),'m:','linewidth',2);
%     title(['Final Result (# Ch:' num2str(Ch) ' in Trial #' num2str(N_trial) ')']);
%     xlabel('Time (sec)');ylabel('Normalized Signal (Ynorm)');
%     set(gca,'fontsize',15,'fontweight','bold');
%     legend([a, b],{'Target','Predicted'},'location','best')
    
    
    % Save image as gif-Animation
    if gif_on
        drawnow
        frame = getframe(1);
        img = frame2im(frame);
        [imind1, cm1] = rgb2ind(img,256);
        imwrite(imind1,cm1,[cur_fig_save '\' fn_sv '.gif'],'WriteMode','append','DelayTime',0.5);
    end
    
    %% Selected
    
    figure(2); set(gcf,'position',[1397.7 212.3 1020.7 390.7]); clf
    plot(t,Y_norm); hold on
    a=plot(t(IX_target)*[1 1],get(gca,'ylim'),'b-','linewidth',2);
    b=plot(t(IND_final(1,Ch))*[1 1],get(gca,'ylim'),'m:','linewidth',2);
    grid on; axis tight;  ylim([-1 1])
    plot(get(gca,'xlim'),mean(Y_norm)*[1 1],'k-')
    title(['Final Result (# Ch:' num2str(Ch) ' in Trial #' num2str(N_trial) ')']);
    set(gca,'fontsize',15,'fontweight','bold');
    legend([a, b],{'Target','Predicted'},'location','best')
    
    % Save image as gif-Animation
    if gif_on>0
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
%         exportgraphics(figure(2),[cur_fig_save '\Final_' fn_sv '.png'],'Resolution',300)
    end
end
ToF_hat{1,N_trial} = t(IND_final);
save([Data_type '.mat']);
