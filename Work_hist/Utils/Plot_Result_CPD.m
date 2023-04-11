
% Plot results from information source #1
% figure(1);
subplot(211); cla
plot(t,Y_norm); hold on; grid on; 
gscatter(t(1:length(idx)),Y_norm(1:length(idx)),idx,color_type,'...'); axis tight;
plot(t(IND_Result(end))*[1 1],get(gca,'ylim'),'k:','linewidth',2);
plot(get(gca,'xlim'),median(Y_norm)*[1 1],'k-','linewidth',0.5);
patch([0 t(ipt) t(ipt) 0],[-1 -1 1 1],'y','FaceAlpha',0.1)
set(gca,'fontsize',15,'fontweight','bold');
xlabel(''); ylabel('Normalized Signal (Ynorm)');
legend off;
title(title_str,'fontsize',17,'fontweight','bold')

subplot(212); cla
plot(t(1:ipt),feature2(1:ipt),'r.'); axis tight; hold on; grid on;
XLIM = get(gca,'Xlim'); YLIM = get(gca,'ylim');
a=plot(t(1:length(y_hat)),y_hat,'k--','linewidth',1);
axis([XLIM YLIM])
plot(t(ipt)*[1 1],get(gca,'ylim'),'k:','linewidth',2);
patch([0 t(ipt) t(ipt) 0],[min(feature2(1:ipt)) min(feature2(1:ipt)) max(feature2(1:ipt)) max(feature2(1:ipt))],'y','FaceAlpha',0.1)
xlabel('Time (sec)'); ylabel('Cumulative Sum of Squares');
title('Segmentation Result using Change Point Detection');
set(gca,'fontsize',15,'fontweight','bold');
legend(a,{['R^2= ' num2str(adjR2(n_step)) ' // RMSE=' num2str(RMSE_resi(n_step))]},'location','best');
% legend([a,b],{['R^2(Piecewise)= ' num2str(R_squared(1))]; ['R^2(All) = ' num2str(R_squared(2))]},'location','best');
% legend([a,b,c],{['R^2(Left)= ' num2str(R_squared(1))]; ['R^2(Right) = ' num2str(R_squared(2))]; ['R^2(Left + Right) = ' num2str(R_squared(3))]},'location','best');


% Save image as gif-Animation
if gif_on == 1
    drawnow
    frame = getframe(1);
    img = frame2im(frame);
    [imind1, cm1] = rgb2ind(img,256);
    imwrite(imind1,cm1,[cur_fig_save '\' fn_sv '.gif'],'Loopcount',5);
    gif_on = gif_on + 1;
    
elseif gif_on == 2
    drawnow
    frame = getframe(1);
    img = frame2im(frame);
    [imind1, cm1] = rgb2ind(img,256);
    imwrite(imind1,cm1,[cur_fig_save '\' fn_sv '.gif'],'WriteMode','append','DelayTime',0.5);
    
end

%     exportgraphics(figure(1),[cur_fig_save '\' fn_sv '.png'],'Resolution',300)
