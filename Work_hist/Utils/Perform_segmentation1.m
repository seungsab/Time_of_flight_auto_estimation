function [IND_Final, idx, ipt, feature1, y_hat, r2adj, IND_Do] = Perform_segmentation1(feature1,ipt)

% [ipt,residual(1)] = findchangepts(feature1,'Statistic','linear','MaxNumChanges',0);
idx(1:ipt(1)) = 0; idx(ipt(1)+1:length(feature1)) = 1;

% Construct Linear Piecewise regression based on result of change point
y_hat = [];
[mdl, y_hat] = Compute_Linear_Piecewise_regression(ipt,feature1);

% Construct Linear regression for all data
mdl0 = fitlm([1:length(feature1)]',feature1);
y_hat(:,end+1) = predict(mdl0,[1:length(feature1)]');

% Compute R-squared values for stopping criterion
r2 = []; r2adj = [];
for i=1:size(y_hat,2)
    SSres=sum((feature1-y_hat(:,i)).^2); % residual sum of squares
    SStot=sum((feature1-mean(feature1)).^2); % total sum of squares
    r2(i)=1-SSres/SStot; % standard rsquared
    r2adj(i) = 1 - SSres/SStot * (length(feature1)-1)/(length(feature1)-2); % adjust for the number of parameters
end

IND_Final = 1:ipt(1); IND_Do = 1;
% residual = [residual; residue];
% Set the Stopping condition based on R-squared values
% if min(r2adj) > 0.99
%     IND_Final = 1:length(feature1); idx = ones(size(idx));
%     IND_Do = 0;
% else
%     IND_Final = 1:ipt(1);
%     IND_Do = 1;
% end
end