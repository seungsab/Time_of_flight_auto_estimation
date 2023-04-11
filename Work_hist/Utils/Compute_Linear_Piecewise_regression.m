function [mdl, y_hat] = Compute_Linear_Piecewise_regression(ipt,feature1)

K = length(ipt); nseg = K+1;
istart = [1; ipt];
istop = [ipt-1; length(feature1)];

y_hat = [];
for s=1:nseg
    ix = (istart(s):istop(s));
    mdl{1,s} = fitlm(ix,feature1(ix,1));
    y_hat(ix,1)= predict(mdl{1,s},ix(:));
end

% p_coeff = []; SE=[]; y_hat = [];
% for s=1:nseg
%     ix = (istart(s):istop(s));
%     mdl{1,s} = fitlm(ix,feature1(ix,1));
%     p_coeff(s,1) = mdl{1,s}.Coefficients{1,1}; % Intercept
%     p_coeff(s,2) = mdl{1,s}.Coefficients{2,1}; % Slope
%     
%     SE(s,:) = diag(sqrt(mdl{1,s}.CoefficientCovariance));
%     y_hat(:,s)= predict(mdl{1,s},[1:length(feature1)]');
% end
end