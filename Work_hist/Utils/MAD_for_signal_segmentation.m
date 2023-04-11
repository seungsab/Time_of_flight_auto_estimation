function [ind_MAD, threshold] = MAD_for_signal_segmentation(Cxx0, alpha)

% Cxx1 = (cumsum(Y(:,2).^2));

% 3) Run Median Absolute Deviation (MAD) for Baseline removal (Cxx1)
C = median(Cxx0); % Median
c=-1/(sqrt(2)*erfcinv(3/2)); % Scaling factor for outlier detection using MAD
threshold(1) = C+alpha*c*median(abs(Cxx0-C)); % Upper Threshold for outlier classification
threshold(2) = C-alpha*c*median(abs(Cxx0-C)); % Lower Threshold for outlier classification

ind_MAD1 = find(Cxx0>threshold(1)); ind_MAD2 = find(Cxx0<threshold(2));
ind_MAD = union(ind_MAD1,ind_MAD2);


end