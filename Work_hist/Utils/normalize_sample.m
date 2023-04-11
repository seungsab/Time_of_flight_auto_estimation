function X = normalize_sample(X0)
X = -1 + 2.*(X0 - min(X0))./(max(X0) - min(X0));
X(isnan(X))=0;
end