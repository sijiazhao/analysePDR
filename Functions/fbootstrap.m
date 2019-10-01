function Y = fbootstrap(X,B)
% X = 3-D matrix subj x cond x time
% B = number of iterations


y = nan(B,size(X,2),size(X,3));


for b = 1:B
    n = size(X,1); %e.g. selecte n samples from the pool and compute the sample average
    idx = datasample(1:n,n,'Replace',true); % draw n subjects with replacement
    
    for cond = 1:size(X,2)
        xsample = squeeze(X(idx,cond,:));
        ysample = nanmean(xsample);
        
        y(b,cond,:) = ysample;
    end
end

Y = squeeze(nanmean(y,1));
end
