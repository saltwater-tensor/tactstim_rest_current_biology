function [cvr_bootstrap, cvr_bootstrap_cell] = covar_bootstrap(data,num_of_bootstraps)

for bstr = 1:num_of_bootstraps
    
    Y = randsample(size(data,1),size(data,1),true);
    cvr_bootstrap(:,:,bstr) = cov(data(Y,:));
    cvr_bootstrap_cell{bstr} = cvr_bootstrap(:,:,bstr);
end
end