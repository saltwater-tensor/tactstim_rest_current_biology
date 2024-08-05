function [p_vals, p_vals_corr] = permtest(distance_vector,Nperms)
p = length(distance_vector);
for np = 1:Nperms+1
    if np == 1
        random_vec(:,np) = distance_vector;
    else
        random_vec(:,np) = distance_vector(randperm(p,p));
    end
end

% test for each p

for j = 1:p
    
    grtr = random_vec(j,2:end) <= random_vec(j,1);
    lssr = random_vec(j,2:end) <= random_vec(j,1);
    p_vals(j,1) = sum(grtr)/Nperms;
    p_vals(j,2) = sum(lssr)/Nperms;
    clear grtr lssr
end
    p_vals_corr(:,1) = mafdr(p_vals(:,1));
    p_vals_corr(:,2) = mafdr(p_vals(:,2));
end