function [ovlp,ovlp_t]= overlap_rate_pw(idx_mat)

[m,n] = size(idx_mat); % m, total neuron number; n, number of index

ovlp   = nan(n,n);
ovlp_t = nan(n,n);

for ii = 1:(n-1)
    
    idx_mat_ii = idx_mat(:,ii);
    
    for jj = (ii+1):n
        
        idx_mat_jj = idx_mat(:,jj);
        
        ovlp(ii,jj) = 2*sum(idx_mat_ii.*idx_mat_jj)/(sum(idx_mat_ii)+sum(idx_mat_jj)+eps);          % actual value
        
        ovlp_t(ii,jj) = 2*sum(idx_mat_ii)*sum(idx_mat_jj)/m/(sum(idx_mat_ii)+sum(idx_mat_jj)+eps);  % chance level
        
    end
    
    
end




end