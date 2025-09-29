function [rval,lag] = xcorrrnd(a,b,nreps)

rval = nan(nreps,1);
lag  = nan(nreps,1);
win  = floor(length(a)/2);

for ii = 1:nreps
    
    [rval_temp,lag_temp] = xcorr(a(:),shuffle(b(:)),'none');
    
    idx_rmax = find(rval_temp==max(rval_temp));
    idx_rmax = idx_rmax((abs(idx_rmax-win-1) == min(abs(idx_rmax-win-1))));
    rval(ii) = rval_temp(idx_rmax);
    lag(ii) = lag_temp(idx_rmax);
    
end

end