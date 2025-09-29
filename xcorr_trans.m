function [r,lag] = xcorr_trans(dff,vel,lagmax)

% dff: TxNxW T-transition number
%            N-cell number
%            W-time window frame number

% vel: TxW   T-transition number
%            W-time window frame number

% win: window for xcorr peaks
%      limits the lag range from -win to win

% dff = dff_lt;
% vel = vel_lt;
% win = 10;

[tn1,cn,wn1] = size(dff);
[tn2,wn2]     = size(vel);

r   = nan(tn1,cn);
lag = nan(tn1,cn);

if (tn1~=tn2)||(wn1~=wn2)
    fprintf('Size of these two vectors are not matching! \n');
    return;
end


for tt = 1:tn1
    
    vel_tt = vel(tt,:);
    
    for n = 1:cn
        
        dff_tc = dff(tt,n,:);
        
        if ~isnan(dff_tc(1))
            
            [rr,ll] = xcorr(vel_tt(:),dff_tc(:),lagmax,'coeff');
            
            idx_rr_max = find(rr == max(rr));
            
            if numel(idx_rr_max)>=1
                
                idx_rr_max = idx_rr_max(1);
                
                if ~isnan(rr(idx_rr_max))
                    
                    r(tt,n)   = rr(idx_rr_max);
                    lag(tt,n) = ll(idx_rr_max);
                    
                end
                
            end
        end
        
    end
    
end

end