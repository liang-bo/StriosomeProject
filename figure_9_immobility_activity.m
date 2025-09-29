%% figure_9_immobility_activity


%% 1.
bins = 0:1:10;
binn = length(bins);
mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];
s = 2;
dff_vbin = nan(mnum,binn-1);
dff_vbin_lt = nan(mnum,binn-1);
dff_vbin_dk = nan(mnum,binn-1);
state_name = {'lt';'drk';};
stnum = numel(state_name);

for m = 1:mnum
    
    dff   = ltdkdata.mouse(m).file(s).dFF;
    dffbin   = ltdkdata.mouse(m).file(s).dFF_bin;
    dff = dff.*dffbin;
    fn = size(dff,2);
    
    if fn<f_lim(s)
        f_use = 1:fn;
    else
        f_use = 1:f_lim(s);
    end
    
    idx_use = ~isnan(dff(:,1));
    cn = sum(idx_use);
    
    dff = dff(idx_use,f_use);
    dff = dffnorm(dff,'max');
    dff_avg = nanmean(dff,1);
    
    vel = ltdkdata.mouse(m).file(s).behavior.velocity(f_use);
    state = ltdkdata.mouse(m).file(s).behavior.state(f_use);
    svar  = zeros(numel(state),numel(state_name));
    
    for ii = 1:stnum
        svar(:,ii) = cell2mat(cellfun(@(x) strcmp(x,state_name{ii}),state,'UniformOutput',false));
    end
    svar = svar>0;
        
    A = 1;
    for b = 2:binn
        
        B = b;
        bin_idx = (vel>bins(A))&(vel<=bins(B));
        bin_idx_lt = bin_idx&svar(:,1);
        bin_idx_dk = bin_idx&svar(:,2);
        
        dff_vbin(m,b-1) = nanmean(dff_avg(bin_idx),2);
        dff_vbin_lt(m,b-1) = nanmean(dff_avg(bin_idx_lt),2);
        dff_vbin_dk(m,b-1) = nanmean(dff_avg(bin_idx_dk),2);
        A = B;
        
    end
    
end


