% figure_8_speed_vs_lightdarkcell


%% number of cells

mnum = numel(ltdkdata.mouse);

for s = 2%:2
    
    cell_perct = zeros(mnum,5);
    
    for m = 1:mnum
        
        dff   = ltdkdata.mouse(m).file(s).dFF(:,10);
        
        cell_num = sum(~isnan(dff));
        spd_idx = ltdkdata.mouse(m).file(s).spd_idx;
        
        lt_dk_idx = ltdkdata.mouse(m).file(s).lt_dk_idx;
        
        spd_idx = spd_idx>0;
        lt_idx = lt_dk_idx==1;
        dk_idx = lt_dk_idx==2;
        
        cell_perct(m,1) = sum(lt_idx&(~spd_idx))/cell_num*100;
        cell_perct(m,2) = sum(lt_idx&spd_idx)/cell_num*100;
        cell_perct(m,3) = sum((~(lt_idx|dk_idx))&spd_idx)/cell_num*100;
        cell_perct(m,4) = sum(dk_idx&spd_idx)/cell_num*100;
        cell_perct(m,5) = sum(dk_idx&(~spd_idx))/cell_num*100;
        fprintf('mouse %d, session %d is done. \n', m, s);
        
    end
    
end

cell_perct_avg = mean(cell_perct,1);

%% speed cell subtypes

mnum = numel(ltdkdata.mouse);

for s = 2%:2
    
    cell_perct_lt_spd   = zeros(mnum,4);
    cell_perct_spd_only = zeros(mnum,4);
    cell_perct_dk_spd   = zeros(mnum,4);
    
    for m = 1:mnum
        
        dff   = ltdkdata.mouse(m).file(s).dFF(:,10);
        
        cell_num = sum(~isnan(dff));
        spd_idx = ltdkdata.mouse(m).file(s).spd_idx;
        
        lt_dk_idx = ltdkdata.mouse(m).file(s).lt_dk_idx;
        
        spd_idx_0 = spd_idx>0;
        lt_idx = lt_dk_idx==1;
        dk_idx = lt_dk_idx==2;
        
        idx_lt_spd = lt_idx&spd_idx_0;
        idx_spd_only = (~(lt_idx|dk_idx))&spd_idx_0;
        idx_dk_spd = dk_idx&spd_idx_0;
        
        
        for ii = 1:4
            
            cell_perct_lt_spd(m,ii)   = sum((spd_idx.*idx_lt_spd) == ii)/(sum(idx_lt_spd)+eps)*100;
            cell_perct_spd_only(m,ii) = sum((spd_idx.*idx_spd_only) == ii)/(sum(idx_spd_only)+eps)*100;
            cell_perct_dk_spd(m,ii) = sum((spd_idx.*idx_dk_spd) == ii)/(sum(idx_dk_spd)+eps)*100;
            
        end
        
        fprintf('mouse %d, session %d is done. \n', m, s);
        
    end
    
end
m_use = abs(sum(cell_perct_lt_spd,2)-100)<1;
lt_spd_avg = mean(cell_perct_lt_spd(m_use,:),1)';

m_use = abs(sum(cell_perct_spd_only,2)-100)<1;
spd_only_avg = mean(cell_perct_spd_only(m_use,:),1)';

m_use = abs(sum(cell_perct_dk_spd,2)-100)<1;
dk_spd_avg = mean(cell_perct_dk_spd(m_use,:),1)';

%%

trans_name = {'elt_lt';'edrk_drk';};
f_lim = [9000;18000];
vel_thr = 1; % cm/s - baseline noise is high here, so i choose a higher threshold
mkprop = {'-ro';'-bo';'-kv';'-kv'};
sname = {'Open field'; 'Light dark box'};
win = 20; % +/- 2s

mnum = numel(ltdkdata.mouse);

dff_lt_avg_pre_all = [];
dff_lt_avg_pst_all = [];
dff_dk_avg_pre_all = [];
dff_dk_avg_pst_all = [];

for s = 2%:2
    
    cell_perct_lt_spd   = zeros(mnum,4);
    cell_perct_spd_only = zeros(mnum,4);
    cell_perct_dk_spd   = zeros(mnum,4);
    
    cell_num_all = zeros(mnum,4);
    
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
        
        idx_use = ~isnan(dff(:,10));
        
        dff = dff(idx_use,f_use);
        dff = dffnorm(dff,'max');
        dff_avg = nanmean(dff,1);
        
        trans = ltdkdata.mouse(m).file(s).behavior.ttype(f_use);
        
        trans_lt = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{1}),trans,'UniformOutput',false)));
        trans_dk = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{2}),trans,'UniformOutput',false)));
        
        dff_lt = extract_from_onset_window(dff,trans_lt,win);
        dff_dk = extract_from_onset_window(dff,trans_dk,win);
        
        dff_lt_avg = squeeze(mean(dff_lt,1));
        dff_dk_avg = squeeze(mean(dff_dk,1));
        
        dff_lt_avg_pre = mean(dff_lt_avg(:,1:14),2);
        dff_lt_avg_pst = mean(dff_lt_avg(:,27:end),2);
        
        dff_dk_avg_pre = mean(dff_dk_avg(:,1:14),2);
        dff_dk_avg_pst = mean(dff_dk_avg(:,27:end),2);
        
        dff_lt_avg_pre_all = [dff_lt_avg_pre_all;dff_lt_avg_pre];
        dff_lt_avg_pst_all = [dff_lt_avg_pst_all;dff_lt_avg_pst];
        dff_dk_avg_pre_all = [dff_dk_avg_pre_all;dff_dk_avg_pre];
        dff_dk_avg_pst_all = [dff_dk_avg_pst_all;dff_dk_avg_pst];
        
        thr = 0.0238; %see below
        dff_lt_avg_pre_all_idx = dff_lt_avg_pre>thr;
        dff_lt_avg_pst_all_idx = dff_lt_avg_pst>thr;
        dff_dk_avg_pre_all_idx = dff_dk_avg_pre>thr;
        dff_dk_avg_pst_all_idx = dff_dk_avg_pst>thr;
        
        cell_num_all(m,1) = sum(dff_lt_avg_pre_all_idx);
        cell_num_all(m,2) = sum(dff_lt_avg_pst_all_idx);
        cell_num_all(m,3) = sum(dff_dk_avg_pre_all_idx);
        cell_num_all(m,4) = sum(dff_dk_avg_pst_all_idx);
        
        spd_idx = ltdkdata.mouse(m).file(s).spd_idx(idx_use);
        lt_dk_idx = ltdkdata.mouse(m).file(s).lt_dk_idx(idx_use);
        
        spd_idx_0 = spd_idx>0;
        lt_idx = lt_dk_idx==1;
        dk_idx = lt_dk_idx==2;
        
        idx_lt_spd = lt_idx&spd_idx_0;
        idx_spd_only = (~(lt_idx|dk_idx))&spd_idx_0;
        idx_dk_spd = dk_idx&spd_idx_0;
        
        idx_12 = zeros(size(spd_idx_0));
        for ii = 1:4
            idx_12 = idx_12 + ((spd_idx.*idx_lt_spd) == ii)*ii;
            idx_12 = idx_12 + ((spd_idx.*idx_spd_only) == ii)*(4+ii);
            idx_12 = idx_12 + ((spd_idx.*idx_dk_spd) == ii)*(8+ii);
        end
        idx_12(idx_12==0) = 13;
        
        pert_all = zeros(13,4);
        
        for ii = 1:13
            pert_all(ii,1) = sum((dff_lt_avg_pre_all_idx.*idx_12)==ii)/(eps+sum(dff_lt_avg_pre_all_idx))*100;
            pert_all(ii,2) = sum((dff_lt_avg_pst_all_idx.*idx_12)==ii)/(eps+sum(dff_lt_avg_pst_all_idx))*100;
            pert_all(ii,3) = sum((dff_dk_avg_pre_all_idx.*idx_12)==ii)/(eps+sum(dff_dk_avg_pre_all_idx))*100;
            pert_all(ii,4) = sum((dff_dk_avg_pst_all_idx.*idx_12)==ii)/(eps+sum(dff_dk_avg_pst_all_idx))*100;
        end
        
        fprintf('mouse %d, session %d is done. \n', m, s);
        
    end
    
end


%%
figure
bins = 0:0.001:0.1;
subplot(2,2,1);histogram(dff_lt_avg_pre_all,bins);
subplot(2,2,2);histogram(dff_lt_avg_pst_all,bins);
subplot(2,2,3);histogram(dff_dk_avg_pre_all,bins);
subplot(2,2,4);histogram(dff_dk_avg_pst_all,bins);

figure;
all = [dff_lt_avg_pre_all;
        dff_lt_avg_pst_all;
        dff_dk_avg_pre_all;
        dff_dk_avg_pst_all;];
histogram(all(all>0),bins);
mean(all(all>0))