%% load the dataset first

%% #1 calculate the cell subcategories

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

for s = 2
    
    cell_perct_lt_spd   = zeros(mnum,4);
    cell_perct_spd_only = zeros(mnum,4);
    cell_perct_dk_spd   = zeros(mnum,4);
    
    cell_num_all = zeros(mnum,4);
    pert_all = zeros(13,4,mnum);
    
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
        
        thr = 0.0238; % determined by session 2
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
        
        for ii = 1:13
            pert_all(ii,1,m) = sum((dff_lt_avg_pre_all_idx.*idx_12)==ii)/(eps+sum(dff_lt_avg_pre_all_idx))*100;
            pert_all(ii,2,m) = sum((dff_lt_avg_pst_all_idx.*idx_12)==ii)/(eps+sum(dff_lt_avg_pst_all_idx))*100;
            pert_all(ii,3,m) = sum((dff_dk_avg_pre_all_idx.*idx_12)==ii)/(eps+sum(dff_dk_avg_pre_all_idx))*100;
            pert_all(ii,4,m) = sum((dff_dk_avg_pst_all_idx.*idx_12)==ii)/(eps+sum(dff_dk_avg_pst_all_idx))*100;
        end
        
        fprintf('mouse %d, session %d is done. \n', m, s);
        
    end
    
end

pert_all_avg = mean(pert_all,3); % average percentage of all 13 categories

%% #2 decide the threshold
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
mean(all(all>0)) % threshold

%% #3 concatenate cell categories all mice
cnum = zeros(10,1);
for m = 1:10
    cnum(m) = sum(ltdkdata.mouse(m).file(s).idx_reg>0);
end
cellnum_all = sum(cnum);

data = cell(cellnum_all,9);
s=2;
p=1;
for m = 1:10
    
    id = ltdkdata.mouse(m).id;
    
    idx_reg = ltdkdata.mouse(m).file(s).idx_reg>0;
    cn = sum(idx_reg);
    spd_idx = ltdkdata.mouse(m).file(s).spd_idx(idx_reg);
    lt_dk_idx = ltdkdata.mouse(m).file(s).lt_dk_idx(idx_reg);
    auc = ltdkdata.mouse(m).file(s).cell.acceleration.auc_comm;
    
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
    
    
    thr = 0.0238; %see below
    dff_lt_avg_pre_all_idx = dff_lt_avg_pre>thr;
    dff_lt_avg_pst_all_idx = dff_lt_avg_pst>thr;
    dff_dk_avg_pre_all_idx = dff_dk_avg_pre>thr;
    dff_dk_avg_pst_all_idx = dff_dk_avg_pst>thr;
    
    for ii = 1:cn
        data{p,1} = id;
        data{p,2} = ii;
        data{p,3} = spd_idx(ii);
        data{p,4} = lt_dk_idx(ii);
        data{p,5} = auc(ii);
        data{p,6} = dff_lt_avg_pre_all_idx(ii);
        data{p,7} = dff_lt_avg_pst_all_idx(ii);
        data{p,8} = dff_dk_avg_pre_all_idx(ii);
        data{p,9} = dff_dk_avg_pst_all_idx(ii);
        p = p + 1;
    end
    
end

xlswrite('Lightdark_data_summary.xlsx',data)

