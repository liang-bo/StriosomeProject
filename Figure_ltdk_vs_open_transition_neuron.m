%% Figure_ltdk_vs_open_transition_neuron

%% 1. calcium activity during transitions - to light and to dark - mouse average

trans_name = {'elt_lt';'edrk_drk';};
               
mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];

mkprop = {'-ro';'-bo';'-kv';'-kv'};
sname = {'Open field'; 'Light dark box'};

win = 20;        % window around transition +/- 2s
win_mid = 11:30; % winddow to judge the act inc neurons +/- 1s

nperm = 100; % random shuffle time
methods = 3; % 1. shuffle transition; 2. shuffle dff; 3. shuffle both

for s = 1:2
    
    for m = 1:mnum
        
        % get the dff
        dff      = ltdkdata.mouse(m).file(s).dFF;
        dffbin   = ltdkdata.mouse(m).file(s).dFF_bin;
        dff = dff.*dffbin;
        fn = size(dff,2);
        
        % only use the first 9000/18000 frames
        if fn<f_lim(s)
            f_use = 1:fn;
        else
            f_use = 1:f_lim(s);
        end
        dff = dff(:,f_use);
        dff = dffnorm(dff,'max');
        
        % get the transition vector
        trans = ltdkdata.mouse(m).file(s).behavior.ttype(f_use);
        trans_lt_vect = transpose(cell2mat(cellfun(@(x) strcmp(x,trans_name{1}),trans,'UniformOutput',false)));
        trans_lt = find(trans_lt_vect);
        
        trans_dk_vect = transpose(cell2mat(cellfun(@(x) strcmp(x,trans_name{2}),trans,'UniformOutput',false)));
        trans_dk = find(trans_dk_vect);
        
        % get the transition activity
        dff_lt = extract_from_onset_window(dff,trans_lt,win);
        dff_dk = extract_from_onset_window(dff,trans_dk,win);
 
        dff_lt_avg = zscore(squeeze(nanmean(dff_lt,1)),0,2);      
        dff_dk_avg = zscore(squeeze(nanmean(dff_dk,1)),0,2);

        tic;
        % get the shuffled data
        switch methods
            
            case 1 % shuffle behavior
                
                dff_lt_s = extract_from_onset_window(dff,find(shuffledmean(trans_lt_vect,nperm)),win);
                dff_dk_s = extract_from_onset_window(dff,find(shuffledmean(trans_dk_vect,nperm)),win);
                
            case 2 % shuffle dff
                
                dffs = shuffledmean(dff,nperm);
                dff_lt_s = extract_from_onset_window(dffs,trans_lt,win);
                dff_dk_s = extract_from_onset_window(dffs,trans_dk,win);
                
            case 3 % shuffle both
                
                dffs = shuffledmean(dff,nperm);
                dff_lt_s = extract_from_onset_window(dffs,find(shuffledmean(trans_lt_vect,nperm)),win);
                dff_dk_s = extract_from_onset_window(dffs,find(shuffledmean(trans_dk_vect,nperm)),win);
                
            otherwise
                
                fprintf('Please use a correct method! \n');
        end
        
        dff_lt_s_avg = zscore(squeeze(nanmean(dff_lt_s,1)),0,2);
        dff_dk_s_avg = zscore(squeeze(nanmean(dff_dk_s,1)),0,2);
        
        % find the neuron with high activity in win_mid (> mean + sigma)
        dff_lt_avg_win   = mean(dff_lt_avg(:,win_mid),2);
        dff_dk_avg_win   = mean(dff_dk_avg(:,win_mid),2);
        
        dff_lt_s_avg_win = mean(dff_lt_s_avg(:,win_mid),2);
        dff_dk_s_avg_win = mean(dff_dk_s_avg(:,win_mid),2);
        
        ltdkdata.mouse(m).file(s).idx(1).value = find_act_inc_neuron(dff_lt_avg_win);
        ltdkdata.mouse(m).file(s).idx(2).value = find_act_inc_neuron(dff_dk_avg_win);
        
        ltdkdata.mouse(m).file(s).idx(3).value = find_act_inc_neuron(dff_lt_s_avg_win);
        ltdkdata.mouse(m).file(s).idx(4).value = find_act_inc_neuron(dff_dk_s_avg_win);
        
        fprintf('mouse %d, session %d is done. run time: %4.2f seconds. \n', m, s, toc);
        
    end
    
end

%% 2. Plot the results

ovlp    = zeros(4,4,mnum);
ovlp_t  = zeros(4,4,mnum);

ovlp_s    = zeros(4,4,mnum);
ovlp_t_s  = zeros(4,4,mnum);

for m = 1:mnum
    
    idx_full = [ltdkdata.mouse(m).file(1).idx(1).value,...% open field: lt transition;
                ltdkdata.mouse(m).file(1).idx(2).value,...% open field: dk transition
                ltdkdata.mouse(m).file(2).idx(1).value,...% light dark: lt transition
                ltdkdata.mouse(m).file(2).idx(2).value,]; % light dark: dk transition    
                                                                
                                                                
                                                                
    
    idx_full_s = [ltdkdata.mouse(m).file(1).idx(3).value,...     % open field: lt transition - shuffle
                  ltdkdata.mouse(m).file(1).idx(4).value,...     % open field: dk transition - shuffle
                  ltdkdata.mouse(m).file(2).idx(3).value,...     % light dark: lt transition - shuffle
                  ltdkdata.mouse(m).file(2).idx(4).value,];      % light dark: dk transition - shuffle

    [ovlp(:,:,m),  ovlp_t(:,:,m)]   = overlap_rate_pw(idx_full);
    [ovlp_s(:,:,m),ovlp_t_s(:,:,m)] = overlap_rate_pw(idx_full_s);
    
end

ovlp_avg     = nanmean(ovlp,3);
ovlp_t_avg   = nanmean(ovlp_t,3);

ovlp_s_avg   = nanmean(ovlp_s,3);
ovlp_t_s_avg = nanmean(ovlp_t_s,3);

crange = [0.0 0.2];
xtlb = {'Open-Cent';'Open-Border';'LightDark-Light';'LightDark-Dark';};
figure;
subplot(2,2,1)
imagesc(ovlp_avg)
title('Observed value')
axis square; caxis(crange); hc = colorbar; ylabel(hc,'Overlap')
set(gca,'xtick',1:4,'xticklabel',xtlb,'XTickLabelRotation',45,...
        'ytick',1:4,'yticklabel',xtlb)

subplot(2,2,2)
imagesc(ovlp_t_avg)
title('Observed value - chance')
axis square; caxis(crange); hc = colorbar; ylabel(hc,'Overlap')
set(gca,'xtick',1:4,'xticklabel',xtlb,'XTickLabelRotation',45,...
        'ytick',1:4,'yticklabel',xtlb)
    
subplot(2,2,3)
imagesc(ovlp_s_avg)
title('Shuffled value')
axis square; caxis(crange); hc = colorbar; ylabel(hc,'Overlap')
set(gca,'xtick',1:4,'xticklabel',xtlb,'XTickLabelRotation',45,...
        'ytick',1:4,'yticklabel',xtlb)
    
subplot(2,2,4)
imagesc(ovlp_t_s_avg)
title('Shuffled value - chance')
axis square; caxis(crange); hc = colorbar; ylabel(hc,'Overlap')
set(gca,'xtick',1:4,'xticklabel',xtlb,'XTickLabelRotation',45,...
        'ytick',1:4,'yticklabel',xtlb)
    