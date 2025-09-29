%% figure_3_locomotion_velocity
% this script is to process the open field data
% the question we'd like to ask here is whether these neurons are related to
% locomotion velocity. If yes, how? 

%% 0. Calculate the smoothed velocity and locomotion distance - 
% DO NOT NEED TO RUN THIS

ltdkdata = add_smooth_vel_dist(ltdkdata);

%% 1. Correlation with locomotion velocity 
% neurons are correlated with locomotion velocity.

s = 1;% open field test
m = 1;
f_use = 1:9000;

figname = 'Figure_3_A';
figure('name', figname, 'unit','normalized', 'position',...
    [0.1 0.1 0.8 0.8],'defaultAxesFontSize',18);
dff   = ltdkdata.mouse(m).file(s).dFF;
dffbin   = ltdkdata.mouse(m).file(s).dFF_bin;
dff = dff.*dffbin;
idx_use = ~isnan(dff(:,1));
dff = dff(idx_use,f_use);
dff = dffnorm(dff,'max');
dff_avg = nanmean(dff,1);

vel = ltdkdata.mouse(m).file(s).behavior.velocity(f_use);

subplot(2,1,1)
imagesc(dff);
ylabel('Cell ID');

subplot(2,1,2)
yyaxis left;
plot(dff_avg);
ylabel('Avg \DeltaF/F (Normalized)');

yyaxis right;
plot(vel);
ylabel('Velocity (cm/s)');
xlabel('Frame number')

saveas(gcf,figname)

%% 2. Correlation coefficients -  individual neuron level vs average-of-all-neuron
figname = 'Figure_3_B_1';
figure('name', figname, 'unit','normalized', 'position',...
    [0.1 0.1 0.8 0.5],'defaultAxesFontSize',18);
sname = {'Open field'; 'Light dark box'};
mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];
rval_g = zeros(mnum,2);
rval_i = zeros(mnum,2);

rbins = -0.2:0.01:0.4;

for s = 1:2
    
    rval_all_neu = [];
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
        
        % calculate the correlation
        rval = zeros(cn,1);
        
        for n = 1:cn
            dff_n = dff(n,:);
            rval(n) = corr(dff_n(:),vel(:));
        end
        
        rval_all_neu0 = [rval_all_neu;rval];
        rval_all_neu  = rval_all_neu0;
        
        rval_g(m,s) = corr(dff_avg(:),vel(:));
        rval_i(m,s) = nanmean(rval);

        fprintf('mouse %d, session %d is done. \n', m, s);
    end
    
    subplot(1,2,s); histogram(rval_all_neu,rbins);
    grid on;
    xlabel('Pearson correlation coefficient');
    ylabel('Neuron number')
    title(sname{s})
end

figname = 'Figure_3_B';
file_dist = [figname,'_population.csv'];
writematrix(rval_g,file_dist)
file_dist = [figname,'_individual.csv'];
writematrix(rval_i,file_dist)

%% 3. Velocity as a function of averaged activity
bins = [0,0.5,1,2,4,8,16,32];
binc = bins(2:end)-1;
binn = length(bins);

mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];

figname = 'Figure_3_C';
figure('name', figname, 'unit','normalized', 'position',...
    [0.1 0.1 0.4 0.6],'defaultAxesFontSize',18);

mkprop = {'-ro';'-b^'};
for s = 1:2
    
    dff_vbin_avg = nan(m,binn-1);
    vel_bin      = nan(m,binn-1);
    
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
        
        dff_vbin = nan(cn,binn-1);
        
        A = 1;
        
        for b = 2:binn
            B = b;
            bin_idx = (vel>bins(A))&(vel<=bins(B));
            A = B;
            
            dff_vbin(:,b-1) = nanmean(dff(:,bin_idx),2);
            vel_bin(m,b-1)  = nanmean(vel(bin_idx));
            
        end
        
        dff_vbin_avg(m,:) = nanmean(dff_vbin,1);
        
        fprintf('mouse %d, session %d is done. \n', m, s);
    end
    
%     dff_vbin_avg = zscore(dff_vbin_avg,0,2);

    xx = nanmean(vel_bin,1);
    yy = nanmean(dff_vbin_avg,1);
    yy_err = std(dff_vbin_avg,0,1,'omitnan')/sqrt(size(dff_vbin_avg,1));
    shadedErrorBar(gca,xx, yy, yy_err, 'lineprops', {mkprop{s},'MarkerFaceColor','w'})
    xlabel('Avg velocity (cm/s)');
    ylabel('Avg \DeltaF/F (normalized)');
    set(gca,'ylim',[0 0.04],'xlim',[0 20])

end
legend('Open field', 'Light dark box')
saveas(gcf,figname)

%% 4. Start/stop neuron activity
figname = 'Figure_3_D';
figure('name', figname, 'unit','normalized', 'position',...
    [0.1 0.1 0.8 0.6],'defaultAxesFontSize',18);

mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];
vel_thr = 1; % cm/s - baseline noise is high here, so i choose a higher threshold
mkprop = {'-ro';'-bo'};
sname = {'Open field'; 'Light dark box'};
win = 20; % +/- 2s

for s = 1:2
    
    dff_onset_m = zeros(m,2*win+1);
    dff_offset_m = zeros(m,2*win+1);
    
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
        
        spk   = ltdkdata.mouse(m).file(s).spk(idx_use,f_use)>0;
        dff = dff(idx_use,f_use);
        dff = dffnorm(dff,'max');
        dff_avg = nanmean(dff,1);
        
        vel = ltdkdata.mouse(m).file(s).behavior.velocity(f_use);
        vel_bin = vel>vel_thr;
        [onset,offset, duration] = find_bin_onset(vel_bin);
        t_idx = duration>10;
        onset = onset(t_idx);
        offset = offset(t_idx);
        
        dff_onset = extract_from_onset_window(dff,onset,win);
        dff_offset = extract_from_onset_window(dff,offset,win);
        
        dff_onset_avg = squeeze(nanmean(dff_onset,1));
        dff_onset_m(m,:) = nanmean(dff_onset_avg,1);
        
        dff_offset_avg = squeeze(nanmean(dff_offset,1));
        dff_offset_m(m,:) = nanmean(dff_offset_avg,1);
        fprintf('mouse %d, session %d is done. \n', m, s);
    end
    
    dff_onset_m  = zscore(dff_onset_m,0,2);
    dff_offset_m = zscore(dff_offset_m,0,2);
    
    subplot(1,2,s)
    xx = (-win:win)/10;
    plot_shadedline_group(xx,dff_onset_m,{mkprop{1},'MarkerFaceColor','w'});hold on;
    xlabel('Time relative to onset (s)'); ylabel('Average \DeltaF/F (zscore)');title(sname{s});
    plot_shadedline_group(xx,dff_offset_m,{mkprop{2},'MarkerFaceColor','w'});
    xlabel('Time relative to onset (s)'); ylabel('Average \DeltaF/F (zscore)');
    legend('Initiation','Termination','location','northwest');
    
end
saveas(gcf, figname)
