%% figure_4_behaviorstate

%% 1. calcium activity during transitions - to light and to dark - mouse average
figname = 'Figure_4_A';
figure('name', figname, 'unit','normalized', 'position',...
    [0.1 0 0.8 0.9],'defaultAxesFontSize',18);

trans_name = {'elt_lt';'edrk_drk';};
               
mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];
vel_thr = 1; % cm/s - baseline noise is high here, so i choose a higher threshold
mkprop = {'-ro';'-bo';'-kv';'-kv'};
sname = {'Open field'; 'Light dark box'};
win = 20; % +/- 2s

for s = 1:2
    
    dff_lt_m = zeros(mnum,2*win+1);
    dff_dk_m = zeros(mnum,2*win+1);
    
    vel_lt_m = zeros(mnum,2*win+1);
    vel_dk_m = zeros(mnum,2*win+1);
    
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
        
        trans = ltdkdata.mouse(m).file(s).behavior.ttype(f_use);
        vel = ltdkdata.mouse(m).file(s).behavior.velocity(f_use)';
        
        trans_lt = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{1}),trans,'UniformOutput',false)));
        trans_dk = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{2}),trans,'UniformOutput',false)));
        
        dff_lt = extract_from_onset_window(dff,trans_lt,win);
        dff_dk = extract_from_onset_window(dff,trans_dk,win);
        
        vel_lt = extract_from_onset_window(vel,trans_lt,win);
        vel_dk = extract_from_onset_window(vel,trans_dk,win);
        
        dff_lt_avg = squeeze(nanmean(dff_lt,1));
        dff_lt_m(m,:) = nanmean(dff_lt_avg,1);
        
        dff_dk_avg = squeeze(nanmean(dff_dk,1));
        dff_dk_m(m,:) = nanmean(dff_dk_avg,1);
        
        vel_lt_m(m,:) = squeeze(nanmean(vel_lt,1));
        vel_dk_m(m,:) = squeeze(nanmean(vel_dk,1));
        
        fprintf('mouse %d, session %d is done. \n', m, s);
    end
    
    dff_lt_m  = zscore(dff_lt_m,0,2);
    dff_dk_m  = zscore(dff_dk_m,0,2);
    
    xx = (-win:win)/10;
    
    hax1 = subplot(2,2,s);
    
    yyaxis(hax1, 'left');
    plot_shadedline_group(hax1,xx,dff_lt_m,{mkprop{1},'MarkerFaceColor','w'});hold on;
    xlabel('Time relative to onset (s)'); ylabel('Average \DeltaF/F (zscore)');
%     legend('Trans to light','Trans to dark','location','northwest');
    yyaxis(hax1, 'right');
    plot_shadedline_group(hax1,xx,vel_lt_m,{mkprop{3},'MarkerFaceColor','w'});
    ylabel('Velocity');title([sname{s},'- Trans to light']);
    
%     legend('Trans to light','Trans to dark','location','northwest');
    
    hax2 = subplot(2,2,s+2);
    yyaxis(hax2, 'left');
%     plot_shadedline_group(hax,xx,dff_lt_m,{mkprop{1},'MarkerFaceColor','w'});hold on;
    plot_shadedline_group(hax2,xx,dff_dk_m,{mkprop{2},'MarkerFaceColor','w'});
    xlabel('Time relative to onset (s)'); ylabel('Average \DeltaF/F (zscore)');
%     legend('Trans to light','Trans to dark','location','northwest');
    
    yyaxis(hax2, 'right');
%     plot_shadedline_group(hax,xx,vel_lt_m,{mkprop{3},'MarkerFaceColor','w'});
    plot_shadedline_group(hax2,xx,vel_dk_m,{mkprop{4},'MarkerFaceColor','w'});
    ylabel('Velocity');title([sname{s},'- Trans to dark']);
%     legend('Trans to light','Trans to dark','location','northwest');
    
end

saveas(gcf, figname)

%% 2. calcium activity during transitions - to light and to dark - single neuron level
% is there any neuron that purely tuned by transition(to light/dark)
% instead of velocity

figname = 'Figure_4_B';
figure('name', figname, 'unit','normalized', 'position',...
    [0.1 0.1 0.8 0.6],'defaultAxesFontSize',18);

trans_name = {'elt_lt';'edrk_drk';};
               
mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];

mkprop = {'-ro';'-bo'};
sname = {'Open field'; 'Light dark box'};
win = 20; % +/- 2s
W = 5; % time window for data smoothing

for s = 1:2
    
    dff_lt_n = [];
    dff_dk_n = [];
    
    for m = 1:mnum
        
        dff   = ltdkdata.mouse(m).file(s).dFF;
        dffbin   = ltdkdata.mouse(m).file(s).dFF_bin;
        dff = dff_smooth(dff,W);
        
%         dff = dff.*dffbin;
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
        
        trans = ltdkdata.mouse(m).file(s).behavior.ttype(f_use);
        
        trans_lt = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{1}),trans,'UniformOutput',false)));
        trans_dk = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{2}),trans,'UniformOutput',false)));
        
        dff_lt = extract_from_onset_window(dff,trans_lt,win);
        dff_dk = extract_from_onset_window(dff,trans_dk,win);
        
        dff_lt_avg = squeeze(mean(dff_lt,1));
        kp_idx = ~isnan(dff_lt_avg(:,1));
        
        dff_lt_n0 = [dff_lt_n;dff_lt_avg(kp_idx,:)];
        dff_lt_n  = dff_lt_n0;
        
        dff_dk_avg = squeeze(nanmean(dff_dk,1));
        kp_idx = ~isnan(dff_dk_avg(:,1));
        dff_dk_n0 = [dff_dk_n;dff_dk_avg(kp_idx,:)];
        dff_dk_n  = dff_dk_n0;

        fprintf('mouse %d, session %d is done. \n', m, s);
    end
    
    dff_lt_n  = zscore(dff_lt_n,0,2);
    dff_dk_n  = zscore(dff_dk_n,0,2);
	
    lt_avg = mean(dff_lt_n(:,win:(win+win/2)),2);
    [~,lt_idx] = sort(lt_avg,'descend');
    
    dk_avg = mean(dff_dk_n(:,win:(win+win/2)),2);
    [~,dk_idx] = sort(dk_avg,'descend');
    
    hax1 = subplot(3,2,s);
    imagesc(dff_lt_n(lt_idx,:))
    ylabel('Cell ID');
    set(gca,'xtick',[])
    
    hax2 = subplot(3,2,s+2);
    imagesc(dff_lt_n(lt_idx,:))
    ylabel('Cell ID');
    set(gca,'xtick',[])
    
    hax3 = subplot(3,2,s+4);
    xx = (-win:win)/10;
    plot_shadedline_group(hax3,xx,dff_lt_n,{mkprop{1},'MarkerFaceColor','w'});hold on;
    plot_shadedline_group(hax3,xx,dff_dk_n,{mkprop{2},'MarkerFaceColor','w'});
    xlabel('Time relative to onset (s)'); 
    ylabel('Average \DeltaF/F (zscore)');
    title(sname{s});
    legend('Trans to light','Trans to dark','location','northwest');
    
end
% saveas(gcf, figname)

%%
figname = 'Figure_4_B';
figure('name', figname, 'unit','normalized', 'position',...
    [0.1 0.1 0.8 0.6],'defaultAxesFontSize',18);

trans_name = {'elt_lt';'edrk_drk';};
               
mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];

mkprop = {'-ro';'-bo'};
sname = {'Open field'; 'Light dark box'};
win = 20; % +/- 2s
W = 5; % time window for data smoothing

for s = 2:2
    
    dff_lt_n = [];
    dff_dk_n = [];
    
    for m = 1:mnum
        
        dff   = ltdkdata.mouse(m).file(s).dFF;
        dffbin   = ltdkdata.mouse(m).file(s).dFF_bin;
        dff = dff_smooth(dff,W);
        
%         dff = dff.*dffbin;
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
        
        trans = ltdkdata.mouse(m).file(s).behavior.ttype(f_use);
        
        trans_lt = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{1}),trans,'UniformOutput',false)));
        trans_dk = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{2}),trans,'UniformOutput',false)));
        
        dff_lt = extract_from_onset_window(dff,trans_lt,win);
        dff_dk = extract_from_onset_window(dff,trans_dk,win);
        
        dff_lt_avg = squeeze(mean(dff_lt,1));
        kp_idx = ~isnan(dff_lt_avg(:,1));
        
        dff_lt_n0 = [dff_lt_n;dff_lt_avg(kp_idx,:)];
        dff_lt_n  = dff_lt_n0;
        
        dff_dk_avg = squeeze(nanmean(dff_dk,1));
        kp_idx = ~isnan(dff_dk_avg(:,1));
        dff_dk_n0 = [dff_dk_n;dff_dk_avg(kp_idx,:)];
        dff_dk_n  = dff_dk_n0;

        fprintf('mouse %d, session %d is done. \n', m, s);
    end
    %%
    dff_lt_n  = zscore(dff_lt_n,0,2);
    dff_dk_n  = zscore(dff_dk_n,0,2);
	
    dff_lt_n_max = max(dff_lt_n,[],2);
    LTN = size(dff_lt_n_max(:),1);
    lt_pk_pos = zeros(LTN,1);
    
    for ii = 1:LTN
       pos_temp = find(dff_lt_n_max(ii)== dff_lt_n(ii,:));
       lt_pk_pos(ii) = pos_temp(1);
    end
    
    dff_dk_n_max = max(dff_dk_n,[],2);
    DKN = size(dff_dk_n_max(:),1);
    lt_dk_pos = zeros(DKN,1);
    
    for ii = 1:DKN
       pos_temp = find(dff_dk_n_max(ii)== dff_dk_n(ii,:));
       lt_dk_pos(ii) = pos_temp(1);
    end
    
    lt_avg = mean(dff_lt_n(:,win:(win+win/2)),2);
    [lt_pk_pos_st,lt_idx] = sort(lt_pk_pos,'ascend');
    
    dk_avg = mean(dff_dk_n(:,win:(win+win/2)),2);
    [dk_pk_pos_st,dk_idx] = sort(lt_dk_pos,'ascend');
    
    TN = size(dff_lt_n,2);
    hax1 = subplot(2,2,1);
    imagesc(dff_lt_n(lt_idx,:))
    line([win+1 win+1],get(gca,'ylim'),'color','r','linestyle','--')
    ylabel('Cell ID');
    set(gca,'xlim',[1 TN],'xtick',[])
    title('Trans Light');
    
    hax2 = subplot(2,2,2);
    imagesc(dff_dk_n(dk_idx,:))
    line([win+1 win+1],get(gca,'ylim'),'color','r','linestyle','--')
    ylabel('Cell ID');
    set(gca,'xtick',[])
    set(gca,'xlim',[1 TN],'xtick',[])
    title('Trans Dark');
    
    bins = 1:(2*win+1);
    hax3 = subplot(2,2,3);
    histogram(lt_pk_pos_st,bins);
    set(gca,'xlim',[1 TN],'ylim',[0 200],'xtick',[])
    line([win+1 win+1],get(gca,'ylim'),'color','r','linestyle','--')
    ylabel('Cell number');
    
    hax4 = subplot(2,2,4);
    histogram(dk_pk_pos_st,bins);
    set(gca,'xlim',[1 TN],'ylim',[0 200],'xtick',[])
    line([win+1 win+1],get(gca,'ylim'),'color','r','linestyle','--')
    ylabel('Cell number');
%     hax3 = subplot(3,2,s+4);
%     xx = (-win:win)/10;
%     plot_shadedline_group(hax3,xx,dff_lt_n,{mkprop{1},'MarkerFaceColor','w'});hold on;
%     plot_shadedline_group(hax3,xx,dff_dk_n,{mkprop{2},'MarkerFaceColor','w'});
%     xlabel('Time relative to onset (s)'); 
%     ylabel('Average \DeltaF/F (zscore)');
%     title(sname{s});
%     legend('Trans to light','Trans to dark','location','northwest');
    
end
% saveas(gcf, figname)

%% 2. - update calcium activity during transitions - to light and to dark - single neuron level
% is there any neuron that purely tuned by transition(to light/dark)
% instead of velocity

figname = 'Figure_4_B';
figure('name', figname, 'unit','normalized', 'position',...
    [0.1 0.1 0.8 0.6],'defaultAxesFontSize',18);

trans_name = {'elt_lt';'edrk_drk';};
               
mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];

mkprop = {'-ro';'-bo'};
sname = {'Open field'; 'Light dark box'};
win = 20; % +/- 2s
W = 5; % time window for data smoothing

for s = 2:2
    
    dff_lt_n = [];
    dff_dk_n = [];
    
    for m = 1:mnum
        
        dff   = ltdkdata.mouse(m).file(s).dFF;
        dffbin   = ltdkdata.mouse(m).file(s).dFF_bin;
        dff = dff_smooth(dff,W);
        
%         dff = dff.*dffbin;
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
        
        trans = ltdkdata.mouse(m).file(s).behavior.ttype(f_use);
        
        trans_lt = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{1}),trans,'UniformOutput',false)));
        trans_dk = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{2}),trans,'UniformOutput',false)));
        
        dff_lt = extract_from_onset_window(dff,trans_lt,win);
        dff_dk = extract_from_onset_window(dff,trans_dk,win);
        
        dff_lt_avg = squeeze(mean(dff_lt,1));
        kp_idx = ~isnan(dff_lt_avg(:,1));
        
        dff_lt_n0 = [dff_lt_n;dff_lt_avg(kp_idx,:)];
        dff_lt_n  = dff_lt_n0;
        
        dff_dk_avg = squeeze(nanmean(dff_dk,1));
        kp_idx = ~isnan(dff_dk_avg(:,1));
        dff_dk_n0 = [dff_dk_n;dff_dk_avg(kp_idx,:)];
        dff_dk_n  = dff_dk_n0;

        fprintf('mouse %d, session %d is done. \n', m, s);
    end
    
%     dff_lt_n  = zscore(dff_lt_n,0,2);
%     dff_dk_n  = zscore(dff_dk_n,0,2);
	
    lt_avg = mean(dff_lt_n(:,win:(win+win/2)),2);
    [~,lt_idx] = sort(lt_avg,'descend');
    
    dk_avg = mean(dff_dk_n(:,win:(win+win/2)),2);
    [~,dk_idx] = sort(dk_avg,'descend');
    
    hax1 = subplot(3,2,s);
    imagesc(dff_lt_n(lt_idx,:))
    ylabel('Cell ID');
    set(gca,'xtick',[])
    
    hax2 = subplot(3,2,s+2);
    imagesc(dff_lt_n(lt_idx,:))
    ylabel('Cell ID');
    set(gca,'xtick',[])
    
    hax3 = subplot(3,2,s+4);
    xx = (-win:win)/10;
    plot_shadedline_group(hax3,xx,dff_lt_n,{mkprop{1},'MarkerFaceColor','w'});hold on;
    plot_shadedline_group(hax3,xx,dff_dk_n,{mkprop{2},'MarkerFaceColor','w'});
    xlabel('Time relative to onset (s)'); 
    ylabel('Average \DeltaF/F (zscore)');
    title(sname{s});
    legend('Trans to light','Trans to dark','location','northwest');
    
end
% saveas(gcf, figname)



