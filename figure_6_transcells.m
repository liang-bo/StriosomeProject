% figure_6_transcells

mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];

% figname = 'Figure_';
% figure('name', figname, 'unit','normalized', 'position',...
%     [0.1 0.1 0.4 0.6],'defaultAxesFontSize',18);

mkprop = {'-ro';'-b^'};

for s = 1%:2
    
    vel_bin      = nan(mnum,binn-1);
    gof_mm = [];
    
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
        

        fprintf('mouse %d, session %d is done. \n', m, s);
        
    end

end

%% 1. light or dark only cell activity in transitions
trans_name = {'elt_lt';'edrk_drk';};
               
mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];
win = 20;

for s = 2%:2
    
    act_consis_lt_lt = [];
    act_consis_lt_dk = [];
    
    act_consis_dk_lt = [];
    act_consis_dk_dk = [];
    
    for m = 1:mnum
        
        dff   = ltdkdata.mouse(m).file(s).dFF;
        dffbin   = ltdkdata.mouse(m).file(s).dFF_bin;
        dff = dff.*dffbin;
        fn = size(dff,2);
        
        spd_idx = ltdkdata.mouse(m).file(s).spd_idx;
        lt_dk_idx = ltdkdata.mouse(m).file(s).lt_dk_idx;
        
        spd_idx = spd_idx>0;
        lt_idx = lt_dk_idx==1;
        dk_idx = lt_dk_idx==2;
        
        idx_lt_only = find(lt_idx&(~spd_idx));
        idx_dk_only = find(dk_idx&(~spd_idx));
        
        lt_num = numel(idx_lt_only);
        dk_num = numel(idx_dk_only);
                
        if fn<f_lim(s)
            f_use = 1:fn;
        else
            f_use = 1:f_lim(s);
        end
        
        trans = ltdkdata.mouse(m).file(s).behavior.ttype(f_use);
        trans_lt = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{1}),trans,'UniformOutput',false)));
        trans_dk = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{2}),trans,'UniformOutput',false)));
        
        dff_lt = extract_from_onset_window(dff,trans_lt,win);
        dff_dk = extract_from_onset_window(dff,trans_dk,win);
        
        trans_lt_num = numel(trans_lt);
        trans_dk_num = numel(trans_dk);
        
        % light cell trans in light
        if trans_lt_num>3
            
            act_mat_lt_lt = zeros(lt_num,trans_lt_num);
            
            for ii = 1:lt_num
                
                for jj = 1:trans_lt_num
                    
                    act_mat_lt_lt(ii,jj) = sum(dff_lt(jj,idx_lt_only(ii),:))>0;
                    
                end
                
            end
            
            act_consis_lt_lt0 = [act_consis_lt_lt;mean(act_mat_lt_lt,2);];
            act_consis_lt_lt = act_consis_lt_lt0;
        end
            
        % light cell trans in dark
        if trans_dk_num>3
            
            act_mat_lt_dk = zeros(lt_num,trans_dk_num);
            
            for ii = 1:lt_num
                
                for jj = 1:trans_dk_num
                    
                    act_mat_lt_dk(ii,jj) = sum(dff_dk(jj,idx_lt_only(ii),:))>0;
                    
                end
                
            end
            
            act_consis_lt_dk0 = [act_consis_lt_dk;mean(act_mat_lt_dk,2);];
            act_consis_lt_dk = act_consis_lt_dk0;
        end
        
        % dark cell trans in light
        if trans_lt_num>3
            
            act_mat_dk_lt = zeros(dk_num,trans_lt_num);
            
            for ii = 1:dk_num
                
                for jj = 1:trans_lt_num
                    
                    act_mat_dk_lt(ii,jj) = sum(dff_lt(jj,idx_dk_only(ii),:))>0;
                    
                end
                
            end
            
            act_consis_dk_lt0 = [act_consis_dk_lt;mean(act_mat_dk_lt,2);];
            act_consis_dk_lt = act_consis_dk_lt0;
        end
            
        % light cell trans in dark
        if trans_dk_num>3
            
            act_mat_dk_dk = zeros(lt_num,trans_dk_num);
            
            for ii = 1:dk_num
                
                for jj = 1:trans_dk_num
                    
                    act_mat_dk_dk(ii,jj) = sum(dff_dk(jj,idx_dk_only(ii),:))>0;
                    
                end
                
            end
            
            act_consis_dk_dk0 = [act_consis_dk_dk;mean(act_mat_dk_dk,2);];
            act_consis_dk_dk = act_consis_dk_dk0;
            
        end
         
    end
    
end

%% Light cell trans to light example

S = act_mat_lt_lt;
[n,k] = size(S);
h = nan*S;
mksz = 12;
rd = [255 0 81]/255;
figure('unit','points','position',[10 10 mksz*k mksz*n])
for i = 1:n
    for j = 1:k
        if S(i,j) ~= 0
           e = 'w'; f = S(i,j)*rd;
        else
           e = 'w'; f = [0.8 0.8 0.8];
        end
        h(i,j) = plot(j,i,'s','MarkerEdgeColor',e,'MarkerFaceColor',f,'MarkerSize',12); hold on
    end
end
set(gca,'xlim',[0 k+0.5],'ylim',[0 n+0.5],...
        'xtick',[1,5:5:k],'ytick',[1,5:5:n])
xlabel('Transition to light');
ylabel('Cell');

%% Light cell trans to dark example

S = act_mat_lt_dk;
[n,k] = size(S);
h = nan*S;
mksz = 12;
bl = [0 128 166]/255;
figure('unit','points','position',[10 10 mksz*k mksz*n])
for i = 1:n
    for j = 1:k
        if S(i,j) ~= 0
           e = 'w'; f = S(i,j)*bl;
        else
           e = 'w'; f = [0.8 0.8 0.8];
        end
        h(i,j) = plot(j,i,'s','MarkerEdgeColor',e,'MarkerFaceColor',f,'MarkerSize',12); hold on
    end
end
set(gca,'xlim',[0 k+0.5],'ylim',[0 n+0.5],...
        'xtick',[1,5:5:k],'ytick',[1,5:5:n])
xlabel('Transition to dark');
ylabel('Cell');

%% Dark cell trans to light example

S = act_mat_dk_lt;
[n,k] = size(S);
h = nan*S;
mksz = 12;
rd = [255 0 81]/255;
figure('unit','points','position',[10 10 mksz*k mksz*n])
for i = 1:n
    for j = 1:k
        if S(i,j) ~= 0
           e = 'w'; f = S(i,j)*rd;
        else
           e = 'w'; f = [0.8 0.8 0.8];
        end
        h(i,j) = plot(j,i,'s','MarkerEdgeColor',e,'MarkerFaceColor',f,'MarkerSize',12); hold on
    end
end
set(gca,'xlim',[0 k+0.5],'ylim',[0 n+0.5],...
        'xtick',[1,5:5:k],'ytick',[1,5:5:n])
xlabel('Transition to light');
ylabel('Cell');

%% Dark cell trans to dark example

S = act_mat_dk_dk;
[n,k] = size(S);
h = nan*S;
mksz = 12;
bl = [0 128 166]/255;
figure('unit','points','position',[10 10 mksz*k mksz*n])
for i = 1:n
    for j = 1:k
        if S(i,j) ~= 0
           e = 'w'; f = S(i,j)*bl;
        else
           e = 'w'; f = [0.8 0.8 0.8];
        end
        h(i,j) = plot(j,i,'s','MarkerEdgeColor',e,'MarkerFaceColor',f,'MarkerSize',12); hold on
    end
end
set(gca,'xlim',[0 k+0.5],'ylim',[0 n+0.5],...
        'xtick',[1,5:5:k],'ytick',[1,5:5:n])
xlabel('Transition to dark');
ylabel('Cell');

%% histogram of four types

figure('Unit','normalized','Position',[0.1 0.1 0.6 0.8]);
bins = 0:10:100;
rd = [255 0 81]/255;
bl = [0 128 166]/255;
subplot(2,2,1)
histogram(act_consis_lt_lt*100,bins,'Normalization','Prob','FaceColor',rd)
xlabel('Consistency (% of transitions)');
ylabel('Probability')
title('Light cell - transition to light');

subplot(2,2,2)
histogram(act_consis_lt_dk*100,bins,'Normalization','Prob','FaceColor',bl)
xlabel('Consistency (% of transitions)');
ylabel('Probability')
title('Light cell - transition to dark');

subplot(2,2,3)
histogram(act_consis_dk_lt*100,bins,'Normalization','Prob','FaceColor',rd)
xlabel('Consistency (% of transitions)');
ylabel('Probability')
title('Dark cell - transition to light');

subplot(2,2,4)
histogram(act_consis_dk_dk*100,bins,'Normalization','Prob','FaceColor',bl)
xlabel('Consistency (% of transitions)');
ylabel('Probability')
title('Dark cell - transition to dark');


%% 1.-update 12/3/20- light or dark only cell activity in transitions
trans_name = {'elt_lt';'edrk_drk';};
               
mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];
win = 20;
%flag = 2;% pre
flag = 1;% post

for s = 1%:2
    
    act_consis_lt_lt = [];
    act_consis_lt_dk = [];
    
    act_consis_dk_lt = [];
    act_consis_dk_dk = [];
    
    for m = 1:mnum
        
        dff   = ltdkdata.mouse(m).file(s).dFF;
        dffbin   = ltdkdata.mouse(m).file(s).dFF_bin;
        dff = dff.*dffbin;
        fn = size(dff,2);
        
        spd_idx = ltdkdata.mouse(m).file(s).spd_idx;
        lt_dk_idx = ltdkdata.mouse(m).file(s).lt_dk_idx;
        
        spd_idx = spd_idx>0;
        lt_idx = lt_dk_idx==1;
        dk_idx = lt_dk_idx==2;
        
        idx_lt_only = find(lt_idx&(~spd_idx));
        idx_dk_only = find(dk_idx&(~spd_idx));
        
        lt_num = numel(idx_lt_only);
        dk_num = numel(idx_dk_only);
                
        if fn<f_lim(s)
            f_use = 1:fn;
        else
            f_use = 1:f_lim(s);
        end
        
        trans = ltdkdata.mouse(m).file(s).behavior.ttype(f_use);
        trans_lt = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{1}),trans,'UniformOutput',false)));
        trans_dk = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{2}),trans,'UniformOutput',false)));
        
        dff_lt = extract_from_onset_window0(dff,trans_lt,win,flag);
        dff_dk = extract_from_onset_window0(dff,trans_dk,win,flag);
        
        trans_lt_num = numel(trans_lt);
        trans_dk_num = numel(trans_dk);
        
        % light cell trans in light
        if trans_lt_num>3
            
            act_mat_lt_lt = zeros(lt_num,trans_lt_num);
            
            for ii = 1:lt_num
                
                for jj = 1:trans_lt_num
                    
                    act_mat_lt_lt(ii,jj) = sum(dff_lt(jj,idx_lt_only(ii),:))>0;
                    
                end
                
            end
            
            act_consis_lt_lt0 = [act_consis_lt_lt;mean(act_mat_lt_lt,2);];
            act_consis_lt_lt = act_consis_lt_lt0;
        end
            
        % light cell trans in dark
        if trans_dk_num>3
            
            act_mat_lt_dk = zeros(lt_num,trans_dk_num);
            
            for ii = 1:lt_num
                
                for jj = 1:trans_dk_num
                    
                    act_mat_lt_dk(ii,jj) = sum(dff_dk(jj,idx_lt_only(ii),:))>0;
                    
                end
                
            end
            
            act_consis_lt_dk0 = [act_consis_lt_dk;mean(act_mat_lt_dk,2);];
            act_consis_lt_dk = act_consis_lt_dk0;
        end
        
        % dark cell trans in light
        if trans_lt_num>3
            
            act_mat_dk_lt = zeros(dk_num,trans_lt_num);
            
            for ii = 1:dk_num
                
                for jj = 1:trans_lt_num
                    
                    act_mat_dk_lt(ii,jj) = sum(dff_lt(jj,idx_dk_only(ii),:))>0;
                    
                end
                
            end
            
            act_consis_dk_lt0 = [act_consis_dk_lt;mean(act_mat_dk_lt,2);];
            act_consis_dk_lt = act_consis_dk_lt0;
        end
            
        % light cell trans in dark
        if trans_dk_num>3
            
            act_mat_dk_dk = zeros(lt_num,trans_dk_num);
            
            for ii = 1:dk_num
                
                for jj = 1:trans_dk_num
                    
                    act_mat_dk_dk(ii,jj) = sum(dff_dk(jj,idx_dk_only(ii),:))>0;
                    
                end
                
            end
            
            act_consis_dk_dk0 = [act_consis_dk_dk;mean(act_mat_dk_dk,2);];
            act_consis_dk_dk = act_consis_dk_dk0;
            
        end
         
    end
    
end

% histogram of four types
win_name = {'Pre+Post-';'Pre-';'Post-'};
figure('Unit','normalized','Position',[0.1 0.1 0.6 0.8]);
bins = 0:10:100;
rd = [255 0 81]/255;
bl = [0 128 166]/255;
subplot(2,2,1)
histogram(act_consis_lt_lt*100,bins,'Normalization','Prob','FaceColor',rd)
xlabel('Consistency (% of transitions)');
ylabel('Probability')
title(['Light cell -',win_name{flag}, 'light']);

subplot(2,2,2)
histogram(act_consis_lt_dk*100,bins,'Normalization','Prob','FaceColor',bl)
xlabel('Consistency (% of transitions)');
ylabel('Probability')
title(['Light cell -',win_name{flag}, 'dark']);

subplot(2,2,3)
histogram(act_consis_dk_lt*100,bins,'Normalization','Prob','FaceColor',rd)
xlabel('Consistency (% of transitions)');
ylabel('Probability')
title(['Dark cell -',win_name{flag}, 'light']);

subplot(2,2,4)
histogram(act_consis_dk_dk*100,bins,'Normalization','Prob','FaceColor',bl)
xlabel('Consistency (% of transitions)');
ylabel('Probability')
title(['Dark cell -',win_name{flag}, 'dark']);

%% 1.-update 12/4/20- light or dark only cell activity in transitions
trans_name = {'elt_lt';'edrk_drk';};
               
mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];
win = 20;

%flag = 2;% pre
for flag = 1% post
    
    for s = 2%:2
        
        act_consis_lt_lt = [];
        act_consis_lt_dk = [];
        
        act_consis_dk_lt = [];
        act_consis_dk_dk = [];
        
        for m = 1:mnum
            
            dff   = ltdkdata.mouse(m).file(s).dFF;
            dffbin   = ltdkdata.mouse(m).file(s).dFF_bin;
            dff = dff.*dffbin;
            fn = size(dff,2);
            
            spd_idx = ltdkdata.mouse(m).file(s).spd_idx;
            lt_dk_idx = ltdkdata.mouse(m).file(s).lt_dk_idx;
            
            spd_idx = spd_idx>0;
            lt_idx = lt_dk_idx==1;
            dk_idx = lt_dk_idx==2;
            
            idx_lt_only = find(lt_idx&(~spd_idx));
            idx_dk_only = find(dk_idx&(~spd_idx));
            
            lt_num = numel(idx_lt_only);
            dk_num = numel(idx_dk_only);
            
            if fn<f_lim(s)
                f_use = 1:fn;
            else
                f_use = 1:f_lim(s);
            end
            
            trans = ltdkdata.mouse(m).file(s).behavior.ttype(f_use);
            trans_lt = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{1}),trans,'UniformOutput',false)));
            trans_dk = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{2}),trans,'UniformOutput',false)));
            
            dff_lt = extract_from_onset_window0(dff,trans_lt,win,flag);
            dff_dk = extract_from_onset_window0(dff,trans_dk,win,flag);
            
            trans_lt_num = numel(trans_lt);
            trans_dk_num = numel(trans_dk);
            
            % light cell trans in light
            if trans_lt_num>3
                
                act_mat_lt_lt = zeros(lt_num,trans_lt_num);
                
                for ii = 1:lt_num
                    
                    for jj = 1:trans_lt_num
                        
                        act_mat_lt_lt(ii,jj) = sum(dff_lt(jj,idx_lt_only(ii),:))>0;
                        
                    end
                    
                end
                
                act_consis_lt_lt0 = [act_consis_lt_lt;mean(act_mat_lt_lt,2);];
                act_consis_lt_lt = act_consis_lt_lt0;
            end
            
            % light cell trans in dark
            if trans_dk_num>3
                
                act_mat_lt_dk = zeros(lt_num,trans_dk_num);
                
                for ii = 1:lt_num
                    
                    for jj = 1:trans_dk_num
                        
                        act_mat_lt_dk(ii,jj) = sum(dff_dk(jj,idx_lt_only(ii),:))>0;
                        
                    end
                    
                end
                
                act_consis_lt_dk0 = [act_consis_lt_dk;mean(act_mat_lt_dk,2);];
                act_consis_lt_dk = act_consis_lt_dk0;
            end
            
            % dark cell trans in light
            if trans_lt_num>3
                
                act_mat_dk_lt = zeros(dk_num,trans_lt_num);
                
                for ii = 1:dk_num
                    
                    for jj = 1:trans_lt_num
                        
                        act_mat_dk_lt(ii,jj) = sum(dff_lt(jj,idx_dk_only(ii),:))>0;
                        
                    end
                    
                end
                
                act_consis_dk_lt0 = [act_consis_dk_lt;mean(act_mat_dk_lt,2);];
                act_consis_dk_lt = act_consis_dk_lt0;
            end
            
            % light cell trans in dark
            if trans_dk_num>3
                
                act_mat_dk_dk = zeros(lt_num,trans_dk_num);
                
                for ii = 1:dk_num
                    
                    for jj = 1:trans_dk_num
                        
                        act_mat_dk_dk(ii,jj) = sum(dff_dk(jj,idx_dk_only(ii),:))>0;
                        
                    end
                    
                end
                
                act_consis_dk_dk0 = [act_consis_dk_dk;mean(act_mat_dk_dk,2);];
                act_consis_dk_dk = act_consis_dk_dk0;
                
            end
            
        end
        
    end
    
    % histogram of four types
    win_name = {'Pre+Post-';'Pre-';'Post-'};
    figure('Unit','normalized','Position',[0.1 0.1 0.6 0.8]);
    bins = 0:5:100;
    rd = [255 0 81]/255;
    bl = [0 128 166]/255;
    subplot(2,2,1)
    histogram(act_consis_lt_lt*100,bins,'Normalization','Prob','FaceColor',rd)
    set(gca,'ylim',[0 1])
    xlabel('Consistency (% of transitions)');
    ylabel('Probability')
    title(['Light cell -',win_name{flag}, 'light']);
    grid on;
    
    subplot(2,2,2)
    histogram(act_consis_lt_dk*100,bins,'Normalization','Prob','FaceColor',bl)
    set(gca,'ylim',[0 1])
    xlabel('Consistency (% of transitions)');
    ylabel('Probability')
    title(['Light cell -',win_name{flag}, 'dark']);
    grid on;
    
    subplot(2,2,3)
    histogram(act_consis_dk_lt*100,bins,'Normalization','Prob','FaceColor',rd)
    set(gca,'ylim',[0 1])
    xlabel('Consistency (% of transitions)');
    ylabel('Probability')
    title(['Dark cell -',win_name{flag}, 'light']);
    grid on;
    
    subplot(2,2,4)
    histogram(act_consis_dk_dk*100,bins,'Normalization','Prob','FaceColor',bl)
    set(gca,'ylim',[0 1])
    xlabel('Consistency (% of transitions)');
    ylabel('Probability')
    title(['Dark cell -',win_name{flag}, 'dark']);
    grid on;
    
end

%% 1.-update 12/4/20- light&spd or dark&spd only cell activity in transitions
trans_name = {'elt_lt';'edrk_drk';};
               
mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];
win = 20;

%flag = 2;% pre
for flag = 1% post
    
    for s = 2%:2
        
        act_consis_lt_lt = [];
        act_consis_lt_dk = [];
        
        act_consis_dk_lt = [];
        act_consis_dk_dk = [];
        
        for m = 1:mnum
            
            dff   = ltdkdata.mouse(m).file(s).dFF;
            dffbin   = ltdkdata.mouse(m).file(s).dFF_bin;
            dff = dff.*dffbin;
            fn = size(dff,2);
            
            spd_idx = ltdkdata.mouse(m).file(s).spd_idx;
            lt_dk_idx = ltdkdata.mouse(m).file(s).lt_dk_idx;
            
            spd_idx = spd_idx>0;
            lt_idx = lt_dk_idx==1;
            dk_idx = lt_dk_idx==2;
            
            idx_lt_only = find(lt_idx&(spd_idx));
            idx_dk_only = find(dk_idx&(spd_idx));
            
            lt_num = numel(idx_lt_only);
            dk_num = numel(idx_dk_only);
            
            if fn<f_lim(s)
                f_use = 1:fn;
            else
                f_use = 1:f_lim(s);
            end
            
            trans = ltdkdata.mouse(m).file(s).behavior.ttype(f_use);
            trans_lt = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{1}),trans,'UniformOutput',false)));
            trans_dk = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{2}),trans,'UniformOutput',false)));
            
            dff_lt = extract_from_onset_window0(dff,trans_lt,win,flag);
            dff_dk = extract_from_onset_window0(dff,trans_dk,win,flag);
            
            trans_lt_num = numel(trans_lt);
            trans_dk_num = numel(trans_dk);
            
            % light cell trans in light
            if trans_lt_num>3
                
                act_mat_lt_lt = zeros(lt_num,trans_lt_num);
                
                for ii = 1:lt_num
                    
                    for jj = 1:trans_lt_num
                        
                        act_mat_lt_lt(ii,jj) = sum(dff_lt(jj,idx_lt_only(ii),:))>0;
                        
                    end
                    
                end
                
                act_consis_lt_lt0 = [act_consis_lt_lt;mean(act_mat_lt_lt,2);];
                act_consis_lt_lt = act_consis_lt_lt0;
            end
            
            % light cell trans in dark
            if trans_dk_num>3
                
                act_mat_lt_dk = zeros(lt_num,trans_dk_num);
                
                for ii = 1:lt_num
                    
                    for jj = 1:trans_dk_num
                        
                        act_mat_lt_dk(ii,jj) = sum(dff_dk(jj,idx_lt_only(ii),:))>0;
                        
                    end
                    
                end
                
                act_consis_lt_dk0 = [act_consis_lt_dk;mean(act_mat_lt_dk,2);];
                act_consis_lt_dk = act_consis_lt_dk0;
            end
            
            % dark cell trans in light
            if trans_lt_num>3
                
                act_mat_dk_lt = zeros(dk_num,trans_lt_num);
                
                for ii = 1:dk_num
                    
                    for jj = 1:trans_lt_num
                        
                        act_mat_dk_lt(ii,jj) = sum(dff_lt(jj,idx_dk_only(ii),:))>0;
                        
                    end
                    
                end
                
                act_consis_dk_lt0 = [act_consis_dk_lt;mean(act_mat_dk_lt,2);];
                act_consis_dk_lt = act_consis_dk_lt0;
            end
            
            % light cell trans in dark
            if trans_dk_num>3
                
                act_mat_dk_dk = zeros(lt_num,trans_dk_num);
                
                for ii = 1:dk_num
                    
                    for jj = 1:trans_dk_num
                        
                        act_mat_dk_dk(ii,jj) = sum(dff_dk(jj,idx_dk_only(ii),:))>0;
                        
                    end
                    
                end
                
                act_consis_dk_dk0 = [act_consis_dk_dk;mean(act_mat_dk_dk,2);];
                act_consis_dk_dk = act_consis_dk_dk0;
                
            end
            
        end
        
    end
    
    % histogram of four types
    win_name = {'Pre+Post-';'Pre-';'Post-'};
    figure('Unit','normalized','Position',[0.1 0.1 0.6 0.8]);
    bins = 0:5:100;
    rd = [255 0 81]/255;
    bl = [0 128 166]/255;
    subplot(2,2,1)
    histogram(act_consis_lt_lt*100,bins,'Normalization','Prob','FaceColor',rd)
    set(gca,'ylim',[0 1])
    xlabel('Consistency (% of transitions)');
    ylabel('Probability')
    title(['Light&speed cell -',win_name{flag}, 'light']);
    grid on;
    
    subplot(2,2,2)
    histogram(act_consis_lt_dk*100,bins,'Normalization','Prob','FaceColor',bl)
    set(gca,'ylim',[0 1])
    xlabel('Consistency (% of transitions)');
    ylabel('Probability')
    title(['Light&speed cell -',win_name{flag}, 'dark']);
    grid on;
    
    subplot(2,2,3)
    histogram(act_consis_dk_lt*100,bins,'Normalization','Prob','FaceColor',rd)
    set(gca,'ylim',[0 1])
    xlabel('Consistency (% of transitions)');
    ylabel('Probability')
    title(['Dark&speed cell -',win_name{flag}, 'light']);
    grid on;
    
    subplot(2,2,4)
    histogram(act_consis_dk_dk*100,bins,'Normalization','Prob','FaceColor',bl)
    set(gca,'ylim',[0 1])
    xlabel('Consistency (% of transitions)');
    ylabel('Probability')
    title(['Dark&speed cell -',win_name{flag}, 'dark']);
    grid on;
    
end

%% speed or other cell activity in transitions
trans_name = {'elt_lt';'edrk_drk';};
               
mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];
win = 20;

for s = 2%:2
    
    act_consis_spd_lt = [];
    act_consis_spd_dk = [];
    
    act_consis_otr_lt = [];
    act_consis_otr_dk = [];
    
    for m = 1:mnum
        
        dff   = ltdkdata.mouse(m).file(s).dFF;
        dffbin   = ltdkdata.mouse(m).file(s).dFF_bin;
        dff = dff.*dffbin;
        fn  = size(dff,2);
        
        spd_idx = ltdkdata.mouse(m).file(s).spd_idx>0;
        lt_dk_idx = ltdkdata.mouse(m).file(s).lt_dk_idx>0;
                
        idx_spd = find(spd_idx); 
        idx_otr = find(~(lt_dk_idx|spd_idx));
        
        spd_num = numel(idx_spd);
        otr_num  = numel(idx_otr);
                
        if fn<f_lim(s)
            f_use = 1:fn;
        else
            f_use = 1:f_lim(s);
        end
        
        trans = ltdkdata.mouse(m).file(s).behavior.ttype(f_use);
        trans_lt = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{1}),trans,'UniformOutput',false)));
        trans_dk = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{2}),trans,'UniformOutput',false)));
        
        dff_lt = extract_from_onset_window(dff,trans_lt,win);
        dff_dk = extract_from_onset_window(dff,trans_dk,win);
        
        trans_lt_num = numel(trans_lt);
        trans_dk_num = numel(trans_dk);
        
        % speed cell trans in light
        if trans_lt_num>3
            
            act_mat_spd_lt = zeros(spd_num,trans_lt_num);
            
            for ii = 1:spd_num
                
                for jj = 1:trans_lt_num
                    
                    act_mat_spd_lt(ii,jj) = sum(dff_lt(jj,idx_spd(ii),:))>0;
                    
                end
                
            end
            
            act_consis_spd_lt0 = [act_consis_spd_lt;mean(act_mat_spd_lt,2);];
            act_consis_spd_lt  = act_consis_spd_lt0;
        end
            
        % spd cell trans in dark
        if trans_dk_num>3
            
            act_mat_spd_dk = zeros(spd_num,trans_dk_num);
            
            for ii = 1:spd_num
                
                for jj = 1:trans_dk_num
                    
                    act_mat_spd_dk(ii,jj) = sum(dff_dk(jj,idx_spd(ii),:))>0;
                    
                end
                
            end
            
            act_consis_spd_dk0 = [act_consis_spd_dk;mean(act_mat_spd_dk,2);];
            act_consis_spd_dk = act_consis_spd_dk0;
        end
        
        % other cell trans in light
        if trans_lt_num>3
            
            act_mat_otr_lt = zeros(otr_num,trans_lt_num);
            
            for ii = 1:otr_num
                
                for jj = 1:trans_lt_num
                    
                    act_mat_otr_lt(ii,jj) = sum(dff_lt(jj,idx_otr(ii),:))>0;
                    
                end
                
            end
            
            act_consis_otr_lt0 = [act_consis_otr_lt;mean(act_mat_otr_lt,2);];
            act_consis_otr_lt = act_consis_otr_lt0;
            
        end
            
        % light cell trans in dark
        if trans_dk_num>3
            
            act_mat_otr_dk = zeros(otr_num,trans_dk_num);
            
            for ii = 1:otr_num
                
                for jj = 1:trans_dk_num
                    
                    act_mat_otr_dk(ii,jj) = sum(dff_dk(jj,idx_otr(ii),:))>0;
                    
                end
                
            end
            
            act_consis_otr_dk0 = [act_consis_otr_dk;mean(act_mat_otr_dk,2);];
            act_consis_otr_dk = act_consis_otr_dk0;
            
        end
         
    end
    
end

%% histogram of four types
figure('Unit','normalized','Position',[0.1 0.1 0.6 0.8]);
bins = 0:10:100;
rd = [255 0 81]/255;
bl = [0 128 166]/255;
subplot(2,2,1)
histogram(act_consis_spd_lt*100,bins,'Normalization','Prob','FaceColor',rd)
xlabel('Consistency (% of transitions)');
ylabel('Probability')
title('Speed cell - transition to light');

subplot(2,2,2)
histogram(act_consis_spd_dk*100,bins,'Normalization','Prob','FaceColor',bl)
xlabel('Consistency (% of transitions)');
ylabel('Probability')
title('Speed cell - transition to dark');

subplot(2,2,3)
histogram(act_consis_otr_lt*100,bins,'Normalization','Prob','FaceColor',rd)
xlabel('Consistency (% of transitions)');
ylabel('Probability')
title('Other cell - transition to light');

subplot(2,2,4)
histogram(act_consis_otr_dk*100,bins,'Normalization','Prob','FaceColor',bl)
xlabel('Consistency (% of transitions)');
ylabel('Probability')
title('Other cell - transition to dark');

%% 2. velocity peak during the transition

trans_name = {'elt_lt';'edrk_drk';};
               
mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];
win = 20;

for s = 2%:2
    
    act_consis_lt_lt = [];
    act_consis_lt_dk = [];
    
    act_consis_dk_lt = [];
    act_consis_dk_dk = [];
    
    for m = 8%:mnum
        
        dff   = ltdkdata.mouse(m).file(s).dFF;
        dffbin   = ltdkdata.mouse(m).file(s).dFF_bin;
        dff = dffnorm(dff,'max');
%         dff = dff.*dffbin;
        fn = size(dff,2);
        act_idx = find(~isnan(dff(:,1)));
        act_num = numel(act_idx);
        
        spd_idx = ltdkdata.mouse(m).file(s).spd_idx;
        lt_dk_idx = ltdkdata.mouse(m).file(s).lt_dk_idx;
        
        spd_idx = spd_idx>0;
        lt_idx = lt_dk_idx==1;
        dk_idx = lt_dk_idx==2;
        
        spd_cell_idx = find(spd_idx);
        idx_lt_only = find(lt_idx&(~spd_idx));
        idx_dk_only = find(dk_idx&(~spd_idx));
        
        lt_num = numel(idx_lt_only);
        dk_num = numel(idx_dk_only);
                
        if fn<f_lim(s)
            f_use = 1:fn;
        else
            f_use = 1:f_lim(s);
        end
        
        vel     = ltdkdata.mouse(m).file(s).behavior.velocity(f_use);
        trans   = ltdkdata.mouse(m).file(s).behavior.ttype(f_use);
        trans_lt = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{1}),trans,'UniformOutput',false)));
        trans_dk = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{2}),trans,'UniformOutput',false)));
        
        dff_lt = extract_from_onset_window(dff,trans_lt,win);
        dff_dk = extract_from_onset_window(dff,trans_dk,win);
        
        vel_lt = extract_from_onset_window(vel',trans_lt,win);
        vel_dk = extract_from_onset_window(vel',trans_dk,win);
        
        % cross-correlation
        trans_lt_num = numel(trans_lt);
        trans_dk_num = numel(trans_dk);
        
        
        if trans_lt_num>3
            
            % single transitions
            act_mat_lt_val = zeros(act_num,trans_lt_num);
            act_mat_lt_lag = zeros(act_num,trans_lt_num);
            
            for jj = 1:trans_lt_num
                
                vel_lt_jj = squeeze(vel_lt(jj,1,:));
                
                for ii = 1:act_num
                    
                    dff_lt_ij = squeeze(dff_lt(jj,act_idx(ii),:));
                    [rval,lag] = xcorr(vel_lt_jj,dff_lt_ij,'none');
                    
                    idx_rmax = find(rval==max(rval));
                    idx_rmax = idx_rmax((abs(idx_rmax-win-1) == min(abs(idx_rmax-win-1))));
                    act_mat_lt_val(ii,jj) = rval(idx_rmax);
                    act_mat_lt_lag(ii,jj) = lag(idx_rmax);
                    
                end
                
            end
            
            %% transition average
            vel_lt_avg = squeeze(nanmean(vel_lt,1));
            vel_lt_avg = (vel_lt_avg-min(vel_lt_avg))./(max(vel_lt_avg)-min(vel_lt_avg));
            
            dff_lt_avg = squeeze(nanmean(dff_lt,1));
            
            act_mat_lt_avg_val = zeros(act_num,1);
            act_mat_lt_avg_lag = zeros(act_num,1);
            figure
            cn = 6;
            for ii = 1:act_num
                
                dff_lt_avg_ii = dff_lt_avg(act_idx(ii),:);
                [rval,lag] = xcorr(vel_lt_avg(:),dff_lt_avg_ii(:),'none');
                
                idx_rmax = find(rval==max(rval));
                idx_rmax = idx_rmax((abs(idx_rmax-win-1) == min(abs(idx_rmax-win-1))));
                act_mat_lt_avg_val(ii) = rval(idx_rmax);
                act_mat_lt_avg_lag(ii) = lag(idx_rmax);
                
                if ii<=cn
                    subplot(3,cn,ii)
                    plot(dff_lt_avg_ii)
                    set(gca,'ylim',[-0.01 0.15],...
                            'xlim',[1 2*win+1],'xtick',[])
                    line([win+1 win+1],[-0.01 0.15],'linestyle','--','color','k')
                    line([0 2*win+1],[0 0],'linestyle','--','color','k')
%                     axis off
                    if ii==1
                        ylabel('\DeltaF/F');
                    end
                    title(['cell ',num2str(ii)])
                    
                    subplot(3,cn,ii+cn)
                    plot(vel_lt_avg,'r')
                    set(gca,'ylim',[-0.01 1],...
                            'xlim',[1 2*win+1],'xtick',[])
                    line([win+1 win+1],[-0.01 1],'linestyle','--','color','k')
                    line([0 2*win+1],[0 0],'linestyle','--','color','k')
                    if ii==1
                        ylabel('Velocity');
                    end
%                     axis off
                    subplot(3,cn,ii+2*cn)
                    plot(rval,'m');%idx_lt_only(jj)
                    set(gca,'ylim',[-0.1 1],...
                        'xlim',[1 4*win+3],'xtick',[])
                    
                    if ii==1
                        ylabel('Xcorr');
                    end
                    line([2*win+1 2*win+1],[0 1],'linestyle','--','color','k')
%                     axis off
                    line([1 4*win+1],[0 0],'linestyle','--','color','k')
                    line([1 4*win+1],[20 20],'linestyle','--','color','k')
                end
            end
            
        end
        
        t_thr = 10;
        idx_use = abs(act_mat_lt_lag)<t_thr;
        act_mat_lt_val(idx_use|(act_mat_lt_val<0)) = 0;
        
        [~,idx_sort] = sort(mean(act_mat_lt_val,2),'descend');
        
        cn = 6;
        trans_lt_num = 8;
        figure
        for ii = 1:trans_lt_num
            for jj = 1:cn
                subplot(cn+1,trans_lt_num,(jj-1)*trans_lt_num+ii)
                plot(squeeze(dff_lt(ii,act_idx(idx_sort(jj)),:)));%idx_lt_only(jj)
                set(gca,'ylim',[-0.2 1])
                line([win+1 win+1],[-0.2 1],'linestyle','--','color','k')
                axis off
            end
            
            subplot(cn+1,trans_lt_num,ii+trans_lt_num*cn)
            plot(squeeze(vel_lt(ii,1,:)),'r','linewidth',1)
            set(gca,'ylim',[-0.05 20])
            line([win+1 win+1],[-0.05 20],'linestyle','--','color','k')
            axis off
        end
        %
        figure
        rdg = [-30 200];
        for ii = 1:trans_lt_num
            
            vel_lt_ii = squeeze(vel_lt(ii,1,:));
            
            for jj = 1:cn
                
                dff_lt_ij = squeeze(dff_lt(ii,act_idx(idx_sort(jj)),:));
                [rval,lag] = xcorr(vel_lt_ii,dff_lt_ij,'none');
                
                subplot(cn,trans_lt_num,(jj-1)*trans_lt_num+ii)
                plot(rval);%idx_lt_only(jj)
                set(gca,'ylim',rdg)
                line([2*win+1 2*win+1],rdg,'linestyle','--','color','k')
                axis off
            end
            
        end

    end
    
end

%% Light cell trans to light example

S = act_mat_lt;
[n,k] = size(S);
h = nan*S;
mksz = 12;
rd = [255 0 81]/255;
figure('unit','points','position',[10 10 mksz*k mksz*n])
for i = 1:n
    for j = 1:k
        if S(i,j) ~= 0
           e = 'w'; f = S(i,j)*rd;
        else
           e = 'w'; f = [0.8 0.8 0.8];
        end
        h(i,j) = plot(j,i,'s','MarkerEdgeColor',e,'MarkerFaceColor',f,'MarkerSize',12); hold on
    end
end
set(gca,'xlim',[0 k+0.5],'ylim',[0 n+0.5],...
        'xtick',[1,5:5:k],'ytick',[1,5:5:n])
xlabel('Transition to light');
ylabel('Cell');

%%
t_thr = 10;
idx_use = abs(act_mat_lt_lag)<t_thr;
act_mat_lt_val(idx_use|(act_mat_lt_val<0)) = 0;

[act_mat_lt_val_sort,idx_sort] = sort(mean(act_mat_lt_val,2));
figure
subplot(1,2,1)
imagesc(act_mat_lt_val(idx_sort,:))
axis xy
xlabel('Transition')
ylabel('Cell');
caxis([0 100])
hc = colorbar;
ylabel(hc,'Cross correlation peak')
subplot(1,2,2)
plot(act_mat_lt_val_sort, 1:length(act_mat_lt_val_sort))
set(gca,'ylim',[1,length(act_mat_lt_val_sort)])
xlabel('Average cross-correlation peak value')

% histogram(act_mat_lt_lag(:))

%% 3. velocity peak during the transition

trans_name = {'elt_lt';'edrk_drk';};
               
mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];
win = 20;

for s = 2%:2
    
    act_consis_lt_lt = [];
    act_consis_lt_dk = [];
    
    act_consis_dk_lt = [];
    act_consis_dk_dk = [];
    
    rval_all_mice = [];
    for m = 1:mnum
        
        dff   = ltdkdata.mouse(m).file(s).dFF;
        dffbin   = ltdkdata.mouse(m).file(s).dFF_bin;
        dff = dffnorm(dff,'max');
%         dff = dff.*dffbin;
        fn = size(dff,2);
        act_idx = find(~isnan(dff(:,1)));
        act_num = numel(act_idx);
        
        spd_idx = ltdkdata.mouse(m).file(s).spd_idx;
        lt_dk_idx = ltdkdata.mouse(m).file(s).lt_dk_idx;
        
        spd_idx = spd_idx>0;
        lt_idx = lt_dk_idx==1;
        dk_idx = lt_dk_idx==2;
        
        spd_cell_idx = find(spd_idx);
        idx_lt_only = find(lt_idx&(~spd_idx));
        idx_dk_only = find(dk_idx&(~spd_idx));
        
        lt_num = numel(idx_lt_only);
        dk_num = numel(idx_dk_only);
                
        if fn<f_lim(s)
            f_use = 1:fn;
        else
            f_use = 1:f_lim(s);
        end
        
        vel     = ltdkdata.mouse(m).file(s).behavior.velocity(f_use);
        trans   = ltdkdata.mouse(m).file(s).behavior.ttype(f_use);
        trans_lt = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{1}),trans,'UniformOutput',false)));
        trans_dk = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{2}),trans,'UniformOutput',false)));
        
        dff_lt = extract_from_onset_window(dff,trans_lt,win);
        dff_dk = extract_from_onset_window(dff,trans_dk,win);
        
        vel_lt = extract_from_onset_window(vel',trans_lt,win);
        vel_dk = extract_from_onset_window(vel',trans_dk,win);
        
        % cross-correlation
        trans_lt_num = numel(trans_lt);
        trans_dk_num = numel(trans_dk);
        
        
        if trans_lt_num>3
            
            % transition average
            vel_lt_avg = squeeze(nanmean(vel_lt,1));
            vel_lt_avg = (vel_lt_avg-min(vel_lt_avg))./(max(vel_lt_avg)-min(vel_lt_avg));
            dff_lt_avg = squeeze(nanmean(dff_lt,1));
            
            act_mat_lt_avg_val = zeros(act_num,1);
            act_mat_lt_avg_lag = zeros(act_num,1);
%             act_mat_lt_avg_val_thr = zeros(act_num,1);
%             nreps = 1000;
            
            for ii = 1:act_num
                
                dff_lt_avg_ii = dff_lt_avg(act_idx(ii),:);
                [rval,lag] = xcorr(vel_lt_avg(:),dff_lt_avg_ii(:),'none');
                
                idx_rmax = find(rval==max(rval));
                idx_rmax = idx_rmax((abs(idx_rmax-win-1) == min(abs(idx_rmax-win-1))));
                act_mat_lt_avg_val(ii) = rval(idx_rmax);
                act_mat_lt_avg_lag(ii) = lag(idx_rmax);
                
%                 [rvaln,lagn] = xcorrrnd(vel_lt_avg(:),dff_lt_avg_ii(:),nreps);
%                 act_mat_lt_avg_val_thr(ii) = quantile(rvaln,0.975);
            end
%             act_mat_lt_avg_sig = act_mat_lt_avg_val>act_mat_lt_avg_val_thr;

            rval_all_mice0 = [rval_all_mice;act_mat_lt_avg_val];
            rval_all_mice  = rval_all_mice0;
            
        end

    end
    
end

%% 4. transition cell

trans_name = {'elt_lt';'edrk_drk';};
               
mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];
win = 20;

for s = 2%:2
    
    for m = 1:mnum
        
        dff      = ltdkdata.mouse(m).file(s).dFF;
        dffbin   = ltdkdata.mouse(m).file(s).dFF_bin;
        dff      = dffnorm(dff,'max');
        
%         dff = dff.*dffbin;
        [cn,fn] = size(dff);
        act_idx = ~isnan(dff(:,1));
        
               
        if fn<f_lim(s)
            f_use = 1:fn;
        else
            f_use = 1:f_lim(s);
        end
        
        vel     = ltdkdata.mouse(m).file(s).behavior.velocity(f_use);
        trans   = ltdkdata.mouse(m).file(s).behavior.ttype(f_use);
        
        % find the transition related activity
        trans_lt = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{1}),trans,'UniformOutput',false)));
        trans_dk = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{2}),trans,'UniformOutput',false)));
        
        dff_lt = extract_from_onset_window(dff,trans_lt,win);
        dff_dk = extract_from_onset_window(dff,trans_dk,win);
        
        vel_lt = extract_from_onset_window(vel',trans_lt,win);
        vel_dk = extract_from_onset_window(vel',trans_dk,win);
        
        trans_lt_num = numel(trans_lt);
        trans_dk_num = numel(trans_dk);
        
        % find the velocity peaks
        pk_thr = 8;
        [pks, locs] = findpeaks(vel,'MinPeakHeight',pk_thr);
        [~,idx_st]  = sort(pks,'descend');
        pks         = pks(idx_st(1:round((trans_lt_num+trans_dk_num)/2)));
        locs        = locs(idx_st(1:round((trans_lt_num+trans_dk_num)/2)));
        pk_num      = numel(locs);
        
        dff_pk = extract_from_onset_window(dff,locs,win);
        vel_pk = extract_from_onset_window(vel',locs,win);
        
        
        % calculate cross-correlation
        % light transition
        cell_xcorr_val_lt = zeros(cn,1);
        cell_xcorr_lag_lt = zeros(cn,1);
        
        if trans_lt_num>3
            
            % transition average
            vel_avg = squeeze(nanmean(vel_lt,1)); % change here
            dff_avg = squeeze(nanmean(dff_lt,1));
            
            vel_avg = (vel_avg-min(vel_avg))./(max(vel_avg)-min(vel_avg));
            
            for ii = 1:cn
                
                if act_idx(ii)
                    
                    dff_avg_ii = dff_avg(ii,:);
                    [rval,lag] = xcorr(vel_avg(:),dff_avg_ii(:),'none');
                    
                    idx_rmax = find(rval==max(rval));
                    idx_rmax = idx_rmax((abs(idx_rmax-win-1) == min(abs(idx_rmax-win-1))));
                    
                    cell_xcorr_lag_lt(ii) = lag(idx_rmax);
                    
                    if abs(cell_xcorr_lag_lt(ii))<(win-5)
                        cell_xcorr_val_lt(ii) = rval(idx_rmax);
                    end
                end
                
            end
            
        end
        
        % dark transition
        cell_xcorr_val_dk = zeros(cn,1);
        cell_xcorr_lag_dk = zeros(cn,1);
        
        if trans_dk_num>3
            
            % transition average
            vel_avg = squeeze(nanmean(vel_dk,1)); % change here
            dff_avg = squeeze(nanmean(dff_dk,1));
            
            vel_avg = (vel_avg-min(vel_avg))./(max(vel_avg)-min(vel_avg));
            
            for ii = 1:cn
                
                if act_idx(ii)
                    
                    dff_avg_ii = dff_avg(ii,:);
                    [rval,lag] = xcorr(vel_avg(:),dff_avg_ii(:),'none');
                    
                    idx_rmax = find(rval==max(rval));
                    idx_rmax = idx_rmax((abs(idx_rmax-win-1) == min(abs(idx_rmax-win-1))));
                    
                    cell_xcorr_lag_dk(ii) = lag(idx_rmax);
                    
                    if abs(cell_xcorr_lag_dk(ii))<(win-5)
                        cell_xcorr_val_dk(ii) = rval(idx_rmax);
                    end
                end
                
            end
            
        end
        
        % peak 
        cell_xcorr_val_pk = zeros(cn,1);
        cell_xcorr_lag_pk = zeros(cn,1);
        
        if trans_lt_num>3
            
            % transition average
            vel_avg = squeeze(nanmean(vel_pk,1)); % change here
            dff_avg = squeeze(nanmean(dff_pk,1));
            
            vel_avg = (vel_avg-min(vel_avg))./(max(vel_avg)-min(vel_avg));
            
            for ii = 1:cn
                
                if act_idx(ii)
                    
                    dff_avg_ii = dff_avg(ii,:);
                    [rval,lag] = xcorr(vel_avg(:),dff_avg_ii(:),'none');
                    
                    idx_rmax = find(rval==max(rval));
                    idx_rmax = idx_rmax((abs(idx_rmax-win-1) == min(abs(idx_rmax-win-1))));
                    
                    cell_xcorr_lag_pk(ii) = lag(idx_rmax);
                    
                    if abs(cell_xcorr_lag_pk(ii))<(win-5)
                        cell_xcorr_val_pk(ii) = rval(idx_rmax);
                    end
                    
                end
                
            end
            
        end

        % decide whether it is a transition cell 
        corr_thr = 1;
        idx_lt = cell_xcorr_val_lt > corr_thr;
        idx_dk = cell_xcorr_val_dk > corr_thr;
        idx_pk = cell_xcorr_val_pk > corr_thr;
        
        idx_lt_trans = idx_lt&(~idx_pk);
        idx_dk_trans = idx_dk&(~idx_pk);
        
        ltdkdata.mouse(m).file(s).idx_lt_trans = idx_lt_trans;
        ltdkdata.mouse(m).file(s).idx_dk_trans = idx_dk_trans;
        
        ltdkdata.mouse(m).file(s).idx_trans.idx_lt = idx_lt;
        ltdkdata.mouse(m).file(s).idx_trans.idx_dk = idx_dk;
        ltdkdata.mouse(m).file(s).idx_trans.idx_pk = idx_pk;
        
    end
    
end

%% numbers

lt_dk_num = zeros(mnum,3);
lttrans_spd_num = zeros(mnum,2);

cell_perct = zeros(mnum,12);

for m = 1:mnum
    
    dff      = ltdkdata.mouse(m).file(s).dFF;
%     act_idx = ~isnan(dff(:,1));
    cell_num = sum(~isnan(dff(:,1)));
    
    spd_idx = ltdkdata.mouse(m).file(s).spd_idx;
    lt_dk_idx = ltdkdata.mouse(m).file(s).lt_dk_idx;
    
    spd_idx = spd_idx>0;
    lt_idx = lt_dk_idx==1;
    dk_idx = lt_dk_idx==2;
    
    idx_lt_trans = ltdkdata.mouse(m).file(s).idx_lt_trans;
    idx_dk_trans = ltdkdata.mouse(m).file(s).idx_dk_trans;
    
    lttrans_spd_num(m,1) = 100*sum(lt_idx&idx_lt_trans)/sum(lt_idx);
    lttrans_spd_num(m,2) = 100*sum(dk_idx&idx_dk_trans)/sum(dk_idx);
    
    lt_dk_num(m,1) = 100*sum(idx_lt_trans)/cell_num;
    lt_dk_num(m,2) = 100*sum(idx_dk_trans)/cell_num;
    
    lt_dk_num(m,3) = 100*sum(idx_lt_trans&idx_dk_trans)/cell_num;
    
    cell_perct(m,1) = sum(idx_lt_trans&(lt_idx&(~spd_idx)))/cell_num*100;
    cell_perct(m,2) = sum((~idx_lt_trans)&(lt_idx&(~spd_idx)))/cell_num*100;
    
    cell_perct(m,3) = sum(idx_lt_trans&(lt_idx&spd_idx))/cell_num*100;
    cell_perct(m,4) = sum((~idx_lt_trans)&(lt_idx&spd_idx))/cell_num*100;
    
    cell_perct(m,5) = sum(idx_lt_trans&((~(lt_idx|dk_idx))&spd_idx))/cell_num*100;
    cell_perct(m,6) = sum((~idx_lt_trans)&((~(lt_idx|dk_idx))&spd_idx))/cell_num*100;
    
    cell_perct(m,7) = sum(idx_lt_trans&(dk_idx&spd_idx))/cell_num*100;
    cell_perct(m,8) = sum((~idx_lt_trans)&(dk_idx&spd_idx))/cell_num*100;
    
    cell_perct(m,9) = sum(idx_lt_trans&(dk_idx&(~spd_idx)))/cell_num*100;
    cell_perct(m,10) = sum((~idx_lt_trans)&(dk_idx&(~spd_idx)))/cell_num*100;
    
    cell_perct(m,11) = sum(idx_lt_trans&(~(lt_idx|dk_idx|spd_idx)))/cell_num*100;
    cell_perct(m,12) = 100-sum(cell_perct(m,1:11),2);
    
end

cell_perct_avg = mean(cell_perct,1)';

%% examples

trans_name = {'elt_lt';'edrk_drk';};
               
mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];
win = 20;

for s = 2%:2
    
    for m = 1%:mnum
        
        dff      = ltdkdata.mouse(m).file(s).dFF;
        dffbin   = ltdkdata.mouse(m).file(s).dFF_bin;
        dff      = dffnorm(dff,'max');
        
%         dff = dff.*dffbin;
        [cn,fn] = size(dff);
        act_idx = ~isnan(dff(:,1));
        
               
        if fn<f_lim(s)
            f_use = 1:fn;
        else
            f_use = 1:f_lim(s);
        end
        
        vel     = ltdkdata.mouse(m).file(s).behavior.velocity(f_use);
        trans   = ltdkdata.mouse(m).file(s).behavior.ttype(f_use);
        
        % find the transition related activity
        trans_lt = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{1}),trans,'UniformOutput',false)));
        trans_dk = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{2}),trans,'UniformOutput',false)));
        
        dff_lt = extract_from_onset_window(dff,trans_lt,win);
        dff_dk = extract_from_onset_window(dff,trans_dk,win);
        
        vel_lt = extract_from_onset_window(vel',trans_lt,win);
        vel_dk = extract_from_onset_window(vel',trans_dk,win);
        
        trans_lt_num = numel(trans_lt);
        trans_dk_num = numel(trans_dk);
        
        % find the velocity peaks
        pk_thr = 8;
        [pks, locs] = findpeaks(vel,'MinPeakHeight',pk_thr);
        [~,idx_st]  = sort(pks,'descend');
        pks         = pks(idx_st(1:round((trans_lt_num+trans_dk_num)/2)));
        locs        = locs(idx_st(1:round((trans_lt_num+trans_dk_num)/2)));
        pk_num      = numel(locs);
        
        dff_pk = extract_from_onset_window(dff,locs,win);
        vel_pk = extract_from_onset_window(vel',locs,win);
        
        
        idx_lt_trans = find(ltdkdata.mouse(m).file(s).idx_lt_trans);
        idx_dk_trans = find(ltdkdata.mouse(m).file(s).idx_dk_trans);
        
        % light trans
        vel_avg = squeeze(nanmean(vel_lt,1)); % change here
        vel_avg = (vel_avg-min(vel_avg))./(max(vel_avg)-min(vel_avg));
        figure
        cn = 6;
        for ii = 1:cn%numel(idx_lt_trans)
            dff_avg_ii = squeeze(nanmean(dff_lt(:,idx_lt_trans(ii),:),1));
            [rval,lag] = xcorr(vel_avg(:),dff_avg_ii(:),'none');
            
            subplot(3,cn,ii)
            plot(dff_avg_ii);hold on;
            axis square
            subplot(3,cn,ii+cn)
            plot(vel_avg);hold on;
            axis square
            subplot(3,cn,ii+cn*2)
            plot(rval); hold on;
            axis square
        end
        
%         vel_avg = squeeze(nanmean(vel_pk,1)); % change here
%         vel_avg = (vel_avg-min(vel_avg))./(max(vel_avg)-min(vel_avg));
%         
%         for ii = 1:cn%numel(idx_lt_trans)
%             dff_avg_ii = squeeze(nanmean(dff_pk(:,idx_lt_trans(ii),:),1));
%             [rval,lag] = xcorr(vel_avg(:),dff_avg_ii(:),'none');
%             
%             subplot(3,cn,ii)
%             plot(dff_avg_ii);hold on;
%             axis square
%             set(gca,'xlim',[1 length(dff_avg_ii)],'xtick',[],...
%                 'ylim',[-0.01 0.15])
%             line([win+1 win+1],get(gca,'ylim'),'linestyle','--','color','k')
%             
%             if ii==1
%                 ylabel('\DeltaF/F');
%                 legend('Light','Peak');
%             end
%             subplot(3,cn,ii+cn)
%             plot(vel_avg);hold on;
%             axis square
%             set(gca,'xlim',[1 length(vel_avg)],'xtick',[],...
%                 'ylim',[-0.1 1.1])
%             line([win+1 win+1],get(gca,'ylim'),'linestyle','--','color','k')
%             if ii==1
%                 ylabel('Velocity');
%                 legend('Light','Peak');
%             end
%             subplot(3,cn,ii+cn*2)
%             plot(rval); hold on;
%             axis square
%             set(gca,'xlim',[1 length(rval)],'xtick',[],...
%                 'ylim',[0 1.5])
%             line([2*win+1 2*win+1],get(gca,'ylim'),'linestyle','--','color','k')
%             if ii==1
%                 ylabel('Cross correlation');
%                 legend('Light','Peak');
%             end
%         end

        
        % dark trans
        vel_avg = squeeze(nanmean(vel_dk,1)); % change here
        vel_avg = (vel_avg-min(vel_avg))./(max(vel_avg)-min(vel_avg));

        cn = 6;
        for ii = 1:cn%numel(idx_lt_trans)
            dff_avg_ii = squeeze(nanmean(dff_dk(:,idx_dk_trans(ii),:),1));
            [rval,lag] = xcorr(vel_avg(:),dff_avg_ii(:),'none');
            
            subplot(3,cn,ii)
            plot(dff_avg_ii);hold on;
            axis square
            subplot(3,cn,ii+cn)
            plot(vel_avg);hold on;
            axis square
            subplot(3,cn,ii+cn*2)
            plot(rval); hold on;
            axis square
        end
        
        vel_avg = squeeze(nanmean(vel_pk,1)); % change here
        vel_avg = (vel_avg-min(vel_avg))./(max(vel_avg)-min(vel_avg));

        for ii = 1:cn%numel(idx_lt_trans)
            dff_avg_ii = squeeze(nanmean(dff_pk(:,idx_dk_trans(ii),:),1));
            [rval,lag] = xcorr(vel_avg(:),dff_avg_ii(:),'none');
            
            subplot(3,cn,ii)
            plot(dff_avg_ii);hold on;
            axis square
            set(gca,'xlim',[1 length(dff_avg_ii)],'xtick',[],...
                'ylim',[-0.01 0.15])
            line([win+1 win+1],get(gca,'ylim'),'linestyle','--','color','k')
            title(['Cell ', num2str(ii)])
            if ii==1
                ylabel('\DeltaF/F');
                legend('Light','Dark','Peak');
            end
            subplot(3,cn,ii+cn)
            plot(vel_avg);hold on;
            axis square
            set(gca,'xlim',[1 length(vel_avg)],'xtick',[],...
                'ylim',[-0.1 1.1])
            line([win+1 win+1],get(gca,'ylim'),'linestyle','--','color','k')
            if ii==1
                ylabel('Velocity');
                legend('Light','Dark','Peak');
            end
            subplot(3,cn,ii+cn*2)
            plot(rval); hold on;
            axis square
            set(gca,'xlim',[1 length(rval)],'xtick',[],...
                'ylim',[0 2])
            line([2*win+1 2*win+1],get(gca,'ylim'),'linestyle','--','color','k')
            line(get(gca,'xlim'),[1 1],'linestyle','--','color','k')
            if ii==1
                ylabel('Cross correlation');
                legend('Light','Dark','Peak');
            end
        end

        
    end
    
end

%%

figure
plot(vel);hold on
plot(locs, vel(locs),'v')
set(gca,'xlim',[0 18000],'xtick',0:2000:18000,'xticklabel',0:200:1800)
ylabel('Velocity (cm/s)');
xlabel('Time (s)');

%%

figure;
plot(idx_lt);hold on;
plot(idx_dk);
plot(idx_pk);



















