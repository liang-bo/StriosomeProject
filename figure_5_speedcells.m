%% load the data
load('lightdark_dataset_Apr-09-2020.mat');

%% the goal of this script is to identify cells whose activity has different 
% types of correlation with locomotion speed

%% 1. Examples
bins = [0:1:20];
binc = bins(2:end)-1;
binn = length(bins);

mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];

% figname = 'Figure_3_C';
% figure('name', figname, 'unit','normalized', 'position',...
%     [0.1 0.1 0.4 0.6],'defaultAxesFontSize',18);

mkprop = {'-ro';'-b^'};

for s = 2%:2
    
    dff_vbin_avg = nan(mnum,binn-1);
    vel_bin      = nan(mnum,binn-1);
    
    for m = 2%:mnum
        
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
        
        figure
        gof_all = zeros(20,1);
        for ii = 1:20
            
            subplot(4,5,ii)
            xx = bins(1:(end-1));xx = xx(:);
            yy = dff_vbin(ii,:);%nanmean(dff_vbin,1);
            yy = yy(:);
            yuse = ~isnan(yy);
            xx = xx(yuse);
            yy = yy(yuse);
            [f,gof]=fit(xx,yy,'poly2');
            gof_all(ii) = gof.adjrsquare;
            
            if gof_all(ii)>0.3
                plot(f,xx,yy,'r.')
            else
                plot(f,xx,yy,'k.')
            end
            title(['Cell ', num2str(ii)])
            xlabel('Speed')
            ylabel('\DeltaF/F')
            legend('off')       
        end
        %
        dff_vbin_avg(m,:) = nanmean(dff_vbin,1);
        fprintf('mouse %d, session %d is done. \n', m, s);
    end
    
%     dff_vbin_avg = zscore(dff_vbin_avg,0,2);

%     xx = nanmean(vel_bin,1);
%     yy = nanmean(dff_vbin_avg,1);
%     yy_err = std(dff_vbin_avg,0,1,'omitnan')/sqrt(size(dff_vbin_avg,1));
%     shadedErrorBar(gca,xx, yy, yy_err, 'lineprops', {mkprop{s},'MarkerFaceColor','w'})
%     xlabel('Avg velocity (cm/s)');
%     ylabel('Avg \DeltaF/F (normalized)');
%     set(gca,'ylim',[0 0.04],'xlim',[0 20])

end
%% 2. Activity vs velocity:
bins = [0:1:20];
binc = bins(2:end)-1;
binn = length(bins);

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
        
        dff_vbin = nan(cn,binn-1);
        
        A = 1;
        
        for b = 2:binn
            B = b;
            bin_idx = (vel>bins(A))&(vel<=bins(B));
            A = B;
            
            dff_vbin(:,b-1) = nanmean(dff(:,bin_idx),2);
            vel_bin(m,b-1)  = nanmean(vel(bin_idx));
            
        end
        
        gof_all = nan(cn,1);
        for ii = 1:cn
            
            xx = bins(1:(end-1));
            xx = xx(:);
            yy = dff_vbin(ii,:);%nanmean(dff_vbin,1);
            yy = yy(:);
            yuse = ~isnan(yy);
            if sum(yuse)>2
                xx = xx(yuse);
                yy = yy(yuse);
                [f,gof]=fit(xx,yy,'poly2');
                gof_all(ii) = gof.adjrsquare;
            end
        end
        
        gof_mm0 = [gof_mm;gof_all];
        gof_mm  = gof_mm0;
        
        fprintf('mouse %d, session %d is done. \n', m, s);
        
    end

end
figure;histogram(gof_mm,20);
%% plot the distribution
bins = [0:1:20];
binc = bins(2:end)-1;
binn = length(bins);

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
        
        dff_vbin = nan(cn,binn-1);
        
        A = 1;
        
        for b = 2:binn
            B = b;
            bin_idx = (vel>bins(A))&(vel<=bins(B));
            A = B;
            
            dff_vbin(:,b-1) = nanmean(dff(:,bin_idx),2);
            vel_bin(m,b-1)  = nanmean(vel(bin_idx));
            
        end
        
        gof_all = nan(cn,1);
        for ii = 1:cn
            
            xx = bins(1:(end-1));
            xx = xx(:);
            yy = dff_vbin(ii,:);%nanmean(dff_vbin,1);
            yy = yy(:);
            yuse = ~isnan(yy);
            if sum(yuse)>2
                xx = xx(yuse);
                yy = yy(yuse);
                [f,gof]=fit(xx,yy,'poly2');
                gof_all(ii) = gof.adjrsquare;
            end
        end
        
        gof_mm0 = [gof_mm;gof_all];
        gof_mm  = gof_mm0;
        
        fprintf('mouse %d, session %d is done. \n', m, s);
        
    end

end
%% distribution of speed

mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];
state_name = {'lt';'drk';};
stnum = numel(state_name);

for s = 2%:2
    
    vel_all = [];
    vel_all_lt = [];
    vel_all_dk = [];
    
    for m = 1:mnum
        
        vel = ltdkdata.mouse(m).file(s).behavior.velocity;
        fn = length(vel);
        if fn<f_lim(s)
            f_use = 1:fn;
        else
            f_use = 1:f_lim(s);
        end
        vel = vel(f_use);
        state = ltdkdata.mouse(m).file(s).behavior.state(f_use);
        svar  = zeros(numel(state),numel(state_name));
        
        for ii = 1:stnum
            svar(:,ii) = ii*cell2mat(cellfun(@(x) strcmp(x,state_name{ii}),state,'UniformOutput',false));
        end
        
        svar_idx = sum(svar,2);
        
        vel_all0 = [vel_all;vel];
        vel_all  = vel_all0;
        
        vel_all_lt0 = [vel_all_lt;vel(svar_idx==1)];
        vel_all_lt  = vel_all_lt0;
        
        vel_all_dk0 = [vel_all_dk;vel(svar_idx==2)];
        vel_all_dk  = vel_all_dk0;
        
        fprintf('mouse %d, session %d is done. \n', m, s);
        
    end

end
%%
figure;
bins = [0:1:18];
% subplot(1,3,1)
histogram(vel_all,bins,'Normalization','Prob'); hold on;
xlabel('Speed (cm/s)')
ylabel('Frame #')
title('Speed distribution')

% subplot(1,3,2)
histogram(vel_all_lt,bins,'Normalization','Prob')
xlabel('Speed (cm/s)')
ylabel('Frame #')
title('Speed distribution')

% subplot(1,3,3)
histogram(vel_all_dk,bins,'Normalization','Prob')
xlabel('Speed (cm/s)')
ylabel('Frame #')
title('Speed distribution')

%% 3. POSITIVE, NEGATIVE OR QUADRITIC

state_name = {'lt';'drk';};
stnum = numel(state_name);

bins = [0:1:18];
binc = bins(2:end)-1;
binn = length(bins);

mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];

methods = 'all';%{'all';'light';'dark'};

for s = 1%:2
    
    vel_bin      = nan(mnum,binn-1);
    
    for m = 1:mnum
        
        dff   = ltdkdata.mouse(m).file(s).dFF;
        dffbin   = ltdkdata.mouse(m).file(s).dFF_bin;
        dff = dff.*dffbin;
        [cn,fn] = size(dff);
        
        if fn<f_lim(s)
            f_use = 1:fn;
        else
            f_use = 1:f_lim(s);
        end
        
        dff = dff(:,f_use);
        dff = dffnorm(dff,'max');
        
        vel = ltdkdata.mouse(m).file(s).behavior.velocity(f_use);
        
        state = ltdkdata.mouse(m).file(s).behavior.state(f_use);
        svar  = zeros(numel(state),numel(state_name));
        for ii = 1:stnum
            svar(:,ii) = ii*cell2mat(cellfun(@(x) strcmp(x,state_name{ii}),state,'UniformOutput',false));
        end
        svar_idx = sum(svar,2);
        
        switch methods
            
            case 'all'
                
            case 'light'
                dff = dff(:,svar_idx==1);
                vel = vel(svar_idx==1);
            case 'dark'
                dff = dff(:,svar_idx==2);
                vel = vel(svar_idx==2);
        end
        
        dff_vbin = nan(cn,binn-1);
        
        A = 1;
        
        for b = 2:binn
            B = b;
            bin_idx = (vel>bins(A))&(vel<=bins(B));
            A = B;
            
            dff_vbin(:,b-1) = nanmean(dff(:,bin_idx),2);
            vel_bin(m,b-1)  = nanmean(vel(bin_idx));
            
        end
        
        spd_idx = zeros(cn,1);
        
        for ii = 1:cn
            
            xx = bins(1:(end-1));
            xx = xx(:);
            yy = dff_vbin(ii,:);%nanmean(dff_vbin,1);
            yy = yy(:);
            yuse = ~isnan(yy);
            if sum(yuse)>2
                xx = xx(yuse);
                yy = yy(yuse);
                
                [rval,pval] = corr(xx,yy,'Type','Spearman');
                
                if pval<0.05/cn
                    if rval>0
                        spd_idx(ii) = 1;
                    else
                        spd_idx(ii) = 2;
                    end
                else
                    [f,gof]=fit(xx,yy,'poly2');
                    gof_rsq = gof.adjrsquare;
                    
                    if gof_rsq>0.3
                        if f.p1>0
                            spd_idx(ii) = 3;
                        else
                            spd_idx(ii) = 4;
                        end
                    end
                    
                end
                
            end
        end
        
        ltdkdata.mouse(m).file(s).spd_idx = spd_idx;
        fprintf('mouse %d, session %d is done. \n', m, s);
        
    end

end

%% number of cells

mnum = numel(ltdkdata.mouse);

for s = 2%:2
    
    cell_perct = zeros(mnum,5);
    
    for m = 1:mnum
        
        dff   = ltdkdata.mouse(m).file(s).dFF(:,10);
        
        cell_num = sum(~isnan(dff));
        spd_idx = ltdkdata.mouse(m).file(s).spd_idx;
        
        cell_perct(m,1) = sum(spd_idx==1)/cell_num*100;
        cell_perct(m,2) = sum(spd_idx==2)/cell_num*100;
        cell_perct(m,3) = sum(spd_idx==3)/cell_num*100;
        cell_perct(m,4) = sum(spd_idx==4)/cell_num*100;
        
        cell_perct(m,5) = 100-sum(spd_idx>0)/cell_num*100;
        fprintf('mouse %d, session %d is done. \n', m, s);
        
    end
    
end


%% examples
state_name = {'lt';'drk';};
stnum = numel(state_name);
bins = [0:1:18];
binc = bins(2:end)-1;
binn = length(bins);

mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];

methods = 'all';%{'all';'light';'dark'};

for s = 2%:2
    
    vel_bin      = nan(mnum,binn-1);
    
    for m = 1%:mnum
        
        dff   = ltdkdata.mouse(m).file(s).dFF;
        dffbin   = ltdkdata.mouse(m).file(s).dFF_bin;
        dff = dff.*dffbin;
        [cn,fn] = size(dff);
        
        if fn<f_lim(s)
            f_use = 1:fn;
        else
            f_use = 1:f_lim(s);
        end
        
        dff = dff(:,f_use);
        dff = dffnorm(dff,'max');
        
        vel = ltdkdata.mouse(m).file(s).behavior.velocity(f_use);
        
        state = ltdkdata.mouse(m).file(s).behavior.state(f_use);
        svar  = zeros(numel(state),numel(state_name));
        for ii = 1:stnum
            svar(:,ii) = ii*cell2mat(cellfun(@(x) strcmp(x,state_name{ii}),state,'UniformOutput',false));
        end
        svar_idx = sum(svar,2);
        
        switch methods
            
            case 'all'
                
            case 'light'
                dff = dff(:,svar_idx==1);
                vel = vel(svar_idx==1);
            case 'dark'
                dff = dff(:,svar_idx==2);
                vel = vel(svar_idx==2);
        end
        
        dff_vbin = nan(cn,binn-1);
        
        A = 1;
        
        for b = 2:binn
            
            B = b;
            bin_idx = (vel>bins(A))&(vel<=bins(B));
            A = B;
            
            dff_vbin(:,b-1) = nanmean(dff(:,bin_idx),2);
            vel_bin(m,b-1)  = nanmean(vel(bin_idx));
            
        end
        
        spd_idx = ltdkdata.mouse(m).file(s).spd_idx;
        figure
        nn = 3;
        for ii = 1:4
            spd_use = find(spd_idx==ii);
            for jj = 1:nn
                subplot(5,nn,(ii-1)*nn+jj)
                xx = bins(1:(end-1));
                xx = xx(:);
                yy = dff_vbin(spd_use(jj),:);%nanmean(dff_vbin,1);
                yy = yy(:);
                if (ii==3)||(ii==4)
                    [f,gof]=fit(xx,yy,'poly2');
                    plot(f,xx,yy,'.')
                    legend off;
                else
                    [f,gof]=fit(xx,yy,'poly1');
                    plot(f,xx,yy,'.')
                    legend off;
                end
                xlabel('Speed bins (cm/s)')
                ylabel('Average activity')
            end
        end
        spd_use = find(spd_idx==0);
        
        for jj = 1:nn
            
            subplot(5,nn,(5-1)*nn+jj)
            xx = bins(1:(end-1));
            xx = xx(:);
            yy = dff_vbin(spd_use(jj),:);%nanmean(dff_vbin,1);
            yy = yy(:);
            [f,gof]=fit(xx,yy,'poly1');
            plot(f,xx,yy,'.')
            legend off;
            xlabel('Speed bins (cm/s)')
            ylabel('Average activity')
        end
        fprintf('mouse %d, session %d is done. \n', m, s);
        
    end

end


%% bin 0-10 cm/s pos or neg

state_name = {'lt';'drk';};
stnum = numel(state_name);

bins = [0:1:10];
binc = bins(2:end)-1;
binn = length(bins);

mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];

methods = 'all';%{'all';'light';'dark'};

for s = 2%:2
    
    vel_bin      = nan(mnum,binn-1);
    
    for m = 1:mnum
        
        dff   = ltdkdata.mouse(m).file(s).dFF;
        dffbin   = ltdkdata.mouse(m).file(s).dFF_bin;
        dff = dff.*dffbin;
        [cn,fn] = size(dff);
        
        if fn<f_lim(s)
            f_use = 1:fn;
        else
            f_use = 1:f_lim(s);
        end
        
        dff = dff(:,f_use);
        dff = dffnorm(dff,'max');
        
        vel = ltdkdata.mouse(m).file(s).behavior.velocity(f_use);
        
        state = ltdkdata.mouse(m).file(s).behavior.state(f_use);
        svar  = zeros(numel(state),numel(state_name));
        for ii = 1:stnum
            svar(:,ii) = ii*cell2mat(cellfun(@(x) strcmp(x,state_name{ii}),state,'UniformOutput',false));
        end
        svar_idx = sum(svar,2);
        
        switch methods
            
            case 'all'
                
            case 'light'
                dff = dff(:,svar_idx==1);
                vel = vel(svar_idx==1);
            case 'dark'
                dff = dff(:,svar_idx==2);
                vel = vel(svar_idx==2);
        end
        
        dff_vbin = nan(cn,binn-1);
        
        A = 1;
        
        for b = 2:binn
            B = b;
            bin_idx = (vel>bins(A))&(vel<=bins(B));
            A = B;
            
            dff_vbin(:,b-1) = nanmean(dff(:,bin_idx),2);
            vel_bin(m,b-1)  = nanmean(vel(bin_idx));
            
        end
        
        spd_idx = zeros(cn,1);
        
        for ii = 1:cn
            
            xx = bins(1:(end-1));
            xx = xx(:);
            yy = dff_vbin(ii,:);%nanmean(dff_vbin,1);
            yy = yy(:);
            yuse = ~isnan(yy);
            if sum(yuse)>2
                xx = xx(yuse);
                yy = yy(yuse);
                
                [rval,pval] = corr(xx,yy,'Type','Spearman');
                
                if pval<0.05/cn
                    if rval>0
                        spd_idx(ii) = 1;
                    else
                        spd_idx(ii) = 2;
                    end
                else
%                     [f,gof]=fit(xx,yy,'poly2');
%                     gof_rsq = gof.adjrsquare;
%                     
%                     if gof_rsq>0.3
%                         if f.p1>0
%                             spd_idx(ii) = 3;
%                         else
%                             spd_idx(ii) = 4;
%                         end
%                     end
                    
                end
                
            end
        end
        
        ltdkdata.mouse(m).file(s).spd_idx_10bin = spd_idx;
        fprintf('mouse %d, session %d is done. \n', m, s);
        
    end

end

%%


mnum = numel(ltdkdata.mouse);

for s = 2%:2
    
    cell_perct = zeros(mnum,3);
    
    for m = 1:mnum
        
        dff   = ltdkdata.mouse(m).file(s).dFF(:,10);
        
        cell_num = sum(~isnan(dff));
        spd_idx = ltdkdata.mouse(m).file(s).spd_idx_10bin;
        
        cell_perct(m,1) = sum(spd_idx==1)/cell_num*100;
        cell_perct(m,2) = sum(spd_idx==2)/cell_num*100;
        cell_perct(m,3) = 100-sum(spd_idx>0)/cell_num*100;
        fprintf('mouse %d, session %d is done. \n', m, s);
        
    end
    
end

cell_percent_avg = mean(cell_perct,1);

%% examples
state_name = {'lt';'drk';};
stnum = numel(state_name);
bins = [0:1:10];
binc = bins(2:end)-1;
binn = length(bins);

mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];

methods = 'all';%{'all';'light';'dark'};

for s = 2%:2
    
    vel_bin      = nan(mnum,binn-1);
    
    for m = 7%:mnum
        
        dff   = ltdkdata.mouse(m).file(s).dFF;
        dffbin   = ltdkdata.mouse(m).file(s).dFF_bin;
        dff = dff.*dffbin;
        [cn,fn] = size(dff);
        
        if fn<f_lim(s)
            f_use = 1:fn;
        else
            f_use = 1:f_lim(s);
        end
        
        dff = dff(:,f_use);
        dff = dffnorm(dff,'max');
        
        vel = ltdkdata.mouse(m).file(s).behavior.velocity(f_use);
        
        state = ltdkdata.mouse(m).file(s).behavior.state(f_use);
        svar  = zeros(numel(state),numel(state_name));
        for ii = 1:stnum
            svar(:,ii) = ii*cell2mat(cellfun(@(x) strcmp(x,state_name{ii}),state,'UniformOutput',false));
        end
        svar_idx = sum(svar,2);
        
        switch methods
            
            case 'all'
                
            case 'light'
                dff = dff(:,svar_idx==1);
                vel = vel(svar_idx==1);
            case 'dark'
                dff = dff(:,svar_idx==2);
                vel = vel(svar_idx==2);
        end
        
        dff_vbin = nan(cn,binn-1);
        
        A = 1;
        
        for b = 2:binn
            
            B = b;
            bin_idx = (vel>bins(A))&(vel<=bins(B));
            A = B;
            
            dff_vbin(:,b-1) = nanmean(dff(:,bin_idx),2);
            vel_bin(m,b-1)  = nanmean(vel(bin_idx));
            
        end
        
        spd_idx = ltdkdata.mouse(m).file(s).spd_idx_10bin;
        figure
        nn = 1;
        for ii = 1:2
            spd_use = find(spd_idx==ii);
            for jj = 1:nn
                subplot(2,nn,(ii-1)*nn+jj)
                xx = bins(1:(end-1));
                xx = xx(:);
                yy = dff_vbin(spd_use(jj),:);%nanmean(dff_vbin,1);
                yy = yy(:);
                if (ii==3)||(ii==4)
                    [f,gof]=fit(xx,yy,'poly2');
                    plot(f,xx,yy,'.')
                    legend off;
                else
                    [f,gof]=fit(xx,yy,'poly1');
                    plot(f,xx,yy,'.')
                    legend off;
                end
                xlabel('Speed bins (cm/s)')
                ylabel('Average activity')
            end
        end
        
        fprintf('mouse %d, session %d is done. \n', m, s);
        
    end

end

