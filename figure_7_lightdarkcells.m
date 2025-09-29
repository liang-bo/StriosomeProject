% figure_7_lightdarkcells

%% identify the cells
mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];

% figname = 'Figure_';
% figure('name', figname, 'unit','normalized', 'position',...
%     [0.1 0.1 0.4 0.6],'defaultAxesFontSize',18);
state_name = {'lt';'drk';};
stnum = numel(state_name);
mkprop = {'-ro';'-b^'};

for s = 1:2
    
%     vel_bin      = nan(mnum,binn-1);
%     gof_mm = [];
    
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

        state = ltdkdata.mouse(m).file(s).behavior.state(f_use);
        svar  = zeros(numel(state),numel(state_name));
        
        for ii = 1:stnum
            svar(:,ii) = cell2mat(cellfun(@(x) strcmp(x,state_name{ii}),state,'UniformOutput',false));
        end
        svar = svar>0;
        
        lt_dk_idx = zeros(cn,1);
        for n = 1:cn
            dff_lt = dff(n,svar(:,1));
            dff_dk = dff(n,svar(:,2));
            if ~isnan(dff_lt(1))
                [P,H,STATS] = ranksum(dff_lt(:),dff_dk(:),'alpha',0.05/cn);
                if H==1
                    if STATS.zval>0
                        lt_dk_idx(n) = 1;
                    else
                        lt_dk_idx(n) = 2;
                    end
                end
            end
        end
        ltdkdata.mouse(m).file(s).lt_dk_idx = lt_dk_idx;
        fprintf('mouse %d, session %d is done. \n', m, s);
        
    end
    
end

%% avg act during light and dark

mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];

% figname = 'Figure_';
% figure('name', figname, 'unit','normalized', 'position',...
%     [0.1 0.1 0.4 0.6],'defaultAxesFontSize',18);
state_name = {'lt';'drk';};
stnum = numel(state_name);
mkprop = {'-ro';'-b^'};

for s = 2%:2
    
    lt_cell_act_lt = [];
    lt_cell_act_dk = [];
    
    dk_cell_act_lt = [];
    dk_cell_act_dk = [];
    
    ot_cell_act_lt = [];
    ot_cell_act_dk = [];
    
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

        state = ltdkdata.mouse(m).file(s).behavior.state(f_use);
        svar  = zeros(numel(state),numel(state_name));
        
        for ii = 1:stnum
            svar(:,ii) = cell2mat(cellfun(@(x) strcmp(x,state_name{ii}),state,'UniformOutput',false));
        end
        svar = svar>0;
        
        lt_dk_idx = ltdkdata.mouse(m).file(s).lt_dk_idx;
        
        lt_cell_act_lt_0 = [lt_cell_act_lt; nanmean(dff(lt_dk_idx==1,svar(:,1)),2)];
        lt_cell_act_lt  = lt_cell_act_lt_0;
        
        lt_cell_act_dk_0 = [lt_cell_act_dk; nanmean(dff(lt_dk_idx==1,svar(:,2)),2)];
        lt_cell_act_dk  = lt_cell_act_dk_0;
        
        dk_cell_act_lt_0 = [dk_cell_act_lt; nanmean(dff(lt_dk_idx==2,svar(:,1)),2)];
        dk_cell_act_lt  = dk_cell_act_lt_0;
        
        dk_cell_act_dk_0 = [dk_cell_act_dk; nanmean(dff(lt_dk_idx==2,svar(:,2)),2)];
        dk_cell_act_dk  = dk_cell_act_dk_0;
        
        ot_cell_act_lt_0 = [ot_cell_act_lt; nanmean(dff(lt_dk_idx==0,svar(:,1)),2)];
        ot_cell_act_lt  = ot_cell_act_lt_0;
        
        ot_cell_act_dk_0 = [ot_cell_act_dk; nanmean(dff(lt_dk_idx==0,svar(:,2)),2)];
        ot_cell_act_dk  = ot_cell_act_dk_0;
        
        fprintf('mouse %d, session %d is done. \n', m, s);
        
    end
    
end

%% example cells

mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];

% figname = 'Figure_';
% figure('name', figname, 'unit','normalized', 'position',...
%     [0.1 0.1 0.4 0.6],'defaultAxesFontSize',18);
state_name = {'lt';'drk';};
stnum = numel(state_name);
mkprop = {'-ro';'-b^'};


for s = 2%:2
    
    cell_perct = zeros(mnum,3);
    
    for m = 1%:mnum
        
        dff   = ltdkdata.mouse(m).file(s).dFF;
        fn = size(dff,2);
        
        if fn<f_lim(s)
            f_use = 1:fn;
        else
            f_use = 1:f_lim(s);
        end
        
        dff = dff(:,f_use);
        dff = dffnorm(dff,'max');
        lt_dk_idx = ltdkdata.mouse(m).file(s).lt_dk_idx;
        state = ltdkdata.mouse(m).file(s).behavior.state(f_use);
        vel = ltdkdata.mouse(m).file(s).behavior.velocity(f_use);
        svar  = zeros(numel(state),numel(state_name));
        
        for ii = 1:stnum
            svar(:,ii) = cell2mat(cellfun(@(x) strcmp(x,state_name{ii}),state,'UniformOutput',false));
        end
        
        lt_cell_idx = find(lt_dk_idx==1);
        dk_cell_idx = find(lt_dk_idx==2);
        spr = 1;
        
        cnn = 10;
        figure%subplot(2,1,1)
        pp = 1;
        for ii = [1,4,6,8]%1:cnn%[5,7,8,10]%
            plot(dff(dk_cell_idx(ii),:)+(pp-1)*spr); hold on;
            pp = pp+1;
        end
%         plot(svar(f_use,1)*(cnn-1));
        
        [x1,y1] = bin2patch(svar(f_use,2));
        ylim = get(gca,'ylim');
        hp = patch('XData',x1,'YData',y1*ylim(2),...
            'facecolor','b',...
            'facealpha',0.2,...
            'edgecolor','none');
        uistack(hp,'bottom')
        set(gca,'ylim',[-0.5,4],'ytick',0:3,'yticklabel',1:4,...
                'xlim',[0 18000],'xtick',0:2000:18000,'xticklabel',0:200:1800)
        xlabel('Time (s)')
        ylabel('Cell');
%         plot(svar(f_use,2)*(cnn-1));
%         plot(vel/15+cnn)
%         figure%subplot(2,1,2)
%         for ii = 1:cnn
%             plot(dff(dk_cell_idx(ii),:)+(ii-1)*spr); hold on;
%         end
%         plot(svar(f_use,2)*(cnn-1));
        fprintf('mouse %d, session %d is done. \n', m, s);
        
    end
    
end

%% numbe of cells

mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];

% figname = 'Figure_';
% figure('name', figname, 'unit','normalized', 'position',...
%     [0.1 0.1 0.4 0.6],'defaultAxesFontSize',18);
state_name = {'lt';'drk';};
stnum = numel(state_name);
mkprop = {'-ro';'-b^'};


for s = 2%:2
    
    cell_perct = zeros(mnum,3);
    
    for m = 1:mnum
        
        dff   = ltdkdata.mouse(m).file(s).dFF(:,10);
        
        cell_num = sum(~isnan(dff));
        lt_dk_idx = ltdkdata.mouse(m).file(s).lt_dk_idx;
        
        cell_perct(m,1) = sum(lt_dk_idx==1)/cell_num*100;
        cell_perct(m,2) = sum(lt_dk_idx==2)/cell_num*100;
        cell_perct(m,3) = 100-sum(lt_dk_idx>0)/cell_num*100;
        fprintf('mouse %d, session %d is done. \n', m, s);
        
    end
    
end
