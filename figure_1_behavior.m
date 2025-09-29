% figure_1_behavior
% basic quantification of the behavior

%% 1. Plot the mouse position for both tests
pixsiz = 38/380; % cm/pix - ??? not 100% sure about this
f_use = 1:3000;
m = 1;

figname = 'Figure_1_A';
figure('name', figname, 'unit','normalized', 'position',...
    [0.1 0.1 0.8 0.6],'defaultAxesFontSize',20);

for s = 1:2
    subplot(1,2,s)
    x = ltdkdata.mouse(m).file(s).behavior.x(f_use);
    y = ltdkdata.mouse(m).file(s).behavior.y(f_use);
    [x,y] = positionsetzero(x,y);
    [vel, dist] = xy2vel(x,y,pixsiz);
    scatter(x,y,20,vel,'filled','MarkerFaceAlpha',0.75);
    hc = colorbar;caxis([0 25])
    ylabel(hc,'Velocity (cm/s)')
    set(gca,'xlim',[0 max(x)],'ylim',[0 max(y)])
    axis square off;
end

%% 2. mouse locomotion distance and time spent in different zone
figname = 'Figure_1_B';
pixsiz = 38/380; % cm/pix - ??? not 100% sure about this    
state_name = {'lt';'drk';'elt';'edrk'};
stnum = numel(state_name);
test_name = {'openfield';'lightdark'};

for s = 1:2
    
    time_state = zeros(mnum, stnum);
    dist_state = zeros(mnum, stnum);
    
    for m = 1:mnum
        
        x = ltdkdata.mouse(m).file(s).behavior.x;
        y = ltdkdata.mouse(m).file(s).behavior.y;
        
        [x,y] = positionsetzero(x,y);
        [vel, dist] = xy2vel(x,y,pixsiz);
        dist = smooth(1:length(dist),dist,10,'rloess'); % you may need this 
        
        % state
        state = ltdkdata.mouse(m).file(s).behavior.state;
        svar  = zeros(numel(state),numel(state_name));
        
        for ii = 1:stnum
            svar(:,ii) = cell2mat(cellfun(@(x) strcmp(x,state_name{ii}),state,'UniformOutput',false));
            dist_state(m,ii) = nanmean(dist(svar(:,ii)==1))*10*60; % cm per minitue
        end
        time_state(m,:) = mean(svar,1)*100;
        fprintf('%s, mouse %d is done! \n ', test_name{s}, m);
    end
    
    time_state_avg = mean(time_state,1)'; 
    
    % save to csv for prism
    file_dist = [figname,'_',test_name{s},'_distance.csv'];
    writematrix(dist_state,file_dist)
    
    file_time = [figname,'_',test_name{s},'_time.csv'];
    writematrix(time_state,file_time)
    
    file_time = [figname,'_',test_name{s},'_time_pie.csv'];
    writematrix(time_state_avg,file_time)
end

%% 3. Transitions - not sure about the transition
figname = 'Figure_1_C';
pixsiz = 38/380; % cm/pix - ??? not 100% sure about this    
trans_name = {'elt_lt';'edrk_drk';};
stnum = numel(trans_name);
test_name = {'openfield';'lightdark'};

for s = 1:2
    
    trans_freq = zeros(mnum, stnum);
%     dist_state = zeros(mnum, stnum);
    
    for m = 1:mnum
        
        
        % state
        trans = ltdkdata.mouse(m).file(s).behavior.ttype;
        svar  = zeros(numel(trans),numel(trans_name));
        
        for ii = 1:stnum
            svar(:,ii) = cell2mat(cellfun(@(x) strcmp(x,trans_name{ii}),trans,'UniformOutput',false));
%             dist_state(m,ii) = nanmean(dist(svar(:,ii)==1))*10*60; % cm per minitue
        end
        trans_freq(m,:) = mean(svar,1)*10*60; % count per min
        fprintf('%s, mouse %d is done! \n ', test_name{s}, m);
    end
    
    % save to csv for prism    
    file_time = [figname,'_',test_name{s},'_trans.csv'];
    writematrix(trans_freq,file_time)

end
