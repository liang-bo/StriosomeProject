% Figure_velocity_peak_xcorr_dff

%% calculate the cross correlation between peak velocity and dff at four
% scenarios: lt_tr, dk_tr, lt, dk 

trans_name = {'elt_lt';'edrk_drk';};
state_name = {'lt';'drk';};

mnum = numel(ltdkdata.mouse);
f_lim = [9000;18000];

win = 20;        % window around transition +/- 2s
lagmax = 10;
se = strel('line',2*win,1);

rval = zeros(mnum,2,4); % average cross correlation value:
                        % dimention: mouse x session x behavior
                        % session: 1-open field; 2-light dark box
                        % behavior: 1-light transition;
                        %           2-dark transition;
                        %           3-in the light
                        %           4-in the dark
                        
vel_thr = 5; % velocity threshold for peaks

for s = 1%:2
    
    for m = 7%:mnum
        
        tic;
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
        
        % transition vector
        trans = ltdkdata.mouse(m).file(s).behavior.ttype(f_use);
        trans_lt_vect = transpose(cell2mat(cellfun(@(x) strcmp(x,trans_name{1}),trans,'UniformOutput',false)));
        trans_dk_vect = transpose(cell2mat(cellfun(@(x) strcmp(x,trans_name{2}),trans,'UniformOutput',false)));
        % extend the transition vector by a +/- win window
        trans_lt_vect = imdilate(trans_lt_vect,se);
        trans_dk_vect = imdilate(trans_dk_vect,se);
        
        % behavior state: light and dark
        state = ltdkdata.mouse(m).file(s).behavior.state(f_use);
        state_lt_vect = cell2mat(cellfun(@(x) strcmp(x,state_name{1}),state,'UniformOutput',false));
        state_dk_vect = cell2mat(cellfun(@(x) strcmp(x,state_name{2}),state,'UniformOutput',false));
        
        % velocity
        vel = transpose(ltdkdata.mouse(m).file(s).behavior.velocity(f_use));
        % find the peaks that above threshold vel_thr
        [pk,locs] = findpeaks(vel,'MinPeakHeight',vel_thr,'MinPeakDistance',2*win);
        vel_pk_vect = false(length(f_use),1);vel_pk_vect(locs) = true;
        
        % peaks belong to each behavior vectors
        trans_lt_pk = find(trans_lt_vect(:).*vel_pk_vect(:));
        trans_dk_pk = find(trans_dk_vect(:).*vel_pk_vect(:));
        state_lt_pk = find(state_lt_vect(:).*vel_pk_vect(:));
        state_dk_pk = find(state_dk_vect(:).*vel_pk_vect(:));
        % remove the overlapping peaks
        vel_pk_all = {  setdiff(trans_lt_pk,state_lt_pk);...
                        setdiff(trans_dk_pk,state_dk_pk);...
                        setdiff(state_lt_pk,trans_lt_pk);...
                        setdiff(state_dk_pk,trans_dk_pk);};
        
        % calculate the average max cross correlation value
        for t = 1:4
            
            vel_pk_t = vel_pk_all{t};
            dff_t  = extract_from_onset_window(dff,vel_pk_t,win);
            vel_t  = extract_from_onset_window(vel,vel_pk_t,win);
            
            [r_t,lag_t] = xcorr_trans(dff_t,vel_t,lagmax);
            rval(m,s,t) = nanmean(r_t(:));
            
        end
        
        fprintf('mouse %d, session %d is done. run time: %4.2f seconds. \n', m, s, toc);
    end
    
    
end


%% plot the transition 

s = 2;
mnum = 10;
f_use = 18000;
trans_name = {'elt_lt';'edrk_drk';};
nreps = 2000;
acc_thr = 0.4;

for m = 2%:mnum
    
    tic;
    dff = ltdkdata.mouse(m).file(s).dFF.*ltdkdata.mouse(m).file(s).dFF_bin;
    n_use = ~isnan(dff(:,1));
    dff = dff(n_use,:);
    dff = dffnorm(dff,'max');
    vel = ltdkdata.mouse(m).file(s).behavior.velocity;
    
%     trans = ltdkdata.mouse(m).file(s).behavior.ttype;
%     trans_lt = cell2mat(cellfun(@(x) strcmp(x,trans_name{1}),trans,'UniformOutput',false));
%     trans_dk = cell2mat(cellfun(@(x) strcmp(x,trans_name{2}),trans,'UniformOutput',false));
    
%     state = ltdkdata.mouse(m).file(s).behavior.state;
%     state_lt = cell2mat(cellfun(@(x) strcmp(x,'lt'),state,'UniformOutput',false))>0;
%     state_dk = cell2mat(cellfun(@(x) strcmp(x,'drk'),state,'UniformOutput',false))>0;
    
    if size(dff,2)>f_use
        dff = dff(:,1:f_use);
        vel = vel(1:f_use);
    end
    acc = diff([vel(1);vel(:);]);
    acc_p = acc>acc_thr;
    acc_n = acc<(-acc_thr);
    
    idx_acc = acc_p+acc_n*2;
    
    rng default % for reproducibility
    A = 1-squareform(pdist(dff,'cosine'));
    A(logical(eye(size(A))))=0;
    A((A<0)|(isnan(A))) = 0;
    k = full(sum(A));
    twom = sum(k);
    B = @(v) A(:,v) - k'*k(v)/twom;
    limit = 10000;
    [S,Q] = genlouvain(B,limit,0);
    Q = Q/twom;
    cmn = max(S);
    
    dff_comm = zeros(cmn,size(dff,2));
    auc_comm = nan(cmn,1);
    auc_p_r  = nan(cmn,1);
    auc_p_l  = nan(cmn,1);
    
    for ii = 1:cmn
        dff_comm(ii,:) = nanmean(dff(S==ii,:),1);
        auc_comm(ii) = roc_bout2(dff_comm(ii,:),idx_acc);
        roc_null = roc_bout2_rand(dff_comm(ii,:),idx_acc,nreps);
        
        auc_p_r(ii) = sum(roc_null>auc_comm(ii))/nreps;
        auc_p_l(ii) = sum(roc_null<auc_comm(ii))/nreps;
    end
    
    
    ltdkdata.mouse(m).file(s).community.S = S;
    ltdkdata.mouse(m).file(s).community.acceleration.auc_comm = auc_comm;
    ltdkdata.mouse(m).file(s).community.acceleration.auc_p_r = auc_p_r;
    ltdkdata.mouse(m).file(s).community.acceleration.auc_p_l = auc_p_l;
    
    
    fprintf('Mouse # %d was done @ %s ! \n',m,datestr(now)); toc;
    
end


