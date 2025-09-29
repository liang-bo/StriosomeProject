
opengl('save', 'software') ;
load('C:\Users\hawessl\Documents\MATLAB\ltdkdata_100219.mat');

%% SECTION ONE
spk_thresh = 0.2; % percent of max spiking signal
smooth_filt_wid = 40; % Width of smoothing. 40=2sec?=max.
movement_thresh = 3; % Velocity threshold for being considered 'moving'
trans_window = 10; % # of frames in either direciton around a transition
mvmt_initiation_window = [-10 10];% # of frames around change from 0 to non0 velocity

nMice = numel(ltdkdata.mouse);

mouse_res = cell(nMice,1);
for xm = 1:nMice
    %idx_reg=1 for cells active both tasks, or 2 for cells only active in That task.
    det_neur_ld = find(ltdkdata.mouse(xm).file(2).idx_reg>0);
    NN = numel(det_neur_ld);
    ca = ltdkdata.mouse(xm).file(2).dFF(det_neur_ld,:).';
    spk = ltdkdata.mouse(xm).file(2).spk(det_neur_ld,:).';
    
    spk_max = max(spk,[],1); %max of col
    
    spk = spk./spk_max; % divide each column by max
    spk(spk<spk_thresh) = 0;
    
    smooth_filt = gausswin(smooth_filt_wid,3); % Gaussian filter to smooth spiking activity
    smooth_filt = smooth_filt(:)./sum(smooth_filt);
    
    spk = conv2(spk, smooth_filt, 'same'); % doing the actual smoothing
    
    ltdkdata.mouse(xm).file(2).filtspk = spk; % save threholded, filtered spike
    
    behavior = ltdkdata.mouse(xm).file(2).behavior;
    preds = table; % new table for predictors
    
    % Light vs. dark
%     preds.ld = categorical(behavior.state); % To use all 4 states
    preds.ld = behavior.state; % To use only light vs. dark
    preds.ld(strcmpi(preds.ld, 'elt')) = {'drk'};
    preds.ld(strcmpi(preds.ld, 'edrk')) = {'lt'};
    preds.ld = categorical(preds.ld);
    
    % Velocity
    preds.velocity = behavior.velocity; % Linear predictor
    
    % moving / still
    preds.moving = behavior.velocity >= movement_thresh; % Produces a logical variable, which will be treated as categorical
    
    % Transitions
    trans_list = behavior.ttype;
    for x = 1:numel(trans_list) % Get rid of "re-transitions"
        if strcmpi(trans_list(x),'edrk_drk')
           next_lt = x + find(strcmpi(behavior.state(x:end), 'lt'), 1, 'first'); % Find next time mouse is in light
           trans_list(x + find(strcmpi(trans_list(x+1:next_lt),'edrk_drk'))) = {'NA'}; % Nullify similar transitions before mouse is next in light
        end
        if strcmpi(trans_list(x),'elt_lt')
           next_lt = x + find(strcmpi(behavior.state(x:end), 'drk'), 1, 'first');
           trans_list(x + find(strcmpi(trans_list(x+1:next_lt),'elt_lt'))) = {'NA'};
        end
    end
    
    % Create windowed predictor variables for transitions
    lt_drk_trans = find(strcmpi(trans_list, 'edrk_drk'));
    drk_lt_trans = find(strcmpi(trans_list, 'elt_lt'));
    preds.trans_l2d = false(height(preds),1);
    preds.trans_d2l = false(height(preds),1);
    
    for x = 1:numel(lt_drk_trans) % Windowing around L->D transitions
        cur_range = (lt_drk_trans(x)-trans_window : lt_drk_trans(x)+trans_window);
        cur_range(cur_range < 1) = [];
        cur_range(cur_range > height(preds)) = [];
        preds.trans_l2d(cur_range) = true;
    end
    for x = 1:numel(drk_lt_trans) % Windowing around D->L transitions
        cur_range = (drk_lt_trans(x)-trans_window : drk_lt_trans(x)+trans_window);
        cur_range(cur_range < 1) = [];
        cur_range(cur_range > height(preds)) = [];
        preds.trans_d2l(cur_range) = true;
    end
    
    
    % Movement initiation
    preds.mvmt_init = false(height(preds),1);
    mvmt_init = find(diff(preds.moving)==1);
    for x = 1:numel(mvmt_init) % Windowing around still-->moving transitions
        cur_range = (mvmt_init(x)+mvmt_initiation_window(1) : mvmt_init(x)+mvmt_initiation_window(2));
        cur_range(cur_range < 1) = [];
        cur_range(cur_range > height(preds)) = [];
        preds.mvmt_init(cur_range) = true;
    end
    
    neur_results = cell(NN,1);
    for xn = 1:NN
        cpreds = preds;
        cpreds.hist = [0; spk(1:end-1,xn)];
        cpreds.Y = spk(:,xn);
        
        % No interactions, 'moving' Or else continuous 'velocity'
        cmdl = fitglm(cpreds, ['Y ~ 1 + ld + velocity + trans_l2d + trans_d2l + hist']); %*********_THE_GLM_*****
        
        coeff = cmdl.Coefficients;
        coeff.sig = coeff.tStat.*(coeff.pValue<(0.05/NN));   % *********************Bonf corrected P-VALUE
        coeff.Properties.RowNames = strrep(coeff.Properties.RowNames, '_1', '');
        coeff.Properties.RowNames = strrep(coeff.Properties.RowNames, 'ld_lt', 'in_light');
%         coeff.Properties.RowNames = regexprep(coeff.Properties.RowNames, '^moving', 'in_dark:moving');
        coeff(strcmpi(coeff.Properties.RowNames, '(Intercept)'),:) = [];
        coeff(strcmpi(coeff.Properties.RowNames, 'hist'),:) = [];
        coeff.Estimate = [];
        coeff.SE = [];
        mCoords = array2table(coeff.tStat.', 'VariableNames', strcat(coeff.Properties.RowNames, '_T'));
        P = array2table(coeff.pValue.', 'VariableNames', strcat(coeff.Properties.RowNames, '_P'));
        SIG = array2table(coeff.sig.', 'VariableNames', strcat(coeff.Properties.RowNames, '_sig'));
        neur_results{xn} = [table(xm, xn, 'VariableNames', {'mouse', 'neuron'}), mCoords, P, SIG];
        % use to determine best model ie model most acurately predicting
        % neural activity.. R2 always goes up with more terms though.
        % just compare models wiht same number of terms.
        neur_results{xn}.r2 = cmdl.Rsquared.Ordinary; 
        
    end
    mouse_res{xm} = cat(1, neur_results{:});
end
mouse_res = cat(1,mouse_res{:});




%% SECTION TWO
% re-write for whatever the predictor names are in 
% mouse_res (results) table

neuron_pattern = repmat('0',height(mouse_res),4); % **************************** fit to number of predictors in GLM****
neuron_pattern(mouse_res.in_light_sig>0,1) = '+';
neuron_pattern(mouse_res.in_light_sig<0,1) = '-';

%***************************************************************************************************
neuron_pattern(mouse_res.velocity_sig>0,2) = '+';%velocity OR moving ***
neuron_pattern(mouse_res.velocity_sig<0,2) = '-';
%  neuron_pattern(mouse_res.moving_sig>0,1) = '+';%velocity OR moving ***
%  neuron_pattern(mouse_res.moving_sig<0,1) = '-';
 
neuron_pattern(mouse_res.trans_l2d_sig>0,3) = '+';
neuron_pattern(mouse_res.trans_l2d_sig<0,3) = '-';
neuron_pattern(mouse_res.trans_d2l_sig>0,4) = '+';
neuron_pattern(mouse_res.trans_d2l_sig<0,4) = '-';
% neuron_pattern(mouse_res.mvmt_init_sig>0,4) = '+';
% neuron_pattern(mouse_res.mvmt_init_sig<0,4) = '-';

mouse_res.pattern = num2cell(neuron_pattern,2);

pat_tab = cell2table(tabulate(mouse_res.pattern), 'VariableNames', {'pattern', 'count', 'percent'});
pat_tab(pat_tab.percent<=1,:) = [];
pat_tab = sortrows(pat_tab, 'count', 'descend');
disp(pat_tab)


%% SECTION 2.1 curve-fitting to relate velocity to activity of individual neurons
% p = polyfit(x,y,n)  ... n=2 being quadratic

close all;

thsh=.001; %threshold for spike activity... unsure of units if any.
DEG=2; %polynomial degree to fit. 2 or 3 work well for many cells. 2 is quadratic 

%pattern order: inLt, Vel, LtD, DtL: 
patTable = mouse_res(strcmpi(mouse_res.pattern, '0+00'),:); 


nMice = numel(ltdkdata.mouse);
for i=1:height(patTable); %which mouse
 figure (i);
 wm=patTable.mouse(i);
 wc=patTable.neuron(i);
 
filtspk=ltdkdata.mouse(wm).file(2).filtspk;
beh=ltdkdata.mouse(wm).file(2).behavior;

idx=(filtspk(:,wc)>thsh);% index as long as frames of exper indicating when tested cell surpasses thrsh
keepSpk=filtspk(idx,wc); %for whichcell you're testing, keep spikes >thresh
keepVel=beh.velocity(idx);

idxBeh=beh.state(idx);
behLt=strcmpi(idxBeh,'lt'); %logical where being in light is true
behDk=strcmpi(idxBeh,'drk'); %logical where being in dark is true

LtSpks=keepSpk(behLt);%The suprathreshCa durring Light
DkSpks=keepSpk(behDk);%The suprathreshCa durring Dark
LtVel=keepVel(behLt);% accompanying Vel in Lt
DkVel=keepVel(behDk);% accompanying Vel in Dk

%  THIS Plots and fits frames in 
%  Light and in Dark separately for each cell:
    lc=[0 1 0];
    dc=[0 0 0];
        scatter(LtVel,LtSpks,10,lc);
        hold on; 
    polfit=polyfit(LtVel,LtSpks,DEG);% p = polyfit(x,y,n)
    polval = polyval(polfit,LtVel);%    y = polyval(p,x)
    plot(LtVel,polval,'--r');
        scatter(DkVel,DkSpks,10,dc);
        %hold on; 
    polfit=polyfit(DkVel,DkSpks,DEG);% p = polyfit(x,y,n)
    polval = polyval(polfit,DkVel);%    y = polyval(p,x)
    plot(DkVel,polval,'--b');
    hold off;
% end  

% %THIS will color-code by frames in light or dark, 
% %but fits all frames together:
% hold off;
%     for i=1:length(keepSpk);
%         c=[0 LDboolLt(i) 0];% BLACK=Dark; Green=Light
%         figure(wc);
%         scatter(keepVel(i),keepSpk(i),10,c);
%         hold on; 
%         polfit=polyfit(keepVel,keepSpk,DEG);% p = polyfit(x,y,n)
%         polval = polyval(polfit,keepVel);%    y = polyval(p,x)
%         plot(keepVel,polval,'--');
%     end

end % end for which mouse

%% SECTION 2.1   graph neural act vs vel per cell, while fitting quadratic.. by GLM-pattern

thsh=.001; %threshold for spike activity... unsure of units if any.
DEG=2; %polynomial degree to fit. 2 or 3 work well for many cells. 2 is quadratic 

for wc=1:3;%1:size(spk,2); % wc=which cell
figure(wc);
idx=(spk(:,wc)>thsh);
keepSpk=spk(idx,wc); %for whichcell you're testing, keep spikes >thresh
keepVel=preds.velocity(idx);
keepLD=preds.ld(idx);
LDboolLt=keepLD=='lt';%LightDarkboolean where Light=True=1
LDlt_idx=LDboolLt;
LDdk_idx=~LDboolLt;
LtSpks=keepSpk(LDlt_idx);%The suprathreshCa durring Light
DkSpks=keepSpk(LDdk_idx);%The suprathreshCa durring Dark
LtVel=keepVel(LDlt_idx);% accompanying Vel in Lt
DkVel=keepVel(LDdk_idx);% accompanying Vel in Dk
%  THIS Plots and fits frames in 
%  Light and in Dark separately for each cell:
    lc=[0 1 0];
    dc=[0 0 0];
        scatter(LtVel,LtSpks,10,lc);
        hold on; 
    polfit=polyfit(LtVel,LtSpks,DEG);% p = polyfit(x,y,n)
    polval = polyval(polfit,LtVel);%    y = polyval(p,x)
    plot(LtVel,polval,'--r');
        scatter(DkVel,DkSpks,10,dc);
        %hold on; 
    polfit=polyfit(DkVel,DkSpks,DEG);% p = polyfit(x,y,n)
    polval = polyval(polfit,DkVel);%    y = polyval(p,x)
    plot(DkVel,polval,'--b');
    hold off;
end

%% SECTION 3  -- use this to find and record coordinates of the edges of the maze.. the corners' x,y positions. To re-orrient to common axies. 
s=2;
for m=1:10
    x=ltdkdata.mouse(m).file(s).behavior.x;
    y=ltdkdata.mouse(m).file(s).behavior.y;
    state=ltdkdata.mouse(m).file(s).behavior.state;
    behavior=ltdkdata.mouse(m).file(s).behavior;
        endlen=length(x);
        seslen=floor(endlen/6); %frames in first session    
        dur=seslen;
        % BUT MOUSE 8 and MOUSE 10 each have only 3 sessions not 6!!!
    
    figure(m)
    subplot(2,2,1)%First of 6 sessions.
        X=behavior.x(1:dur);
        Y=behavior.y(1:dur);
        gscatter(X,Y,state(1:dur))
    subplot(2,2,2)%4th of 6 sessions
        P=behavior.x(4*dur:5*dur);
        Q=behavior.y(4*dur:5*dur);
        gscatter(P,Q,state(4*dur:5*dur))
    subplot(2,2,3)%6th of 6 sessions
        R=behavior.x(4*dur:6*dur);
        N=behavior.y(4*dur:6*dur);
        gscatter(R,N,state(4*dur:6*dur))
%      subplot(2,2,4)
%         P=behavior.x(4*dur:5*dur);
%         Q=behavior.y(4*dur:5*dur);
%         gscatter(P,Q,state(4*dur:5*dur))
end


%% REmap mouse xy coordinates to make compatible across mice and sessions -SECTION 4
% [remapped_mouse_coords] = knn_graph(ref_pixel_xy_coords, new_ref_coords, mouse_xy_coords)
mCoords=table();
mCoords.coordcell=cell(10,4); %order of 4 cols: coords first 3ses, coords last 3ses, target coords, flag for whether or not to flip L&D.
mCoords.coordcell{1,1}= [380,760;380,496;751,496];
mCoords.coordcell{2,1}= [367,651;376,392;732,392];
mCoords.coordcell{3,1}= [428,732;426,472;787,472];
mCoords.coordcell{4,1}= [397,689;386,427;751,427];
mCoords.coordcell{5,1}= [395,840;395,584;761,584];
mCoords.coordcell{6,1}= [376,648;376,390;736,390];
mCoords.coordcell{7,1}= [551,231;326,231;326,600];
mCoords.coordcell{8,1}= [1192,355;1192,616;829,616];
mCoords.coordcell{9,1}= [569,126;833,126;833,497];
mCoords.coordcell{10,1}= [1200,348;1200,618;816,618];

mCoords.coordcell{1,2}= [758,580;758,844;610,844];
mCoords.coordcell{2,2}= [762,592;762,850;392,850];
mCoords.coordcell{3,2}= [785,620;793,878;430,878];
mCoords.coordcell{4,2}= [710,593;710,851;417,851];
mCoords.coordcell{5,2}= [748,433;748,686;378,686];
mCoords.coordcell{6,2}= [777,589;777,850;398,850];
mCoords.coordcell{7,2}= [551,231;326,231;326,600];
mCoords.coordcell{8,2}= [1192,355;1192,616;829,616];
mCoords.coordcell{9,2}= [470,502;207,502;207,123];
mCoords.coordcell{10,2}= [1200,348;1200,618;816,618];

for i=1:10
mCoords.coordcell{i,3}=[0,1;0,0;1,0]; %new_ref_coords
mCoords.coordcell{i,4}= 1 %to flip or not to flip. 1 is flag to flip them.
end
mCoords.coordcell{7,4}= 0
mCoords.coordcell{8,4}= 0
mCoords.coordcell{10,4}= 0 %mouse 7,8,10 dont flip the light side
%%   REMAPPING -SECTION 5

s=2;
for m=1:10;
    
%       firstcnt=0;
%       for cnt=1:length(ltdkdata.mouse(m).file(s).behavior.snum);
%           if ltdkdata.mouse(m).file(s).behavior.snum(cnt)>=4;
%               firstcnt=firstcnt+1;
%           end
%       end
%       endFirsthalf=cnt-firstcnt
       
        new_ref_coords=mCoords.coordcell{m,3}; %[0,1;0,0;1,0];
        ref_pixel_xy_coordsA=mCoords.coordcell{m,1}; % first half of session.
        ref_pixel_xy_coordsB=mCoords.coordcell{m,2}; % last half of session.
        
      %sessions 1,2,3 of 6.
            idxA=ltdkdata.mouse(m).file(s).behavior.snum<=3;
            dataA=ltdkdata.mouse(m).file(s).behavior(idxA,:);%(1:endFirsthalf,:);
        mouse_xy_coordsA=horzcat(dataA.x,dataA.y);
        [remapped_mouse_coordsA] = knn_graph(ref_pixel_xy_coordsA, new_ref_coords, mouse_xy_coordsA);
        
       % sessions 4,5,6 of 6.
            idxB=ltdkdata.mouse(m).file(s).behavior.snum>=4;
            dataB=ltdkdata.mouse(m).file(s).behavior(idxB,:);%(endFirsthalf+1:end,:);
        mouse_xy_coordsB=horzcat(dataB.x,dataB.y);
        [remapped_mouse_coordsB] = knn_graph(ref_pixel_xy_coordsB, new_ref_coords, mouse_xy_coordsB);
            %all but mouse 7,8,10 need light side flipped in sessions 4,5,6.
            flipLightmice=[1,2,3,4,5,6,9];
            if ismember(m,flipLightmice)
                remapped_mouse_coordsB(:,1)=1-remapped_mouse_coordsB(:,1);%horizontally flip by subtracting x coords from 1.
            end
 
                 
    
    ltdkdata.mouse(m).file(s).behavior.newcoords=vertcat(remapped_mouse_coordsA,remapped_mouse_coordsB);
    %mouse 8 frames 15308:end are in new location/ not mapped.
    %mouse 10 frames 8975:end are in new location/ not mapped.
end


    %mouse 8 frames 15308:end are in new location/ not mapped.
    %mouse 10 frames 8975:end are in new location/ not mapped.
    ltdkdata.mouse(8).file(2).behavior.newcoords(15308:end,:)=NaN;
    ltdkdata.mouse(10).file(2).behavior.newcoords(8975:end,:)=NaN;
    
    
%%  SECTION 5.5 -- SHOWS REMAP NOT WORKING GREAT.. by flagging the points outside a 0-1 by 0-1 common goal axiz.
s=2;
for m=1:10
    indx=(ltdkdata.mouse(m).file(s).behavior.newcoords(:,1)>=1.10)|(ltdkdata.mouse(m).file(s).behavior.newcoords(:,1)<=-1.10);
    indy=(ltdkdata.mouse(m).file(s).behavior.newcoords(:,2)>=1.1)|(ltdkdata.mouse(m).file(s).behavior.newcoords(:,2)<=-1.10);

sum(indx)
sum(indy)

end


%% GRAPHS Poisiton during suprathresh activation for neurons fitting particular GLM pattern -- SECTION 6a - single cells

%*********************************Alex method of dff "spk" used in this
%section but not in section 6b below****************************

%  ALL maze coordinates Should be overlapping, with Light to one side and
%  Dark to the other.

s=2; %2=LDbox
cellcount=0;
L_actperframe=nan(300,10);%300 is more than numcells in any mouse I think. m= num mice.
D_actperframe=nan(300,10);
CorSpkToColstate=nan(300,10);
CorSpkToLightstate=nan(300,10);
CorSpkToDarkstate=nan(300,10);
for whichcell=1:(height(mouse_res))
    % 4 pattern entries :  light, vel, l2d, d2l
    if strcmp((mouse_res.pattern{whichcell}),'-000' )
            m=mouse_res.mouse(whichcell);  
            neurid=mouse_res.neuron(whichcell);
            % cellcount = cellcount+1;
            % disp ([m,neurid,cellcount])


            det_neur_ld = find(ltdkdata.mouse(m).file(2).idx_reg>0);
            NN = numel(det_neur_ld);
            ca = ltdkdata.mouse(m).file(2).dFF(det_neur_ld,:).';
            spk = ltdkdata.mouse(m).file(2).spk(det_neur_ld,:);

            spk_max = max(spk,[],1); %max of col

            spk = spk./spk_max; % divide each column by max
            spk(spk<spk_thresh) = 0;

            smooth_filt = gausswin(smooth_filt_wid,3); % Gaussian filter to smooth spiking activity
            smooth_filt = smooth_filt(:)./sum(smooth_filt);

            spk = conv2(spk, smooth_filt, 'same'); % doing the actual smoothing
            % *********************************Alex method of dff


        f_lim = [9000;18000];
        fn=size(spk,2); % *********************************Alex method of dff

        if fn<f_lim(s);
            len=fn;
        else
            len=f_lim(s);
        end;
        f_use=1:len; %frames to use

        idx_use = ~isnan(spk(:,1));% boolean column of 1s where ~isnan
        spk = spk(idx_use,f_use);% *********************************Alex method of dff

        state=ltdkdata.mouse(m).file(s).behavior.state;
        velocity=ltdkdata.mouse(m).file(s).behavior.velocity;
        behavior=ltdkdata.mouse(m).file(s).behavior;

        %thresh=(  nanmean(spk(neurid,:))+ (nanstd(spk(neurid,:)))  ); %threshold normalized per neuron  % was neur while neur=1:cn above 
        thresh=(max(spk(neurid,:))*0.2); % 20% of max spk is whats used in glm I think
        X=nan(len,1);%cn); %initiate mouse location matrix
            Y=nan(len,1);%cn);
            colstate=nan(len,1);
            fulldff=nan(len,1); %len is frames to use
            for i = 1:len %len is frames to use
                if spk(neurid,i)>thresh; %whereever neur dff breaches thresh at given frame, i     %was neur not neurid
                     X(i,neurid)=behavior.newcoords(i,1); %was behavior.x(i) ..record mouse position in location matrix     %was neur not neurid
                     Y(i,neurid)=behavior.newcoords(i,2); %was behavior.y(i)      %was neur not neurid
                     fulldff(i,neurid)=spk(neurid,i);% ***************spk is Alex method of dff
                     colstate(i)=strcmpi(state(i), 'lt');
                end
            end


        figure(whichcell);colormap('winter');
        scatter3(X(:,neurid),Y(:,neurid),fulldff(:,neurid),5,colstate==1);% Blue=0=false=dark, Yellow/Green=1=true=light
        ylim([-.2 1.2]);
        xlim([-.2 1.2]);
        hold on
        
%         q=((height(mouse_res))+(whichcell));
%         figure(q);colormap('winter');
%         scatter(fulldff(:,neurid),velocity(1:len,:),5,colstate==1);% Blue=0=false=dark, Yellow/Green=1=true=light
% %         ylim([-.2 1.2]);
% %         xlim([-.2 1.2]);
%         hold on

        Lightstate=colstate==1; 
        Darkstate=colstate==0; 
        % 
        % nansum(colstate)
        % nansum (Lightstate)
        % nansum(Darkstate)

        L_actperframe(neurid,m)=(nansum(spk(neurid,Lightstate)))/(nansum(Lightstate));
        D_actperframe(neurid,m)=(nansum(spk(neurid,Darkstate)))/(nansum(Darkstate));
        CorSpkToLightstate(neurid,m)=corr(spk(neurid,:)',Lightstate(:),'rows','complete');
        CorSpkToDarkstate(neurid,m)=corr(spk(neurid,:)',Darkstate(:),'rows','complete');

end

end
hold off

AvgactpframeinLt=nanmean(L_actperframe,1)
AvgactpframeinDk=nanmean(D_actperframe,1)
AvgCorSpkToLightstate=nanmean(CorSpkToLightstate,1)
AvgCorSpkToDarkstate=nanmean(CorSpkToDarkstate,1)



%% Alex wrote this% Identify +000 and -000 neurons in mouse_res
mouse_res.L_only = strcmpi(mouse_res.pattern, '+000');
mouse_res.D_only = strcmpi(mouse_res.pattern, '-000');
 
% Define smoothing filter
smooth_filt = gausswin(smooth_filt_wid,3);
smooth_filt = smooth_filt(:)./sum(smooth_filt);
 
% Get list of mice
mlist = unique(mouse_res.mouse);
 
% Initialize result variables for this calculation
mouse_res.LD_corr = nan(height(mouse_res),1);
mouse_res.LD_mean_diff = nan(height(mouse_res),1);
 
for x = 1:numel(mlist)
    curmouse = mlist(x);
    
    % Get logical index for rows of mouse_res corresponding to current
    % mouse
    mouseindex = mouse_res.mouse==curmouse;
    
    % Obtain spiking activity for this mouse
    det_neur_ld = find(ltdkdata.mouse(curmouse).file(2).idx_reg>0);
    spk = ltdkdata.mouse(curmouse).file(2).spk(det_neur_ld,:).';
    state=ltdkdata.mouse(curmouse).file(2).behavior.state;
    
    lnth=size(state,1);
    colstate=nan(lnth,1);
    for i =1:lnth
        colstate(i)=strcmpi(state(i), 'lt');
    end
    
    % Reorder columns of spk to match the order of rows for this mouse in
    % mouse_res
    spk = spk(:,mouse_res.neuron(mouseindex));
    
    % Process spiking activity data
    spk_max = max(spk,[],1); %max of col
    spk = spk./spk_max; % divide each column by max
    spk(spk<spk_thresh) = 0;
    spk = conv2(spk, smooth_filt, 'same'); % doing the actual smoothing
      
    
    % Calculate correlation and mean diff for L vs D
    mouse_res.LD_corr(mouseindex) = cellfun(@(X) corr(X, colstate(:), 'rows', 'complete'), num2cell(spk,1));
    mouse_res.LD_mean_diff(mouseindex) = cellfun(@(X) mean(X(colstate==1))-mean(X(colstate==0)), num2cell(spk,1));
end
 
 clf % Plot results for both values (mean diff and correlation)
subplot(1,2,1) % mean diff
[~,hbins1] = histcounts(mouse_res.LD_mean_diff,20);
histogram(mouse_res.LD_mean_diff(mouse_res.L_only), hbins1, 'Normalization', 'probability')% Light only in Blue
hold on
histogram(mouse_res.LD_mean_diff(mouse_res.D_only), hbins1, 'Normalization', 'probability')% Dark only in Red
 
subplot(1,2,2) % correlation
[~,hbins2] = histcounts(mouse_res.LD_corr,20);
histogram(mouse_res.LD_corr(mouse_res.L_only), hbins2,'Normalization', 'probability')% Light only in Blue
hold on
histogram(mouse_res.LD_corr(mouse_res.D_only), hbins2, 'Normalization', 'probability') % Dark only in Red


%% GRAPHS Poisiton during suprathresh activation for neurons fitting particular GLM patter -- SECTION 6b - combining those fitting a pattern
%  ALL maze coordinates Should be overlapping, with Light to one side and
%  Dark to the other.


% I need help identifying the sessions where I gave poor coordinates/ where
% they're not being remapped well.



s=2; %2=LDbox
cellcount=0;
for whichcell=1:(height(mouse_res))
    % 4 pattern entries :  light, vel, l2d, d2l
    if mouse_res.pattern{whichcell}=='+000'


m=mouse_res.mouse(whichcell);  
neurid=mouse_res.neuron(whichcell);
cellcount = cellcount+1;
disp ([m,neurid,cellcount])


f_lim = [9000;18000];
dff   = ltdkdata.mouse(m).file(s).dFF; 
fn = size(dff,2);

if fn<f_lim(s);
    len=fn;
else
    len=f_lim(s);
end;
f_use=1:len;

idx_use = ~isnan(dff(:,1));% boolean column of 1s where ~isnan
dff = dff(idx_use,f_use);  %dff to use
dff = dffnorm(dff,'max');  %normalized
%halflen=floor(len/2); %because I switched maze location after 3rd of 6 5min trials
state=ltdkdata.mouse(m).file(s).behavior.state;
velocity=ltdkdata.mouse(m).file(s).behavior.velocity;
behavior=ltdkdata.mouse(m).file(s).behavior;

thresh=(  mean(dff(neurid,:))+ (2*std(dff(neurid,:)))  ); %threshold normalized per neuron  % was neur while neur=1:cn above 
    X=nan(len,1);%cn); %initiate mouse location matrix
    Y=nan(len,1);%cn);
    colstate=nan(len,1);
    fulldff=nan(len,1);
    for i = 1:len; 
        if dff(neurid,i)>thresh; %whereever neur dff breaches thresh at given frame, i     %was neur not neurid
             X(i,neurid)=behavior.newcoords(i,1); %was behavior.x(i) ..record mouse position in location matrix     %was neur not neurid
             Y(i,neurid)=behavior.newcoords(i,2); %was behavior.y(i)      %was neur not neurid
             fulldff(i,neurid)=dff(neurid,i);     %was neur not neurid
             colstate(i)=strcmpi(state(i), 'lt');
        end
    end

figure(2);
scatter3(X(:,neurid),Y(:,neurid),fulldff(:,neurid),5,colstate)% Blue=0=false=dark, Yellow=1=true=light
ylim([-.2 1.2])
xlim([-.2 1.2])
hold on
end

end
hold off







%% SECTION 7
% PLOT CERTAIN NEURONS IN GRIN LENS FOV
%  in_light; Vel/Moving; L2D; D2L
load('C:\Users\hawessl\Documents\MATLAB\cellcoords.mat');
res_subtable = mouse_res(strcmpi(mouse_res.pattern, '00-0'),:);


newtable=res_subtable{:,{'mouse','neuron'}};
for i=1:length(newtable)
            mi=newtable(i,1); %mouse number
            ci=newtable(i,2); %cell number within mouse
    figure(mi);
     plot((cellcoords.mouse(mi).file(2).row_inds(ci)),(cellcoords.mouse(mi).file(2).col_inds(ci)),'bo','MarkerSize',7);
     xlim([0,400])
     ylim([0,400])
     hold on

           
end



%% SECTION 8  -uses pattern named in section 3, makes raster plot of those cells at transitions
close all;
ci=[];
mi=[];

spk_thresh = 0.2; % percent of max spiking signal
smooth_filt_wid = 40; % Width of smoothing
movement_thresh = 3; % Velocity threshold for being considered 'moving'
trans_window = 10; % # of frames in either direciton around a transition
mvmt_initiation_window = [-10 10];% # of frames around change from 0 to non0 velocity

nMice = numel(ltdkdata.mouse);
f_lim = [9000;18000];
newlt=[];
newdk=[];
newltci=[];
newdkci=[];
dffnrm_lt=[];
dffnrm_dk=[];

s=2;% LDbox
for i=1:length(newtable)
    mi=(newtable(i,1))
    ci=(newtable(i,2))
    
    det_neur_ld = find(ltdkdata.mouse(mi).file(s).idx_reg>0); %find present cells
    NN = numel(det_neur_ld);
    ca = ltdkdata.mouse(mi).file(s).dFF(det_neur_ld,:).';
    spk = ltdkdata.mouse(mi).file(s).spk(det_neur_ld,:).';
    
    spk_max = max(spk,[],1); %max of col
    
    spk = spk./spk_max; % divide each column by max
    spk(spk<spk_thresh) = 0;
    
    smooth_filt = gausswin(smooth_filt_wid,3); % Gaussian filter to smooth spiking activity
    smooth_filt = smooth_filt(:)./sum(smooth_filt);
    
    spk = conv2(spk, smooth_filt, 'same'); % doing the actual smoothing
    
                   %df = ltdkdata.mouse(mi).file(s).dFF;
                   fn = size(spk,1);

                   if fn<f_lim(s)  
                        f_use = 1:fn;
                   else
                        f_use = 1:f_lim(s); 
                   end
    spkuse=spk(f_use',:);
    trans_name = {'elt_lt';'edrk_drk';};
    within_name = {'elt_lt';'edrk_drk';};
    win = 20; % +/- 2s
    
    behavior = ltdkdata.mouse(mi).file(s).behavior;
     trans = behavior.ttype(f_use); 
     wthn = behavior.state(f_use);
     
                    trans_lt = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{1}),trans,'UniformOutput',false)));
                    trans_dk = find(cell2mat(cellfun(@(x) strcmp(x,trans_name{2}),trans,'UniformOutput',false)));
                    in_lt = find(cell2mat(cellfun(@(x) strcmp(x,within_name{1}),wthn,'UniformOutput',false)));
                    in_dk = find(cell2mat(cellfun(@(x) strcmp(x,within_name{2}),wthn,'UniformOutput',false)));
   
                    
                    
                    dffnrm_lt = extract_from_onset_window(spkuse',trans_lt,win);
                    dffnrm_ltavg = squeeze(mean(dffnrm_lt,1));%for each cell, avg over num of transitions made
                                     
                    dffnrm_dk = extract_from_onset_window(spkuse',trans_dk,win);
                    dffnrm_dkavg = squeeze(mean(dffnrm_dk,1));%for each cell, avg over num of transitions made

                    newlt=vertcat(newlt,dffnrm_ltavg);
                    newdk=vertcat(newdk,dffnrm_dkavg);
                    
                    
                    newltci=vertcat(newltci,dffnrm_ltavg(ci,:));
                    newdkci=vertcat(newdkci,dffnrm_dkavg(ci,:));
        
    
end

midnewlt=[];
midnewdk=[];
                    for i=1:(size(newltci,1))
                        midnewlt(i)=mean(newltci(i,10:30));
                    end 
                    [B,I]=sort(midnewlt,'descend');
                    newltci=newltci(I,:);
                    
                    for i=1:(size(newdkci,1))
                        midnewdk(i)=mean(newdkci(i,10:30));
                    end 
                    [C,J]=sort(midnewdk,'descend');
                    newdkci=newdkci(J,:);

                
figure(1);imagesc(newltci); % should be pattern-identified neuorons at transitions into Light
title('trans into Light');
figure(2);imagesc(newdkci);% should be pattern-identified neuorons at transitions into Dark
title('trans into Dark');
figure(3);



  
                    
                    







%%

[vel_r_pearson,vel_p_pearson] = corr(spk, preds.velocity, 'type', 'pearson');
[vel_r_spearman,vel_p_spearman] = corr(spk, preds.velocity, 'type', 'spearman'); % correlate Ca activity to velocity. n.s. overall.
[tl2d_r_pearson,tl2d_p_pearson] = corr(spk, preds.trans_l2d, 'type', 'pearson');

figure(1); histogram(vel_r_spearman)
figure(2); histogram(mouse_res.velocity_T)

figure(3); scatter(mouse_res.velocity_T, mouse_res.in_light_T)



%% Save all neurons' coordinates in cell map/ FOV
% 
% for m=1:10
%     for s=1:2
%         
%      A=ltdkdata.mouse(m).file(s).A;
%     [~,max_ind]=max(reshape(A,160000,[]),[],1); %reshapes to get out of 3D storage
%     [row_inds,col_inds]=ind2sub(size(A),max_ind);
%     cellcoords.mouse(m).file(s).row_inds=row_inds;
%     cellcoords.mouse(m).file(s).col_inds=col_inds;
%     % Save cellcoords.mat
% 
%     end
% end

%




%% Section 2.from SLH_MouseMaps..  plots calculations from sec1.SLHMouseMaps
%Plot cells active while mouse in more or less spatially restricted location

%cells active while mouse spatial locatiopn least wide-spread across maze
numlowmeand=sum(sortd<=(mean(sortd)-std(sortd)));% 
lowmeandidx=ltdkdata.mouse(m).file(s).meanD(1:numlowmeand,1);

%cells active while mouse spatial location most wide-spread across maze
fulllen=length(ltdkdata.mouse(m).file(s).meanD);
numhighmeand=sum(sortd>=(mean(sortd)+std(sortd)));% 
highmeandidx=ltdkdata.mouse(m).file(s).meanD((fulllen-numhighmeand):end,1);


%***********************
% *** CHANGE next two lines to either "low" or "high" versions ***
% low=low distance ie cell acts in more spatially compact area. high
% opposite.
%***********************
for idx=1:numlowmeand; % numlowmeand; OR numhighmeand; 
    neur=lowmeandidx(idx);%lowmeandidx(idx); OR highmeandidx(idx);
    
    
    figure; %one fig per neuron
    for i = 1:halflen; %frames first three sessions.. before moving maze. 
        if dff(neur,i)>thresh;
                 X=behavior.x(i);
                 Y=behavior.y(i);
                    if strcmp((state(i)),'lt');
                        scatter(X,Y,5,'b');hold on;
                    elseif strcmp((state(i)),'drk');
                        scatter(X,Y,5,'k');hold on;
                    elseif strcmp((state(i)),'elt');
                        scatter(X,Y,5,'r');hold on;
                    elseif strcmp((state(i)),'edrk');
                        scatter(X,Y,5,'g');hold on;
                    end
            
        end
    end
end





