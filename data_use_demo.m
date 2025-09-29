%% load the data
load('CNMF_DATA_REG.mat')

%% Plot calcium traces

m = 2; % mouse number
f = 2; % 2- light dark test; 1- open field test

sp = 1.5;
N = 5;% smoothing window 500 ms
isplotsmooth = 1;% plot the smoothed traces

figure;
for ii = 1:2:20
    
    dff_ii = ltdkdata.mouse(m).file(f).dFF(ii,:);
    dff_ii_smooth = smooth(dff_ii,N);
    if isplotsmooth
        plot(dff_ii_smooth+(ii-1)*sp);
    else
        plot(dff_ii+(ii-1)*sp); 
    end
    hold on;
    
end

%% Raster plot 

figure;
subplot(2,1,1)
dff = ltdkdata.mouse(m).file(f).dFF;
imagesc(dff);colorbar;caxis([0 10])
title('all neurons')
subplot(2,1,2)
dff_nonnan = ~isnan(sum(dff,2));
imagesc(dff(dff_nonnan,:));colorbar;caxis([0 10])
title('active neurons only')

%% Plot the activity only on dark state
state_name = 'drk'; % the name you defined as dark in the state vector
state = ltdkdata.mouse(m).file(f).behavior.state;  % extract the behavior state
% trans = ltdkdata.mouse(m).file(f).behavior.trans; % similarly you can get the transition vector
idx_dk = cell2mat(cellfun(@(x) strcmp(x,state_name),state,'UniformOutput',false)); % compare the vecor with the name 'drk', 
                                                                                   % return a binary variable representing the state of dark
                                                                                   % 1- dark; 0- other states

figure;
subplot(2,1,1)
imagesc(dff(dff_nonnan,:));colorbar;caxis([0 10]);% plot active neurons only
axis xy; % reverse the y axis for visualization
hold on;
plot(idx_dk*80,'r') % superpose the dark state binary vector on the raster plot;
                    % scale to 80 for visualization 
title('all frames')
subplot(2,1,2)
dff_dk = dff(dff_nonnan,idx_dk); % calcium traces only at dark frames
imagesc(dff_dk);
colorbar;caxis([0 10])%color bar and color scale
axis xy;
title('dark state frames')