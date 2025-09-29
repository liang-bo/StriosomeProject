%% figure_2_calcium_dynamics
% basic quantification of calcium signals

%% 1. Plot the example traces with spatial footprint superposed on a 
% correlation image
figname = 'Figure_2_A';
figure('name',figname,...
        'unit','normalized', 'position',...
        [0 0.1 1 0.5],'defaultAxesFontSize',14);

m = 3; % which mouse to show
s = 2;  % session number: 1- openfield 2-lightdark box
pix2um = 2.75; %micron per pixel
rangeij = [50 350];
thr = 0.4;

% plot the correlation image
subplot(1,3,1)
corr_img = double(ltdkdata.mouse(m).file(s).Cn);
corr_img = uint8((corr_img/max(corr_img(:)))*255);% normalize to 8-bit gray scale
imshow(corr_img,[5 150]);hold on;
axis square off; colorbar;
L = 30; H = 25; scale = 100;
line([L L+scale/pix2um],[H H],'color','w','linewidth', 2) % scale bar represent 100 micron

% superpose the cell footprint and plot the calcium traces
idx = 10:2:58;
if ~exist('cbrewer.m','file');addpath('cbrewer');end
clrmap =cbrewer('qual', 'Dark2', numel(idx));

A = ltdkdata.mouse(m).file(s).A;
plotfootprint(gca,A(:,:,idx), thr, clrmap)
set(gca,'xlim',rangeij,'ylim',rangeij)
title('Calcium image (correlation image)');

subplot(1,3,2:3)
dff = ltdkdata.mouse(m).file(s).dFF;
spr = 2;
fplot = 1:3000;
for ii = 1:numel(idx)
    plot(dff(idx(ii),fplot)+(ii-1)*spr, 'color',clrmap(ii,:));hold on;
end
xlabel('Frame number')
ylabel('Cell number')
set(gca,'ylim',[-0.5 ii*spr+2.5],...
        'ytick',0:spr:spr*(ii-1),...
        'yticklabel',1:ii)
title('Calcium trace');
saveas(gcf,figname)

%% 1. Plot the example traces with spatial footprint superposed on a 
% correlation image
figname = 'Figure_2_A';
figure('name',figname,...
        'unit','normalized', 'position',...
        [0 0.1 1 0.5],'defaultAxesFontSize',14);

m = 2; % which mouse to show
s = 2;  % session number: 1- openfield 2-lightdark box
pix2um = 2.75; %micron per pixel
rangeij = [50 350];
thr = 0.4;
fplot = 1:3000;

% plot the correlation image
subplot(2,3,[1,4])
corr_img = double(ltdkdata.mouse(m).file(s).Cn);
corr_img = uint8((corr_img/max(corr_img(:)))*255);% normalize to 8-bit gray scale
imshow(corr_img,[5 150]);hold on;
axis square off; colorbar;
L = 30; H = 25; scale = 100;
line([L L+scale/pix2um],[H H],'color','w','linewidth', 2) % scale bar represent 100 micron

id_use = ltdkdata.mouse(m).file(s).idx_reg>0;
dff = ltdkdata.mouse(m).file(s).dFF(id_use,:);
vel = ltdkdata.mouse(m).file(s).behavior.velocity;
dff_avg = nanmean(dff(:,fplot),2);
[~,sidx] = sort(dff_avg,'descend');

% idx = 10:2:58;
idx = sidx(1:15);
% superpose the cell footprint and plot the calcium traces
if ~exist('cbrewer.m','file');addpath('cbrewer');end
clrmap =cbrewer('qual', 'Dark2', numel(idx));

A = ltdkdata.mouse(m).file(s).A;
plotfootprint(gca,A(:,:,idx), thr, clrmap)
set(gca,'xlim',rangeij,'ylim',rangeij)
title('Calcium image (correlation image)');

subplot(2,3,2:3)

spr = 8;

for ii = 1:numel(idx)
    plot(dff(idx(ii),fplot)+(ii-1)*spr, 'color',clrmap(ii,:));hold on;
end
xlabel('Frame number')
ylabel('Cell number')
set(gca,'ylim',[-2 ii*spr+2.5],...
        'ytick',0:spr:spr*(ii-1),...
        'yticklabel',1:ii)
title('Calcium trace');

subplot(2,3,5:6)
plot(vel(fplot),'r')
xlabel('Frame number')
ylabel('Velocity (cm/s)')
% saveas(gcf,figname)

%% 2. Neuron numbers detected in each test and the percentage of commonly 
% detected neurons in each test
mnum = 10;
figname = 'Figure_2_B';
neu_num = zeros(mnum,2);
for m = 1:mnum   
    for s = 1:2        
       neu_num(m,s) = size(ltdkdata.mouse(m).file(s).A,3); 
    end
    
end

pval = ranksum(neu_num(:,1),neu_num(:,2)); % non-para t-test

% plot the results
neu_num_mean = mean(neu_num,1);
neu_num_se   = std(neu_num,0,1)/sqrt(size(neu_num,1));
xx = 1:length(neu_num_mean);

figure('name', figname, 'unit','normalized', 'position',...
    [0.1 0.1 0.25 0.8],'defaultAxesFontSize',18);

hb = bar(xx, neu_num_mean);hold on;
hb.EdgeColor = 'none';

he = errorbar(xx, neu_num_mean,neu_num_se,'CapSize',0, 'LineStyle','none','LineWidth',2);
he.YNegativeDelta = [];
he.Color   = hb.FaceColor;

set(gca,'xticklabel',{'Open field';'Light dark'})
xtickangle(45)
ylabel('Neuron number');
saveas(gcf,figname)

%% 3. The averaged calcium activities (calcium traces and inferred spikes) 
% under different behavioral states - light/dark or central/peripheral
figname = 'Figure_2_C';
figure('name', figname, 'unit','normalized', 'position',...
    [0.1 0.1 0.8 0.6],'defaultAxesFontSize',18);

bins_dff = linspace(0,0.05,20);
bins_spk = linspace(0,0.5,20);
f_use = 9000;% first 3 sessions

for s = 2:-1:1
    
   dff_all = [];
   spk_all = [];
   
   for m = 1:mnum
       
       dff = ltdkdata.mouse(m).file(s).dFF.*ltdkdata.mouse(m).file(s).dFF_bin;
       n_use = ~isnan(dff(:,1));
       dff = dff(n_use,:);
       dff = dffnorm(dff,'max');
       
       if size(dff,2)>f_use
           dff = dff(:,1:f_use);
       end
       
       dff_all0 = [dff_all;mean(dff,2)];
       dff_all = dff_all0;
       
       spk = ltdkdata.mouse(m).file(s).spk(n_use,:)>0;
       
       if size(spk,2)>f_use
           spk = spk(:,1:f_use);
       end
       
       spk_all0 = [spk_all;mean(spk,2)*10];
       spk_all = spk_all0;
   end
   fprintf('total neuron number = %d from %d mice \n', numel(spk_all), mnum)
   hax1 = subplot(1,2,1); histogram(dff_all,bins_dff);hold on;axis square;
   xlabel('        Avg calcium activity\newline (Normalized \DeltaF/F per second)',...
          'HorizontalAlignment','center');
   ylabel('Neuron number');
   hax2 = subplot(1,2,2); histogram(spk_all,bins_spk);hold on;axis square;
   xlabel('Avg inferred spikes\newline (count per second)');
   ylabel('Neuron number');
   
end

legend(hax1,{'Light dark';'Open field'})
legend(hax2,{'Light dark';'Open field'})
saveas(gcf,figname)
