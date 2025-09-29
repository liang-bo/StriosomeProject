function plotfootprint(hax,A, thr, cmap)

windowSize = 3;
kernel = ones(windowSize) / windowSize ^ 2;

for jj = 1:size(A,3)
    A_temp = medfilt2(A(:,:,jj),[3,3]);
    A_temp(A_temp<thr*max(A_temp(:))) = 0;
    A_temp = conv2(single(A_temp), kernel, 'same');
    
    BW = bwareafilt(A_temp>0.95,1);
    BW2 = bwboundaries(BW);
    if ~isempty(BW2)
        for ii = 1:length(BW2)
            BW2{ii} = fliplr(BW2{ii});
            plot(hax, BW2{ii}(:,1),BW2{ii}(:,2),'Color',cmap(jj,:), 'linewidth', 1);
        end
    end
    hold on;
end

end