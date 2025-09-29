function dff = dff_smooth(dff,N)

parfor n = 1:size(dff,1)   
   dff(n,:) = fastsmooth(dff(n,:),N);
end

end