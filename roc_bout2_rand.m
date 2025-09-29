function roc_null = roc_bout2_rand(A,B,nrep)

roc_null = nan(nrep,1);

for ii = 1:nrep
    roc_null(ii) = roc_bout2(A,cell2mat(shuffle(bout_count(B))));
end

end