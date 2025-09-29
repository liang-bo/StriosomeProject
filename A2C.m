function C = A2C(A)

[ih,iw,cn] = size(A);
C = zeros(cn,2);

for ii = 1:cn
    A_TEMP = squeeze(A(:,:,ii));
    [idx1,idx2] = ind2sub([ih,iw],find(A_TEMP>0));
    C(ii,1) = mean(idx1);
    C(ii,2) = mean(idx2);
end


end