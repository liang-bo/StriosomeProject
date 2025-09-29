function Y = shuffledmean(X,nperm)
% shuffle X along the second dimension
[n,m] = size(X);
Y = zeros(n,m);
X0 = zeros(n,m);

for ii = 1:nperm
    
    [~,index] = sort(rand(n,m),2);  
    
    for jj = 1:n; X0(jj,:) = X(jj,index(jj,:));end
    
    Y = Y + X0/nperm;
    
end


end