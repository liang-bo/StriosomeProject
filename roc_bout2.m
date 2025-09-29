function auc = roc_bout2(A,B)
A = A(:);
B = B(:);
A_ON = A(B==1);
A_OF = A(B==2);

A_max = max(A);
bins = linspace(0, A_max,100);
prob_on = histcounts(A_ON,bins,'Normalization','cdf');hold on
prob_of = histcounts(A_OF,bins,'Normalization','cdf');
auc = (trapz(prob_on,prob_of)-0.5)*2;

end