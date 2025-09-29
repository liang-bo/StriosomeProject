function [du,onset,offset] = ca_boutduration2(bin)
% Calculate the behavior duration

[onset,offset] = be_findonset(bin);
onset = onset(:);offset = offset(:);

onset_num = length(onset);
du = zeros(onset_num,1);

for s = 1:onset_num
    du(s) = offset(s)-onset(s);
end

end