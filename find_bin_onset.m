function [onset,offset,duration] = find_bin_onset(bin)

set = bin(2:end)-bin(1:(end-1));

onset = find(set>0)+1;
offset = find(set<0);

if bin(1)==1
    onset = [1;onset];
end

onset_num = length(onset);
offset_num = length(offset);

if onset_num>offset_num
    offset = [offset;length(bin)];
end

duration = offset-onset;
end