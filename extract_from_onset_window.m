function dff_onset = extract_from_onset_window(dff,onset,win)

[cn,fn] = size(dff);
osn = length(onset);

dff_onset = nan(osn,cn,2*win+1);

for n = 1:osn
    
    % for the onset
    A = onset(n)-win;
    B = onset(n)+win;
    
    tempmat = zeros(cn,2*win+1);
    
    if A<1
        A=1; tempmat(:,(end-length(A:B)+1):end) = dff(:,A:B);
    elseif B>fn
        B = fn;
        tempmat(:,1:length(A:B)) = dff(:,A:B);
    else
        tempmat = dff(:,A:B);
    end
    
    dff_onset(n,:,:) = tempmat;
 
end

end