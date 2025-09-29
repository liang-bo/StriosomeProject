function dff_onset = extract_from_onset_window0(dff,onset,win,flag)

[cn,fn] = size(dff);
osn = length(onset);
switch flag
    case 1
        dff_onset = nan(osn,cn,2*win+1);
    case 2
        dff_onset = nan(osn,cn,win+1);
    case 3
        dff_onset = nan(osn,cn,win+1);
    otherwise
        fprintf('something is wrong!');
end


for n = 1:osn
    
    % for the onset
    switch flag
        case 1
            A = onset(n)-win;
            B = onset(n)+win;
            tempmat = zeros(cn,2*win+1);
        case 2
            A = onset(n)-win;
            B = onset(n);
            tempmat = zeros(cn,win+1);
        case 3
            A = onset(n);
            B = onset(n)+win;
            tempmat = zeros(cn,win+1);
        otherwise
            fprintf('something is wrong!');
    end
 
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