function bloc = bin2bloc2(bin)

[du,onset,offset] = ca_boutduration(bin);

if (bin(1)+bin(end))== 1
    
    dnum = 2*numel(du);
    bloc = cell(dnum,1);
    
    if bin(1)==1
        
        onset = [onset;6001];
        
        for n = 1:dnum/2
            
            bloc{(n-1)*2+1} = true(du(n),1);
            bloc{(n-1)*2+2} = false(onset(n+1)-offset(n),1);
            
        end
        
    else
        
        offset = [1;offset];
        
        for n = 1:dnum/2
            
            bloc{(n-1)*2+1} = false(onset(n)-offset(n),1);
            bloc{(n-1)*2+2} = true(du(n),1);
            
        end
        
    end
    
elseif (bin(1)+bin(end))== 0
    
    dnum = 2*numel(du)+1;
    bloc = cell(dnum,1);
    
    offset = [1;offset];
    
    for n = 1:(dnum-1)/2
        
        bloc{(n-1)*2+1} = false(onset(n)-offset(n),1);
        bloc{(n-1)*2+2} = true(du(n),1);
        
    end
    
    bloc{dnum} = false(6001-offset(n+1),1);
    
else
    
    dnum = 2*numel(du)-1;
    bloc = cell(dnum,1);
    
    %     onset = [onset;6000];
    for n = 1:((dnum-1)/2)
        
        bloc{(n-1)*2+1} = true(du(n),1);
        bloc{(n-1)*2+2} = false(onset(n+1)-offset(n),1);
        
    end
    
    bloc{dnum} = true(6001-onset(n+1),1);
    
end

end