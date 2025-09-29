function bin_cell = bout_count(bin)

bn      = length(bin);
sign    = bin(1);
count   = 1;
bnum    = 1;
bin_cell = cell(1,1);

for ii = 2:bn
    
    
    if bin(ii) == sign
        
        if ii==bn
            count = count+1;
            bin_cell{bnum} = sign*ones(count,1);
        else
            count = count+1;
        end
        
    else
        
        bin_cell{bnum} = sign*ones(count,1);
        
        count = 1;
        bnum = bnum+1;
        
    end
    
    sign = bin(ii);
    
end

bin_cell = bin_cell(:);

end