function dffn = dffnorm(dff,method)
% normalize the calcium traces

switch method
    
    case 'max'
        dffmax = max(dff,[],2);
        dffmax = repmat(dffmax,1,size(dff,2));
        dffn = dff./(dffmax+eps);
        
    case 'zscore'
        dffn = zscore(dff,0,2);
        
    otherwise
        fprintf('This is not an existing method!')
end

end