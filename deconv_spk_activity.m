% deconv_spk_activity
% need caiman toolbox
sn = 3;
mnum = numel(ltdkdata.mouse);

for m = 1:mnum
    
    for f = 1:numel(ltdkdata.mouse(m).file)
        
        tic;
        dff   = ltdkdata.mouse(m).file(f).dFF;
        ltdkdata.mouse(m).file(f).spk = dff2spk(dff,sn);
        fprintf('Mouse #%d, session #%d was done!\n', m, f); toc;
        
    end
    
end

