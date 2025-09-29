function spk = dff2spk(dff,sn)

tau_d = 40; 
tau_r = 1; 
nMax = 20; 
pars = [tau_d, tau_r]; 
kernel = create_kernel('exp2', pars, nMax); 
kernel.bound_pars = false; 

[cn,fn] = size(dff);
spk = nan(cn,fn);

parfor ii = 1:cn
    [~,spk(ii,:),~,~] = deconvCa(dff(ii,:), kernel, sn, false, false);
end

end