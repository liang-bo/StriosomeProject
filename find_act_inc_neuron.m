function idx = find_act_inc_neuron(dff_avg)

idx = dff_avg > (nanmean(dff_avg)+std(dff_avg,'omitnan'));

end