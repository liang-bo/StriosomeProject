function hln = plot_shadedline_group(hax,x,y,lnprop)

xx = nanmean(x,1);
yy = nanmean(y,1);
yy_err = std(y,0,1,'omitnan')/sqrt(size(y,1));
hln = shadedErrorBar(hax, xx, yy, yy_err, 'lineprops', lnprop);

end