function dsp3_plot_psd_new(varargin)

defaults = dsp3.get_common_plot_defaults( dsp3.get_common_make_defaults() );
defaults.is_choice = true;
defaults.subset_good_days = true;

params = dsp3.parsestruct( defaults, varargin );

conf = params.config;

is_choice = params.is_choice;
% 
trial_type_label = ternary( is_choice, 'choice', 'cued' );
% event_name = ternary( is_choice, 'targAcq-150-cc', 'targOn-150-cc' );

event_names = { 'targAcq-150-cc', 'targOn-150-cc' };

psd_p = dsp3.get_intermediate_dir( 'original_summarized_psd', conf );
full_psd_p = shared_utils.io.fullfiles( psd_p, event_names);
psd_mats = shared_utils.io.findmat( full_psd_p );

[psd, psd_labs, freqs, t] = bfw.load_time_frequency_measure( psd_mats ...
  , 'get_labels_func', @(x) x.labels ...
);

choice_mask = find( psd_labs, {'targAcq', 'choice'} );
cued_mask = find( psd_labs, {'targOn', 'cued'} );
mask = union( choice_mask, cued_mask );
mask = findnone( psd_labs, 'errors', mask );

%%

if ( params.subset_good_days )
  mask = fcat.mask( psd_labs, mask ...
    , @findnone, dsp2.process.format.get_bad_days() ...
    , @findnone, 'day__05272017' ...
  );
end

% use_mask = fcat.mask( psd_labs, choice_mask ...
%   , @findnone, dsp2.process.format.get_bad_days() ...
%   , @findnone, 'day__05272017' ...
% );
% 
% regs = findall( psd_labs, 'regions', use_mask );
% n = cellfun( @(x) numel(findall(psd_labs, {'channels','days','sites'}, x)), regs );
%%

psd = rowref( psd, mask );
keep( psd_labs, mask );

prune( dsp3.add_context_labels(psd_labs) );

setcat( psd_labs, 'contexts', 'context__received', find(psd_labs, {'self', 'both', 'targOn', 'cued'}) );
setcat( psd_labs, 'contexts', 'context__forgone', find(psd_labs, {'other', 'none', 'targOn', 'cued'}) );
prune( psd_labs );

%%

proanti_each = { 'days', 'sites', 'channels', 'regions', 'trialtypes' };
[proanti, proanti_labs] = dsp3.pro_v_anti( psd, psd_labs', proanti_each );

%%

f_ind = mask_gele( freqs, 10, 80 );
t_ind = mask_gele( t, -500, 500 );

plt_freqs = freqs(f_ind);
plt_t = t(t_ind);

pl = plotlabeled.make_spectrogram( round(plt_freqs), plt_t );

mask = fcat.mask( proanti_labs ...
  , @find, {trial_type_label, 'acc'} ...
);

pltdat = proanti(mask, f_ind, t_ind);
pltdat(isinf(pltdat)) = nan;

pltlabs = prune( proanti_labs(mask) );

nan_rows = any( any(isnan(pltdat), 3), 2 );
pltdat = pltdat(~nan_rows, :, :);
keep( pltlabs, find(~nan_rows) );

axs = pl.imagesc( pltdat, pltlabs, {'regions', 'outcomes', 'trialtypes'} );

shared_utils.plot.tseries_xticks( axs, plt_t );
shared_utils.plot.fseries_yticks( axs, round(flip(plt_freqs)), 5 );

% shared_utils.plot.set_clims( axs, [-2.6, 2.1] );
shared_utils.plot.hold( axs, 'on' );
shared_utils.plot.add_vertical_lines( axs, find(plt_t == 0) ); %#ok

%%

save_p = char( dsp3.plotp({'psd_lines', dsp3.datedir}) );
do_save = params.do_save;

is_pro_anti = true;
is_pro_minus_anti = true;
per_outcome = false;

bands = dsp3.get_bands( 'map' );

[band_dat, band_labs] = dsp3.get_band_means( psd, psd_labs', freqs, bands );

band_dat(isinf(band_dat)) = nan;
nan_rows = any( any(isnan(band_dat), 3), 2 );
band_dat = band_dat(~nan_rows, :, :);
keep( band_labs, find(~nan_rows) );

proanti_each = { 'days', 'sites', 'channels', 'regions', 'trialtypes', 'bands' };
[proanti, proanti_labs] = dsp3.pro_v_anti( band_dat, band_labs', proanti_each );

if ( is_pro_minus_anti )
  [proanti, proanti_labs] = dsp3.pro_minus_anti( proanti, proanti_labs', proanti_each );
end

fcats = { 'contexts', 'bands', 'trialtypes' };

if ( is_pro_minus_anti )
  gcats = { 'outcomes', 'regions' };
  pcats = { 'trialtypes', 'bands' };
else
  gcats = { 'outcomes', 'trialtypes' };
  pcats = { 'bands', 'regions' };
end

if ( is_pro_anti )
  pltdat = proanti;
  pltlabs = proanti_labs';
else
  pltdat = band_dat;
  pltlabs = band_labs';
end

if ( ~per_outcome )
  prune( pltlabs );
  collapsecat( pltlabs, {'contexts', 'outcomes'} );
end

t_ind = mask_gele( t, -350, 350 );

pl = plotlabeled.make_common();
pl.x = t(t_ind);
pl.add_errors = true;
pl.smooth_func = @(x) smooth(x, 3);
pl.add_smoothing = true;
% pl.y_lims = [ -2.5e-7, 2.5e-7 ];

mask = fcat.mask( pltlabs...
  , @find, {'new_gamma', 'beta'} ...
  , @findnone, 'errors' ...
);

fig_I = findall_or_one( pltlabs, fcats, mask );

for i = 1:numel(fig_I)
  pltdat_ = pltdat(fig_I{i}, t_ind);
  pltlabs_ = prune( pltlabs(fig_I{i}) );
  
  [tot_ns, ns, n_I, n_C] = figure_s4_n_calculation( pltlabs_' );

  axs = pl.lines( pltdat_, pltlabs_, gcats, pcats );

  if ( do_save )
    shared_utils.plot.fullscreen( gcf );
    dsp3.req_savefig( gcf, save_p, pltlabs_, [pcats, fcats] );
  end
end

%%

pre0 = t >= -300 & t < 0;
post0 = t >= 0 & t <= 300;

pre_data = nanmean( pltdat(mask, pre0), 2 );
post_data = nanmean( pltdat(mask, post0), 2 );
period_labels = prune( pltlabs(mask) );

addcat( period_labels, 'period' );
repset( period_labels, 'period', {'before', 'after'} );

periods = [ pre_data; post_data ];

figs = dsp3.multi_plot( @bar, periods, period_labels ...
  , 'bands', {'contexts', 'bands'}, {'period', 'trialtypes'}, {'outcomes', 'bands', 'regions'} ...
);
%%

save_p = char( dsp3.plotp({'psd_hists', dsp3.datedir}) );
remove_outliers = false;
n_devs = 3;
prefix = ternary( remove_outliers, sprintf('%d_outliers_removed', n_devs), 'all_data' );
do_save = params.do_save;

mask = fcat.mask( pltlabs...
  , @find, {'gamma', 'beta'} ...
  , @findnone, 'errors' ...
);

site_spec = setdiff( dsp3_ct.site_specificity(), {'unit_uuid', 'outcomes'} );
site_spec = union( site_spec, 'bands' );
[site_labs, site_I] = keepeach( pltlabs', site_spec, mask );
site_means = bfw.row_nanmean( pltdat, site_I );

peaks = max( site_means, [], 2 );

if ( remove_outliers )
  check_I = findall( site_labs, {'trialtypes', 'bands', 'regions'} );
  check_means = bfw.row_nanmean( peaks, check_I );
  check_devs = rowop( peaks, check_I, @(x) nanstd(x, [], 1) );

  for i = 1:numel(check_I)
    subset_peaks = peaks(check_I{i});
    oob = subset_peaks < (check_means(i) - check_devs(i)*n_devs) | ...
      subset_peaks > (check_means(i) + check_devs(i)*n_devs);
    check_I{i}(oob) = [];
  end

  to_keep = vertcat( check_I{:} );
else
  to_keep = rowmask( peaks );
end

[figs, axs, labs] = dsp3.multi_plot( @hist, peaks, site_labs ...
  , {'bands'}, {'trialtypes', 'bands', 'regions'} ...
  , 'plot_func_inputs', {100} ...
  , 'mask', to_keep ...
);

ks_outs = dsp3.kstest2( peaks, site_labs ...
  , {'bands', 'regions'}, 'choice', 'cued' ...
  , 'mask', to_keep ...
);

if ( do_save )
  stat_save_p = fullfile( save_p, ternary(remove_outliers, 'outliers_removed', 'with_outliers') );
  
  shared_utils.plot.fullscreen( figs );
  dsp3.req_savefigs( figs, save_p, labs, {'trialtypes', 'bands', 'regions'}, prefix );
  dsp3.save_kstest_outputs( ks_outs, stat_save_p  );
end

end

function [tot_n, ns, n_I, n_C] = figure_s4_n_calculation(labels)

mask = fcat.mask( labels );

[n_I, n_C] = findall( labels, {'bands', 'regions', 'outcomes'}, mask );
ns = cellfun( @(x) numel(findall(labels, {'sites', 'channels', 'days'}, x)), n_I );
tot_n = sum( ns );

end

