function plot_sfcoh_by_look_type(coh, coh_labs, freqs, t, look_outs)

if ( nargin == 0 )
  conf = dsp3.set_dataroot( '/Volumes/My Passport/NICK/Chang Lab 2016/dsp3/' );
  [coh, coh_labs, freqs, t] = dsp3_sfq.load_per_day_sfcoh( conf );

  look_outs = dsp3_find_iti_looks( ...
      'config', conf ...
    , 'require_fixation', false ...
    , 'look_back', -3.3 ...
    , 'is_parallel', true ...
  );
else
  conf = dsp3.config.load();
end

labels = dsp3_add_iti_first_look_labels( look_outs.labels', look_outs, 0.15 );

%%

new_labs = dsp3_ct.add_first_look_labels_to_sfcoh( coh_labs', labels );
dsp3_ct.add_site_labels( new_labs );

setcat( new_labs, 'duration', 'long_enough__true', find(new_labs, 'no_look') );

%%

require_all_sites = false;
keep_prop = 0.2;

site_mask = fcat.mask( new_labs ...
  , @find, {'choice', 'selected-site'} ...
  , @findnone, 'errors' ...
);

site_spec = dsp3_ct.site_specificity();
look_spec = csunion( site_spec, {'duration', 'looks_to'} );

site_I = findall( new_labs, site_spec, site_mask );
all_to_keep = {};
for i = 1:numel(site_I)
  no_look_ind = find( new_labs, 'no_look', site_I{i} );
  rest_ind = setdiff( site_I{i}, no_look_ind );
  num_to_keep = ceil( numel(no_look_ind) * keep_prop );
  to_keep = sort( no_look_ind(randperm(numel(no_look_ind), num_to_keep)) );
  all_to_keep{end+1, 1} = sort( [to_keep; rest_ind] );
end

site_mask = vertcat( all_to_keep{:} );

if ( require_all_sites )
  [site_coh, site_labs] = dsp3_ct.site_meaned_sfcoh( coh, new_labs', site_mask );
  
  [~, missing_removed] = dsp3_ct.remove_missing_sites( site_coh, site_labs' );
  present_sites = combs( missing_removed, 'sites' );
  present_ind = find( site_labs, present_sites );

  site_mask = find( new_labs, present_sites, site_mask );
end

[site_coh, site_labs] = dsp3_ct.site_meaned_sfcoh( coh, new_labs', site_mask, look_spec );

%%

tmp_coh = site_coh;
tmp_labs = site_labs';

plot_p = char( dsp3.plotp({'sfcoh_by_gaze', dsp3.datedir}) );

pro_v_anti = false;
pro_minus_anti = false;
per_outcome = false;

if ( ~per_outcome )
  collapsecat( tmp_labs, 'outcomes' );
end

minus_no_look = true;

% clims = [-0.02, 0.02];

proanti_each = setdiff( look_spec, 'outcomes' );
look_each = setdiff( look_spec, {'looks_to', 'outcomes'} );

if ( minus_no_look )
  [bottle_coh, bottle_labs] = ...
    dsp3.sbop( tmp_coh, tmp_labs', look_each, 'bottle', 'no_look', @minus, @(x) nanmean(x, 1) );
  setcat( bottle_labs, 'looks_to', 'bottle - no_look' );
  [monk_coh, monk_labs] = ...
    dsp3.sbop( tmp_coh, tmp_labs', look_each, 'monkey', 'no_look', @minus, @(x) nanmean(x, 1) );
  setcat( monk_labs, 'looks_to', 'monkey - no_look' );
  
  [tmp_coh, tmp_labs] = appendpair( bottle_coh, bottle_labs', monk_coh, monk_labs' );
end

if ( pro_v_anti )
  [tmp_coh, tmp_labs] = dsp3.pro_v_anti( tmp_coh, tmp_labs', proanti_each );
end
if ( pro_minus_anti )
  [tmp_coh, tmp_labs] = dsp3.pro_minus_anti( tmp_coh, tmp_labs', proanti_each );
end

if ( ~per_outcome )
  [tmp_labs, look_I] = keepeach( tmp_labs', proanti_each );
  tmp_coh = bfw.row_nanmean( tmp_coh, look_I );
end

%%  save

do_save_coh = false;

save_mask = fcat.mask( tmp_labs ...
  , @find, {'choice', 'selected-site'} ...
  , @findnone, 'errors' ...
  , @find, 'long_enough__true' ...
);

save_p = fullfile( dsp3.dataroot(conf), 'data', 'sfcoh', 'gaze' );
save_labs = gather( prune(tmp_labs(save_mask)) );
save_coh = tmp_coh(save_mask, :, :);

if ( do_save_coh )
  save( fullfile(save_p, 'gaze_sf_coherence.mat') ...
    , 'save_coh', 'save_labs', 'freqs', 't', '-v7.3' );
end

%%

do_save = false;
clims = [];

f_ind = freqs >= 10 & freqs <= 100;
t_ind = t >= -300 & t <= 300;

plt_f = freqs(f_ind);
plt_t = t(t_ind);

fig_cats = { 'trialtypes', 'regions' };
pcats = { 'outcomes', 'regions', 'trialtypes', 'looks_to' };

if ( ~per_outcome )
  pcats = setdiff( pcats, 'outcomes' );
end

formats = { 'epsc', 'png', 'fig', 'svg' };

plt_mask = fcat.mask( tmp_labs ...
  , @find, {'long_enough__true'} ...
);

fig_I = findall_or_one( tmp_labs, fig_cats, plt_mask );

store_labs = cell( size(fig_I) );
store_axs = cell( size(fig_I) );
figs = gobjects( size(fig_I) );

pl = plotlabeled.make_spectrogram( plt_f, plt_t );
pl.sort_combinations = true;
pl.add_smoothing = true;

for i = 1:numel(fig_I)
  fig = figure(i);
  pl.fig = fig;
  
  plt_coh = tmp_coh(fig_I{i}, f_ind, t_ind);
  plt_labs = prune( tmp_labs(fig_I{i}) );

  axs = pl.imagesc( plt_coh, plt_labs, pcats );
  shared_utils.plot.hold( axs, 'on' );
  
  shared_utils.plot.tseries_xticks( axs, plt_t );
  shared_utils.plot.fseries_yticks( axs, round(flip(plt_f)), 5 );
  shared_utils.plot.add_vertical_lines( axs, find(plt_t == 0) );
  
%   shared_utils.plot.set_clims( axs, [-15e-3, 15e-3] );
  
  store_axs{i} = axs(:);
  store_labs{i} = plt_labs;
  figs(i) = fig;
end

axs = vertcat( store_axs{:} );

if ( isempty(clims) )
  shared_utils.plot.match_clims( axs );
else
  shared_utils.plot.set_clims( axs, clims );
end

if ( do_save )
  for i = 1:numel(fig_I)
    spectra_p = fullfile( plot_p, 'spectra' );
    shared_utils.plot.fullscreen( figs(i) );
    dsp3.req_savefig( figs(i), spectra_p, store_labs{i}, pcats, '', formats ); 
  end
end

%%  stat frequency windows

analysis_p = char( dsp3.analysisp({'sfcoh_by_gaze', dsp3.datedir}) );
analysis_p = fullfile( analysis_p, 'compare_within_band' );

f_ind = freqs >= 10 & freqs <= 80;
t_ind = t >= 0 & t <= 150;

plt_mask = fcat.mask( tmp_labs ...
  , @find, {'long_enough__true'} ...
);

stat_coh = tmp_coh(plt_mask, f_ind, t_ind);
stat_labs = prune( tmp_labs(plt_mask) );

stat_coh = nanmean( stat_coh, 3 );
[stat_coh, stat_labs] = dsp3.get_band_means( stat_coh, stat_labs, freqs(f_ind), dsp3.get_bands('map') );

rs_outs = dsp3.ranksum( stat_coh, stat_labs, {'outcomes', 'trialtypes', 'regions', 'bands'} ...
  , 'monkey - no_look', 'bottle - no_look' ...
  , 'mask', find(stat_labs, {'beta', 'new_gamma'}) ...
);

ind1 = find( rs_outs.rs_labels, {'new_gamma', 'acc_bla'} );
ind2 = find( rs_outs.rs_labels, {'beta', 'bla_acc'} );

outs = rs_outs;
[outs.rs_tables, outs.rs_labels] = indexpair( outs.rs_tables, outs.rs_labels', [ind1; ind2] );

dsp3.save_ranksum_outputs( outs, analysis_p );

%%  lines

over_freqs = [ true ];
% smoothings = [ false ];
smoothings = [ true, false ];

plt_combs = dsp3.numel_combvec( over_freqs, smoothings );

do_save = false;

for idx = 1:size(plt_combs, 2)
  
use_coh = tmp_coh;
use_labs = tmp_labs';

prefix = '';
ylims = [];
over_freq = over_freqs(plt_combs(1, idx));
is_smoothed = smoothings(plt_combs(2, idx));

f_ind = freqs >= 10 & freqs <= 80;

if ( over_freq )
  t_ind = t >= 0 & t <= 150;
else
  t_ind = t >= -300 & t <= 300;
  
  [use_coh, use_labs] = dsp3.get_band_means( use_coh, use_labs', freqs, dsp3.get_bands('map') );
end

fig_cats = { 'trialtypes', 'outcomes' };
gcats = { 'looks_to' };
pcats = { 'trialtypes', 'regions', 'outcomes' };

if ( ~per_outcome )
  fig_cats = setdiff( fig_cats, 'outcomes' );
  pcats = setdiff( pcats, 'outcomes' );
end

if ( ~over_freq )
  fig_cats{end+1} = 'bands';
  pcats{end+1} = 'bands';
end

formats = { 'epsc', 'png', 'fig', 'svg' };

plt_mask = fcat.mask( use_labs ...
  , @find, {'long_enough__true'} ...
);

if ( ~over_freq )
  plt_mask = find( use_labs, {'beta', 'new_gamma'}, plt_mask );
end

fig_I = findall_or_one( use_labs, fig_cats, plt_mask );

store_labs = cell( size(fig_I) );
store_axs = cell( size(fig_I) );
figs = gobjects( size(fig_I) );

for i = 1:numel(fig_I)
  pl = plotlabeled.make_common();
  pl.smooth_func = @(x) smoothdata(x, 'smoothingfactor', 0.5 );
  pl.add_smoothing = is_smoothed;
  
  if ( ~isempty(ylims) )
    pl.y_lims = ylims;
  end
  
  fig = figure(i);
  pl.fig = fig;
  
  plt_labs = prune( use_labs(fig_I{i}) );
  
  if ( over_freq )
    plt_coh = use_coh(fig_I{i}, f_ind, t_ind);
    plt_coh = squeeze( nanmean(plt_coh, 3) );
    x = freqs(f_ind);
  else
    plt_coh = use_coh(fig_I{i}, t_ind);
    x = t(t_ind);
  end
  
  pl.x = x;

  [axs, hs, inds] = pl.lines( plt_coh, plt_labs, gcats, pcats );
%   shared_utils.plot.hold( axs, 'on' );

  figure_3ik_n_calculation( plt_labs' );
  figure_3ik_stat_time_window( plt_coh, plt_labs', freqs(f_ind) );

  dsp3.compare_series( axs, inds, plt_coh, @ranksum ...
    , 'x', x ...
  );

  store_axs{i} = axs(:);
  store_labs{i} = plt_labs;
  figs(i) = fig;
end

axs = vertcat( store_axs{:} );

if ( isempty(ylims) )
  shared_utils.plot.match_ylims( axs );
else
  shared_utils.plot.set_ylims( axs, ylims );
end

smooth_prefix = ternary( is_smoothed, 'smoothed', 'nonsmoothed' );
freq_prefix = ternary( over_freq, 'overfreq', 'overtime' );

if ( do_save )
  use_prefix = sprintf( '%s-%s-%s', smooth_prefix, freq_prefix, prefix );
  
  for i = 1:numel(fig_I)
    line_p = fullfile( plot_p, 'lines' );
    shared_utils.plot.fullscreen( figs(i) );
    dsp3.req_savefig( figs(i), line_p, store_labs{i}, pcats, use_prefix, formats ); 
  end
end

end

%%  boxes / lines

do_save = false;
prefix = 'lines__';

over_time = true;

if ( ~over_time )
  t_ind = t >= 0 & t <= 150;
else
  t_ind = t >= -300 & t <= 300;
end

plt_t = t(t_ind);

ylims = [];

fig_cats = { 'trialtypes', 'outcomes', 'bands' };
gcats = { 'looks_to' };
pcats = { 'trialtypes', 'regions', 'outcomes', 'bands' };

if ( ~per_outcome )
  fig_cats = setdiff( fig_cats, 'outcomes' );
  pcats = setdiff( pcats, 'outcomes' );
end

formats = { 'epsc', 'png', 'fig', 'svg' };

use_coh = tmp_coh;
use_labs = tmp_labs';

[use_coh, use_labs] = dsp3.get_band_means( use_coh, use_labs', freqs, dsp3.get_bands('map'), @nanmedian );

plt_mask = fcat.mask( use_labs ...
  , @find, {'long_enough__true'} ...
  , @find, {'beta'} ...
);

fig_I = findall_or_one( use_labs, fig_cats, plt_mask );

store_labs = cell( size(fig_I) );
store_axs = cell( size(fig_I) );
figs = gobjects( size(fig_I) );

for i = 1:numel(fig_I)
  pl = plotlabeled.make_common();
  pl.summary_func = @plotlabeled.nanmedian;
  
  fig = figure(i);
  pl.fig = fig;
  
  plt_coh = use_coh(fig_I{i}, t_ind);
  plt_labs = prune( use_labs(fig_I{i}) );

  if ( ~over_time )
    plt_coh = nanmean( plt_coh, 2 );
  end

  if ( over_time )
    pl.x = plt_t;
    axs = pl.lines( plt_coh, plt_labs, gcats, pcats );
  else
    axs = pl.boxplot( plt_coh, plt_labs, gcats, pcats );
  end

  shared_utils.plot.hold( axs, 'on' );
  
  store_axs{i} = axs(:);
  store_labs{i} = plt_labs;
  figs(i) = fig;
end

axs = vertcat( store_axs{:} );

if ( isempty(ylims) )
  shared_utils.plot.match_ylims( axs );
else
  shared_utils.plot.set_ylims( axs, ylims );
end

if ( do_save )
  for i = 1:numel(fig_I)
    line_p = fullfile( plot_p, 'boxes' );
    shared_utils.plot.fullscreen( figs(i) );
    dsp3.req_savefig( figs(i), line_p, store_labs{i}, pcats, prefix, formats ); 
  end
end

end

function figure_3ik_stat_time_window(coh, labels, freqs)

assert( numel(freqs) == size(coh, 2) );
[band_ranges, band_names] = dsp3.some_bands( {'new_gamma', 'beta'} );

[banddat, bandlabs] = dsp3.get_band_means( coh, labels', freqs, band_ranges, band_names );

a = 'bottle - no_look';
b = 'monkey - no_look';

acc_mask = fcat.mask( bandlabs ...
  , @find, {'acc_bla', 'new_gamma'} ...
);

bla_mask = fcat.mask( bandlabs ...
  , @find, {'bla_acc', 'beta'} ...
);

mask = union( acc_mask, bla_mask );

rs_outs = dsp3.ranksum( banddat, bandlabs', {'bands', 'regions'}, a, b ...
  , 'mask', mask ...
);

end

function figure_3ik_n_calculation(labels)

[I, C] = findall( labels, {'regions', 'looks_to'} );
n = cellfun( @numel, I );

end