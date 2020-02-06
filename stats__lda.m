function stats__lda(varargin)

import shared_utils.cell.percell;

defaults = dsp3.get_behav_stats_defaults();
defaults.specificity = 'sites';
defaults.is_cued = false;
defaults.underlying_measure = 'coherence';
defaults.analysis_type = 'lda';
defaults.smooth_func = @(x) smooth( x, 5 );
defaults.plot_spec = { 'days', 'channels', 'regions' };
defaults.plot_per_site = true;
defaults.time_window = [-250, 0];
defaults.freq_window = [ 15, 25 ];
defaults.xlims = [];
defaults.over_frequency = true;
defaults.keep_one_per_site = false;
defaults.stats_time_window = [ -250, 0 ];
defaults.take_mean_time_window_for_stats = false;
defaults.freq_roi_name = '';
defaults.line_evaluation = 'ttest2';
defaults.pro_v_anti_lines_are = 'p';
defaults.pro_v_anti_ylims = [];
defaults.lines_v_null_ylims = [];
defaults.make_figs = true;
defaults.site_smooth_func = @(x) x;
defaults.is_high_resolution_sf_coherence = false;

params = dsp3.parsestruct( defaults, varargin );

drug_type = params.drug_type;
bsd = params.base_subdir;
plot_spec = params.plot_spec;
is_cued = params.is_cued;
freq_roi_name = params.freq_roi_name;

underlying_measure = validatestring( params.underlying_measure ...
  , {'coherence', 'sfcoherence'}, mfilename, 'underlying_measure' );
analysis_type = validatestring( params.analysis_type ...
  , {'lda', 'rf'}, mfilename, 'analysis_type' );

conf = dsp3.config.load();

switch ( underlying_measure )
  case 'coherence'
    % date_dir = '061118';  % per day
    % date_dir = '061218';  % across days
    % date_dir = '072418';  % per site
    date_dir = get_data_dir( params.specificity );
  case 'sfcoherence'
    date_dir = get_sf_data_dir( params.specificity, analysis_type, is_cued ...
      , params.is_high_resolution_sf_coherence );
end

lda_dir = fullfile( conf.PATHS.dsp2_analyses, analysis_type, date_dir );

lda = get_messy_lda_data( lda_dir );

if ( strcmp(drug_type, 'nondrug') ), lda('drugs') = '<drugs>'; end
if ( ~strcmp(params.specificity, 'sites') && contains_fields(lda, 'channels') )
  lda('channels') = '<channels>'; 
end

if ( ischar(freq_roi_name) )
  freq_p = freq_roi_name;
else
  freq_p = strjoin( freq_roi_name, '_' );
end

path_components = { analysis_type, dsp3.datedir, bsd, underlying_measure, drug_type, freq_p };

params.plotp = char( dsp3.plotp(path_components) );
params.analysisp = char( dsp3.analysisp(path_components) );

%
%
%

ldalabs = setdisp( fcat.from(lda.labels), 'short' );
ldadat = lda.data * 100;  % convert to percent
freqs = lda.frequencies;

if ( ~params.is_high_resolution_sf_coherence )
  time = -500:50:500;
else
  time = -300:5:300;
end

if ( strcmp(underlying_measure, 'sfcoherence') && ~hascat(ldalabs, 'channels') )
  renamecat( ldalabs, 'sites', 'channels' );
end

%  reshape such that each permutation, currently the 3-dimension, is concatenated
%   along the first dimension, reducing the 3d-array to a matrix

if ( params.over_frequency )
  ts = params.time_window;

  t_ind = time >= ts(1) & time <= ts(2);

  t_meaned = squeeze( nanmean(ldadat(:, :, t_ind, :), 3) );
  
  x_values = freqs;
else  
  fs = params.freq_window;
  
  t_meaned = [];
  
  if ( ~iscell(fs) )
    fs = { fs };
  end
  
  assert( numel(fs) == numel(freq_roi_name), 'Rois must match band names.' ); 
  
  for i = 1:numel(fs)
    f_ind = freqs >= fs{i}(1) & freqs <= fs{i}(2);
  
    tmp_t_meaned = squeeze( nanmean(ldadat(:, f_ind, :, :), 2) );
    t_meaned = [ t_meaned; tmp_t_meaned ];
  
    x_values = time;
  end
  
  addcat( ldalabs, 'band' );
  repset( ldalabs, 'band', freq_roi_name );
end

nrows = size( t_meaned, 1 );
ncols = size( t_meaned, 2 );
niters = size( t_meaned, 3 );

newdat = zeros( niters*nrows, ncols );
newlabs = repmat( ldalabs', niters );

stp = 1;

for i = 1:niters
  newdat(stp:stp+nrows-1, :) = squeeze( t_meaned(:, :, i) );
  stp = stp + nrows;
end

if ( params.keep_one_per_site )
  newdat = keep_one_per_site( newdat, newlabs );
  newdat = keep_valid_rows( newdat, newlabs );
end

if ( ~params.plot_per_site )
  collapsecat( newlabs, {'channels', 'days'} );
end

% Stats for pro minus anti lines.
% pro_minus_anti_stat( newdat, newlabs', x_values, params );

I = findall( newlabs, plot_spec );

for i = 1:numel(I)
  shared_utils.general.progress( i, numel(I), mfilename );
  
  if ( params.make_figs )
    compare_lines( newdat, newlabs', x_values, I{i}, params );
  end
  
  handle_stats( newdat, newlabs', x_values, I{i}, params );
  
  if ( params.make_figs )
    pro_v_anti_lines( newdat, newlabs', x_values, I{i}, params );
  end
  
  if ( params.make_figs )
    pro_minus_anti_lines( newdat, newlabs', x_values, I{i}, params );
  end
end

end

function [dat, labs] = keep_valid_rows(dat, labs)

is_missing = all( dat == 0, 2 );

dat(is_missing, :) = [];
keep( labs, find(~is_missing) );

end

function [dat, labs] = keep_one_per_site(dat, labs)

I = findall( labs, {'measure', 'drugs', 'days', 'channels', 'regions' ...
  , 'trialtypes', 'contexts'} );

unqs = unique( cellfun(@numel, I) );

assert( numel(unqs) == 1 );

to_keep = cellfun( @(x) x(1), I );

dat = dat(to_keep, :);
keep( labs, to_keep );

end

function handle_stats(tdata, labels, x_values, mask, params)

mask = findnot( labels, {'targAcq', 'cued'}, mask );

test_each = { 'trialtypes', 'drugs' ...
  , 'administration', 'contexts', 'days', 'channels', 'regions', 'band' };

test_each = dsp3.nonun_or_all( labels, test_each );
describe_each = csunion( test_each, 'measure' );

[newlabs, I] = keepeach( labels', test_each, mask );
addcat( newlabs, 'time_window' );

t_window = params.stats_time_window;
t_ind = find( x_values >= t_window(1) & x_values <= t_window(2) );
subset_data = tdata(:, t_ind);

if ( params.take_mean_time_window_for_stats )
  subset_data = nanmean( subset_data, 2 );
  t_ind = t_ind(1);
end

mean_subset_data = nanmean( subset_data, 2 );

stat_labels = fcat();
all_ps = [];

for i = 1:numel(I)
  real_ind = find( labels, 'real_percent', I{i} );
  shuffled_ind = find( labels, 'shuffled_percent', I{i} );
  
  assert( numel(real_ind) == 100 && numel(shuffled_ind) == 100 );
  
  ps = zeros( numel(t_ind), 1 );
  
  for j = 1:numel(t_ind)
    ct_ind = t_ind(j);
    
    if ( strcmp(params.line_evaluation, 'ttest2') )
      [~, ps(j)] = ttest2( subset_data(shuffled_ind, j), subset_data(real_ind, j) );
    elseif ( strcmp(params.line_evaluation, 'signrank') )
      ps(j) = signrank( subset_data(shuffled_ind, j), subset_data(real_ind, j) );
    else
      ps(j) = nnz( subset_data(shuffled_ind, j) > subset_data(real_ind, j) ) / 100;
    end
    
    time_lab = sprintf( '%dms', x_values(ct_ind) );
    
    setcat( newlabs, 'time_window', time_lab );
    append( stat_labels, newlabs, i );
  end
  
  ps = dsp3.fdr( ps );
  all_ps = [ all_ps; ps ];
end

[t, rc] = tabular( stat_labels, 'time_window', {'contexts', 'regions', 'band'} );
p_tbl = fcat.table( cellrefs(all_ps, t), rc{:} );

[d_tbl, ~, d_labels] = ...
  dsp3.descriptive_table( mean_subset_data, labels', describe_each, [], mask );

context_comparison_outs = dsp3.ranksum( mean_subset_data, labels' ...
  , setdiff(describe_each, 'contexts'), 'othernone', 'selfboth', 'mask', mask );

% context_comparison_outs = dsp3.ttest2( mean_subset_data, labels' ...
%   , setdiff(describe_each, 'contexts'), 'othernone', 'selfboth', 'mask', mask );

if ( params.do_save )
  stat_prefix = sprintf( '%s_p_values', params.line_evaluation );
  
  dsp3.savetbl( p_tbl, params.analysisp, stat_labels, test_each, stat_prefix );
  dsp3.savetbl( d_tbl, params.analysisp, d_labels, describe_each, 'descriptives_' );
  
  comparison_p = fullfile( params.analysisp, 'compare_contexts' );
  
%   dsp3.save_ttest2_outputs( context_comparison_outs, comparison_p );
  dsp3.save_ranksum_outputs( context_comparison_outs, comparison_p );
end

end

function pro_minus_anti_stat(tdata, labels, x_values, params)

t_ind = x_values >= params.stats_time_window(1) & x_values <= params.stats_time_window(2);
t_data = nanmean( tdata(:, t_ind), 2 );

mask = union( find(labels, {'targAcq', 'choice'}), find(labels, {'targOn', 'cued'}) );
mask = fcat.mask( labels, mask ...
  , @find, 'real_percent' ...
);

each = { 'trialtypes', 'administration', 'days', 'channels', 'regions', 'measure', 'band' };

[pm_data, pm_labels] = dsp3.sbop( t_data, labels', each, 'othernone', 'selfboth', @minus, @(x) nanmean(x, 1) );

end

function pro_minus_anti_lines(tdata, labels, x_values, mask, params)

mask = fcat.mask( labels, mask ...
  , @findnot, {'targAcq', 'cued'} ...
  , @find, 'real_percent' ...
);

compare_each = { 'trialtypes', 'drugs' ...
  , 'administration', 'days', 'channels', 'regions', 'measure', 'band' };

I = findall( labels, compare_each, mask );

pm_labs = fcat();
pm_dat = [];

for i = 1:numel(I)
  on_ind = find( labels, 'othernone', I{i} );
  sb_ind = find( labels, 'selfboth', I{i} );
  
  pro_minus_anti = tdata(on_ind, :) - tdata(sb_ind, :);
  
  pm_dat = [ pm_dat; pro_minus_anti ];
  
  lab_copy = labels';
  setcat( lab_copy, 'contexts', 'pro_minus_anti' );
  append( pm_labs, lab_copy, on_ind );
end

% addsetcat( pm_labs, 'band', params.freq_roi_name );

assert_ispair( pm_dat, pm_labs );

gcats = { 'measure' };
pcats = { 'band', 'trialtypes', 'drugs' ...
  , 'administration', 'contexts', 'days', 'channels', 'regions' };

pl = plotlabeled.make_common();
pl.x = x_values;

axs = pl.lines( pm_dat, pm_labs, gcats, pcats );
shared_utils.plot.hold( axs, 'on' );

if ( numel(axs) == 1 )

  if ( ~isempty(params.pro_v_anti_ylims) )
    shared_utils.plot.set_ylims( axs, params.pro_v_anti_ylims );
  end

  ps = zeros( numel(x_values), 1 );
  for i = 1:numel(x_values)
    ps(i) = signrank( pm_dat(:, i) );
  end
  ps = dsp3.fdr( ps );

  add_stars( axs, x_values, ps );
end

if ( ~isempty(params.xlims) )
  shared_utils.plot.set_xlims( axs, params.xlims );
end

shared_utils.plot.add_horizontal_lines( axs, 0 );
shared_utils.plot.add_vertical_lines( axs, 50 );

if ( ~isempty(params.pro_v_anti_ylims) )
  shared_utils.plot.set_ylims( axs, params.pro_v_anti_ylims );
end

if ( params.do_save )
  use_plot_p = fullfile( params.plotp, 'pro_minus_anti_lines' );
  use_p_tbl_p = fullfile( params.analysisp, 'pro_minus_anti_lines' );
  
  dsp3.req_savefig( gcf, use_plot_p, pm_labs, pcats );
  
  x_v = arrayfun( @num2str, x_values, 'un', 0 );
  rc = fcat.strjoin( combs(pm_labs, pcats), ' | ' );
  tbl = fcat.table( ps(:)', rc, x_v );
  
  dsp3.req_writetable( tbl, use_p_tbl_p, pm_labs, pcats );
end

%%

stats_t_window = params.stats_time_window;
stats_t_ind = x_values >= stats_t_window(1) & x_values <= stats_t_window(2);

mean_stats = nanmean( pm_dat(:, stats_t_ind), 2 );

sr_spec = setdiff( compare_each, {'band', 'regions'} );
sr_labs = add_target_region_band_labels( pm_labs' );

try
  compare_beta_gamma_signrank_outs = dsp3.signrank2( mean_stats, sr_labs', sr_spec, 'beta_region', 'gamma_region' );
  compare_to_0_outs = dsp3.signrank1( mean_stats, sr_labs', compare_each ...
    , 'mask', findor(sr_labs, {'beta_region', 'gamma_region'}) );

  if ( params.do_save )
    prefix = sprintf( '%d_%dms__', stats_t_window(1), stats_t_window(2) );

    beta_gamma_p = fullfile( params.analysisp, 'compare_beta_gamma' );
    to_0_p = fullfile( params.analysisp, 'compare_to_0' );

    dsp3.save_signrank1_outputs( compare_beta_gamma_signrank_outs, beta_gamma_p, {}, prefix );
    dsp3.save_signrank1_outputs( compare_to_0_outs, to_0_p, {}, prefix );
  end
catch err
  warning( err.message );
end

end

function labels = add_target_region_band_labels(labels)

beta_blas_accf = find( labels, {'beta', 'bla_spike_acc_field'} );
gamma_accs_blaf = find( labels, {'new_gamma', 'acc_spike_bla_field'} );

if ( isempty(beta_blas_accf) || isempty(gamma_accs_blaf) )
  warning( 'No rows matched target region-band' );
end

addcat( labels, 'band_region' );
setcat( labels, 'band_region', 'beta_region', beta_blas_accf );
setcat( labels, 'band_region', 'gamma_region', gamma_accs_blaf );

end

function add_stars(ax, xs, ps)

y_max = max( get(ax, 'ylim') );

threshs = [0.05, 0.001, 0.0001];
colors = { 'y', 'g', 'r' };

for i = 1:numel(ps)
  
  if ( ps(i) > threshs(1) ), continue; end
  
  ind = 2;
  
  while ( ind <= numel(threshs) && ps(i) < threshs(ind) )
    ind = ind + 1;
  end
  
  plot( xs(i), y_max, sprintf('%s*', colors{ind-1}) );
end

end

function pro_v_anti_lines(tdata, labels, x_values, mask, params)

% context_comparison_outs = dsp3.signrank2( mean_subset_data, labels' ...
%   , setdiff(describe_each, 'contexts'), 'othernone', 'selfboth', 'mask', mask );

mask = findnot( labels, {'targAcq', 'cued'}, mask );

compare_each = { 'trialtypes', 'drugs' ...
  , 'administration', 'days', 'channels', 'regions', 'measure' };

[compare_labs, compare_I] = keepeach( labels', compare_each ...
  , find(labels, 'real_percent', mask) );

ps = zeros( numel(compare_I), 1 );
errs = nan( size(ps) );

comparison_kind = validatestring( params.pro_v_anti_lines_are ...
  , {'p', 'pro_minus_anti'} );

for i = 1:numel(compare_I)
  selfboth_ind = find( labels, 'selfboth', compare_I{i} );
  othernone_ind = find( labels, 'othernone', compare_I{i} );
  
  for j = 1:numel(x_values)
    sb_data = tdata(selfboth_ind, j);
    on_data = tdata(othernone_ind, j);
    
    if ( strcmp(comparison_kind, 'p' ) )
      ps(i, j) = signrank( sb_data, on_data );
    else
      assert( strcmp(comparison_kind, 'pro_minus_anti') );
      
      pro_minus_anti = on_data - sb_data;
      
      ps(i, j) = nanmean( pro_minus_anti );
      errs(i, j) = plotlabeled.nansem( pro_minus_anti );
    end
  end
  
  if ( strcmp(comparison_kind, 'p') )
    ps(i, :) = dsp3.fdr( ps(i, :) );
  end
end

setcat( compare_labs, 'measure', 'p value' );
% addsetcat( compare_labs, 'band', params.freq_roi_name );

gcats = { 'measure' };
pcats = { 'band', 'trialtypes', 'drugs' ...
  , 'administration', 'contexts', 'days', 'channels', 'regions' };

pl = plotlabeled.make_common();
pl.x = x_values;

axs = pl.lines( ps, compare_labs, gcats, pcats );

if ( ~isempty(params.xlims) )
  shared_utils.plot.set_xlims( axs, params.xlims );
end

if ( ~isempty(params.pro_v_anti_ylims) )
  shared_utils.plot.set_ylims( axs, params.pro_v_anti_ylims );
end

if ( strcmp(comparison_kind, 'p') )
  arrayfun( @(x) set(x, 'yscale', 'log'), axs );
else
  shared_utils.plot.hold( axs, 'on' );
  shared_utils.plot.add_horizontal_lines( axs, 0 );
end

if ( params.do_save )
  use_plot_p = fullfile( params.plotp, sprintf('%s_lines', comparison_kind) );
  
  dsp3.req_savefig( gcf, use_plot_p, compare_labs, pcats );
end

end

function compare_lines( tdata, labels, x_values, basemask, params )

F = figure(1);
clf( F );
set( F, 'defaultLegendAutoUpdate', 'off' );

mask = findnot( labels, {'targAcq', 'cued'}, basemask );

[threshs, sort_ind] = sort( [0.05, 0.001, 0.0001], 'descend' );
colors = { 'y', 'g', 'r' };
colors = colors( sort_ind );

assert( numel(colors) == numel(threshs) );

gcats = { 'measure' };
pcats = { 'trialtypes', 'drugs' ...
  , 'administration', 'contexts', 'days', 'channels', 'regions', 'band' };

pcats = dsp3.nonun_or_all( labels, pcats );

[newlabs, p_i, p_c] = keepeach( labels', pcats, mask );
plabs = fcat.strjoin( p_c, [], ' | ' );

shp = plotlabeled.get_subplot_shape( numel(p_i) );

all_ps = cell( size(p_i) );
axs = gobjects( size(all_ps) );

sfunc = dsp3.field_or_default( params, 'smooth_func', @(x) x );
do_save = dsp3.field_or_default( params, 'do_save', false );

for i = 1:numel(p_i)
  ax = subplot( shp(1), shp(2), i );
  
  hold( ax, 'on' );
  
  [g_i, g_c] = findall( labels, gcats, p_i{i} );
  glabs = fcat.strjoin( g_c, [], ' | ' );
  
  assert( numel(g_i) == 2, 'Expected 2 outcomes; got %d', numel(g_i) );
  
  first = rowref( tdata, g_i{1} );
  sec = rowref( tdata, g_i{2} );
  
  n_freqs = size( first, 2 );
  ps = zeros( 1, n_freqs );
  
  for j = 1:n_freqs  
    if ( strcmp(params.line_evaluation, 'ttest2') )
      [~, ps(j)] = ttest2( dimref(first, j, 2), dimref(sec, j, 2) );
    elseif ( strcmp(params.line_evaluation, 'signrank') )
      ps(j) = signrank( dimref(first, j, 2), dimref(sec, j, 2) );
    else
      assert( strcmp(params.line_evaluation, 'perm') ...
        , 'Unrecognized line_evaluation "%s".', params.line_evaluation );
      
      is_real = cellfun( @(x) ~isempty(strfind(x, 'real')), g_c );
      
      real_dat = tdata(g_i{is_real}, :);
      shuff_dat = tdata(g_i{~is_real}, :);
      
      ps(j) = nnz( real_dat(:, j) < shuff_dat(:, j) ) / rows( shuff_dat );
    end
  end
  
  all_ps{i} = dsp3.fdr( ps );
  
  first = rowop( first, arrayfun(@identity, 1:rows(first), 'un', 0), params.site_smooth_func );
  sec = rowop( sec, arrayfun(@identity, 1:rows(first), 'un', 0), params.site_smooth_func );
  
  mean1 = plotlabeled.nanmean( first );
  mean2 = plotlabeled.nanmean( sec );
  errs1 = plotlabeled.nansem( first );
  errs2 = plotlabeled.nansem( sec );
  
  h1 = plot( ax, x_values, sfunc(mean1) );
  h2 = plot( ax, x_values, sfunc(mean2) );
  
  ops = { @plus, @minus };
  
  for j = 1:numel(ops)
    h3 = plot( ax, x_values, ops{j}(sfunc(mean1), sfunc(errs1)) );
    h4 = plot( ax, x_values, ops{j}(sfunc(mean2), sfunc(errs2)) );

    set( h3, 'color', get(h1, 'color') );
    set( h4, 'color', get(h2, 'color') );
    set( h3, 'linewidth', get(h1, 'linewidth')/2 );
    set( h4, 'linewidth', get(h1, 'linewidth')/2 );
  end

  lines = [ h1; h2 ];
  
  legend( lines, strrep(glabs, '_', ' ') );
  title( ax, strrep(plabs{i}, '_', ' ') );
  
  axs(i) = ax;
end

shared_utils.plot.hold( axs );
shared_utils.plot.match_xlims( axs );
shared_utils.plot.match_ylims( axs );

if ( ~isempty(params.xlims) )
  shared_utils.plot.set_xlims( axs, params.xlims );
end

if ( ~isempty(params.lines_v_null_ylims) )
  shared_utils.plot.set_ylims( axs, params.lines_v_null_ylims );
end

% arrayfun( @(x) set(x, 'ylim', [-0.15, 0.15]), axs );

markersize = 8;

% add stars
for i = 1:numel(axs)
  ax = axs(i);
  lims = get( ax, 'ylim' );
  
  for j = 1:numel(threshs)
    inds = find( all_ps{i} < threshs(j) );
    colorspec = sprintf( '%s*', colors{j} );
    
    for k = 1:numel(inds)
      plot( ax, x_values(inds(k)), lims(2), colorspec, 'markersize', markersize );
    end
  end
end

shared_utils.plot.add_horizontal_lines( axs, 50 );

if ( do_save )
  use_plot_p = fullfile( params.plotp, 'lines_vs_null' );
  use_p_value_p = fullfile( params.analysisp, 'lines_vs_null'  );
  
  prefix = 'pro_anti_coh';
  shared_utils.io.require_dir( use_plot_p );
  
  fname = dsp3.fname( newlabs, csunion(pcats, 'channels') );
  fname = dsp3.prefix( prefix, fname );

  dsp3.savefig( gcf, fullfile(use_plot_p, fname) );
  
  
  ps = vertcat( all_ps{:} );
  cols = arrayfun( @num2str, x_values, 'un', 0 );
  
  p_tbl = fcat.table( ps, plabs, cols );
  dsp3.req_writetable( p_tbl, use_p_value_p, newlabs, pcats );
end

end

function date_dir = get_sf_data_dir(spec, analysis_type, is_cued, is_high_resolution)

if ( is_high_resolution )
  date_dir = '082019';
  return
end

switch ( spec )
  case 'contexts'
    if ( strcmp(analysis_type, 'lda') )
      if ( is_cued )
        date_dir = '121918';
      else
        date_dir = '121018';
      end
    else
      assert( strcmp(analysis_type, 'rf') );
      date_dir = '121718';
    end
  case 'sites'
    assert( strcmp(analysis_type, 'lda') );
    date_dir = '121318';
  otherwise
    error( 'Unrecognized specificty "%s".', spec );
end

% date_dir = '120618';

end

function date_dir = get_data_dir(spec)

switch ( spec )
  case 'contexts'
    date_dir = '061218';
  case 'days'
    date_dir = '061118';
  case 'sites'
    date_dir = '072618';
%     date_dir = '072418';
  otherwise
    error( 'Unrecognized specificty "%s".', spec );
end
end


