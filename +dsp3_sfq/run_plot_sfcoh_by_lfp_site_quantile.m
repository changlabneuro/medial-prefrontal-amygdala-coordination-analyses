function run_plot_sfcoh_by_lfp_site_quantile(varargin)

defaults = dsp3.get_common_plot_defaults( dsp3.get_common_make_defaults() );
defaults.config = dsp3.set_dataroot( '/Volumes/My Passport/NICK/Chang Lab 2016/dsp3/' );
defaults.time_window_mean = true;
defaults.coh = [];
defaults.coh_labs = fcat();
defaults.freqs = [];
defaults.t = [];
defaults.psd = [];
defaults.psd_labs = fcat();
defaults.include_unit_uuid_specificity = false;
defaults.line_y_lims = [];

params = dsp3.parsestruct( defaults, varargin );

conf = params.config;
coh = params.coh;

coh_t_min = 0;
coh_t_max = 150;

if ( isempty(coh) )
  % Load coherence
  [coh, coh_labs, freqs, t] = dsp3_sfq.load_per_day_sfcoh( conf );
  [coh, coh_labs] = dsp3_sfq.band_meaned_data( coh, coh_labs', freqs );
  
  if ( params.time_window_mean )
    coh = nanmean( coh(:, mask_gele(t, coh_t_min, coh_t_max)), 2 );
  end

  dsp3_sfq.add_spike_lfp_region_labels( coh_labs );
  
  % Load power
  [psd, psd_labs, psd_freqs, psd_t] = dsp3_sfq.load_summarized_psd();
  psd = psd(:, mask_gele(psd_freqs, 0, 100), :);
  assert( size(psd, 2) == numel(freqs) && size(psd, 3) == numel(t) );

  psd_mask = fcat.mask( psd_labs ...
    , @find, {'choice', 'pre'} ...
    , @findnone, 'errors' ...
  );

  psd_each = { 'trialtypes', 'regions', 'channels', 'days' };
  [psd, psd_labs] = dsp3_sfq.band_meaned_data( psd, psd_labs', freqs, psd_each, psd_mask );

  psd_t_min = 0;
  psd_t_max = 150;

  psd = nanmean( psd(:, mask_gele(t, psd_t_min, psd_t_max)), 2 );
else
  coh_labs = params.coh_labs';
  freqs = params.freqs;
  t = params.t;
  psd = params.psd;
  psd_labs = params.psd_labs';
  
  assert_ispair( coh, coh_labs );
  assert_ispair( psd, psd_labs );
end

[quant_labs, quant_mask, quants_of] = make_quantile_labels( psd, psd_labs' );

%%

is_pro_minus_anti = true;
do_save = params.do_save;

save_components = { 'sfcoh_by_quantile', 'by_lfp_site', dsp3.datedir };
plot_p = char( dsp3.plotp(save_components) );
analysis_p = char( dsp3.analysisp(save_components) );

to_label = addcat( coh_labs', 'quantile' );

[quant_I, quant_C] = findall( quant_labs, {'regions', 'bands'}, quant_mask );

for i = 1:numel(quant_I)
  region = quant_C{1, i};
  band = quant_C{2, i};
  
  region_search_str = sprintf( 'lfp_%s', region );
  
  [site_I, site_C] = findall( quant_labs, quants_of, quant_I{i} );
  
  for j = 1:numel(site_I)
    quant_name = combs( quant_labs, 'quantile', site_I{j} );
    assert( numel(quant_name) == 1 );
    
    match_site_ind = find( to_label, [site_C(:, j)', region_search_str] );
    
    if ( isempty(match_site_ind) )
      continue;
    end
    
    setcat( to_label, 'quantile', quant_name, match_site_ind );
  end
  
  proanti_mask = fcat.mask( to_label ...
    , @find, {'beta', 'new_gamma', region_search_str} ...
    , @find, 'selected-site' ...
  );

  mask_a = intersect( proanti_mask, find(to_label, {'acc_bla', 'new_gamma'}) );
  mask_b = intersect( proanti_mask, find(to_label, {'bla_acc', 'beta'}) );
  proanti_mask = union( mask_a, mask_b );

  proanti_spec = { 'trialtypes', 'bands', 'channels', 'regions', 'days' };
  
  if ( params.include_unit_uuid_specificity )
    proanti_spec{end+1} = 'unit_uuid';
  end

  [proanti_dat, proanti_labs] = dsp3.pro_v_anti( coh, to_label', proanti_spec, proanti_mask );

  if ( is_pro_minus_anti )
    [proanti_dat, proanti_labs] = dsp3.pro_minus_anti( proanti_dat, proanti_labs', proanti_spec );
  end
  
  if ( params.time_window_mean )
  %   box_plot_and_anova( proanti_dat, proanti_labs', band, region, plot_p, analysis_p, do_save );
    scatter_plot_and_stats( proanti_dat, proanti_labs', band, region, plot_p, analysis_p, do_save );
  else
    lines_over_time( proanti_dat, t, proanti_labs', band, region, plot_p, analysis_p, do_save, params );
  end
end

end

function lines_over_time(proanti_dat, t, proanti_labs, band, region, plot_p, analysis_p, do_save, params)

%%
pl = plotlabeled.make_common();
pl.x = t;
pl.smooth_func = @(x) smoothdata(x, 'SmoothingFactor', 0.75);
pl.add_smoothing = true;
  
gcats = { 'outcomes', 'quantile' };
pcats = { 'regions', 'bands' };

pltdat = proanti_dat;
pltlabs = proanti_labs';

axs = pl.lines( pltdat, pltlabs, gcats, pcats );
shared_utils.plot.set_xlims( axs, [-300, 300] );

if ( ~isempty(params.line_y_lims) )
  shared_utils.plot.set_ylims( axs, params.line_y_lims );
end

window_means = nanmean( proanti_dat(:, t >= 0 & t <= 150), 2 );
no_nans = find( ~isnan(window_means) );
quantiles = fcat.parse( cellstr(pltlabs, 'quantile'), 'quantile_' );

[scatter_I, scatter_C] = findall( pltlabs, pcats, no_nans );
tbls = cell( numel(scatter_I), 1 );

correlation_types = { 'spearman', 'pearson' };

corr_outs = dsp3.corr( quantiles, window_means, pltlabs', pcats ...
  , 'mask', no_nans ...
  , 'corr_inputs', {'rows', 'complete', 'type', 'spearman'} ...
);

for idx = 1:numel(correlation_types)
  for i = 1:numel(scatter_I)
    ind = scatter_I{i};

    [rho, p] = corr( quantiles(ind), window_means(ind), 'rows', 'complete', 'type', correlation_types{idx} );
    tbls{i} = table( rho, p );
    tbls{i}.Properties.RowNames = fcat.strjoin( scatter_C(:, i), ' | ');
  end
  
  prefix = sprintf( 'lines__quantiles_of_%s_%s_field__', band, region );
  tbl_prefix = sprintf( '%s__%s', correlation_types{idx}, prefix );

  if ( do_save )
    shared_utils.plot.fullscreen( gcf );
    dsp3.req_savefig( gcf, plot_p, pltlabs, [gcats, pcats], prefix );

    for i = 1:numel(tbls)
      tbl_labs = prune( pltlabs(scatter_I{i}) );

      dsp3.req_writetable( tbls{i}, analysis_p, tbl_labs, pcats, tbl_prefix );
    end
  end
end

end

function scatter_plot_and_stats(proanti_dat, proanti_labs, band, region, plot_p, analysis_p, do_save)

%%

pl = plotlabeled.make_common();
  
gcats = { 'outcomes' };
pcats = { 'regions', 'bands' };

no_nans = find( ~isnan(proanti_dat) );

pltdat = proanti_dat(no_nans);
pltlabs = keep( proanti_labs', no_nans );

% [ns, tot_n, n_I, n_C] = figure_s3b_n_calculation( pltlabs' );
[ns, tot_n, n_I, n_C] = figure_s3b_n_calculation( proanti_labs' );

quantiles = fcat.parse( cellstr(proanti_labs, 'quantile', no_nans), 'quantile_' );

[axs, ids] = pl.scatter( quantiles, pltdat, pltlabs, gcats, pcats );
plotlabeled.scatter_addcorr( ids, quantiles, pltdat );

if ( do_save )
  prefix = sprintf( 'quantiles_of_%s_%s_field__', band, region );
  shared_utils.plot.fullscreen( gcf );
  dsp3.req_savefig( gcf, plot_p, pltlabs, [gcats, pcats], prefix );
%   dsp3.save_anova_outputs( anova_outs, analysis_p, [gcats, pcats], prefix );
end

end

function [ns, tot_n, I, C] = figure_s3b_n_calculation(labels)

%%

bands = combs( labels, 'bands' );
[I, C] = findall( labels, {'regions', 'bands', 'quantile'}, find(labels, bands(1)) );

[~, sort_ind] = sort( fcat.parse(C(3, :), 'quantile_') );
I = I(sort_ind);
C = C(:, sort_ind);

ns = cellfun( @(x) numel(findall(labels, {'unit_uuid', 'channels', 'days'}, x)), I );

% ns = cellfun( @numel, I );
tot_n = sum( ns );

end

function box_plot_and_anova(proanti_dat, proanti_labs, band, region, plot_p, analysis_p, do_save)

pl = plotlabeled.make_common();
  
xcats = { 'quantile' };
gcats = { 'outcomes' };
pcats = { 'regions', 'bands' };

pltdat = proanti_dat;
pltlabs = proanti_labs';

axs = pl.boxplot( pltdat, pltlabs, xcats, [gcats, pcats] );
ylabel( axs(1), 'Spike-field coherence' );

anova_spec = { 'trialtypes', 'regions', 'bands', 'outcomes' };
anova_outs = dsp3.anova1( pltdat, pltlabs', anova_spec, 'quantile' );

if ( do_save )
  prefix = sprintf( 'quantiles_of_%s_%s_field__', band, region );
  shared_utils.plot.fullscreen( gcf );
  dsp3.req_savefig( gcf, plot_p, pltlabs, [gcats, pcats], prefix );
  dsp3.save_anova_outputs( anova_outs, analysis_p, [gcats, pcats], prefix );
end

end

function [quant_labs, quant_mask, quants_of] = make_quantile_labels(psd, psd_labs)

num_tiles = 3;

% separate quantiles for each region + band
quants_each = { 'regions', 'bands' }; 

% sites
quants_of = { 'channels', 'days' };

% only choice + pre
quant_mask = fcat.mask( psd_labs ...
  , @find, {'choice', 'pre', 'beta', 'new_gamma'} ...
);

quant_labs = dsp3_sfq.quantiles_each( psd, psd_labs, num_tiles, quants_each, quants_of, quant_mask );

end


