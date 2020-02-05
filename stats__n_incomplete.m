function stats__n_incomplete(varargin)

defaults = dsp3.get_behav_stats_defaults();
params = dsp3.parsestruct( defaults, varargin );

conf = params.config;
drug_type = params.drug_type;
bs = params.base_subdir;

if ( isempty(params.consolidated) )
  consolidated = dsp3.get_consolidated_data( conf );
else
  consolidated = params.consolidated;
end

path_components = { 'behavior', dsp3.datedir, bs, drug_type, 'p_incompleted_trials' };

params.plot_p = char( dsp3.plotp(path_components, conf) );

%%

labs = fcat.from( consolidated.trial_data.labels );

event_labels = fcat.from( consolidated.events.labels );
event_key = consolidated.event_key;
events = consolidated.events.data;

is_choice_error = dsp3.is_choice_error( events, event_labels, event_key );

dsp3.label_error_types( labs, is_choice_error );

%%

[subsetlabs] = dsp3.get_subset( labs', drug_type );
keep( subsetlabs, findnone(subsetlabs, params.remove) );

prune( subsetlabs );

%%

error_inds = find( subsetlabs, 'errors' );

addsetcat( subsetlabs, 'completed_trial', 'complete' );
setcat( subsetlabs, 'completed_trial', 'incomplete', error_inds );

prune( subsetlabs );

%%

proportion_spec = { 'days', 'contexts', 'administration' };
props_of = 'completed_trial';

prop_mask = findnone( subsetlabs, 'init_error' );

[n_complete_props, proportion_labels, prop_I] = ...
  proportions_of( subsetlabs, proportion_spec, props_of, prop_mask );

%%

analysis_complete( subsetlabs' );

plot_cued_and_choice_together( n_complete_props, proportion_labels', params );

plot_choice( n_complete_props, proportion_labels', params );
plot_cued( n_complete_props, proportion_labels', params );


% plot_cue_and_choice_together( n_complete_props, proportion_labels', params );

end

function plot_cued_and_choice_together(n_complete_props, proportion_labels, params)

%%
xcats = { 'contexts' };
gcats = { 'completed_trial' };
pcats = {};

mask = fcat.mask( proportion_labels ...
  , @find, 'complete' ...
  , @findnot, 'errors' ...
);

pl = plotlabeled.make_common();
pl.add_points = true;
pl.marker_size = 10;
pl.x_order = { 'context__self', 'context__both', 'context__other', 'none' };
pl.y_lims = [ 0, 1.1 ];
pl.point_jitter = 0;
pl.marker_type = 'o';
pl.x_tick_rotation = 0;

pltdat = n_complete_props(mask);
pltlabs = prune( proportion_labels(mask) );

axs = pl.bar( pltdat, pltlabs, xcats, gcats, pcats );

figure_1e_n_calculation( prune(pltlabs') );

if ( params.do_save )
  pltcats = unique( cshorzcat(gcats, pcats) );
  shared_utils.plot.fullscreen( gcf );
  dsp3.req_savefig( gcf, params.plot_p, prune(proportion_labels(mask)), pltcats, 'not-jittered' );
end

end

function figure_1e_n_calculation(pltlabs)

[I, C] = findall( pltlabs, {'contexts', 'outcomes'} );
n = unique( cellfun(@numel, I) );

end

function plot_cued(n_complete_props, proportion_labels, params)
%%

xcats = { 'contexts' };
gcats = { 'completed_trial' };
pcats = { 'trialtypes' };

mask = fcat.mask( proportion_labels ...
  , @find, 'complete' ...
  , @findnot, 'errors' ...
  , @find, 'cued' ...
);

pl = plotlabeled.make_common();
pl.add_points = true;
pl.marker_size = 10;
pl.x_order = { 'context__self', 'context__both', 'context__other', 'none' };
pl.y_lims = [ 0.4, 1.1 ];
pl.point_jitter = 0.8;
pl.marker_type = 'o';
pl.x_tick_rotation = 0;

pltdat = n_complete_props(mask);
pltlabs = prune( proportion_labels(mask) );

axs = pl.bar( pltdat, pltlabs, xcats, gcats, pcats );

anova_outs = dsp3.anova1( pltdat, pltlabs', {}, 'contexts' ...
  , 'remove_nonsignificant_comparisons', false ...
);

% axs = pl.boxplot( pltdat, pltlabs, gcats, pcats );

if ( params.do_save )
  pltcats = unique( cshorzcat(gcats, pcats) );
  shared_utils.plot.fullscreen( gcf );
  dsp3.req_savefig( gcf, params.plot_p, prune(proportion_labels(mask)), pltcats );
end

end

function analysis_complete(subsetlabs)

%%  N complete per sessions

monk_spec = { 'monkeys' };

complete_each = { 'days' };

mask = fcat.mask( subsetlabs ...
  , @findnone, 'errors' ...
);

[complete_labs, complete_I] = keepeach( subsetlabs', complete_each, mask );

trial_counts = cellfun( @numel, complete_I );

[mean_labs, mean_I] = keepeach_or_one( complete_labs', monk_spec );

trial_table = dsp3.descriptive_table( trial_counts, complete_labs', monk_spec );

end

function plot_choice(n_complete_props, proportion_labels, params)
%%

xcats = { 'contexts' };
gcats = { 'completed_trial' };
pcats = { 'trialtypes' };

mask = fcat.mask( proportion_labels ...
  , @find, 'complete' ...
  , @findnot, 'errors' ...
  , @find, 'choice' ...
);

pl = plotlabeled.make_common();
pl.add_points = true;
pl.y_lims = [ 0.4, 1.1 ];
pl.x_tick_rotation = 0;

pltdat = n_complete_props(mask);
pltlabs = proportion_labels(mask);

axs = pl.bar( pltdat, pltlabs, xcats, gcats, pcats );

% axs = pl.boxplot( pltdat, pltlabs, [xcats, gcats], pcats );


if ( params.do_save )
  pltcats = unique( cshorzcat(gcats, pcats) );
  shared_utils.plot.fullscreen( gcf );
  dsp3.req_savefig( gcf, params.plot_p, prune(proportion_labels(mask)), pltcats );
end

end

function plot_cue_and_choice_together(n_complete, proportion_labels, params)

%%

pltlabels = proportion_labels';

cue_outs = { 'self', 'both', 'other', 'none' };
ctx_outs = cellfun( @(x) ['context__', x], cue_outs, 'un', 0 );
replace_outs = cellfun( @(x) [x, '-cue'], cue_outs, 'un', 0 );

cellfun( @(x, y) replace(pltlabels, x, y), ctx_outs, replace_outs, 'un', 0 );
replace( pltlabels, 'othernone', 'other/bottle' );
replace( pltlabels, 'selfboth', 'self/both' );
replace( pltlabels, 'none-cue', 'bottle-cue' );

mask = fcat.mask( pltlabels ...
  , @find, 'incomplete' ...
);

pl = plotlabeled.make_common();
pl.y_lims = [ 0, 1 ];
pl.x_tick_rotation = 0;
pl.x_order = { 'other/bottle', 'self/both', 'self-cue', 'other-cue', 'both-cue' };

fcats = {};
xcats = { 'contexts' };
gcats = { 'completed_trial' };
pcats = {};

f_I = findall_or_one( pltlabels, fcats, mask );

figs = gobjects( numel(f_I), 1 );
axs = gobjects;
for i = 1:numel(f_I)
  figs(i) = figure(i);
  clf( figs(i) );
  pl.fig = figs(i);
  
  ax = pl.bar( n_complete_props(f_I{i}), pltlabels(f_I{i}), xcats, gcats, pcats );
  axs = [ axs, ax(:)' ];
  
  shared_utils.plot.ylabel( ax, '% Incompleted trials' );
end

if ( params.do_save )
  pltcats = unique( cshorzcat(fcats, xcats, gcats, pcats) );
  
  for i = 1:numel(f_I)
    dsp3.req_savefig( figs(i), params.plot_p, prune(pltlabels(f_I{i})), pltcats );
  end
end

end