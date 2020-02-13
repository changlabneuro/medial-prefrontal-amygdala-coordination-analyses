function [spike_trials_sua_all, include, region_label, name, day_cnt] = get_spike_trial

%% align spikes to different epochs of the task

% load single unit activity data
out = load('dictator_game_SUAdata_pre');
% add dependence
addpath(genpath('/Users/chuchengchi/Documents/MATLAB/repositories/dsp3'))
addpath(genpath('/Users/chuchengchi/Documents/MATLAB/repositories/shared_utils'))
addpath(genpath('/Users/chuchengchi/Documents/MATLAB/repositories/dsp'))
addpath(genpath('/Users/chuchengchi/Documents/MATLAB/repositories/global'))

% consolidated = dsp3.get_consolidated_data();
consolidated = dsp3.load_fixed_consolidated_data();
% event time
fxON = consolidated.event_key('fixOn');
cueON = consolidated.event_key('cueOn');
targON = consolidated.event_key('targOn');
targACQ = consolidated.event_key('targAcq');
rwON = consolidated.event_key('rwdOn');

% free or forced choice
choice_tag = 'choice'; % cued

% get region label and list for removing units don't have enough trials 
region_label = {};
cnt = 0;
include = [];
name = {};
day_cnt = {}; % idx within each day
cnt_s = 1;
for iday = 1:numel(out.all_spike_time)  
    events = out.all_event_time{iday}.event;
    events_choice = events.only(choice_tag);
    nU = numel(out.all_spike_time{iday}.data);
    day_cnt{iday} = cnt_s:(cnt_s+nU-1);
    cnt_s = (cnt_s+nU-1)+1;
    for i_u = 1:nU     
%         i_u;
        cnt = cnt+1;
        name{cnt} = [events('days') iday i_u];
        spike = out.all_spike_time{iday}.data{i_u}.data;     
        % bin spike, finding the units have no spikes in 3/4 of all trials
        mx = max(max(events_choice.data));
        t = 0:mx;
        C = histc(spike,t);
        if sum(C > 0.5) < length(t)*(1/4)
           include = [include cnt];
        end  
        region_label{cnt} = out.all_spike_time{iday}.data{i_u}.name{1};
    end
end

% get spike count for trials
cnt = 0;
spike_trials_sua_all = {};
for iday = 1:numel(out.all_spike_time)
    % event 
    events = out.all_event_time{iday}.event;
    events_choice = events.only(choice_tag);
%     events_choice = events_choice.only('pre'); % 
    nU = numel(out.all_spike_time{iday}.data);
    for i_u = 1:nU
         spike = out.all_spike_time{iday}.data{i_u}.data;                   
         spk_s_trials_a = get_spks_trials(spike,events_choice,cueON,[-0.55 0.55]);
         spk_s_trials_b = get_spks_trials(spike,events_choice,targON,[-0.55 0.55]);   
         [spk_s_trials_c, Lb] = get_spks_trials(spike,events_choice,targACQ,[-0.55 0.55]);
         spk_s_trials_d = get_spks_trials(spike,events_choice,rwON,[-0.55 0.55]);
         cnt = cnt + 1;
         spike_trials_sua_all.cueon{cnt} = spk_s_trials_a;
         spike_trials_sua_all.targon{cnt} = spk_s_trials_b;
         spike_trials_sua_all.targacq{cnt} = spk_s_trials_c;
         spike_trials_sua_all.rwon{cnt} = spk_s_trials_d;
         spike_trials_sua_all.L{cnt} = Lb;
    end     
end


function [spk_s_trials, L] = get_spks_trials(spike,events_choice,event_label,trng)
% cond
cs = {'self','both','other','none'};
L = {};
spk_s_trials = {};
for i_c = 1:numel(cs)
    events_cond = events_choice.only(cs{i_c});
    events_data = events_cond.data;
    L_v = {};
    high_ind = find(where(events_cond, 'high'));
    medium_ind = find(where(events_cond, 'medium'));
    low_ind = find(where(events_cond, 'low'));
    
    for r = 1:size(events_data,1)
        if ismember(r,high_ind)
           L_v{r} = 'h';
        elseif ismember(r,medium_ind)
           L_v{r} = 'm';
        elseif ismember(r,low_ind)   
           L_v{r} = 'l';
        else
           error('not found') 
        end    
    end
    L{i_c} = L_v;
    
    event_t = events_data(:,event_label);     
    for i_t = 1:length(event_t)        
        wins = [event_t(i_t)+trng(1) event_t(i_t)+trng(2)];              
        idc = find(spike>=wins(1) & spike<wins(2)); 
        spk_time_trial = spike(idc)-event_t(i_t);
        spk_s_trials{i_c}{i_t} = spk_time_trial; 
    end    
end
    
