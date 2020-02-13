

% this script was used to compute the spike field coherence results in the figure 3 of NN paper.
% Specialized medial prefrontal-amygdala coordination in other-regarding
% decision preference, Copyright 2020 chengchi chu

% add chronux toolbox
addpath(genpath('/Users/chengchichu/Desktop/BRAINanalysis/chronux_2_12/spectral_analysis'))

% data roots and load data
root = '/Volumes/172.28.142.132/data/cc_dictator/mua/';
addpath('/Users/chengchichu/Desktop/dictator_analysis/coher/mat')
out = load('dictator_game_SUAdata_pre');
m = load('sua_lfp_daymatch.mat');  
z = load('days_targacq.mat');

% using blocks/sessions without drug administration
drug = 0;
if drug == 1
   session = 'pre';
end   

% which epoch
epoch = 'choice';
event_key = 4; 
region_spike = 'bla';
region_lfp = 'acc';


cnt2 = 0; % how many spike unit
cnt3 = 0; % how many field unit
coher_data_all = {};
phase_data_all = {};
% main loop begins
for iday = 1:numel(out.all_spike_time)
    iday
       
    % spike and event
    dataName_spike = out.all_spike_time{iday}.filename;
    data_spike = out.all_spike_time{iday}.data;
    data_event = out.all_event_time{iday}.event;
     
    % spikes in that region
    rm = [];
    for i_u = 1:numel(data_spike)
        r_name = data_spike{i_u}.name{1};
        if strcmp(r_name,region_spike) ~= 1
           rm = [rm i_u];
        end    
    end
    data_spike(rm) = [];
    if isempty(data_spike)
    continue
    end
      
    % lfp
    dataName_lfp = strcat(root,'lfp_',z.days{m.match_idx(iday)},'_targacq');
    data_lfp = D.lfp_all{m.match_idx(iday)};
    region1_lfp = data_lfp.lfp.only(region_lfp);
      
    % choice only
    region1_lfp_choice = region1_lfp.only(epoch);

    % pre only
    if drug == 0
        if region1_lfp_choice.contains('oxytocin') || region1_lfp_choice.contains('saline')
           if region1_lfp_choice.contains('pre')
              region1_lfp_choice = region1_lfp_choice.only('pre');
        
           end   
        end   
    else 
        % do it pre and post 
        if region1_lfp_choice.contains(session) == 0
        continue
        else
           region1_lfp_choice = region1_lfp_choice.only(session);   
        end
    end
    
    % block match
    % make sure they have same block number 
    b1 = region1_lfp_choice('blocks');
    b2 = data_event('blocks');
    if numel(b1) ~= numel(b2)
       if size(b1,1) > size(b2,1)
          minb = b2;
       else
          minb = b1; 
       end
       region1_lfp_choice=region1_lfp_choice.only(minb);      
       data_event=data_event.only(minb);
    end       
        
    % method
    method = 'spike-field';
    cnt2 = cnt2 + numel(data_spike);
    cnt3 = cnt3 + numel(region1_lfp_choice('channels'));
    [coher_data, ~, phase_data, info] = coher_analysis_sua(region1_lfp_choice,[],method,data_spike,data_event,iday,[],epoch,event_key);
          coher_data_all{iday} = coher_data;
%         phase_data_all{iday} = phase_data;

end

spike_unit_n = cnt2
field_unit_n = cnt3

