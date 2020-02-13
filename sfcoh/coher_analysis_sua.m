function [coher_data, f, phase_data, info] = coher_analysis_sua(rs1,rs2,method,data_spike,data_event,iday,pair_idc,epoch,event_key)

% spike field coherence analysis

% condition
cs = {'self','both','other','none'};

% parameters
% sliding window
t2 = -0.5:0.05:0.5;
params.tapers = [1.5 2];
params.Fs = [1000];
dt = 1/[params.Fs];
slide_win_size = 150;
step_win_size = 50;

info = {};
if strcmp(method, 'spike-field') 
   
   % lfp
   region1_lfp_choice = rs1;    
   
   % event
   event_label = data_event; 
   event_label = event_label.only(epoch);
   
   coher_data = {};     
   for i_c = 1:numel(cs)         
       event_cond = event_label.only(cs{i_c}); 
       event = event_cond.data;
       targacq = event_key; 
       event_ts = event(:,targacq);   
       nCh_lfp = region1_lfp_choice('channels'); 
       cnt = 0; % combination nLFP x nSUA
       for nL = 1:numel(nCh_lfp)
           lfp_per_ = region1_lfp_choice.only(nCh_lfp{nL}); % lfp per channel         
           lfp_cond_ = lfp_per_.only(cs{i_c}); % lfp per condition  
           lfp_cond_data_full = lfp_cond_.data;
           lfp_cond_data_full = lfp_cond_data_full'; % sample x trials
                                      
           % begin spike loop, under one lfp ch       
           for i_u = 1:numel(data_spike)
               spike_time_per_unit = data_spike{i_u}.data;                      
               Cs = {}; % coher across time
               Ps = {};
               for nt = 1:length(t2)                 
                   % spike part
                   spike_cond_data = [];                  
                   for i_trial = 1:length(event_ts)        
                       t_s = t2(nt)+event_ts(i_trial);
                       t_e = t2(nt)+[slide_win_size/1000-dt]+event_ts(i_trial);
                       spike_idx = find(spike_time_per_unit> t_s & spike_time_per_unit<= t_e);
                       spike_time_per_trial = spike_time_per_unit(spike_idx);             
                       c_spike_time = [spike_time_per_trial - t_s]; % converted spike time, align to start of sliding window     
                       spike_cond_data(i_trial).time = c_spike_time';     
                   end                         
                   % lfp part
                   ws = step_win_size*(nt-1)+1; % sample window for lfp
                   we = ws+slide_win_size-1;
                   lfp_cond_data = lfp_cond_data_full(ws:we,:); 
                   % coher fun
                   [C,phi,~,~,~,f] = coherencycpt(lfp_cond_data,spike_cond_data,params);
                   Cs{nt} = C; % coher across time
                   Ps{nt} = phi;
               end     
               
               cnt = cnt+1; % count for nLFP x nSUA
               info{cnt} = {nCh_lfp{nL} data_spike{i_u}.name num2str(i_u)};
               coher_data{i_c}{cnt}.C = Cs;      
               phase_data{i_c}{cnt}.P = Ps;
           end        
       end       
    end
end

