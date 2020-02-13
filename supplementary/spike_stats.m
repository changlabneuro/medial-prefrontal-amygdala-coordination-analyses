
% run get_spike_trial first, set choice_tag = "choice" or "cued" 

% anova reward vs cond
region = 'bla';
epoch = 'cued'; % cued, choice
fname = 'targon'; % targon, targacq
data = spike_trials_sua_all.(fname); 
tt = 1:numel(data); % total unit number, ACC:1~253 BLA:254~343
dbase = spike_trials_sua_all.cueon;
label = spike_trials_sua_all.L; % rw magnitude
data_rw = spike_trials_sua_all.rwon; 

% bla or acc
acc_idx = find(tt<=253);
if strcmp(region,'acc')
   data = data(1:length(acc_idx));
   label = label(1:length(acc_idx));
   dbase = dbase(1:length(acc_idx));
   data_rw = data_rw(1:length(acc_idx));
else
   data = data(length(acc_idx)+1:length(tt));
   label = label(length(acc_idx)+1:length(tt));
   dbase = dbase(length(acc_idx)+1:length(tt));
   data_rw = data_rw(length(acc_idx)+1:length(tt));
end    

%% 
t_rng = [-0.55 0.55]; % for spike cnt
tH = -0.5:0.01:0.5; % general time axis
tw_idx = find(tH>=-0.3 & tH<=0.3);
pM = zeros(numel(data),length(tH)); % test  recieved vs forgone
pM_anti = pM; % test  self vs both (baseline subtracted)
pM_pro = pM;  % test  other vs none (baseline subtracted)
pM_outcome = []; 
pM_outcome.self = pM; % test self vs baseline
pM_outcome.both = pM; % test both vs baseline
pM_outcome.other = pM; % test other vs baseline
pM_outcome.none = pM; % test none vs baseline
pM_outcome.Or = pM; % any of condition sig from baseline
p_self_other = []; % how many units sig between self vs other in either 
% choice or rw period
spk_choice = {};
spk_rw = {};
received = [];
forgone = [];
% latency 
self_L = []; both_L = self_L; other_L = self_L; none_L = self_L; 
p_self_all = []; p_both_all = []; p_other_all = []; p_none_all = [];
for i_u = 1:numel(data)   
    i_u
    yM  = []; % matrix [selfcnt bothcnt]
    B = {}; % matrix [selfcntBaseline bothcntBaseline]
    B2 = {};
    yM2 = []; % matrix [othercnt nonecnt]
    ANTIM = {}; 
    PROM = {};
    for i_c = 1:2 % self and both
        [T, ~, cnt, ~] = estimate_rate_bins(data{i_u}{i_c}, [100 10]/1000, t_rng(1), t_rng(2));   
        yM = [yM; cnt];
        ANTIM{i_c} = cnt;       
        [~, ~, bcnt, ~] = estimate_rate_bins(dbase{i_u}{i_c}, [100 10]/1000, -0.15, 0); 
        B{i_c} = mean(bcnt,2);
    end
    for i_c = 3:4 % other and none
        [~, ~, cnt2, ~] = estimate_rate_bins(data{i_u}{i_c}, [100 10]/1000, t_rng(1), t_rng(2));   
        size(cnt2);
        yM2 = [yM2; cnt2];
        PROM{i_c-2} = cnt2;
        [~, ~, bcnt2, ~] = estimate_rate_bins(dbase{i_u}{i_c}, [100 10]/1000, -0.15, 0); 
        B2{i_c-2} = mean(bcnt2,2);
    end    
    
    % choice by outcome 
    yM3 = {}; yM3_a = {};
    for i_c = 1:4
        [~, ~, cnt3, ~] = estimate_rate_bins(data{i_u}{i_c}, [100 10]/1000, t_rng(1), t_rng(2));   
        yM3{i_c} = cnt3;
        [~, ~, cnt3a, ~] = estimate_rate_bins(data{i_u}{i_c}, 0.15, 0, 0.15);   
        yM3_a{i_c} = cnt3a;
        spk_choice{i_u} = yM3_a;
    end
  
    % rw by outcome
    yM4 = {};
    for i_c = 1:4
        [~, ~, cnt4, ~] = estimate_rate_bins(data_rw{i_u}{i_c}, 0.4, 0.05, 0.45);   
        yM4{i_c} = cnt4;
        spk_rw{i_u} = yM4;
    end

    % choice and rw sig distiguish between self vs other
    t_1 = find(tH>=0 & tH<0.15); % for choice
    t_2 = find(tH>=0.05 & tH<0.45); % for reward
    p_self_other_choice = ranksum(yM3_a{1}, yM3_a{3});
    p_self_other_rw = ranksum(yM4{1}, yM4{3});
    p_self_other(i_u) = p_self_other_choice<0.05 | p_self_other_rw<0.05;
    
    % 
    p_self = [];
    p_both = [];
    p_other = [];
    p_none = [];
    ps = [];
    ps_anti = [];
    ps_pro = [];
%     ps_ctx = [];
    l_self = zeros(1,101); l_both = l_self; l_other = l_self; l_none = l_self;
    for t = 1:length(T)
        ps(t) = ranksum(yM(:,t),yM2(:,t));
        ps_anti(t) = ranksum(ANTIM{1}(:,t)-B{1},ANTIM{2}(:,t)-B{2});
        ps_pro(t) = ranksum(PROM{1}(:,t)-B2{1},PROM{2}(:,t)-B2{2});
        p_self(t) = ranksum(yM3{1}(:,t),B{1});
        p_both(t) = ranksum(yM3{2}(:,t),B{2});
        p_other(t) = ranksum(yM3{3}(:,t),B2{1});
        p_none(t) = ranksum(yM3{4}(:,t),B2{2});
    end 
    
    p_self(find(isnan(p_self))) = 1; % if that bin has no spikes, ranksum return NaN
    p_both(find(isnan(p_both))) = 1;
    p_other(find(isnan(p_other))) = 1;
    p_none(find(isnan(p_none))) = 1;
    
    p_self_all = [p_self_all; p_self(tw_idx)];
    p_both_all = [p_both_all; p_both(tw_idx)];
    p_other_all = [p_other_all; p_other(tw_idx)];
    p_none_all = [p_none_all; p_none(tw_idx)];
    
    pM(i_u,:) = ps<0.05;
    pM_anti(i_u,:) = ps_anti<0.05;
    pM_pro(i_u,:) = ps_pro<0.05;
    pM_outcome.self(i_u,:) = p_self<0.05;
    pM_outcome.both(i_u,:) = p_both<0.05;
    pM_outcome.other(i_u,:) = p_other<0.05;
    pM_outcome.none(i_u,:) = p_none<0.05;

    pOr = (p_self<0.05)+(p_both<0.05)+(p_other<0.05)+(p_none<0.05);
    
    pM_outcome.Or(i_u,:) = pOr>0;
    
    % find latency
%     tw_idx = find(tH>=-0.3 & tH<=0.3);
    p_self_ = p_self(tw_idx);
    p_both_ = p_both(tw_idx);
    p_other_ = p_other(tw_idx);
    p_none_ = p_none(tw_idx);
    self_latency_idx = find_consective_bins(find(p_self_<0.05),3); % at least 3 consecutive bins
    both_latency_idx = find_consective_bins(find(p_both_<0.05),3);
    other_latency_idx = find_consective_bins(find(p_other_<0.05),3);
    none_latency_idx = find_consective_bins(find(p_none_<0.05),3);
   
    if isempty(self_latency_idx) == 0
       self_L(i_u) = self_latency_idx(1);
    else
       self_L(i_u) = NaN;
    end    
    if isempty(both_latency_idx) == 0
       both_L(i_u) = both_latency_idx(1);
    else
       both_L(i_u) = NaN;
    end    
    if isempty(other_latency_idx) == 0
       other_L(i_u) = other_latency_idx(1);
    else
       other_L(i_u) = NaN;
    end  
    if isempty(none_latency_idx) == 0
       none_L(i_u) = none_latency_idx(1);
    else
       none_L(i_u) = NaN;
    end  
end      

%% 

figure(8),clf
plot(T,smooth(100*mean(pM)))
ylabel('percent significant')
xlabel(fname)
lim = max(100*mean(pM));
ylim([0 lim])

name = strcat(region,'_',epoch,'_',fname,'_difference anti-pro allunits - base smooth')
% pirntfigs('/eps',8,name)

figure(9), clf
h = plot(T,smooth(100*mean(pM_anti))), hold on
h2 = plot(T,smooth(100*mean(pM_pro)))
ylabel('percent significant')
legend('anti','pro')
xlabel(fname)

name = strcat(region,'_',epoch,'_',fname,'_difference anti&&pro allunits - base smooth')
% pirntfigs('/eps',9,name)

figure(10), clf
h1 = plot(T,smooth(100*mean(pM_outcome.self),10),'r'); hold on
h2 = plot(T,smooth(100*mean(pM_outcome.both),10),'b'); hold on
h3 = plot(T,smooth(100*mean(pM_outcome.other),10),'g'); hold on
h4 = plot(T,smooth(100*mean(pM_outcome.none),10),'m');
legend('self','both','other','none')
ylabel('percent significant')
xlabel(fname)
name = strcat(region,'_',epoch,'_',fname,'_difference outcome allunits compare to baseline smooth10')
% pirntfigs('/eps',10,name)

figure(11), clf
h = plot(T,smooth(100*mean(pM_outcome.Or),10),'m');
legend('ANY')
ylabel('percent significant')
xlabel(fname)
name = strcat(region,'_',epoch,'_',fname,'_any outcome sig allunits compare to baseline smooth10')
% pirntfigs('/eps',11,name)


