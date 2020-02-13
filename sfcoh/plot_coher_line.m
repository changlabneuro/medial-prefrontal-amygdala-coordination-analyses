
%% plot results from spike field coherence analysis, either in time or frequency dimension across pairs

side1 = load('/Volumes/ccDrive3/oxy_coher_data/sfc_pre_bla_spike_all_withsame.mat');
addpath('/Users/chuchengchi/Desktop/dictator_analysis/coher/mat');
coher_data_all = side1.coher_data_all;
plot_sfc;

side1_pro = pro_data;
side1_anti = anti_data;
side1_ap = a_p_data;

side2 = load('/Volumes/ccDrive3/oxy_coher_data/sfc_pre_acc_spike_all_withsame.mat');
coher_data_all = side2.coher_data_all;
plot_sfc;

side2_pro = pro_data;
side2_anti = anti_data;
side2_ap = a_p_data;

%% 

dir1_data = side1_ap;

dir2_data = side2_ap;

%% plot coher line by time

% frequency axis 
fh = load('f_handle.mat');
tt = -0.5:0.05:0.5;
t_idx = find(tt>=0 & tt<=0.15);
beta_fr_rng = [15 25];
betafr_dx = find(fh.f>beta_fr_rng(1) & fh.f<beta_fr_rng(2));
gamma_fr_rng = [45 70];
gammafr_dx = find(fh.f>gamma_fr_rng(1) & fh.f<gamma_fr_rng(2));
f_w = gammafr_dx; % or betafr_dx

% get data
sfc1 = []; sfc2 = [];
for t = 1:length(tt)
    xs_dir1 = [];
    for ipair = 1:numel(dir1_data)
        smoothed_data = imgaussfilt(dir1_data{ipair},2);
        x = mean(smoothed_data(f_w,t),1);
        xs_dir1 = [xs_dir1 x];
    end
    sfc1 = [sfc1; xs_dir1];
    xs_dir2 = [];
    for ipair = 1:numel(dir2_data)
        smoothed_data = imgaussfilt(dir2_data{ipair},2);
        x = mean(smoothed_data(f_w,t),1);
        xs_dir2 = [xs_dir2 x];
    end
    sfc2 = [sfc2; xs_dir2];       
end    

% plot
figure(1),clf
patcherrorbar(tt,sfc1');
patcherrorbar(tt,sfc2');
xlim([-0.3 0.3])
p1 = []; p2 = []; p3 = [];
for t = 1:length(tt)
    p1(t) = signrank(sfc1(t,:));
    p2(t) = signrank(sfc2(t,:));
    p3(t) = ranksum(sfc1(t,:),sfc2(t,:));
end
% put on * for significance
v1 = ones(1,length(tt));
v1(find(p1>0.05)) = NaN;
% the position for putting *
y1 = 0;
plot(tt,y1*ones(1,length(tt)).*v1,'*r');

mean_val = mean(mean(sfc1(t_idx,:)));
mean_val2 = mean(mean(sfc2(t_idx,:)));
% stats between 2 conditions
ranksum(mean(sfc1(t_idx,:),1),mean(sfc2(t_idx,:),1))
 
%% plot coher line by frequency

ff_idx = find(fh.f<105);

sfc1 = []; sfc2 = [];
xs_dir1 = [];
for ipair = 1:numel(dir1_data)
    smoothed_data = imgaussfilt(dir1_data{ipair},2);
    x = mean(smoothed_data(:,t_idx),2);
    xs_dir1 = [xs_dir1 x];
end
sfc1 = [sfc1; xs_dir1];
xs_dir2 = [];
for ipair = 1:numel(dir2_data)
    smoothed_data = imgaussfilt(dir2_data{ipair},2);
    x = mean(smoothed_data(:,t_idx),2);
    xs_dir2 = [xs_dir2 x];
end
sfc2 = [sfc2; xs_dir2];   
    
% plot
figure(1),clf
h1 = patcherrorbar(fh.f(ff_idx),sfc1(ff_idx,:)');
xlim([10 80])
legend('ACCspikeBLAfield')
h2 = patcherrorbar(fh.f(ff_idx),sfc2(ff_idx,:)');
xlim([10 80])
legend('BLAspikeBLAfield')

beta_fr_rng = [15 25];
betafr_dx = find(fh.f>beta_fr_rng(1) & fh.f<beta_fr_rng(2));
p = signrank(mean(sfc1(betafr_dx,:),1))

gamma_fr_rng = [45 70];
gammafr_dx = find(fh.f>gamma_fr_rng(1) & fh.f<gamma_fr_rng(2));
p = signrank(mean(sfc1(gammafr_dx,:),1))

beta_fr_rng = [15 25];
betafr_dx = find(fh.f>beta_fr_rng(1) & fh.f<beta_fr_rng(2));
p = signrank(mean(sfc2(betafr_dx,:),1))

gamma_fr_rng = [45 70];
gammafr_dx = find(fh.f>gamma_fr_rng(1) & fh.f<gamma_fr_rng(2));
p = signrank(mean(sfc2(gammafr_dx,:),1))


