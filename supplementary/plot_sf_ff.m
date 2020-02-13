

% plot spike-field and field-field coherence

%% load field-field data
nf_ff = load('/Users/chengchichu/Downloads/combined2.mat'); % free choice

self_ff = nf_ff.cont.only('self');
both_ff = nf_ff.cont.only('both');
other_ff = nf_ff.cont.only('other');
none_ff = nf_ff.cont.only('none');

side_ff_anti = self_ff.data - both_ff.data;
side_ff_pro = other_ff.data - none_ff.data;
side_ff_anti2 = none_ff.data - other_ff.data;
side_ff_pro2 = both_ff.data - self_ff.data;
side_ff_ap = side_ff_pro - side_ff_anti;

%% load spike field data
side1 = load('/Volumes/ccDrive3/oxy_coher_data/sfc_pre_acc_spike_all_withsame.mat');
coher_data_all = side1.coher_data_all;
addpath('/Users/chengchichu/Desktop/dictator_analysis/coher/mat')
plot_sfc;

side1_pro = pro_data;
side1_anti = anti_data;
side1_ap = a_p_data;

side2 = load('/Volumes/ccDrive3/oxy_coher_data/sfc_pre_bla_spike_all_withsame.mat');
coher_data_all = side2.coher_data_all;
plot_sfc;

side2_pro = pro_data;
side2_anti = anti_data;
side2_ap = a_p_data;

%% common parameter
fh = load('/Users/chengchichu/Desktop/dictator_analysis/coher/mat/f_handle.mat');

beta_fr_rng = [15 25];
betafr_dx = find(fh.f>beta_fr_rng(1) & fh.f<beta_fr_rng(2));

gamma_fr_rng = [45 70];
gammafr_dx = find(fh.f>gamma_fr_rng(1) & fh.f<gamma_fr_rng(2));

ff_idx = find(fh.f<105);
tt = -0.5:0.05:0.5;
t_idx = find(tt>=0 & tt<=0.15);

%% plot by frequency

%s-f
dir1_data = side1_ap;
dir2_data = side2_ap;

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

%f-f
xs_dirff = [];
for ipair = 1:size(side_ff_ap,1)
  smoothed_data = imgaussfilt(squeeze(side_ff_ap(ipair,:,:)),2);
  x = mean(smoothed_data(:,t_idx),2);
  xs_dirff = [xs_dirff x];
end    
ff3 = xs_dirff;

%% plot together

figure(1),clf
h1 = patcherrorbar(fh.f(ff_idx),sfc1(ff_idx,:)');
h2 = patcherrorbar(fh.f(ff_idx),sfc2(ff_idx,:)');
h3 = patcherrorbar(fh.f(ff_idx),ff3(ff_idx,:)');
legend([h1 h2 h3], 'ACC>BLA','BLA>ACC','FF')
xlim([10 80])

% add stats
p_sf_ff1 = []; p_sf_ff2 = [];
for i = 1:length(fh.f(ff_idx))
    p_sf_ff1(i) = ranksum(sfc1(i,:), ff3(i,:));
    p_sf_ff2(i) = ranksum(sfc2(i,:), ff3(i,:));
end

%% by time

f_w = betafr_dx; % or gammafr_dx

%s-f
dir1_data = side1_ap;
dir2_data = side2_ap;

sfc1 = []; sfc2 = [];
xs_dir1 = [];
for ipair = 1:numel(dir1_data)
    smoothed_data = imgaussfilt(dir1_data{ipair},2);
    x = mean(smoothed_data(f_w,:),1);
    xs_dir1 = [xs_dir1; x];
end
sfc1 = [sfc1; xs_dir1];  

xs_dir2 = [];
for ipair = 1:numel(dir2_data)
    smoothed_data = imgaussfilt(dir2_data{ipair},2);   
    x = mean(smoothed_data(f_w,:),1);
    xs_dir2 = [xs_dir2; x];
end
sfc2 = [sfc2; xs_dir2];    

%f-f
xs_dirff = [];
for ipair = 1:size(side_ff_ap,1)
    smoothed_data = imgaussfilt(squeeze(side_ff_ap(ipair,:,:)),2);
    x = mean(smoothed_data(f_w,:),1);
    xs_dirff = [xs_dirff; x];
end    
ff3 = xs_dirff;

