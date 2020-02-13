
%% plot results from spike field coherence analysis

% first to load the data or run nnfigs3_sfc.m

h = load('/Users/chengchichu/Desktop/dictator_analysis/coher/mat/f_handle');
f_resolution = length(h.f);
t2 = -0.5:0.05:0.5;
t2s = length(t2); % axis for computing
t3 = -0.3:0.05:0.3; % axis for plotting
label_flip = 0; 
% 
self_all = zeros(f_resolution,t2s);
both_all = self_all;
other_all = self_all;
none_all = self_all;
anti_all = self_all;
pro_all = self_all;
a_p_all = self_all;
anti_data = {};
pro_data = {};
a_p_data = {};
cnt = 0;
ctr2 = [];
ctr3 = {};
day_cnt = [];

% main loop begins
for iday = 1:numel(coher_data_all)
    iday
    if isempty(coher_data_all{iday})
    continue
    end
    % remove day 0527, which is post drug session
    if iday == 33 || iday == 67
    continue
    end   
    pairs = {};
    rmp = []; % remove pairs with too few spikes to reliably estimate sfc
    for i_c = 1:4
        pairdata = coher_data_all{iday}{i_c};
        for i_p = 1:numel(pairdata) 
            coherm = cell2mat(cellfun(@nanmean,pairdata{i_p}.C,repmat({2},[1 t2s]),'UniformOutput',false));         
            if sum(sum(isnan(coherm)==1)) > 1
               rmp = [rmp i_p];         
            end  
            pairs{i_c}{i_p} = coherm;
        end   
    end
    ctr2 = [ctr2 numel(pairdata)];
    rmp = unique(rmp); 
    ctr3{iday} = setdiff(1:numel(pairdata),rmp);
    
    for i_c = 1:4 % 
        pairs{i_c}(rmp) = [];
    end
    remain_size = numel(pairdata)-length(rmp);
    day_cnt = [day_cnt remain_size];
   
    for i_p = 1:remain_size
       
        cnt = cnt + 1;
        
        self = pairs{1}{i_p};
        both = pairs{2}{i_p};
        other = pairs{3}{i_p};
        none = pairs{4}{i_p};     
        
        if label_flip == 0
            anti = pairs{1}{i_p} - pairs{2}{i_p};
            pro = pairs{3}{i_p} - pairs{4}{i_p};
        else    
            pro = pairs{2}{i_p} - pairs{1}{i_p};
            anti = pairs{4}{i_p} - pairs{3}{i_p};
        end
        
        a_p = pro - anti;
        
        self_all = self_all + self;
        both_all = both_all + both;
        other_all = other_all + other;
        none_all = none_all + none;
        anti_all = anti_all + anti;
        pro_all = pro_all + pro;
        a_p_all = a_p_all + a_p;
        a_p_data{cnt} = a_p;
        anti_data{cnt} = anti;
        pro_data{cnt} = pro;
    end    
end

% averaging across pairs and smoothing
self_all = imgaussfilt(self_all/cnt,2);
both_all = imgaussfilt(both_all/cnt,2);
other_all = imgaussfilt(other_all/cnt,2);
none_all = imgaussfilt(none_all/cnt,2);
anti_all = imgaussfilt(anti_all/cnt,2);
pro_all = imgaussfilt(pro_all/cnt,2);
a_p_all = imgaussfilt(a_p_all/cnt,2);

% crop the range for showing the data
f_idc = find(h.f<105); % we plot anything under 101hz
t_idc = find(t2>=-0.3 & t2<=0.3);
self_all = self_all(f_idc,t_idc);
both_all = both_all(f_idc,t_idc);
other_all = other_all(f_idc,t_idc);
none_all = none_all(f_idc,t_idc);
anti_all = anti_all(f_idc,t_idc);
pro_all = pro_all(f_idc,t_idc);
a_p_all = a_p_all(f_idc,t_idc);

fh = h.f(f_idc);
figure(1),clf,set(1,'Position',[-1331 319 985 356])
subplot(221)
imagesc(t3,fh,self_all,[0.675 0.695]),set(gca,'YDir','normal'), colormap(jet), colorbar
title('self')
ylim([10 80])
subplot(222)
imagesc(t3,fh,both_all,[0.675 0.695]),set(gca,'YDir','normal'), colormap(jet), colorbar
title('both')
ylim([10 80])
subplot(223)
imagesc(t3,fh,other_all,[0.675 0.695]),set(gca,'YDir','normal'), colormap(jet), colorbar
title('other')
ylim([10 80])
subplot(224)
imagesc(t3,fh,none_all,[0.675 0.695]),set(gca,'YDir','normal'), colormap(jet), colorbar
title('none')
ylim([10 80])

figure(2),clf
subplot(211)
imagesc(t3,fh,anti_all,[-0.016 0.016]),set(gca,'YDir','normal'), colormap(jet), colorbar
title('anti')
ylim([10 80])
subplot(212)
imagesc(t3,fh,pro_all,[-0.016 0.016]),set(gca,'YDir','normal'), colormap(jet), colorbar
title('pro')
ylim([10 80])

figure(3),clf
imagesc(t3,fh,a_p_all,[-0.016 0.016]),set(gca,'YDir','normal'), h = colormap(jet), colorbar
title('pro-anti')
ylim([10 80])

h = colorbar;
ylabel(h, 'coherence difference')

