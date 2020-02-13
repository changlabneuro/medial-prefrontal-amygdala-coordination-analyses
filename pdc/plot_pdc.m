
% load saved data
% load('/Volumes/ccdisk2/isDiranalysisData/sivPDCabs12data');
% ivPDCs = ivgoPDCs;

%% 

tt = -0.5:0.05:0.5;
% time window for stats, 150ms after choice
t_idx = find(tt>=0 & tt<=0.15); 

sbPDC_line_dir1 = [];
sbPDC_line_dir2 = [];
onPDC_line_dir1 = [];
onPDC_line_dir2 = [];

selfPDC_ = []; bothPDC_ = []; otherPDC_ = []; nonePDC_ = []; sbPDC_ = []; onPDC_ = [];
for d_t = 1:length(tt)
    cnt = 0;
    selfPDC = zeros(2,2,Nf);
    bothPDC = selfPDC;
    otherPDC = selfPDC;
    nonePDC = selfPDC; sbdir1 = []; sbdir2 = []; ondir1 = []; ondir2 = [];
    for iday = 1:numel(ivPDCs)
        iday;
        for i_p = 1:numel(ivPDCs{iday}{1})
            i_p;
            cnt = cnt+1;
            selfPDC = selfPDC+abs(ivPDCs{iday}{1}{i_p}{d_t});
            bothPDC = bothPDC+abs(ivPDCs{iday}{2}{i_p}{d_t});
            otherPDC = otherPDC+abs(ivPDCs{iday}{3}{i_p}{d_t});
            nonePDC = nonePDC+abs(ivPDCs{iday}{4}{i_p}{d_t});
            sb_p = abs(ivPDCs{iday}{1}{i_p}{d_t})-abs(ivPDCs{iday}{2}{i_p}{d_t});
            on_p = abs(ivPDCs{iday}{3}{i_p}{d_t})-abs(ivPDCs{iday}{4}{i_p}{d_t});
            sbdir1 = [sbdir1 squeeze(sb_p(1,2,:))];
            sbdir2 = [sbdir2 squeeze(sb_p(2,1,:))];
            ondir1 = [ondir1 squeeze(on_p(1,2,:))];
            ondir2 = [ondir2 squeeze(on_p(2,1,:))];
        end
    end  
    sbPDC_line_dir1 = cat(3,sbPDC_line_dir1,sbdir1);
    sbPDC_line_dir2 = cat(3,sbPDC_line_dir2,sbdir2);
    onPDC_line_dir1 = cat(3,onPDC_line_dir1,ondir1);
    onPDC_line_dir2 = cat(3,onPDC_line_dir2,ondir2);
    selfPDCavg = selfPDC/cnt;
    bothPDCavg = bothPDC/cnt;
    otherPDCavg = otherPDC/cnt;
    nonePDCavg = nonePDC/cnt;
    sbPDC = selfPDCavg-bothPDCavg;
    onPDC = otherPDCavg-nonePDCavg;
    selfPDC_ = cat(4,selfPDC_,selfPDCavg);
    bothPDC_= cat(4,bothPDC_,bothPDCavg);
    otherPDC = cat(4,otherPDC_,otherPDCavg);
    nonePDC = cat(4,nonePDC_,nonePDCavg);
    sbPDC_ = cat(4,sbPDC_,sbPDC);
    onPDC_ = cat(4,onPDC_,onPDC);
end

% plot the results
figure(1),clf
f_size = 2; % 2D gaussin smoothing
clims = [-0.005 0.003];
subplot(211)
imagesc(tt,1:100,imgaussfilt(squeeze(sbPDC_(1,2,:,:)),f_size),clims)
set(gca,'YDir','normal'), colormap(jet), colorbar, ylim([10 80]),xlim([-0.3 0.3])
title('anti BLA>ACC')
subplot(212)
imagesc(tt,1:100,imgaussfilt(squeeze(sbPDC_(2,1,:,:)),f_size),clims)
set(gca,'YDir','normal'), colormap(jet), colorbar, ylim([10 80]),xlim([-0.3 0.3])
title('anti ACC>BLA')
figure(2),clf
subplot(211)
imagesc(tt,1:100,imgaussfilt(squeeze(onPDC_(1,2,:,:)),f_size),clims)
set(gca,'YDir','normal'), colormap(jet), colorbar, ylim([10 80]),xlim([-0.3 0.3])
title('pro BLA>ACC')
subplot(212)
imagesc(tt,1:100,imgaussfilt(squeeze(onPDC_(2,1,:,:)),f_size),clims)
set(gca,'YDir','normal'), colormap(jet), colorbar, ylim([10 80]),xlim([-0.3 0.3])
title('pro ACC>BLA')

