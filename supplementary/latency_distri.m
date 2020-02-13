

% run spike_stats.m to get latency distri
% 

new_t = -0.3:0.01:0.3;
bins = 1:5:length(new_t);
cnt = histc(self_L,bins);
cnt2 = histc(both_L,bins);
cnt3 = histc(other_L,bins);
cnt4 = histc(none_L,bins);

figure(1),clf
subplot(231), bar(new_t(bins),cnt), title('self latency')
subplot(232), bar(new_t(bins),cnt2), title('both latency')
subplot(233), bar(new_t(bins),cnt3), title('other latency')
subplot(234), bar(new_t(bins),cnt4), title('none latency')
subplot(235)
plot(new_t(bins),cumsum(cnt)/sum(cnt),'r'); hold on
plot(new_t(bins),cumsum(cnt2)/sum(cnt2),'g'); hold on
plot(new_t(bins),cumsum(cnt3)/sum(cnt3),'b'); hold on
plot(new_t(bins),cumsum(cnt4)/sum(cnt4),'k'); hold on
legend('self','both','other','none')
