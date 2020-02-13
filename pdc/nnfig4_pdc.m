
% this script was used to compute the partial directed coherence results in the figure 4 of NN paper.
% Specialized medial prefrontal-amygdala coordination in other-regarding
% decision preference, Copyright 2020 chengchi chu

% I adopted the code from Omidvarnia, A. H. et al., 2013 which is a
% modified version of classical PDC (generalized orthogonal PDC)
% 

clear all
% add_depends
addpath(genpath('/gpfs/milgram/home/cc2586/global'))
addpath(genpath('/gpfs/milgram/home/cc2586/dsp'))
addpath('/gpfs/milgram/home/cc2586/GPDC/artool')

% load lfp data
Data = load('/gpfs/milgram/project/chang/CHANG_LAB/cc2586/cued_dt_LFP.mat');

epoch = 'cued'; % 'choice'
varname = 'var'; % 'lfp'
name = 'cued_sivPDCabs12data.mat'; % 

outcomes = {'self', 'both', 'other', 'none'};

% selecting subsets of pairs from all pairs/avoid using redundant pairs
v = load('pairs');
% labels for recording sessions
z = load('days_targacq');

% use fixed model order, coumpte a range of model order by using arfit.m 
model_order = 12;
Fs = 1000;
Nf = 100;
Fmax = 100;
slide_t = -0.5:0.05:0.5;
t_ = -0.5:0.001:0.65-0.001;

% main loop
ivgoPDCs = {};
for iday = 1:numel(z.days)
    iday
    data_lfp_choice = Data.lfp_all{iday}.(varname).only(epoch); 
    acc_lfp =  data_lfp_choice.only('acc');
    bla_lfp =  data_lfp_choice.only('bla');
    
    % use pre-drug administration epoch
    if acc_lfp.contains('oxytocin') || acc_lfp.contains('saline')
       if acc_lfp.contains('pre')
          acc_lfp = acc_lfp.only('pre');
       end
    end
    if bla_lfp.contains('oxytocin') || bla_lfp.contains('saline')
       if bla_lfp.contains('pre')
          bla_lfp = bla_lfp.only('pre');
       end
    end       

    pairs = v.pairs.channels{iday};
    gout = {};
    for i_c = 1:4
        acc_cond_lfp = acc_lfp.only(outcomes{i_c});
        bla_cond_lfp = bla_lfp.only(outcomes{i_c});
        pdc_result = {}; gpdc_result = {};
        for i_p = 1:size(pairs,1)   
            acc_cond_lfp_per = acc_cond_lfp.only(pairs{i_p,2});
            bla_cond_lfp_per = bla_cond_lfp.only(pairs{i_p,1});
            for d_t = 1:length(slide_t)              
                t_idx = find(t_>= slide_t(d_t) & t_< slide_t(d_t)+0.15); % 150ms sliding time window    
                acc_ = acc_cond_lfp_per.data(:,t_idx); %  
                bla_ = bla_cond_lfp_per.data(:,t_idx);
                channel_data = cat(3,acc_,bla_);

                gpdc_avg = zeros(2,2,Nf); % nch = 2
                for i_t = 1:size(channel_data,1)
                    y = squeeze(channel_data(i_t,:,:));
                    [w, A, C, sbc, fpe, th] = arfit(y, model_order, model_order);  
    %                 [~,p_opt] = min(sbc)
                    p_opt = model_order;
                    [GPDC,OPDC,PDC,GOPDC,S] = PDC_dDTF_imag_cc(A,C,p_opt,Fs,Fmax,Nf);
                    gpdc_avg = gpdc_avg + abs(GOPDC);
                end
                gpdc_avg = gpdc_avg/size(channel_data,1); % divide by how many trials
                gpdc_result{i_p}{d_t} = gpdc_avg;
            end
        end
        gout{i_c} = gpdc_result;
    end
    ivgoPDCs{iday} = gout;
end

save_name = strcat('/gpfs/milgram/project/chang/CHANG_LAB/cc2586/',name);
save(save_name,'ivgoPDCs')
