clc;
close all;
clearvars;

%%%%%%%%%%%%%%%%%%%%%%%%%% SELECT NILEARN TIMESERIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0418/timeseries/
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0417/timeseries/
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0416/timeseries/
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0407_denoised/timeseries_hipass/
cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0405_denoised_GS/timeseries/
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_timeseries_bandpass/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load atlas labels
[~,labels] = readvars('atlas_labels_harvardoxford_cort-maxprob-thr25-2mm.csv');
labels = labels(2:49,:);
% Load Nilearn time-series per condition
load('ctrl_premanip_neut_timeseries.mat');
load('ctrl_premanip_heat_timeseries.mat');
load('ctrl_posmanip_neut_timeseries.mat');
load('ctrl_posmanip_heat_timeseries.mat');
load('medt_premanip_neut_timeseries.mat');
load('medt_premanip_heat_timeseries.mat');
load('medt_posmanip_neut_timeseries.mat');
load('medt_posmanip_heat_timeseries.mat');
% Permute dimensions (now: ROIs, timepts, sub[1-20])
ctrl_premanip_neut_corr_timeseries = permute(ctrl_premanip_neut_corr_timeseries,[3 2 1]);
ctrl_premanip_heat_corr_timeseries = permute(ctrl_premanip_heat_corr_timeseries,[3 2 1]);
ctrl_posmanip_neut_corr_timeseries = permute(ctrl_posmanip_neut_corr_timeseries,[3 2 1]);
ctrl_posmanip_heat_corr_timeseries = permute(ctrl_posmanip_heat_corr_timeseries,[3 2 1]);
medt_premanip_neut_corr_timeseries = permute(medt_premanip_neut_corr_timeseries,[3 2 1]);
medt_premanip_heat_corr_timeseries = permute(medt_premanip_heat_corr_timeseries,[3 2 1]);
medt_posmanip_neut_corr_timeseries = permute(medt_posmanip_neut_corr_timeseries,[3 2 1]);
medt_posmanip_heat_corr_timeseries = permute(medt_posmanip_heat_corr_timeseries,[3 2 1]);
% Separate into runs
ctrl_premanip_neut_1_timeseries = ctrl_premanip_neut_corr_timeseries(1:48,1:147,1:20);
ctrl_premanip_neut_2_timeseries = ctrl_premanip_neut_corr_timeseries(1:48,1:147,21:40);
clear ctrl_premanip_neut_corr_timeseries;
ctrl_premanip_heat_1_timeseries = ctrl_premanip_heat_corr_timeseries(1:48,1:147,1:20);
ctrl_premanip_heat_2_timeseries = ctrl_premanip_heat_corr_timeseries(1:48,1:147,21:40);
clear ctrl_premanip_heat_corr_timeseries;
ctrl_posmanip_neut_3_timeseries = ctrl_posmanip_neut_corr_timeseries(1:48,1:147,1:20);
ctrl_posmanip_neut_4_timeseries = ctrl_posmanip_neut_corr_timeseries(1:48,1:147,21:40);
clear ctrl_posmanip_neut_corr_timeseries;
ctrl_posmanip_heat_3_timeseries = ctrl_posmanip_heat_corr_timeseries(1:48,1:147,1:20);
ctrl_posmanip_heat_4_timeseries = ctrl_posmanip_heat_corr_timeseries(1:48,1:147,21:40);
clear ctrl_posmanip_heat_corr_timeseries;
medt_premanip_neut_1_timeseries = medt_premanip_neut_corr_timeseries(1:48,1:147,1:20);
medt_premanip_neut_2_timeseries = medt_premanip_neut_corr_timeseries(1:48,1:147,21:40);
clear medt_premanip_neut_corr_timeseries;
medt_premanip_heat_1_timeseries = medt_premanip_heat_corr_timeseries(1:48,1:147,1:20);
medt_premanip_heat_2_timeseries = medt_premanip_heat_corr_timeseries(1:48,1:147,21:40);
clear medt_premanip_heat_corr_timeseries;
medt_posmanip_neut_3_timeseries = medt_posmanip_neut_corr_timeseries(1:48,1:147,1:20);
medt_posmanip_neut_4_timeseries = medt_posmanip_neut_corr_timeseries(1:48,1:147,21:40);
clear medt_posmanip_neut_corr_timeseries;
medt_posmanip_heat_3_timeseries = medt_posmanip_heat_corr_timeseries(1:48,1:147,1:20);
medt_posmanip_heat_4_timeseries = medt_posmanip_heat_corr_timeseries(1:48,1:147,21:40);
clear medt_posmanip_heat_corr_timeseries;

%% Load VAS Pain scores
% get VASint+unp
[id	v6_int_h1 v6_unpl_h1 v6_int_h2 v6_unpl_h2 v6_int_h3 v6_unpl_h3 v6_int_h4 v6_unpl_h4] = readvars('pain_ratings.xlsx');
% split into ctrl/medt 1/run
vas_int_h1_ctrl = (v6_int_h1([1:20],:));
vas_int_h1_medt = (v6_int_h1([21:40],:));
vas_int_h2_ctrl = (v6_int_h2([1:20],:));
vas_int_h2_medt = (v6_int_h2([21:40],:));
vas_int_h3_ctrl = (v6_int_h3([1:20],:));
vas_int_h3_medt = (v6_int_h3([21:40],:));
vas_int_h4_ctrl = (v6_int_h4([1:20],:));
vas_int_h4_medt = (v6_int_h4([21:40],:));
vas_unp_h1_ctrl = (v6_unpl_h1([1:20],:));
vas_unp_h1_medt = (v6_unpl_h1([21:40],:));
vas_unp_h2_ctrl = (v6_unpl_h2([1:20],:));
vas_unp_h2_medt = (v6_unpl_h2([21:40],:));
vas_unp_h3_ctrl = (v6_unpl_h3([1:20],:));
vas_unp_h3_medt = (v6_unpl_h3([21:40],:));
vas_unp_h4_ctrl = (v6_unpl_h4([1:20],:));
vas_unp_h4_medt = (v6_unpl_h4([21:40],:));
clear v6*
% change in pain post-pre (mean of h3,h4 - mean of h1,h2)
vas_chg_int_ctrl = ((vas_int_h3_ctrl + vas_int_h4_ctrl)./2) - ((vas_int_h1_ctrl + vas_int_h2_ctrl)./2);
vas_chg_int_medt = ((vas_int_h3_medt + vas_int_h4_medt)./2) - ((vas_int_h1_medt + vas_int_h2_medt)./2);
vas_chg_unp_ctrl = ((vas_unp_h3_ctrl + vas_unp_h4_ctrl)./2) - ((vas_unp_h1_ctrl + vas_unp_h2_ctrl)./2);
vas_chg_unp_medt = ((vas_unp_h3_medt + vas_unp_h4_medt)./2) - ((vas_unp_h1_medt + vas_unp_h2_medt)./2);
% permute
vas_chg_int_ctrl = permute(vas_chg_int_ctrl,[2,1]);
vas_chg_int_medt = permute(vas_chg_int_medt,[2,1]);
vas_chg_unp_ctrl = permute(vas_chg_unp_ctrl,[2,1]);
vas_chg_unp_medt = permute(vas_chg_unp_medt,[2,1]);
clear vas_unp* vas_int*
%% Design Butterworth Filter (data bandpass filtered in Nilearn so unnecessary)

% Band-pass filter (0.08 - 0.12 Hz) using 6th order Butterworth filter (hi
% freq range bc Ries et al 2018, Wu et al 2008)

% % plot time series sample (1 roi, 1 subject)
% t = 1:1:147;
% roi = ctrl_posmanip_heat_1_timeseries(1,:,1);
% plot(t,roi)
% 
% % design a 6th order butterworth bandpass filter (0.08 - 0.12 Hz)
% % method 1
% [A,B,C,D] = butter(3,[0.08 0.12]/15); 
% sos = ss2sos(A,B,C,D);
% % method 2
% bpFilt = designfilt('bandpassiir', 'FilterOrder',6, ...
%     'HalfPowerFrequency1', 0.08, 'HalfPowerFrequency2', 0.12, ...
%     'DesignMethod', 'butter',...
%     'SampleRate',30);
% % compare freq responses
% fvt = fvtool(sos,bpFilt,'Fs',30);
% legend(fvt,'butter','designfilt')
% 
% % test filter and plot 
% dataIn = roi;
% dataOut = filter(bpFilt,dataIn);
% plot(t,dataOut)

%% Compute Shannon entropy for each run using wentropy (wavelet toolbox), take neg & log10
% S = wentropy(roi,'shannon');
% S = -S;
% logS = log10(S)
% S2 = log10(-wentropy(roi,'shannon'))

% preassign
ctrl_pos_h3_S = zeros(48,20);ctrl_pos_h4_S = zeros(48,20);
ctrl_pos_n3_S = zeros(48,20);ctrl_pos_n4_S = zeros(48,20);
ctrl_pre_h1_S = zeros(48,20);ctrl_pre_h2_S = zeros(48,20);
ctrl_pre_n1_S = zeros(48,20);ctrl_pre_n2_S = zeros(48,20);
medt_pos_h3_S = zeros(48,20);medt_pos_h4_S = zeros(48,20);
medt_pos_n3_S = zeros(48,20);medt_pos_n4_S = zeros(48,20);
medt_pre_h1_S = zeros(48,20);medt_pre_h2_S = zeros(48,20);
medt_pre_n1_S = zeros(48,20);medt_pre_n2_S = zeros(48,20);

for sub = 1:20
    for roi = 1:48
        ctrl_pos_h3_S(roi,sub) = log10(-wentropy(ctrl_posmanip_heat_3_timeseries(roi,:,sub),'shannon'));
        ctrl_pos_n3_S(roi,sub) = log10(-wentropy(ctrl_posmanip_neut_3_timeseries(roi,:,sub),'shannon'));
        ctrl_pre_h1_S(roi,sub) = log10(-wentropy(ctrl_premanip_heat_1_timeseries(roi,:,sub),'shannon'));
        ctrl_pre_n1_S(roi,sub) = log10(-wentropy(ctrl_premanip_neut_1_timeseries(roi,:,sub),'shannon'));
        ctrl_pos_h4_S(roi,sub) = log10(-wentropy(ctrl_posmanip_heat_4_timeseries(roi,:,sub),'shannon'));
        ctrl_pos_n4_S(roi,sub) = log10(-wentropy(ctrl_posmanip_neut_4_timeseries(roi,:,sub),'shannon'));
        ctrl_pre_h2_S(roi,sub) = log10(-wentropy(ctrl_premanip_heat_2_timeseries(roi,:,sub),'shannon'));
        ctrl_pre_n2_S(roi,sub) = log10(-wentropy(ctrl_premanip_neut_2_timeseries(roi,:,sub),'shannon'));
        medt_pos_h3_S(roi,sub) = log10(-wentropy(medt_posmanip_heat_3_timeseries(roi,:,sub),'shannon'));
        medt_pos_n3_S(roi,sub) = log10(-wentropy(medt_posmanip_neut_3_timeseries(roi,:,sub),'shannon'));
        medt_pre_h1_S(roi,sub) = log10(-wentropy(medt_premanip_heat_1_timeseries(roi,:,sub),'shannon'));
        medt_pre_n1_S(roi,sub) = log10(-wentropy(medt_premanip_neut_1_timeseries(roi,:,sub),'shannon'));
        medt_pos_h4_S(roi,sub) = log10(-wentropy(medt_posmanip_heat_4_timeseries(roi,:,sub),'shannon'));
        medt_pos_n4_S(roi,sub) = log10(-wentropy(medt_posmanip_neut_4_timeseries(roi,:,sub),'shannon'));
        medt_pre_h2_S(roi,sub) = log10(-wentropy(medt_premanip_heat_2_timeseries(roi,:,sub),'shannon'));
        medt_pre_n2_S(roi,sub) = log10(-wentropy(medt_premanip_neut_2_timeseries(roi,:,sub),'shannon'));
    end
end
clear sub roi *timeseries

%% dS (post-pre) per ROI per sub
S_chg_ctrl_heat = ((ctrl_pos_h3_S + ctrl_pos_h4_S)./2) - ((ctrl_pre_h1_S + ctrl_pre_h2_S)./2);
S_chg_ctrl_neut = ((ctrl_pos_n3_S + ctrl_pos_n4_S)./2) - ((ctrl_pre_n1_S + ctrl_pre_n2_S)./2);
S_chg_medt_heat = ((medt_pos_h3_S + medt_pos_h4_S)./2) - ((medt_pre_h1_S + medt_pre_h2_S)./2);
S_chg_medt_neut = ((medt_pos_n3_S + medt_pos_n4_S)./2) - ((medt_pre_n1_S + medt_pre_n2_S)./2);
clear medt* ctrl*

% SE for dS per ROI
SE_dS_ctrl_heat = std(S_chg_ctrl_heat,0,2) ./ sqrt(20);
SE_dS_ctrl_neut = std(S_chg_ctrl_neut,0,2) ./ sqrt(20);
SE_dS_medt_heat = std(S_chg_medt_heat,0,2) ./ sqrt(20);
SE_dS_medt_neut = std(S_chg_medt_neut,0,2) ./ sqrt(20);

%% plot change in local entropy per ROI averaged over all sub per condition 

ymin = -.1; ymax = .1;

figure
tiledlayout(2,2)
sgtitle('Change in ROI logS, Post-Pre, Group Level')

nexttile
y = mean(S_chg_ctrl_heat,2);
b = bar(y, 'FaceColor','g');
hold on
er = errorbar(y,SE_dS_ctrl_heat); er.Color = [0 0 0]; er.LineStyle = 'none';
hold off
b.FaceColor='flat';
b.CData(y<0,:) = repmat([1 0 0],sum(y<0),1);
title('Controls, Heat Runs'); axis([0 48 ymin ymax]); ylabel('\Delta log(S)'); %yticks(0:0.2:1);
xticks(1:48); xticklabels(labels); xtickangle(45);

nexttile
y = mean(S_chg_ctrl_neut,2);
b = bar(y, 'FaceColor','g');
hold on
er = errorbar(y,SE_dS_ctrl_neut); er.Color = [0 0 0]; er.LineStyle = 'none';
hold off
b.FaceColor='flat';
b.CData(y<0,:) = repmat([1 0 0],sum(y<0),1);
title('Controls, Neutral Runs'); axis([0 48 ymin ymax]); ylabel('\Delta log(S)'); %yticks(0:0.2:1);
xticks(1:48); xticklabels(labels); xtickangle(45);

nexttile
y = mean(S_chg_medt_heat,2);
b = bar(y, 'FaceColor','g');
hold on
er = errorbar(y,SE_dS_medt_heat); er.Color = [0 0 0]; er.LineStyle = 'none';
hold off
b.FaceColor='flat';
b.CData(y<0,:) = repmat([1 0 0],sum(y<0),1);
title('Meditators, Heat Runs'); axis([0 48 ymin ymax]); ylabel('\Delta log(S)'); %yticks(0:0.2:1);
xticks(1:48); xticklabels(labels); xtickangle(45);

nexttile
y = mean(S_chg_medt_neut,2);
b = bar(y, 'FaceColor','g');
hold on
er = errorbar(y,SE_dS_medt_neut); er.Color = [0 0 0]; er.LineStyle = 'none';
hold off
b.FaceColor='flat';
b.CData(y<0,:) = repmat([1 0 0],sum(y<0),1);
title('Meditators, Neutral Runs'); axis([0 48 ymin ymax]); ylabel('\Delta log(S)'); %yticks(0:0.2:1);
xticks(1:48); xticklabels(labels); xtickangle(45);

clear y*
%% Subset of ROIs implicated in pain?

%% Changes in local roi entropy vs changes in VAS pain scores by Condition (next: by run & compute R all together per ROI)
% Correlation analysis

corr_ctrl_int_neut = zeros(48,2);
corr_ctrl_int_heat = zeros(48,2);
corr_ctrl_unp_neut = zeros(48,2);
corr_ctrl_unp_heat = zeros(48,2);
corr_medt_int_neut = zeros(48,2);
corr_medt_int_heat = zeros(48,2);
corr_medt_unp_neut = zeros(48,2);
corr_medt_unp_heat = zeros(48,2);

for roi = 1:48
    % controls, pain int vs S in neut/heat
    [R,P] = corrcoef(vas_chg_int_ctrl,S_chg_ctrl_neut(roi,:));
    corr_ctrl_int_neut(roi,1) = R(1,2);
    corr_ctrl_int_neut(roi,2) = P(1,2);
    [R,P] = corrcoef(vas_chg_int_ctrl,S_chg_ctrl_heat(roi,:));
    corr_ctrl_int_heat(roi,1) = R(1,2);
    corr_ctrl_int_heat(roi,2) = P(1,2);
    % controls, pain unp vs S in neut/heat
    [R,P] = corrcoef(vas_chg_unp_ctrl,S_chg_ctrl_neut(roi,:));
    corr_ctrl_unp_neut(roi,1) = R(1,2);
    corr_ctrl_unp_neut(roi,2) = P(1,2);
    [R,P] = corrcoef(vas_chg_unp_ctrl,S_chg_ctrl_heat(roi,:));
    corr_ctrl_unp_heat(roi,1) = R(1,2);
    corr_ctrl_unp_heat(roi,2) = P(1,2);
    % meditators, pain int vs S in neut/heat
    [R,P] = corrcoef(vas_chg_int_medt,S_chg_medt_neut(roi,:));
    corr_medt_int_neut(roi,1) = R(1,2);
    corr_medt_int_neut(roi,2) = P(1,2);
    [R,P] = corrcoef(vas_chg_int_medt,S_chg_medt_heat(roi,:));
    corr_medt_int_heat(roi,1) = R(1,2);
    corr_medt_int_heat(roi,2) = P(1,2);
    % meditators, pain unp vs S in neut/heat
    [R,P] = corrcoef(vas_chg_unp_medt,S_chg_medt_neut(roi,:));
    corr_medt_unp_neut(roi,1) = R(1,2);
    corr_medt_unp_neut(roi,2) = P(1,2);
    [R,P] = corrcoef(vas_chg_unp_medt,S_chg_medt_heat(roi,:));
    corr_medt_unp_heat(roi,1) = R(1,2);
    corr_medt_unp_heat(roi,2) = P(1,2);
end
clear roi P R

% Plot Correlation Values per Condition
x = 1:48;
l = categorical(labels);
figure
tiledlayout(4,2)
sgtitle('Change in VAS Pain vs Change in ROI logS - Correlation R')
nexttile
bar(x,corr_ctrl_int_heat(:,1),'b')
title('Controls, Intensity, S(heat)'); axis([0 48 -.4 .5]); ylabel('r');
nexttile
bar(x,corr_ctrl_int_neut(:,1),'b')
title('Controls, Intensity, S(neut)'); axis([0 48 -.4 .5]); ylabel('r');
nexttile
bar(x,corr_ctrl_unp_heat(:,1),'b')
title('Controls, Unpleasantness, S(heat)'); axis([0 48 -.4 .5]); ylabel('r');
nexttile
bar(x,corr_ctrl_unp_neut(:,1),'b')
title('Controls, Unpleasantness, S(neut)'); axis([0 48 -.4 .5]); ylabel('r');
nexttile
bar(x,corr_medt_int_heat(:,1),'r')
title('Meditators, Intensity, S(heat)'); axis([0 48 -.4 .5]); ylabel('r');
nexttile
bar(x,corr_medt_int_neut(:,1),'r')
title('Meditators, Intensity, S(neut)'); axis([0 48 -.4 .5]); ylabel('r');
nexttile
bar(x,corr_medt_unp_heat(:,1),'r')
title('Meditators, Unpleasantness, S(heat)'); axis([0 48 -.4 .5]); ylabel('r');
xticks(1:48); xticklabels(labels); xtickangle(45);
nexttile
bar(x,corr_medt_unp_neut(:,1),'r')
title('Meditators, Unpleasantness, S(neut)'); axis([0 48 -.4 .5]); ylabel('r');
xticks(1:48); xticklabels(labels); xtickangle(45);

% Plot Correlation p-values per condition
pthresh = 0.0549;
figure
tiledlayout(4,2)
sgtitle('Change in VAS Pain vs Change in ROI logS - Corr P Values')

nexttile
y=corr_ctrl_int_heat(:,2);
b = bar(x,y, 'FaceColor','#d3d3d3');
b.FaceColor='flat';
b.CData(y<pthresh,:) = repmat([1 0 0],sum(y<pthresh),1);
title('Controls, Intensity, S(heat)'); axis([0 48 0 1]); ylabel('p'); yticks(0:0.2:1);
nexttile
y=corr_ctrl_int_neut(:,2);
b = bar(x,y, 'FaceColor','#d3d3d3');
b.FaceColor='flat';
b.CData(y<pthresh,:) = repmat([1 0 0],sum(y<pthresh),1);
title('Controls, Intensity, S(neut)'); axis([0 48 0 1]); ylabel('p'); yticks(0:0.2:1);
nexttile
y=corr_ctrl_unp_heat(:,2);
b = bar(x,y, 'FaceColor','#d3d3d3');
b.FaceColor='flat';
b.CData(y<pthresh,:) = repmat([1 0 0],sum(y<pthresh),1);
title('Controls, Unpleasantness, S(heat)'); axis([0 48 0 1]); ylabel('p'); yticks(0:0.2:1);
nexttile
y=corr_ctrl_unp_neut(:,2);
b = bar(x,y, 'FaceColor','#d3d3d3');
b.FaceColor='flat';
b.CData(y<pthresh,:) = repmat([1 0 0],sum(y<pthresh),1);
title('Controls, Unpleasantness, S(neut)'); axis([0 48 0 1]); ylabel('p'); yticks(0:0.2:1);
nexttile
y=corr_medt_int_heat(:,2);
b = bar(x,y, 'FaceColor','#d3d3d3');
b.FaceColor='flat';
b.CData(y<pthresh,:) = repmat([1 0 0],sum(y<pthresh),1);
title('Meditators, Intensity, S(heat)'); axis([0 48 0 1]); ylabel('p'); yticks(0:0.2:1);
nexttile
y=corr_medt_int_neut(:,2);
b = bar(x,y, 'FaceColor','#d3d3d3');
b.FaceColor='flat';
b.CData(y<pthresh,:) = repmat([1 0 0],sum(y<pthresh),1);
title('Meditators, Intensity, S(neut)'); axis([0 48 0 1]); ylabel('p'); yticks(0:0.2:1);
nexttile
y=corr_medt_unp_heat(:,2);
b = bar(x,y, 'FaceColor','#d3d3d3');
b.FaceColor='flat';
b.CData(y<pthresh,:) = repmat([1 0 0],sum(y<pthresh),1);
title('Meditators, Unpleasantness, S(heat)'); axis([0 48 0 1]); ylabel('p'); yticks(0:0.2:1);
xticks(1:48); xticklabels(labels); xtickangle(45);
nexttile
y=corr_medt_unp_neut(:,2);
b = bar(x,y, 'FaceColor','#d3d3d3');
b.FaceColor='flat';
b.CData(y<pthresh,:) = repmat([1 0 0],sum(y<pthresh),1);
title('Meditators, Unpleasantness, S(neut)'); axis([0 48 0 1]); ylabel('p'); yticks(0:0.2:1);
xticks(1:48); xticklabels(labels); xtickangle(45);
clear x y l pthresh b

%% Changes in local roi entropy vs changes in VAS pain scores FOR ONE ROI (all sub, all cond)
% Correlation analysis

% One ROI
roi = 30; % select ROI number (31 = precuneous cortex)

% Change in log(S) for ROI for all sub (x)
dS_ctrl_h_roi1 = S_chg_ctrl_heat(roi,:);
dS_ctrl_n_roi1 = S_chg_ctrl_neut(roi,:);
dS_medt_h_roi1 = S_chg_medt_heat(roi,:);
dS_medt_n_roi1 = S_chg_medt_neut(roi,:); 

% Corr Plot for Chosen ROI
figure
tiledlayout(2,2)
nROI = cell2mat(labels(roi)); str = sprintf('Change in log(S) vs Change in VAS Pain for %s ROI',nROI);
sgtitle(str);
% heat, intensity
nexttile
scatter(dS_ctrl_h_roi1,vas_chg_int_ctrl)
hold on
scatter(dS_medt_h_roi1,vas_chg_int_medt)
hold off
title('S During Heat Runs, VAS Intensity'); xlabel('\Delta log(S)'); ylabel('\Delta VAS Intensity'); legend('controls','meditators');
% neut, intensity
nexttile
scatter(dS_ctrl_n_roi1,vas_chg_int_ctrl)
hold on
scatter(dS_medt_n_roi1,vas_chg_int_medt)
hold off
title('S During Neutral Runs, VAS Intensity'); xlabel('\Delta log(S)'); ylabel('\Delta VAS Intensity'); legend('controls','meditators');
% heat, unpleasantness
nexttile
scatter(dS_ctrl_h_roi1,vas_chg_unp_ctrl)
hold on
scatter(dS_medt_h_roi1,vas_chg_unp_medt)
hold off
title('S During Heat Runs, VAS Unpleasantness'); xlabel('\Delta log(S)'); ylabel('\Delta VAS unpleasantness'); legend('controls','meditators');
% neut, unpleasantness
nexttile
scatter(dS_ctrl_n_roi1,vas_chg_unp_ctrl)
hold on
scatter(dS_medt_n_roi1,vas_chg_unp_medt)
hold off
title('S During Neutral Runs, VAS Unpleasantness'); xlabel('\Delta log(S)'); ylabel('\Delta VAS unpleasantness'); legend('controls','meditators');

clear str roi nROI
%% Student t-tests to compare mean local entropy changes in ROIs between
% meditation and control groups (uncorrected & Bonferroni corr.)

