clc;
close all;
clearvars;

%%%%%%%%%%%%%%%% section of ENTROPY PER sub V2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This version takes matrices per individual sub, averages over 2 subs and then
% computes ddist, S per sub to make sub & group comparisons.
% Threshold_robustness_2 and adjmatrix.m average 2 sub matrices 1st, then
% compute ddist, S on the averaged sub-condition matrix.
% ---> this sub-version generates fitlm plots for visualization (change in
% S vs change in VAS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHOOSE NILEARN MATRIX VERSION
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0418/
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0417/
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0416/
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0407_denoised/
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0405_denoised_GS/
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0404_denoised/
cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0321
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
% SET THRESH
thresh = 0.35;
%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET 1 IF LOG(S), 0 NO LOG
log_S = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath '/Volumes/Seagate_Desktop_Drive/mfc/code'

% load
load('ctrl_premanip_neut_matrices.mat');
load('ctrl_premanip_heat_matrices.mat');
load('ctrl_posmanip_neut_matrices.mat');
load('ctrl_posmanip_heat_matrices.mat');
load('medt_premanip_neut_matrices.mat');
load('medt_premanip_heat_matrices.mat');
load('medt_posmanip_neut_matrices.mat');
load('medt_posmanip_heat_matrices.mat');

% permute
ctrl_pre_neut_corr_matrices = permute(ctrl_premanip_neut_corr_matrices,[3 2 1]);
ctrl_pre_heat_corr_matrices = permute(ctrl_premanip_heat_corr_matrices,[3 2 1]);
ctrl_pos_neut_corr_matrices = permute(ctrl_posmanip_neut_corr_matrices,[3 2 1]);
ctrl_pos_heat_corr_matrices = permute(ctrl_posmanip_heat_corr_matrices,[3 2 1]);
medt_pre_neut_corr_matrices = permute(medt_premanip_neut_corr_matrices,[3 2 1]);
medt_pre_heat_corr_matrices = permute(medt_premanip_heat_corr_matrices,[3 2 1]);
medt_pos_neut_corr_matrices = permute(medt_posmanip_neut_corr_matrices,[3 2 1]);
medt_pos_heat_corr_matrices = permute(medt_posmanip_heat_corr_matrices,[3 2 1]);

% get # rois
nrois = size(ctrl_pre_neut_corr_matrices);
nrois = nrois(1);

% clear imports
clear *posmanip*
clear *premanip*

%% Flatten/Average 2 subs/condition into 1 matrix

% preallocate sub 
sub_ctrl_premanip_heat = zeros(nrois,nrois,20);
sub_ctrl_posmanip_heat = zeros(nrois,nrois,20);
sub_ctrl_premanip_neut = zeros(nrois,nrois,20);
sub_ctrl_posmanip_neut = zeros(nrois,nrois,20);
sub_medt_premanip_heat = zeros(nrois,nrois,20);
sub_medt_posmanip_heat = zeros(nrois,nrois,20);
sub_medt_premanip_neut = zeros(nrois,nrois,20);
sub_medt_posmanip_neut = zeros(nrois,nrois,20);

% Create condition matrices w/1 matrix per subject (i.e. flatten 2 into 1)
% Controls (20 sub, 160 scans)
% premanip-heat
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,ctrl_pre_heat_corr_matrices(:,:,i),ctrl_pre_heat_corr_matrices(:,:,i+1));
    sub_ctrl_premanip_heat(:,:,i) = mean(sub,3);
    sub_ctrl_premanip_heat(:,:,i) = sub_ctrl_premanip_heat(:,:,i) - diag(diag(sub_ctrl_premanip_heat(:,:,i)));
end
sub_ctrl_premanip_heat = sub_ctrl_premanip_heat(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% posmanip-heat
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,ctrl_pos_heat_corr_matrices(:,:,i),ctrl_pos_heat_corr_matrices(:,:,i+1));
    sub_ctrl_posmanip_heat(:,:,i) = mean(sub,3);
    sub_ctrl_posmanip_heat(:,:,i) = sub_ctrl_posmanip_heat(:,:,i) - diag(diag(sub_ctrl_posmanip_heat(:,:,i)));
end
sub_ctrl_posmanip_heat = sub_ctrl_posmanip_heat(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% premanip-neut
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,ctrl_pre_neut_corr_matrices(:,:,i),ctrl_pre_neut_corr_matrices(:,:,i+1));
    sub_ctrl_premanip_neut(:,:,i) = mean(sub,3);
    sub_ctrl_premanip_neut(:,:,i) = sub_ctrl_premanip_neut(:,:,i) - diag(diag(sub_ctrl_premanip_neut(:,:,i)));
end
sub_ctrl_premanip_neut = sub_ctrl_premanip_neut(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% posmanip-neut
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,ctrl_pos_neut_corr_matrices(:,:,i),ctrl_pos_neut_corr_matrices(:,:,i+1));
    sub_ctrl_posmanip_neut(:,:,i) = mean(sub,3);
    sub_ctrl_posmanip_neut(:,:,i) = sub_ctrl_posmanip_neut(:,:,i) - diag(diag(sub_ctrl_posmanip_neut(:,:,i)));
end
sub_ctrl_posmanip_neut = sub_ctrl_posmanip_neut(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% Meditators (20 sub, 160 scans)
% premanip-heat
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,medt_pre_heat_corr_matrices(:,:,i),medt_pre_heat_corr_matrices(:,:,i+1));
    sub_medt_premanip_heat(:,:,i) = mean(sub,3);
    sub_medt_premanip_heat(:,:,i) = sub_medt_premanip_heat(:,:,i) - diag(diag(sub_medt_premanip_heat(:,:,i)));
end
sub_medt_premanip_heat = sub_medt_premanip_heat(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% posmanip-heat
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,medt_pos_heat_corr_matrices(:,:,i),medt_pos_heat_corr_matrices(:,:,i+1));
    sub_medt_posmanip_heat(:,:,i) = mean(sub,3);
    sub_medt_posmanip_heat(:,:,i) = sub_medt_posmanip_heat(:,:,i) - diag(diag(sub_medt_posmanip_heat(:,:,i)));
end
sub_medt_posmanip_heat = sub_medt_posmanip_heat(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% premanip-neut
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,medt_pre_neut_corr_matrices(:,:,i),medt_pre_neut_corr_matrices(:,:,i+1));
    sub_medt_premanip_neut(:,:,i) = mean(sub,3);
    sub_medt_premanip_neut(:,:,i) = sub_medt_premanip_neut(:,:,i) - diag(diag(sub_medt_premanip_neut(:,:,i)));
end
sub_medt_premanip_neut = sub_medt_premanip_neut(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% posmanip-neut
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,medt_pos_neut_corr_matrices(:,:,i),medt_pos_neut_corr_matrices(:,:,i+1));
    sub_medt_posmanip_neut(:,:,i) = mean(sub,3);
    sub_medt_posmanip_neut(:,:,i) = sub_medt_posmanip_neut(:,:,i) - diag(diag(sub_medt_posmanip_neut(:,:,i)));
end
sub_medt_posmanip_neut = sub_medt_posmanip_neut(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);

clear sub i 
clear ctrl* medt*

%% Diag-zeros, threshold, binarize, compute entropy, ddist per sub:

% Controls, pre, heat
sub_bin_abs_ctrl_pre_heat = zeros(nrois,nrois,20);
sub_ddist_ctrl_pre_heat = zeros(1,nrois);
sub_nS_ctrl_pre_heat = zeros(1,20);
for i = 1:20
    unthresh_matrix = sub_ctrl_premanip_heat(:,:,i); % sub_ctrl_pre_heat_corr_matrices
    [w_bin_abs, ddist, nS] = proc_corr_matrix(unthresh_matrix, thresh);
    sub_bin_abs_ctrl_pre_heat(:,:,i) = w_bin_abs;
    sub_ddist_ctrl_pre_heat(i,:) = ddist;
    sub_nS_ctrl_pre_heat(1,i) = nS;
end
% Controls, pos, heat
sub_bin_abs_ctrl_pos_heat = zeros(nrois,nrois,20);
sub_ddist_ctrl_pos_heat = zeros(1,nrois);
sub_nS_ctrl_pos_heat = zeros(1,20);
for i = 1:20
    unthresh_matrix = sub_ctrl_posmanip_heat(:,:,i);
    [w_bin_abs, ddist, nS] = proc_corr_matrix(unthresh_matrix, thresh);
    sub_bin_abs_ctrl_pos_heat(:,:,i) = w_bin_abs;
    sub_ddist_ctrl_pos_heat(i,:) = ddist;
    sub_nS_ctrl_pos_heat(1,i) = nS;
end
% Controls, pre, neut
sub_bin_abs_ctrl_pre_neut = zeros(nrois,nrois,20);
sub_ddist_ctrl_pre_neut = zeros(1,nrois);
sub_nS_ctrl_pre_neut = zeros(1,20);
for i = 1:20
    unthresh_matrix = sub_ctrl_premanip_neut(:,:,i);
    [w_bin_abs, ddist, nS] = proc_corr_matrix(unthresh_matrix, thresh);
    sub_bin_abs_ctrl_pre_neut(:,:,i) = w_bin_abs;
    sub_ddist_ctrl_pre_neut(i,:) = ddist;
    sub_nS_ctrl_pre_neut(1,i) = nS;
end
% Controls, pos, neut
sub_bin_abs_ctrl_pos_neut = zeros(nrois,nrois,20);
sub_ddist_ctrl_pos_neut = zeros(1,nrois);
sub_nS_ctrl_pos_neut = zeros(1,20);
for i = 1:20
    unthresh_matrix = sub_ctrl_posmanip_neut(:,:,i);
    [w_bin_abs, ddist, nS] = proc_corr_matrix(unthresh_matrix, thresh);
    sub_bin_abs_ctrl_pos_neut(:,:,i) = w_bin_abs;
    sub_ddist_ctrl_pos_neut(i,:) = ddist;
    sub_nS_ctrl_pos_neut(1,i) = nS;
end
% Meditators, pre, heat
sub_bin_abs_medt_pre_heat = zeros(nrois,nrois,20);
sub_ddist_medt_pre_heat = zeros(1,nrois);
sub_nS_medt_pre_heat = zeros(1,20);
for i = 1:20
    unthresh_matrix = sub_medt_premanip_heat(:,:,i);
    [w_bin_abs, ddist, nS] = proc_corr_matrix(unthresh_matrix, thresh);
    sub_bin_abs_medt_pre_heat(:,:,i) = w_bin_abs;
    sub_ddist_medt_pre_heat(i,:) = ddist;
    sub_nS_medt_pre_heat(1,i) = nS;
end
% Meditators, pos, heat
sub_bin_abs_medt_pos_heat = zeros(nrois,nrois,20);
sub_ddist_medt_pos_heat = zeros(1,nrois);
sub_nS_medt_pos_heat = zeros(1,20);
for i = 1:20
    unthresh_matrix = sub_medt_posmanip_heat(:,:,i);
    [w_bin_abs, ddist, nS] = proc_corr_matrix(unthresh_matrix, thresh);
    sub_bin_abs_medt_pos_heat(:,:,i) = w_bin_abs;
    sub_ddist_medt_pos_heat(i,:) = ddist;
    sub_nS_medt_pos_heat(1,i) = nS;
end
% Meditators, pre, neut
sub_bin_abs_medt_pre_neut = zeros(nrois,nrois,20);
sub_ddist_medt_pre_neut = zeros(1,nrois);
sub_nS_medt_pre_neut = zeros(1,20);
for i = 1:20
    unthresh_matrix = sub_medt_premanip_neut(:,:,i);
    [w_bin_abs, ddist, nS] = proc_corr_matrix(unthresh_matrix, thresh);
    sub_bin_abs_medt_pre_neut(:,:,i) = w_bin_abs;
    sub_ddist_medt_pre_neut(i,:) = ddist;
    sub_nS_medt_pre_neut(1,i) = nS;
end
% Meditators, pos, neut
sub_bin_abs_medt_pos_neut = zeros(nrois,nrois,20);
sub_ddist_medt_pos_neut = zeros(1,nrois);
sub_nS_medt_pos_neut = zeros(1,20);
for i = 1:20
    unthresh_matrix = sub_medt_posmanip_neut(:,:,i);
    [w_bin_abs, ddist, nS] = proc_corr_matrix(unthresh_matrix, thresh);
    sub_bin_abs_medt_pos_neut(:,:,i) = w_bin_abs;
    sub_ddist_medt_pos_neut(i,:) = ddist;
    sub_nS_medt_pos_neut(1,i) = nS;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% change in entropy post-pre
nS_chg_ctrl_heat = sub_nS_ctrl_pos_heat - sub_nS_ctrl_pre_heat;
nS_chg_ctrl_neut = sub_nS_ctrl_pos_neut - sub_nS_ctrl_pre_neut;
nS_chg_medt_heat = sub_nS_medt_pos_heat - sub_nS_medt_pre_heat;
nS_chg_medt_neut = sub_nS_medt_pos_neut - sub_nS_medt_pre_neut;

% SE of change in S
SE_dS_ctrl_heat = std(nS_chg_ctrl_heat) / sqrt(20);
SE_dS_ctrl_neut = std(nS_chg_ctrl_neut) / sqrt(20);
SE_dS_medt_heat = std(nS_chg_medt_heat) / sqrt(20);
SE_dS_medt_neut = std(nS_chg_medt_neut) / sqrt(20);

%% RANOVA S or log(s)

% NEUTRAL
% Construct datatable with the 4 clusters of data and name the variables.
datatable = table(sub_nS_ctrl_pre_neut',sub_nS_ctrl_pos_neut',sub_nS_medt_pre_neut',sub_nS_medt_pos_neut');
datatable.Properties.VariableNames = {'pre_controls','pre_meditators','post_controls','post_meditators'};
% Since we have more than one repeated-measures factor, we set up a table
% to indicate the levels on each factor for each of the different variables.
WithinStructure = table([1 1 2 2]',[1 2 1 2]','VariableNames',{'Treatment','Group'});
% The 4 rows of the WithinStructure table correspond to 4 different cols
% 'pre_A','pre_B','post_A','post_B' in data table. Each 'pre_A','pre_B','post_A','post_B' column is coded as 1/2 on the Treatment factor
% and as 1/2 on the Group factor. Now pass the WithinStructure table to fitrm so that it knows how the different
% columns correspond to the different levels of the repeated-measures factors.
rm = fitrm(datatable, 'pre_controls,pre_meditators,post_controls,post_meditators~1','WithinDesign',WithinStructure);
% Finally, specify the repeated-measures factors again when you call ranova & display results:
disp('Repeated Measures ANOVA: log(S), Neutral, 2 by 2:');
ranovatable = ranova(rm,'WithinModel','Treatment*Group')
% disp('Coefficients:'); 
% rm.Coefficients
% disp('Covariance:'); 
% rm.Covariance
% Clear temp vars.
clear datatable WithinStructure rm ranovatable

% HEAT
% Construct datatable with the 4 clusters of data and name the variables.
datatable = table(sub_nS_ctrl_pre_heat',sub_nS_ctrl_pos_heat',sub_nS_medt_pre_heat',sub_nS_medt_pos_heat');
datatable.Properties.VariableNames = {'pre_controls','pre_meditators','post_controls','post_meditators'};
% Since we have more than one repeated-measures factor, we set up a table
% to indicate the levels on each factor for each of the different variables.
WithinStructure = table([1 1 2 2]',[1 2 1 2]','VariableNames',{'Treatment','Group'});
% The 4 rows of the WithinStructure table correspond to 4 different cols
% 'pre_A','pre_B','post_A','post_B' in data table. Each 'pre_A','pre_B','post_A','post_B' column is coded as 1/2 on the Treatment factor
% and as 1/2 on the Group factor. Now pass the WithinStructure table to fitrm so that it knows how the different
% columns correspond to the different levels of the repeated-measures factors.
rm = fitrm(datatable, 'pre_controls,pre_meditators,post_controls,post_meditators~1','WithinDesign',WithinStructure);
% Finally, specify the repeated-measures factors again when you call ranova & display results:
disp('Repeated Measures ANOVA: log(S), Heat, 2 by 2:');
ranovatable = ranova(rm,'WithinModel','Treatment*Group')
% disp('Coefficients:'); 
% rm.Coefficients
% disp('Covariance:'); 
% rm.Covariance
% Clear temp vars.
clear datatable WithinStructure rm ranovatable


% T-tests per subject for sig diff in S pre vs post

disp('T-tests for Difference in S Pre vs Post');

[h,p] = ttest(sub_nS_ctrl_pre_neut, sub_nS_ctrl_pos_neut); % Controls, Neutral
if h == 0
    h = '*not* rejected';
else
    h = '***rejected***';
end
d = sprintf('\nControls,   Neutral - Null H (no mean difference) is %s with p = %.3f', h, p); disp(d); 

[h,p] = ttest(sub_nS_ctrl_pre_heat, sub_nS_ctrl_pos_heat); % Controls, Heat
if h == 0
    h = '*not* rejected';
else
    h = '***rejected***';
end
d = sprintf('Controls,   Heat    - Null H (no mean difference) is %s with p = %.3f', h, p); disp(d); 

[h,p] = ttest(sub_nS_medt_pre_neut, sub_nS_medt_pos_neut); % Meditators, Neutral
if h == 0
    h = '*not* rejected';
else
    h = '***rejected***';
end
d = sprintf('Meditators, Neutral - Null H (no mean difference) is %s with p = %.3f', h, p); disp(d); 
 
[h,p] = ttest(sub_nS_medt_pre_heat, sub_nS_medt_pos_heat); % Meditators, Heat
if h == 0
    h = '*not* rejected';
else
    h = '***rejected***';
end
d = sprintf('Meditators, Heat    - Null H (no mean difference) is %s with p = %.3f', h, p); disp(d); 

clear h p d

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VAS PAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get VASint+unp
[id	v6_int_h1 v6_unpl_h1 v6_int_h2 v6_unpl_h2 v6_int_h3 v6_unpl_h3 v6_int_h4 v6_unpl_h4] = readvars('pain_ratings.xlsx');

% split into ctrl/medt 1/sub
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
clear v6_int_h1 v6_unpl_h1 v6_int_h2 v6_unpl_h2 v6_int_h3 v6_unpl_h3 v6_int_h4 v6_unpl_h4

% vas by condition (e.g. h1 h2 = h12, h3 h4 = h34)
vas_int_h12_ctrl = (vas_int_h1_ctrl + vas_int_h2_ctrl) ./ 2;
vas_int_h34_ctrl = (vas_int_h3_ctrl + vas_int_h4_ctrl) ./ 2;
vas_int_h12_medt = (vas_int_h1_medt + vas_int_h2_medt) ./ 2;
vas_int_h34_medt = (vas_int_h3_medt + vas_int_h4_medt) ./ 2;
vas_unp_h12_ctrl = (vas_unp_h1_ctrl + vas_unp_h2_ctrl) ./ 2;
vas_unp_h34_ctrl = (vas_unp_h3_ctrl + vas_unp_h4_ctrl) ./ 2;
vas_unp_h12_medt = (vas_unp_h1_medt + vas_unp_h2_medt) ./ 2;
vas_unp_h34_medt = (vas_unp_h3_medt + vas_unp_h4_medt) ./ 2;

% change in pain post-pre (mean of h3,h4 - mean of h1,h2)
vas_chg_int_ctrl = ((vas_int_h3_ctrl + vas_int_h4_ctrl)./2) - ((vas_int_h1_ctrl + vas_int_h2_ctrl)./2);
vas_chg_int_medt = ((vas_int_h3_medt + vas_int_h4_medt)./2) - ((vas_int_h1_medt + vas_int_h2_medt)./2);
vas_chg_unp_ctrl = ((vas_unp_h3_ctrl + vas_unp_h4_ctrl)./2) - ((vas_unp_h1_ctrl + vas_unp_h2_ctrl)./2);
vas_chg_unp_medt = ((vas_unp_h3_medt + vas_unp_h4_medt)./2) - ((vas_unp_h1_medt + vas_unp_h2_medt)./2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Heat subs - regress chgS vs chg VAS
%fprintf('\nChange in Entropy during Heat vs Change in VAS Pain during Heat\n');
%fprintf('\nVAS INT, Controls\n'); 
dS_heat_dInt_ctrl = fitlm(vas_chg_int_ctrl,nS_chg_ctrl_heat);
%fprintf('\nVAS INT, Meditators\n'); 
dS_heat_dInt_medt = fitlm(vas_chg_int_medt,nS_chg_medt_heat);
%fprintf('\nVAS UNP, Controls\n'); 
dS_heat_dUnp_ctrl = fitlm(vas_chg_unp_ctrl,nS_chg_ctrl_heat);
%fprintf('\nVAS UNP, Meditators\n'); 
dS_heat_dUnp_medt = fitlm(vas_chg_unp_medt,nS_chg_medt_heat);


%% plot regressions

x_tick = [-4:1:3]; y_tick = [-.5:.2:.6];

figure; tiledlayout(2,2); 
str = sprintf('Change in Entropy (Heat Runs) vs Change in VAS Pain at %.2f Correlation Threshold',thresh);
sgtitle(str)

nexttile
plot(dS_heat_dInt_ctrl); title({'Controls, VAS Pain Int';''},'FontUnits','points','FontWeight','normal','FontSize',13,'FontName','Calibri'); 
xlabel('Change in pain intensity'); ylabel('Change in S'); axis([-3 3 -.5 .6]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_heat_dInt_ctrl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_heat_dInt_ctrl.Rsquared.Ordinary)];
annotation('textbox',[.20 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

nexttile
plot(dS_heat_dInt_medt); title({'Meditators, VAS Pain Int';''},'FontUnits','points','FontWeight','normal','FontSize',13,'FontName','Calibri'); 
xlabel('Change in pain intensity'); ylabel('Change in S'); axis([-4 3 -.5 .6]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_heat_dInt_medt.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_heat_dInt_medt.Rsquared.Ordinary)];
annotation('textbox',[.66 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

nexttile
plot(dS_heat_dUnp_ctrl); title({'Controls, VAS Pain Unp';''},'FontUnits','points','FontWeight','normal','FontSize',13,'FontName','Calibri');   
xlabel('Change in pain unpleasantness'); ylabel('Change in S'); axis([-3 3 -.5 .6]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_heat_dUnp_ctrl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_heat_dUnp_ctrl.Rsquared.Ordinary)];
annotation('textbox',[.20 0.44 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

nexttile
plot(dS_heat_dUnp_medt); title({'Meditators, VAS Pain Unp';''},'FontUnits','points','FontWeight','normal','FontSize',13,'FontName','Calibri'); 
xlabel('Change in pain unpleasantness'); ylabel('Change in S'); axis([-4 3 -.5 .6]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_heat_dUnp_medt.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_heat_dUnp_medt.Rsquared.Ordinary)];
annotation('textbox',[.66 0.44 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

set(gcf,'Position',[100,100,550,400]);
print -depsc2 lm_dS_dVAS_heat.eps  % export

% Neut subs - regress chg S vs chg VAS
%fprintf('\nChange in Entropy during Neut vs Change in VAS Pain during Heat\n');
%fprintf('\nVAS INT, Controls\n'); 
dS_neut_dInt_ctrl = fitlm(vas_chg_int_ctrl,nS_chg_ctrl_neut);
%fprintf('\nVAS INT, Meditators\n'); 
dS_neut_dInt_medt = fitlm(vas_chg_int_medt,nS_chg_medt_neut);
%fprintf('\nVAS UNP, Controls\n'); 
dS_neut_dUnp_ctrl = fitlm(vas_chg_unp_ctrl,nS_chg_ctrl_neut);
%fprintf('\nVAS UNP, Meditators\n'); 
dS_neut_dUnp_medt = fitlm(vas_chg_unp_medt,nS_chg_medt_neut);

% plot regressions
figure; tiledlayout(2,2); 
str = sprintf('Change in Entropy (Neutral Runs) vs Change in VAS Pain at %.2f Correlation Threshold',thresh);
sgtitle(str)

nexttile
plot(dS_neut_dInt_ctrl); title({'Controls, VAS Pain Int';''},'FontUnits','points','FontWeight','normal','FontSize',13,'FontName','Calibri'); 
xlabel('Change in pain intensity'); ylabel('Change in S'); axis([-3 3 -.5 .6]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_neut_dInt_ctrl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_neut_dInt_ctrl.Rsquared.Ordinary)];
annotation('textbox',[.20 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

nexttile
plot(dS_neut_dInt_medt); title({'Meditators, VAS Pain Int';''},'FontUnits','points','FontWeight','normal','FontSize',13,'FontName','Calibri'); 
xlabel('Change in pain intensity'); ylabel('Change in S'); axis([-4 3 -.5 .6]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_neut_dInt_medt.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_neut_dInt_medt.Rsquared.Ordinary)];
annotation('textbox',[.66 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

nexttile
plot(dS_neut_dUnp_ctrl); title({'Controls, VAS Pain Unp';''},'FontUnits','points','FontWeight','normal','FontSize',13,'FontName','Calibri');   
xlabel('Change in pain unpleasantness'); ylabel('Change in S'); axis([-3 3 -.5 .6]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_neut_dUnp_ctrl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_neut_dUnp_ctrl.Rsquared.Ordinary)];
annotation('textbox',[.20 0.44 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

nexttile
plot(dS_neut_dUnp_medt); title({'Meditators, VAS Pain Unp';''},'FontUnits','points','FontWeight','normal','FontSize',13,'FontName','Calibri'); 
xlabel('Change in pain unpleasantness'); ylabel('Change in S'); axis([-4 3 -.5 .6]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_neut_dUnp_medt.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_neut_dUnp_medt.Rsquared.Ordinary)];
annotation('textbox',[.66 0.44 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

set(gcf,'Position',[100,100,550,400]);
print -depsc2 lm_dS_dVAS_neut.eps  % export

%% functions
% takes unthresholded corr matrix, makes diag zeroes, and returns 
% thresholded binarized matrix, its degree distribution, and its entropy
function [w_bin_abs, ddist, nS] = proc_corr_matrix(unthresh_matrix, thresh)
   unthresh_matrix = unthresh_matrix - diag(diag(unthresh_matrix));
   abs_thresh_matrix = threshold_absolute(unthresh_matrix, thresh);
   w_bin_abs = weight_conversion(abs_thresh_matrix, 'binarize');
   ddist = degrees_und(w_bin_abs);
   
   entropy_bins = 120;
   f = figure('visible','off'); %initialize fig invisible
   hist1 = histogram(ddist,entropy_bins);
   hist1.Normalization = 'probability';
   binvalues = hist1.Values;
   nS = sum(-(binvalues(binvalues>0).*(log(binvalues(binvalues>0)))));
   clear f; clear hist1;
end