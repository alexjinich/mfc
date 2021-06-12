

clc;
close all;
clearvars;

%% Choose version of matrices:

cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0321
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0404_denoised/
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0405_denoised_GS/
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0407_denoised/
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0416/
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0417/
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0418/

%% test robustness of >S during meditation to different abs-threshold levels 

%loads
load('ctrl_premanip_neut_matrices.mat');
load('ctrl_premanip_heat_matrices.mat');
load('ctrl_posmanip_neut_matrices.mat');
load('ctrl_posmanip_heat_matrices.mat');
load('medt_premanip_neut_matrices.mat');
load('medt_premanip_heat_matrices.mat');
load('medt_posmanip_neut_matrices.mat');
load('medt_posmanip_heat_matrices.mat');

% permute dimensions
ctrl_premanip_neut_corr_matrices = permute(ctrl_premanip_neut_corr_matrices,[3 2 1]);
ctrl_premanip_heat_corr_matrices = permute(ctrl_premanip_heat_corr_matrices,[3 2 1]);
ctrl_posmanip_neut_corr_matrices = permute(ctrl_posmanip_neut_corr_matrices,[3 2 1]);
ctrl_posmanip_heat_corr_matrices = permute(ctrl_posmanip_heat_corr_matrices,[3 2 1]);
medt_premanip_neut_corr_matrices = permute(medt_premanip_neut_corr_matrices,[3 2 1]);
medt_premanip_heat_corr_matrices = permute(medt_premanip_heat_corr_matrices,[3 2 1]);
medt_posmanip_neut_corr_matrices = permute(medt_posmanip_neut_corr_matrices,[3 2 1]);
medt_posmanip_heat_corr_matrices = permute(medt_posmanip_heat_corr_matrices,[3 2 1]);

% get # rois
nrois = size(ctrl_premanip_heat_corr_matrices);
nrois = nrois(1);

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
    sub = cat(3,ctrl_premanip_heat_corr_matrices(:,:,i),ctrl_premanip_heat_corr_matrices(:,:,i+1));
    sub_ctrl_premanip_heat(:,:,i) = mean(sub,3);
    sub_ctrl_premanip_heat(:,:,i) = sub_ctrl_premanip_heat(:,:,i) - diag(diag(sub_ctrl_premanip_heat(:,:,i)));
end
sub_ctrl_premanip_heat = sub_ctrl_premanip_heat(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% posmanip-heat
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,ctrl_posmanip_heat_corr_matrices(:,:,i),ctrl_posmanip_heat_corr_matrices(:,:,i+1));
    sub_ctrl_posmanip_heat(:,:,i) = mean(sub,3);
    sub_ctrl_posmanip_heat(:,:,i) = sub_ctrl_posmanip_heat(:,:,i) - diag(diag(sub_ctrl_posmanip_heat(:,:,i)));
end
sub_ctrl_posmanip_heat = sub_ctrl_posmanip_heat(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% premanip-neut
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,ctrl_premanip_neut_corr_matrices(:,:,i),ctrl_premanip_neut_corr_matrices(:,:,i+1));
    sub_ctrl_premanip_neut(:,:,i) = mean(sub,3);
    sub_ctrl_premanip_neut(:,:,i) = sub_ctrl_premanip_neut(:,:,i) - diag(diag(sub_ctrl_premanip_neut(:,:,i)));
end
sub_ctrl_premanip_neut = sub_ctrl_premanip_neut(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% posmanip-neut
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,ctrl_posmanip_neut_corr_matrices(:,:,i),ctrl_posmanip_neut_corr_matrices(:,:,i+1));
    sub_ctrl_posmanip_neut(:,:,i) = mean(sub,3);
    sub_ctrl_posmanip_neut(:,:,i) = sub_ctrl_posmanip_neut(:,:,i) - diag(diag(sub_ctrl_posmanip_neut(:,:,i)));
end
sub_ctrl_posmanip_neut = sub_ctrl_posmanip_neut(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% Meditators (20 sub, 160 scans)
% premanip-heat
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,medt_premanip_heat_corr_matrices(:,:,i),medt_premanip_heat_corr_matrices(:,:,i+1));
    sub_medt_premanip_heat(:,:,i) = mean(sub,3);
    sub_medt_premanip_heat(:,:,i) = sub_medt_premanip_heat(:,:,i) - diag(diag(sub_medt_premanip_heat(:,:,i)));
end
sub_medt_premanip_heat = sub_medt_premanip_heat(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% posmanip-heat
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,medt_posmanip_heat_corr_matrices(:,:,i),medt_posmanip_heat_corr_matrices(:,:,i+1));
    sub_medt_posmanip_heat(:,:,i) = mean(sub,3);
    sub_medt_posmanip_heat(:,:,i) = sub_medt_posmanip_heat(:,:,i) - diag(diag(sub_medt_posmanip_heat(:,:,i)));
end
sub_medt_posmanip_heat = sub_medt_posmanip_heat(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% premanip-neut
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,medt_premanip_neut_corr_matrices(:,:,i),medt_premanip_neut_corr_matrices(:,:,i+1));
    sub_medt_premanip_neut(:,:,i) = mean(sub,3);
    sub_medt_premanip_neut(:,:,i) = sub_medt_premanip_neut(:,:,i) - diag(diag(sub_medt_premanip_neut(:,:,i)));
end
sub_medt_premanip_neut = sub_medt_premanip_neut(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% posmanip-neut
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,medt_posmanip_neut_corr_matrices(:,:,i),medt_posmanip_neut_corr_matrices(:,:,i+1));
    sub_medt_posmanip_neut(:,:,i) = mean(sub,3);
    sub_medt_posmanip_neut(:,:,i) = sub_medt_posmanip_neut(:,:,i) - diag(diag(sub_medt_posmanip_neut(:,:,i)));
end
sub_medt_posmanip_neut = sub_medt_posmanip_neut(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);

clear sub i 

%% DDist and S per subject for all corr thresh's

thresh_vals = [0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90];
num_thresh_vals = length(thresh_vals);
results = zeros(1,num_thresh_vals);

% preallocate
sub_ddist_ctrl_premanip_heat = zeros(20,nrois,num_thresh_vals);
sub_nS_ctrl_premanip_heat = zeros(num_thresh_vals,20);
sub_ddist_ctrl_posmanip_heat = zeros(20,nrois,num_thresh_vals);
sub_nS_ctrl_posmanip_heat = zeros(num_thresh_vals,20);
sub_ddist_ctrl_premanip_neut = zeros(20,nrois,num_thresh_vals);
sub_nS_ctrl_premanip_neut = zeros(num_thresh_vals,20);
sub_ddist_ctrl_posmanip_neut = zeros(20,nrois,num_thresh_vals);
sub_nS_ctrl_posmanip_neut = zeros(num_thresh_vals,20);
sub_ddist_medt_premanip_heat = zeros(20,nrois,num_thresh_vals);
sub_nS_medt_premanip_heat = zeros(num_thresh_vals,20);
sub_ddist_medt_posmanip_heat = zeros(20,nrois,num_thresh_vals);
sub_nS_medt_posmanip_heat = zeros(num_thresh_vals,20);
sub_ddist_medt_premanip_neut = zeros(20,nrois,num_thresh_vals);
sub_nS_medt_premanip_neut = zeros(num_thresh_vals,20);
sub_ddist_medt_posmanip_neut = zeros(20,nrois,num_thresh_vals);
sub_nS_medt_posmanip_neut = zeros(num_thresh_vals,20);

for thresh = 1:num_thresh_vals
    
    % Get subject-level thresholded matrices, ddist, and nS
    % controls, premanip, heat
    sub_ddist = zeros(20,nrois);
    sub_nS = zeros(1,20);
    for i = 1:20
        unthresh_matrix = sub_ctrl_premanip_heat(:,:,i);
        [~, ddist, nS] = preproc_corr_matrix(unthresh_matrix, thresh_vals(thresh));
        sub_ddist(i,:) = ddist;
        sub_nS(i) = nS;
    end
    sub_ddist_ctrl_premanip_heat(:,:,thresh) = sub_ddist(:,:);
    sub_nS_ctrl_premanip_heat(thresh,:) = sub_nS(:,:);
    clear sub_ddist ddist nS sub_nS unthresh_matrix; 

    % controls, posmanip, heat
    sub_ddist = zeros(20,nrois);
    sub_nS = zeros(1,20);
    for i = 1:20
        unthresh_matrix = sub_ctrl_posmanip_heat(:,:,i);
        [~, ddist, nS] = preproc_corr_matrix(unthresh_matrix, thresh_vals(thresh));
        sub_ddist(i,:) = ddist;
        sub_nS(i) = nS;
    end
    sub_ddist_ctrl_posmanip_heat(:,:,thresh) = sub_ddist(:,:);
    sub_nS_ctrl_posmanip_heat(thresh,:) = sub_nS(:,:);
    clear sub_ddist ddist nS sub_nS unthresh_matrix; 

    % controls, premanip, neut
    sub_ddist = zeros(20,nrois);
    sub_nS = zeros(1,20);
    for i = 1:20
        unthresh_matrix = sub_ctrl_premanip_neut(:,:,i);
        [~, ddist, nS] = preproc_corr_matrix(unthresh_matrix, thresh_vals(thresh));
        sub_ddist(i,:) = ddist;
        sub_nS(i) = nS;
    end
    sub_ddist_ctrl_premanip_neut(:,:,thresh) = sub_ddist(:,:);
    sub_nS_ctrl_premanip_neut(thresh,:) = sub_nS(:,:);
    clear sub_ddist ddist nS sub_nS unthresh_matrix; 

    % controls, posmanip, neut
    sub_ddist = zeros(20,nrois);
    sub_nS = zeros(1,20);
    for i = 1:20
        unthresh_matrix = sub_ctrl_posmanip_neut(:,:,i);
        [~, ddist, nS] = preproc_corr_matrix(unthresh_matrix, thresh_vals(thresh));
        sub_ddist(i,:) = ddist;
        sub_nS(i) = nS;
    end
    sub_ddist_ctrl_posmanip_neut(:,:,thresh) = sub_ddist(:,:);
    sub_nS_ctrl_posmanip_neut(thresh,:) = sub_nS(:,:);
    clear sub_ddist ddist nS sub_nS unthresh_matrix; 
 
    % meditators, premanip, heat
    sub_ddist = zeros(20,nrois);
    sub_nS = zeros(1,20);
    for i = 1:20
        unthresh_matrix = sub_medt_premanip_heat(:,:,i);
        [~, ddist, nS] = preproc_corr_matrix(unthresh_matrix, thresh_vals(thresh));
        sub_ddist(i,:) = ddist;
        sub_nS(i) = nS;
    end
    sub_ddist_medt_premanip_heat(:,:,thresh) = sub_ddist(:,:);
    sub_nS_medt_premanip_heat(thresh,:) = sub_nS(:,:);
    clear sub_ddist ddist nS sub_nS unthresh_matrix; 
 
    % meditators, posmanip, heat
    sub_ddist = zeros(20,nrois);
    sub_nS = zeros(1,20);
    for i = 1:20
        unthresh_matrix = sub_medt_posmanip_heat(:,:,i);
        [~, ddist, nS] = preproc_corr_matrix(unthresh_matrix, thresh_vals(thresh));
        sub_ddist(i,:) = ddist;
        sub_nS(i) = nS;
    end
    sub_ddist_medt_posmanip_heat(:,:,thresh) = sub_ddist(:,:);
    sub_nS_medt_posmanip_heat(thresh,:) = sub_nS(:,:);
    clear sub_ddist ddist nS sub_nS unthresh_matrix; 
 
    % meditators, premanip, neut
    sub_ddist = zeros(20,nrois);
    sub_nS = zeros(1,20);
    for i = 1:20
        unthresh_matrix = sub_medt_premanip_neut(:,:,i);
        [~, ddist, nS] = preproc_corr_matrix(unthresh_matrix, thresh_vals(thresh));
        sub_ddist(i,:) = ddist;
        sub_nS(i) = nS;
    end
    sub_ddist_medt_premanip_neut(:,:,thresh) = sub_ddist(:,:);
    sub_nS_medt_premanip_neut(thresh,:) = sub_nS(:,:);
    clear sub_ddist ddist nS sub_nS unthresh_matrix; 
 
    % meditators, posmanip, neut
    sub_ddist = zeros(20,nrois);
    sub_nS = zeros(1,20);
    for i = 1:20
        unthresh_matrix = sub_medt_posmanip_neut(:,:,i);
        [~, ddist, nS] = preproc_corr_matrix(unthresh_matrix, thresh_vals(thresh));
        sub_ddist(i,:) = ddist;
        sub_nS(i) = nS;
    end
    sub_ddist_medt_posmanip_neut(:,:,thresh) = sub_ddist(:,:);
    sub_nS_medt_posmanip_neut(thresh,:) = sub_nS(:,:);
    clear sub_ddist ddist nS sub_nS unthresh_matrix;
    
end

%% S and chgS (with SEs)

% mean S by group (all thresh)
S_medt_pre_neut = mean(sub_nS_medt_premanip_neut,2);
S_medt_pos_neut = mean(sub_nS_medt_posmanip_neut,2);
S_medt_pre_heat = mean(sub_nS_medt_premanip_heat,2);
S_medt_pos_heat = mean(sub_nS_medt_posmanip_heat,2);
S_ctrl_pre_neut = mean(sub_nS_ctrl_premanip_neut,2);
S_ctrl_pos_neut = mean(sub_nS_ctrl_posmanip_neut,2);
S_ctrl_pre_heat = mean(sub_nS_ctrl_premanip_heat,2);
S_ctrl_pos_heat = mean(sub_nS_ctrl_posmanip_heat,2);

% SE of S by group (all thresh)
S_SE_medt_pre_neut = std(sub_nS_medt_premanip_neut,0,2) / sqrt(20);
S_SE_medt_pos_neut = std(sub_nS_medt_posmanip_neut,0,2) / sqrt(20);
S_SE_medt_pre_heat = std(sub_nS_medt_premanip_heat,0,2) / sqrt(20);
S_SE_medt_pos_heat = std(sub_nS_medt_posmanip_heat,0,2) / sqrt(20);
S_SE_ctrl_pre_neut = std(sub_nS_ctrl_premanip_neut,0,2) / sqrt(20);
S_SE_ctrl_pos_neut = std(sub_nS_ctrl_posmanip_neut,0,2) / sqrt(20);
S_SE_ctrl_pre_heat = std(sub_nS_ctrl_premanip_heat,0,2) / sqrt(20);
S_SE_ctrl_pos_heat = std(sub_nS_ctrl_posmanip_heat,0,2) / sqrt(20);

h_ctrl_neut = zeros(1,17);
p_ctrl_neut = zeros(1,17);
h_ctrl_heat = zeros(1,17);
p_ctrl_heat = zeros(1,17);
h_medt_neut = zeros(1,17);
p_medt_neut = zeros(1,17);
h_medt_heat = zeros(1,17);
p_medt_heat = zeros(1,17);

for thresh = 1:num_thresh_vals
    % paired t-tests
    [h_ctrl_neut(thresh),p_ctrl_neut(thresh)] = ttest(sub_nS_ctrl_premanip_neut(thresh,:), sub_nS_ctrl_posmanip_neut(thresh,:));
    [h_ctrl_heat(thresh),p_ctrl_heat(thresh)] = ttest(sub_nS_ctrl_premanip_heat(thresh,:), sub_nS_ctrl_posmanip_heat(thresh,:));
    [h_medt_neut(thresh),p_medt_neut(thresh)] = ttest(sub_nS_medt_premanip_neut(thresh,:), sub_nS_medt_posmanip_neut(thresh,:));
    [h_medt_heat(thresh),p_medt_heat(thresh)] = ttest(sub_nS_medt_premanip_heat(thresh,:), sub_nS_medt_posmanip_heat(thresh,:));
end

% post-pre (chg) at sub level
chgS_ctrl_neut = sub_nS_ctrl_posmanip_neut - sub_nS_ctrl_premanip_neut;
chgS_ctrl_heat = sub_nS_ctrl_posmanip_heat - sub_nS_ctrl_premanip_heat;
chgS_medt_neut = sub_nS_medt_posmanip_neut - sub_nS_medt_premanip_neut;
chgS_medt_heat = sub_nS_medt_posmanip_heat - sub_nS_medt_premanip_heat;

% mean chgS per thresh
mean_chgS_ctrl_neut = mean(chgS_ctrl_neut,2);
mean_chgS_ctrl_heat = mean(chgS_ctrl_heat,2);
mean_chgS_medt_neut = mean(chgS_medt_neut,2);
mean_chgS_medt_heat = mean(chgS_medt_heat,2);

% SE of chgS at thresh level
chgS_SE_ctrl_neut = std(chgS_ctrl_neut,0,2) / sqrt(20);
chgS_SE_ctrl_heat = std(chgS_ctrl_heat,0,2) / sqrt(20);
chgS_SE_medt_neut = std(chgS_medt_neut,0,2) / sqrt(20);
chgS_SE_medt_heat = std(chgS_medt_heat,0,2) / sqrt(20);

%% VAS Pain

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

%% FitLM chgVAS - chgS, all thresh

dS_heat_dInt_ctrl = {};
dS_heat_dInt_medt = {};
dS_heat_dUnp_ctrl = {};
dS_heat_dUnp_medt = {};
dS_neut_dInt_ctrl = {};
dS_neut_dInt_medt = {};
dS_neut_dUnp_ctrl = {};
dS_neut_dUnp_medt = {};

for thresh = 1:num_thresh_vals
    
    % Heat
    mdl = fitlm(vas_chg_int_ctrl, chgS_ctrl_heat(thresh,:));
    dS_heat_dInt_ctrl{thresh} = mdl;
    mdl = fitlm(vas_chg_int_medt, chgS_medt_heat(thresh,:));
    dS_heat_dInt_medt{thresh} = mdl;
    mdl = fitlm(vas_chg_unp_ctrl, chgS_ctrl_heat(thresh,:));
    dS_heat_dUnp_ctrl{thresh} = mdl;
    mdl = fitlm(vas_chg_unp_medt, chgS_medt_heat(thresh,:));
    dS_heat_dUnp_medt{thresh} = mdl;
    
    % Neut
    mdl = fitlm(vas_chg_int_ctrl, chgS_ctrl_neut(thresh,:));
    dS_neut_dInt_ctrl{thresh} = mdl;
    mdl = fitlm(vas_chg_int_medt, chgS_medt_neut(thresh,:));
    dS_neut_dInt_medt{thresh} = mdl;
    mdl = fitlm(vas_chg_unp_ctrl, chgS_ctrl_neut(thresh,:));
    dS_neut_dUnp_ctrl{thresh} = mdl;
    mdl = fitlm(vas_chg_unp_medt, chgS_medt_neut(thresh,:));
    dS_neut_dUnp_medt{thresh} = mdl;

    % p-vals 
    dS_heat_dInt_ctrl_p(thresh) = dS_heat_dInt_ctrl{thresh}.Coefficients.pValue(2);
    dS_heat_dInt_medt_p(thresh) = dS_heat_dInt_medt{thresh}.Coefficients.pValue(2);
    dS_heat_dUnp_ctrl_p(thresh) = dS_heat_dUnp_ctrl{thresh}.Coefficients.pValue(2);
    dS_heat_dUnp_medt_p(thresh) = dS_heat_dUnp_medt{thresh}.Coefficients.pValue(2);
    dS_neut_dInt_ctrl_p(thresh) = dS_neut_dInt_ctrl{thresh}.Coefficients.pValue(2);
    dS_neut_dInt_medt_p(thresh) = dS_neut_dInt_medt{thresh}.Coefficients.pValue(2);
    dS_neut_dUnp_ctrl_p(thresh) = dS_neut_dUnp_ctrl{thresh}.Coefficients.pValue(2);
    dS_neut_dUnp_medt_p(thresh) = dS_neut_dUnp_medt{thresh}.Coefficients.pValue(2);

    % coefficients
    dS_heat_dInt_ctrl_c(thresh) = dS_heat_dInt_ctrl{thresh}.Coefficients.Estimate(2);
    dS_heat_dInt_medt_c(thresh) = dS_heat_dInt_medt{thresh}.Coefficients.Estimate(2);
    dS_heat_dUnp_ctrl_c(thresh) = dS_heat_dUnp_ctrl{thresh}.Coefficients.Estimate(2);
    dS_heat_dUnp_medt_c(thresh) = dS_heat_dUnp_medt{thresh}.Coefficients.Estimate(2);
    dS_neut_dInt_ctrl_c(thresh) = dS_neut_dInt_ctrl{thresh}.Coefficients.Estimate(2);
    dS_neut_dInt_medt_c(thresh) = dS_neut_dInt_medt{thresh}.Coefficients.Estimate(2);
    dS_neut_dUnp_ctrl_c(thresh) = dS_neut_dUnp_ctrl{thresh}.Coefficients.Estimate(2);
    dS_neut_dUnp_medt_c(thresh) = dS_neut_dUnp_medt{thresh}.Coefficients.Estimate(2);

end

% cat ctrl and medt p-vals and coefficients
fitlm_p_ctrl = cat(2,dS_heat_dInt_ctrl_p,dS_heat_dUnp_ctrl_p,dS_neut_dInt_ctrl_p,dS_neut_dUnp_ctrl_p);
fitlm_p_medt = cat(2,dS_heat_dInt_medt_p,dS_heat_dUnp_medt_p,dS_neut_dInt_medt_p,dS_neut_dUnp_medt_p);
fitlm_c_ctrl = cat(2,dS_heat_dInt_ctrl_c,dS_heat_dUnp_ctrl_c,dS_neut_dInt_ctrl_c,dS_neut_dUnp_ctrl_c);
fitlm_c_medt = cat(2,dS_heat_dInt_medt_c,dS_heat_dUnp_medt_c,dS_neut_dInt_medt_c,dS_neut_dUnp_medt_c);

% cat ctrl and medt p-vals and coefficients - only neutral
fitlm_p_ctrl_neut = cat(2,dS_neut_dInt_ctrl_p,dS_neut_dUnp_ctrl_p);
fitlm_p_medt_neut = cat(2,dS_neut_dInt_medt_p,dS_neut_dUnp_medt_p);
fitlm_c_ctrl_neut = cat(2,dS_neut_dInt_ctrl_c,dS_neut_dUnp_ctrl_c);
fitlm_c_medt_neut = cat(2,dS_neut_dInt_medt_c,dS_neut_dUnp_medt_c);

% cat ctrl and medt p-vals and coefficients - only heat
fitlm_p_ctrl_heat = cat(2,dS_heat_dInt_ctrl_p,dS_heat_dUnp_ctrl_p);
fitlm_p_medt_heat = cat(2,dS_heat_dInt_medt_p,dS_heat_dUnp_medt_p);
fitlm_c_ctrl_heat = cat(2,dS_heat_dInt_ctrl_c,dS_heat_dUnp_ctrl_c);
fitlm_c_medt_heat = cat(2,dS_heat_dInt_medt_c,dS_heat_dUnp_medt_c);

%% Q-Q p-value plot and p-value histograms (all conditions)

figure
tiledlayout(2,6)

nexttile([2,2])
qqplot(fitlm_p_ctrl,fitlm_p_medt); title('Q-Q Plot: p-values','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
xlabel('p-values (controls)'); ylabel('p-values (meditators)'); ylim([-.2 0.5]); xlim([0 0.5]); xticks(0:.1:.5); yticks(-.2:.1:.5);

nexttile([1,2])
histogram(fitlm_p_ctrl,40,'FaceColor','#EDB120'); 
title('p-Values','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
xl = xline(0.05,'--b','p = 0.05','LineWidth',1,'LabelHorizontalAlignment','left');
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
xlabel('p-value'); ylabel('frequency');
ylim([0 14])
legend('Controls')

nexttile([1,2])
histogram(fitlm_c_ctrl,20,'FaceColor','#EDB120'); 
title('Coefficients','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
xl = xline(0,'--b','zero','LineWidth',1,'LabelHorizontalAlignment','left');
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
xlabel('Coefficient (slope)'); ylabel('frequency');
ylim([0 10]); xlim([-.15 .11]);
legend('Controls')

nexttile([1,2])
histogram(fitlm_p_medt,40,'FaceColor','#A2142F');
%title('Meditators','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
xl = xline(0.05,'--b','p = 0.05','LineWidth',1,'LabelHorizontalAlignment','left');
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
xlabel('p-value'); ylabel('frequency');
ylim([0 14])
legend('Meditators')

nexttile([1,2])
histogram(fitlm_c_medt,20,'FaceColor','#A2142F');
%title('Meditators','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
xl = xline(0,'--b','zero','LineWidth',1,'LabelHorizontalAlignment','left');
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
xlabel('Coefficient (slope)'); ylabel('frequency');
ylim([0 10]); xlim([-.15 .11]);
legend('Meditators')

set(gcf,'Position',[100,100,1400,450]); %xl = xlabel('d'); xl.FontSize = 20;

print -depsc2 qqplot_pvalhist.eps  % export

% %% Overlayed p-value histograms
% 
% figure
% tiledlayout(1,2)
% 
% nexttile
% qqplot(fitlm_p_ctrl,fitlm_p_medt); title('Q-Q Plot: Linear Model p-values');
% xlabel('controls'); ylabel('meditators');
% 
% nexttile
% histogram(fitlm_p_medt,20,'FaceColor','#EDB120'); 
% hold on
% histogram(fitlm_p_ctrl,20,'FaceColor','#7E2F8E'); title('p-vals');
% xl = xline(0.05,'--b','p = 0.05','LineWidth',2,'LabelHorizontalAlignment','left');
% xl.LabelVerticalAlignment = 'middle';
% xl.LabelHorizontalAlignment = 'center';
% xl.Annotation.LegendInformation.IconDisplayStyle = 'off'; % skip line in legend
% legend('meditation','control');
% hold off
% xlabel('p-value'); ylabel('frequency');

%% Q-Q p-value plot and p-value histograms (neutral-only)

figure
tiledlayout(2,6)

nexttile([2,2])
qqplot(fitlm_p_ctrl_neut,fitlm_p_medt_neut); title('Q-Q Plot: p-values (Neutral Series Only)','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
xlabel('p-values (controls)'); ylabel('p-values (meditators)'); ylim([-.2 1]); xlim([0 1]); xticks(0:.1:1); yticks(-.2:.1:1);

nexttile([1,2])
histogram(fitlm_p_ctrl_neut,20,'FaceColor','#EDB120'); 
title('p-Values','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
xl = xline(0.05,'--b','p = 0.05','LineWidth',1,'LabelHorizontalAlignment','left');
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
xlabel('p-value'); ylabel('frequency');
ylim([0 14]); xlim([-0.05 1.05]);
legend('Controls')

nexttile([1,2])
histogram(fitlm_c_ctrl_neut,20,'FaceColor','#EDB120'); 
title('Coefficients','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
xl = xline(0,'--b','zero','LineWidth',1,'LabelHorizontalAlignment','left');
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
xlabel('Coefficient (slope)'); ylabel('frequency');
ylim([0 10]); xlim([-.15 .11]);
legend('Controls')

nexttile([1,2])
histogram(fitlm_p_medt_neut,20,'FaceColor','#A2142F');
%title('Meditators','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
xl = xline(0.05,'--b','p = 0.05','LineWidth',1,'LabelHorizontalAlignment','left');
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
xlabel('p-value'); ylabel('frequency');
ylim([0 14]); xlim([-0.05 1.05]);
legend('Meditators')

nexttile([1,2])
histogram(fitlm_c_medt_neut,20,'FaceColor','#A2142F');
%title('Meditators','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
xl = xline(0,'--b','zero','LineWidth',1,'LabelHorizontalAlignment','left');
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
xlabel('Coefficient (slope)'); ylabel('frequency');
ylim([0 10]); xlim([-.15 .11]);
legend('Meditators')

set(gcf,'Position',[100,100,1400,450]); %xl = xlabel('d'); xl.FontSize = 20;

print -depsc2 qqplot_pvalhist_neut.eps  % export

%% Q-Q p-value plot and p-value histograms (heat-only)

figure
tiledlayout(2,6)

nexttile([2,2])
qqplot(fitlm_p_ctrl_heat,fitlm_p_medt_heat); title('Q-Q Plot: p-values (Heat Series Only)','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
xlabel('p-values (controls)'); ylabel('p-values (meditators)'); ylim([-.2 1]); xlim([0 1]); xticks(0:.1:1); yticks(-.2:.1:1);

nexttile([1,2])
histogram(fitlm_p_ctrl_heat,20,'FaceColor','#EDB120'); 
title('p-Values','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
xl = xline(0.05,'--b','p = 0.05','LineWidth',1,'LabelHorizontalAlignment','left');
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
xlabel('p-value'); ylabel('frequency');
ylim([0 14]); xlim([-0.05 .95]);
legend('Controls')

nexttile([1,2])
histogram(fitlm_c_ctrl_heat,20,'FaceColor','#EDB120'); 
title('Coefficients','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
xl = xline(0,'--b','zero','LineWidth',1,'LabelHorizontalAlignment','left');
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
xlabel('Coefficient (slope)'); ylabel('frequency');
ylim([0 10]); xlim([-.15 .11]);
legend('Controls')

nexttile([1,2])
histogram(fitlm_p_medt_heat,20,'FaceColor','#A2142F');
%title('Meditators','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
xl = xline(0.05,'--b','p = 0.05','LineWidth',1,'LabelHorizontalAlignment','left');
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
xlabel('p-value'); ylabel('frequency');
ylim([0 14]); xlim([-.05 .95]);
legend('Meditators')

nexttile([1,2])
histogram(fitlm_c_medt_heat,20,'FaceColor','#A2142F');
%title('Meditators','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
xl = xline(0,'--b','zero','LineWidth',1,'LabelHorizontalAlignment','left');
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
xlabel('Coefficient (slope)'); ylabel('frequency');
ylim([0 10]); xlim([-.15 .11]);
legend('Meditators')

set(gcf,'Position',[100,100,1400,450]); %xl = xlabel('d'); xl.FontSize = 20;

print -depsc2 qqplot_pvalhist_heat.eps  % export