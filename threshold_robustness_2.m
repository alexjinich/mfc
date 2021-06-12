% Computes entropy (and ddist) at the individual subject level. Then
% aggregates them into group level and plots this for all correlation
% thresholds. Unlike version 1, which flattens matrices into group level
% first and computes entropy later (more sensitive to ddist ave). V2
% preferred. Pending: compute per run, express as log(S) w/ log(0)=0.

clc;
close all;
clearvars;

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

%% S and deltaS by group

nS_medt_posmanip_heat = zeros(1,num_thresh_vals);
nS_medt_premanip_heat = zeros(1,num_thresh_vals);
nS_medt_posmanip_neut = zeros(1,num_thresh_vals);
nS_medt_premanip_neut = zeros(1,num_thresh_vals);
nS_ctrl_posmanip_heat = zeros(1,num_thresh_vals);
nS_ctrl_premanip_heat = zeros(1,num_thresh_vals);
nS_ctrl_posmanip_neut = zeros(1,num_thresh_vals);
nS_ctrl_premanip_neut = zeros(1,num_thresh_vals);

chg_nS_medt_heat = zeros(1,num_thresh_vals);
chg_nS_medt_neut = zeros(1,num_thresh_vals);
chg_nS_ctrl_heat = zeros(1,num_thresh_vals);
chg_nS_ctrl_neut = zeros(1,num_thresh_vals);

for thresh = 1:num_thresh_vals
    % mean S by group
    nS_medt_posmanip_heat(thresh) = mean(sub_nS_medt_posmanip_heat(thresh,:));
    nS_medt_premanip_heat(thresh) = mean(sub_nS_medt_premanip_heat(thresh,:));
    nS_medt_posmanip_neut(thresh) = mean(sub_nS_medt_posmanip_neut(thresh,:));
    nS_medt_premanip_neut(thresh) = mean(sub_nS_medt_premanip_neut(thresh,:));
    nS_ctrl_posmanip_heat(thresh) = mean(sub_nS_ctrl_posmanip_heat(thresh,:));
    nS_ctrl_premanip_heat(thresh) = mean(sub_nS_ctrl_premanip_heat(thresh,:));
    nS_ctrl_posmanip_neut(thresh) = mean(sub_nS_ctrl_posmanip_neut(thresh,:));
    nS_ctrl_premanip_neut(thresh) = mean(sub_nS_ctrl_premanip_neut(thresh,:));

    % delta S (post - pre) 
    chg_nS_medt_heat(thresh) = nS_medt_posmanip_heat(thresh) - nS_medt_premanip_heat(thresh);
    chg_nS_medt_neut(thresh) = nS_medt_posmanip_neut(thresh) - nS_medt_premanip_neut(thresh);
    chg_nS_ctrl_heat(thresh) = nS_ctrl_posmanip_heat(thresh) - nS_ctrl_premanip_heat(thresh);
    chg_nS_ctrl_neut(thresh) = nS_ctrl_posmanip_neut(thresh) - nS_ctrl_premanip_neut(thresh);
end
% mean change per group (over all thresholds)
mean_chg_nS_medt_heat = mean(chg_nS_medt_heat);
mean_chg_nS_medt_neut = mean(chg_nS_medt_neut);
mean_chg_nS_ctrl_heat = mean(chg_nS_ctrl_heat);
mean_chg_nS_ctrl_neut = mean(chg_nS_ctrl_neut);

%% S and chgS (v2, with SEs)

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

% %% Plot: Change in S by r Threshold (Group) - v2
% 
% x = categorical(["Controls Neutral" "Controls Heat" "Meditators Neutral" "Meditators Heat"]);
% x = reordercats(x,{'Controls Neutral' 'Meditators Neutral' 'Controls Heat' 'Meditators Heat'});
% 
% y = [mean_chgS_ctrl_neut(1) mean_chgS_medt_neut(1) mean_chgS_ctrl_heat(1) mean_chgS_medt_heat(1); ...
%      mean_chgS_ctrl_neut(2) mean_chgS_medt_neut(2) mean_chgS_ctrl_heat(2) mean_chgS_medt_heat(2); ...
%      mean_chgS_ctrl_neut(3) mean_chgS_medt_neut(3) mean_chgS_ctrl_heat(3) mean_chgS_medt_heat(3); ...
%      mean_chgS_ctrl_neut(4) mean_chgS_medt_neut(4) mean_chgS_ctrl_heat(4) mean_chgS_medt_heat(4); ...
%      mean_chgS_ctrl_neut(5) mean_chgS_medt_neut(5) mean_chgS_ctrl_heat(5) mean_chgS_medt_heat(5); ...
%      mean_chgS_ctrl_neut(6) mean_chgS_medt_neut(6) mean_chgS_ctrl_heat(6) mean_chgS_medt_heat(6); ...
%      mean_chgS_ctrl_neut(7) mean_chgS_medt_neut(7) mean_chgS_ctrl_heat(7) mean_chgS_medt_heat(7); ...
%      mean_chgS_ctrl_neut(8) mean_chgS_medt_neut(8) mean_chgS_ctrl_heat(8) mean_chgS_medt_heat(8); ...
%      mean_chgS_ctrl_neut(9) mean_chgS_medt_neut(9) mean_chgS_ctrl_heat(9) mean_chgS_medt_heat(9); ...
%      mean_chgS_ctrl_neut(10) mean_chgS_medt_neut(10) mean_chgS_ctrl_heat(10) mean_chgS_medt_heat(10); ...
%      mean_chgS_ctrl_neut(11) mean_chgS_medt_neut(11) mean_chgS_ctrl_heat(11) mean_chgS_medt_heat(11); ...
%      mean_chgS_ctrl_neut(12) mean_chgS_medt_neut(12) mean_chgS_ctrl_heat(12) mean_chgS_medt_heat(12); ...
%      mean_chgS_ctrl_neut(13) mean_chgS_medt_neut(13) mean_chgS_ctrl_heat(13) mean_chgS_medt_heat(13); ...
%      mean_chgS_ctrl_neut(14) mean_chgS_medt_neut(14) mean_chgS_ctrl_heat(14) mean_chgS_medt_heat(14); ...
%      mean_chgS_ctrl_neut(15) mean_chgS_medt_neut(15) mean_chgS_ctrl_heat(15) mean_chgS_medt_heat(15); ...
%      mean_chgS_ctrl_neut(16) mean_chgS_medt_neut(16) mean_chgS_ctrl_heat(16) mean_chgS_medt_heat(16); ...
%      mean_chgS_ctrl_neut(17) mean_chgS_medt_neut(17) mean_chgS_ctrl_heat(17) mean_chgS_medt_heat(17)];
% 
%  err = [chgS_SE_ctrl_neut(1) chgS_SE_ctrl_heat(1) chgS_SE_medt_neut(1) chgS_SE_medt_heat(1); ...
%         chgS_SE_ctrl_neut(2) chgS_SE_ctrl_heat(2) chgS_SE_medt_neut(2) chgS_SE_medt_heat(2); ...
%         chgS_SE_ctrl_neut(3) chgS_SE_ctrl_heat(3) chgS_SE_medt_neut(3) chgS_SE_medt_heat(3); ...
%         chgS_SE_ctrl_neut(4) chgS_SE_ctrl_heat(4) chgS_SE_medt_neut(4) chgS_SE_medt_heat(4); ...
%         chgS_SE_ctrl_neut(5) chgS_SE_ctrl_heat(5) chgS_SE_medt_neut(5) chgS_SE_medt_heat(5); ...
%         chgS_SE_ctrl_neut(6) chgS_SE_ctrl_heat(6) chgS_SE_medt_neut(6) chgS_SE_medt_heat(6); ...
%         chgS_SE_ctrl_neut(7) chgS_SE_ctrl_heat(7) chgS_SE_medt_neut(7) chgS_SE_medt_heat(7); ...
%         chgS_SE_ctrl_neut(8) chgS_SE_ctrl_heat(8) chgS_SE_medt_neut(8) chgS_SE_medt_heat(8); ...
%         chgS_SE_ctrl_neut(9) chgS_SE_ctrl_heat(9) chgS_SE_medt_neut(9) chgS_SE_medt_heat(9); ... 
%         chgS_SE_ctrl_neut(10) chgS_SE_ctrl_heat(10) chgS_SE_medt_neut(10) chgS_SE_medt_heat(10); ...
%         chgS_SE_ctrl_neut(11) chgS_SE_ctrl_heat(11) chgS_SE_medt_neut(11) chgS_SE_medt_heat(11); ...
%         chgS_SE_ctrl_neut(12) chgS_SE_ctrl_heat(12) chgS_SE_medt_neut(12) chgS_SE_medt_heat(12); ...
%         chgS_SE_ctrl_neut(13) chgS_SE_ctrl_heat(13) chgS_SE_medt_neut(13) chgS_SE_medt_heat(13); ...
%         chgS_SE_ctrl_neut(14) chgS_SE_ctrl_heat(14) chgS_SE_medt_neut(14) chgS_SE_medt_heat(14); ...
%         chgS_SE_ctrl_neut(15) chgS_SE_ctrl_heat(15) chgS_SE_medt_neut(15) chgS_SE_medt_heat(15); ...
%         chgS_SE_ctrl_neut(16) chgS_SE_ctrl_heat(16) chgS_SE_medt_neut(16) chgS_SE_medt_heat(16); ...
%         chgS_SE_ctrl_neut(17) chgS_SE_ctrl_heat(17) chgS_SE_medt_neut(17) chgS_SE_medt_heat(17)];
% 
% f=figure; 
% b=bar(x,y); 
% hold on
% %er = errorbar(x,y,err); er.Color = 'k'; er.LineStyle = 'none';
% hold off
% c = parula(17);
% newcolors = c;
% colororder(newcolors)
% yline(mean_chg_nS_medt_heat,':','medt heat','LabelHorizontalAlignment','center','LabelVerticalAlignment','top');
% yline(mean_chg_nS_medt_neut,':','medt neut','LabelHorizontalAlignment','center','LabelVerticalAlignment','top');
% yline(mean_chg_nS_ctrl_heat,':','ctrl heat','LabelHorizontalAlignment','center','LabelVerticalAlignment','top');
% yline(mean_chg_nS_ctrl_neut,':','ctrl neut','LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
% lgd = legend('0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50',...
%     '0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','Location','eastoutside','Orientation','horizontal'); 
% lgd.NumColumns = 2; lgd.Title.String = 'Correlation Threshold Value';
% sgtitle('Pre-to-Post Change in Shannon Entropy, Group-level'); ylabel('\DeltaS','Rotation',0');
% set(gcf,'Position',[100,100,1300,700]);


%% Plot: Change in S by r Threshold (Group) - v1

x = categorical(["Controls Neutral" "Controls Heat" "Meditators Neutral" "Meditators Heat"]);
x = reordercats(x,{'Controls Neutral' 'Meditators Neutral' 'Controls Heat' 'Meditators Heat'});
y = [chg_nS_ctrl_neut(1) chg_nS_ctrl_heat(1) chg_nS_medt_neut(1) chg_nS_medt_heat(1); ...
    chg_nS_ctrl_neut(2) chg_nS_ctrl_heat(2) chg_nS_medt_neut(2) chg_nS_medt_heat(2); ...
    chg_nS_ctrl_neut(3) chg_nS_ctrl_heat(3) chg_nS_medt_neut(3) chg_nS_medt_heat(3); ...
    chg_nS_ctrl_neut(4) chg_nS_ctrl_heat(4) chg_nS_medt_neut(4) chg_nS_medt_heat(4); ...
    chg_nS_ctrl_neut(5) chg_nS_ctrl_heat(5) chg_nS_medt_neut(5) chg_nS_medt_heat(5); ...
    chg_nS_ctrl_neut(6) chg_nS_ctrl_heat(6) chg_nS_medt_neut(6) chg_nS_medt_heat(6); ...
    chg_nS_ctrl_neut(7) chg_nS_ctrl_heat(7) chg_nS_medt_neut(7) chg_nS_medt_heat(7); ...
    chg_nS_ctrl_neut(8) chg_nS_ctrl_heat(8) chg_nS_medt_neut(8) chg_nS_medt_heat(8); ...
    chg_nS_ctrl_neut(9) chg_nS_ctrl_heat(9) chg_nS_medt_neut(9) chg_nS_medt_heat(9); ... 
    chg_nS_ctrl_neut(10) chg_nS_ctrl_heat(10) chg_nS_medt_neut(10) chg_nS_medt_heat(10); ...
    chg_nS_ctrl_neut(11) chg_nS_ctrl_heat(11) chg_nS_medt_neut(11) chg_nS_medt_heat(11); ...
    chg_nS_ctrl_neut(12) chg_nS_ctrl_heat(12) chg_nS_medt_neut(12) chg_nS_medt_heat(12); ...
    chg_nS_ctrl_neut(13) chg_nS_ctrl_heat(13) chg_nS_medt_neut(13) chg_nS_medt_heat(13); ...
    chg_nS_ctrl_neut(14) chg_nS_ctrl_heat(14) chg_nS_medt_neut(14) chg_nS_medt_heat(14); ...
    chg_nS_ctrl_neut(15) chg_nS_ctrl_heat(15) chg_nS_medt_neut(15) chg_nS_medt_heat(15); ...
    chg_nS_ctrl_neut(16) chg_nS_ctrl_heat(16) chg_nS_medt_neut(16) chg_nS_medt_heat(16); ...
    chg_nS_ctrl_neut(17) chg_nS_ctrl_heat(17) chg_nS_medt_neut(17) chg_nS_medt_heat(17)];
f=figure; b=bar(x,y); 
c = parula(17);
newcolors = c;
colororder(newcolors)
yline(mean_chg_nS_medt_heat,':','mean (meditation, heat)','LabelHorizontalAlignment','center','LabelVerticalAlignment','top');
yline(mean_chg_nS_medt_neut,':','mean (meditation, neutral)','LabelHorizontalAlignment','center','LabelVerticalAlignment','top');
yline(mean_chg_nS_ctrl_heat,':','mean (control, heat)','LabelHorizontalAlignment','center','LabelVerticalAlignment','top');
yline(mean_chg_nS_ctrl_neut,':','mean (control, neutral)','LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
lgd = legend('0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50',...
    '0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','Location','eastoutside','Orientation','horizontal'); 
lgd.NumColumns = 2; lgd.Title.String = 'Correlation Threshold';
ylim([-0.07 0.38]);
sgtitle('Change in Shannon Entropy (Post-Pre), Group-level'); ylabel('\DeltaS','Rotation',0','FontUnits','points','FontWeight','normal','FontSize',14);
set(gcf,'Position',[100,100,1100,400]); %xl = xlabel('d'); xl.FontSize = 20;

print -depsc2 chgS_allthresh.eps  % export


%% Plot: S (Pre & Post) for all thresholds - all groups, no t-test plots

x = 0:16;
y = 3 * ones(1,17);
thresh_labels = {'0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90'};

y1 = 1 * ones(1,17);

figure
str = sprintf('Shannon Entropy Pre- & Post-Manipulation, All Thresholds');
sgtitle(str);
tiledlayout(8,8)

nexttile([4,4])
s1 = errorbar(x,S_ctrl_pre_neut',S_SE_ctrl_pre_neut');
s1.Marker = 's'; s1.MarkerSize = 5; s1.MarkerEdgeColor = 'blue'; s1.MarkerFaceColor = 'blue';
hold on
s2 = errorbar(x,S_ctrl_pos_neut',S_SE_ctrl_pos_neut');
s2.Marker = 's'; s2.MarkerSize = 5; s2.MarkerEdgeColor = 'red'; s2.MarkerFaceColor = 'red';
stem(x,y,':','MarkerEdgeColor','w')  % draw vertical ref lines
title('Controls, Neutral')
ylabel('S     ','Rotation',0'); legend('Pre','Post'); 
xticks(0:1:16); xticklabels(thresh_labels); xlabel('Threshold');

nexttile([4,4])
s1 = errorbar(x,S_medt_pre_neut',S_SE_medt_pre_neut');
s1.Marker = 's'; s1.MarkerSize = 5; s1.MarkerEdgeColor = 'blue'; s1.MarkerFaceColor = 'blue';
hold on
s2 = errorbar(x,S_medt_pos_neut',S_SE_medt_pos_neut');
s2.Marker = 's'; s2.MarkerSize = 5; s2.MarkerEdgeColor = 'red'; s2.MarkerFaceColor = 'red';
stem(x,y,':','MarkerEdgeColor','w')  % draw vertical ref lines
title('Meditators, Neutral')
ylabel('S     ','Rotation',0'); legend('Pre','Post'); 
xticks(0:1:16); xticklabels(thresh_labels); xlabel('Threshold');

nexttile([4,4])
s1 = errorbar(x,S_ctrl_pre_heat',S_SE_ctrl_pre_heat');
s1.Marker = 's'; s1.MarkerSize = 5; s1.MarkerEdgeColor = 'blue'; s1.MarkerFaceColor = 'blue';
hold on
s2 = errorbar(x,S_ctrl_pos_heat',S_SE_ctrl_pos_heat');
s2.Marker = 's'; s2.MarkerSize = 5; s2.MarkerEdgeColor = 'red'; s2.MarkerFaceColor = 'red';
stem(x,y,':','MarkerEdgeColor','w')  % draw vertical ref lines
title('Controls, Heat')
ylabel('S     ','Rotation',0'); legend('Pre','Post'); 
xticks(0:1:16); xticklabels(thresh_labels); xlabel('Threshold');

nexttile([4,4])
s1 = errorbar(x,S_medt_pre_heat',S_SE_medt_pre_heat');
s1.Marker = 's'; s1.MarkerSize = 5; s1.MarkerEdgeColor = 'blue'; s1.MarkerFaceColor = 'blue';
hold on
s2 = errorbar(x,S_medt_pos_heat',S_SE_medt_pos_heat');
s2.Marker = 's'; s2.MarkerSize = 5; s2.MarkerEdgeColor = 'red'; s2.MarkerFaceColor = 'red';
stem(x,y,':','MarkerEdgeColor','w')  % draw vertical ref lines
title('Meditators, Heat')
ylabel('S     ','Rotation',0'); legend('Pre','Post'); 
xticks(0:1:16); xticklabels(thresh_labels); xlabel('Threshold');

%% Plot: S (Pre & Post) for all thresholds - Indiv plots with t-tests

x = 0:16;
y = 3 * ones(1,17);
thresh_labels = {'0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90'};
y1 = 1 * ones(1,17);

% Controls Neutral
figure
tiledlayout(5,4)

nexttile([4,4]) 
s1 = errorbar(x,S_ctrl_pre_neut',S_SE_ctrl_pre_neut');
s1.Marker = 's'; s1.MarkerSize = 5; s1.MarkerEdgeColor = 'blue'; s1.MarkerFaceColor = 'blue';
hold on
s2 = errorbar(x,S_ctrl_pos_neut',S_SE_ctrl_pos_neut');
s2.Marker = 's'; s2.MarkerSize = 5; s2.MarkerEdgeColor = 'red'; s2.MarkerFaceColor = 'red';
stem(x,y,':','MarkerEdgeColor','w')  % draw vertical ref lines
title('Controls, Neutral','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
ylabel('S     ','Rotation',0'); legend('Pre','Post'); 
xticks(0:1:16); set(gca,'xtick',[])

nexttile([1,4]) % t-test labels
stem(x,y1,':','MarkerEdgeColor','w');
lblc = sprintfc('%.3f',p_ctrl_neut); set(gca,'ytick',[])
t = text(x, y1-.5, lblc, 'HorizontalAlignment','center', 'VerticalAlignment','middle','Rotation',90);
lblc_p = find(p_ctrl_neut < 0.05); % find p<0.05 indices
for i = 1:length(lblc_p)
    t(lblc_p(i)).Color = 'blue'; t(lblc_p(i)).FontWeight = 'bold';
end
xticks(0:1:16); xticklabels(thresh_labels); xlabel('Threshold'); ylabel('t-test   ','Rotation',0,'HorizontalAlignment','right');

print -depsc2 S_allthresh_ctrl_neut.eps  % export

% Meditators Neutral
figure
tiledlayout(5,4)

nexttile([4,4]) 
s1 = errorbar(x,S_medt_pre_neut',S_SE_medt_pre_neut');
s1.Marker = 's'; s1.MarkerSize = 5; s1.MarkerEdgeColor = 'blue'; s1.MarkerFaceColor = 'blue';
hold on
s2 = errorbar(x,S_medt_pos_neut',S_SE_medt_pos_neut');
s2.Marker = 's'; s2.MarkerSize = 5; s2.MarkerEdgeColor = 'red'; s2.MarkerFaceColor = 'red';
stem(x,y,':','MarkerEdgeColor','w')  % draw vertical ref lines
title('Meditators, Neutral','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
ylabel('S     ','Rotation',0'); legend('Pre','Post'); 
xticks(0:1:16); set(gca,'xtick',[])

nexttile([1,4]) % t-test labels
stem(x,y1,':','MarkerEdgeColor','w');
lblc = sprintfc('%.3f',p_medt_neut); set(gca,'ytick',[])
t = text(x, y1-.5, lblc, 'HorizontalAlignment','center', 'VerticalAlignment','middle','Rotation',90);
lblc_p = find(p_medt_neut < 0.05); % find p<0.05 indices
for i = 1:length(lblc_p)
    t(lblc_p(i)).Color = 'blue'; t(lblc_p(i)).FontWeight = 'bold';
end
xticks(0:1:16); xticklabels(thresh_labels); xlabel('Threshold'); ylabel('t-test   ','Rotation',0,'HorizontalAlignment','right');

print -depsc2 S_allthresh_medt_neut.eps  % export


% Controls Heat
figure
tiledlayout(5,4)

nexttile([4,4])
s1 = errorbar(x,S_ctrl_pre_heat',S_SE_ctrl_pre_heat');
s1.Marker = 's'; s1.MarkerSize = 5; s1.MarkerEdgeColor = 'blue'; s1.MarkerFaceColor = 'blue';
hold on
s2 = errorbar(x,S_ctrl_pos_heat',S_SE_ctrl_pos_heat');
s2.Marker = 's'; s2.MarkerSize = 5; s2.MarkerEdgeColor = 'red'; s2.MarkerFaceColor = 'red';
stem(x,y,':','MarkerEdgeColor','w')  % draw vertical ref lines
title('Controls, Heat','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
ylabel('S     ','Rotation',0'); legend('Pre','Post'); 
xticks(0:1:16); set(gca,'xtick',[])

nexttile([1,4]) % t-test labels
stem(x,y1,':','MarkerEdgeColor','w');
lblc = sprintfc('%.3f',p_ctrl_heat); set(gca,'ytick',[])
t = text(x, y1-.5, lblc, 'HorizontalAlignment','center', 'VerticalAlignment','middle','Rotation',90);
lblc_p = find(p_ctrl_heat < 0.05); % find p<0.05 indices
for i = 1:length(lblc_p)
    t(lblc_p(i)).Color = 'blue'; t(lblc_p(i)).FontWeight = 'bold';
end
xticks(0:1:16); xticklabels(thresh_labels); xlabel('Threshold'); ylabel('t-test   ','Rotation',0,'HorizontalAlignment','right');

print -depsc2 S_allthresh_ctrl_heat.eps  % export


% Meditators Heat
figure
tiledlayout(5,4)

nexttile([4,4])
s1 = errorbar(x,S_medt_pre_heat',S_SE_medt_pre_heat');
s1.Marker = 's'; s1.MarkerSize = 5; s1.MarkerEdgeColor = 'blue'; s1.MarkerFaceColor = 'blue';
hold on
s2 = errorbar(x,S_medt_pos_heat',S_SE_medt_pos_heat');
s2.Marker = 's'; s2.MarkerSize = 5; s2.MarkerEdgeColor = 'red'; s2.MarkerFaceColor = 'red';
stem(x,y,':','MarkerEdgeColor','w')  % draw vertical ref lines
title('Meditators, Heat','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
ylabel('S     ','Rotation',0'); legend('Pre','Post'); 
xticks(0:1:16); set(gca,'xtick',[])

nexttile([1,4]) % t-test labels
stem(x,y1,':','MarkerEdgeColor','w');
lblc = sprintfc('%.3f',p_medt_heat); set(gca,'ytick',[])
t = text(x, y1-.5, lblc, 'HorizontalAlignment','center', 'VerticalAlignment','middle','Rotation',90);
lblc_p = find(p_medt_heat < 0.05); % find p<0.05 indices
for i = 1:length(lblc_p)
    t(lblc_p(i)).Color = 'blue'; t(lblc_p(i)).FontWeight = 'bold';
end
xticks(0:1:16); xticklabels(thresh_labels); xlabel('Threshold'); ylabel('t-test   ','Rotation',0,'HorizontalAlignment','right');

print -depsc2 S_allthresh_medt_heat.eps  % export
