clc;
close all;
clearvars;

% load correlation matrices from nilearn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHOOSE NILEARN MATRICES VERSION
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0418/
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0417/
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0416/
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0407_denoised/
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0405_denoised_GS/
%cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0404_denoised/
cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0321
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%
% CHOOSE THRESH
thresh_prop = 0.3;
thresh_abs = 0.35;
%%%%%%%%%%%%%%%%%%%%

addpath '/Volumes/Seagate_Desktop_Drive/mfc/code'
% load('ctrl_corr_matrices.mat');
% load('medt_corr_matrices.mat');
% load('ctrl_heat_corr_matrices.mat');
% load('ctrl_neut_corr_matrices.mat');
% load('medt_heat_corr_matrices.mat');
% load('medt_neut_corr_matrices.mat');
load('ctrl_premanip_neut_matrices.mat');
load('ctrl_premanip_heat_matrices.mat');
load('ctrl_posmanip_neut_matrices.mat');
load('ctrl_posmanip_heat_matrices.mat');
load('medt_premanip_neut_matrices.mat');
load('medt_premanip_heat_matrices.mat');
load('medt_posmanip_neut_matrices.mat');
load('medt_posmanip_heat_matrices.mat');

% Load atlas labels
[~,labels] = readvars('atlas_labels_harvardoxford_cort-maxprob-thr25-2mm.csv');
labels = labels(2:49,:);

% permute dimensions

% ctrl_corr_matrices = permute(ctrl_corr_matrices,[3 2 1]);
% medt_corr_matrices = permute(medt_corr_matrices,[3 2 1]);
% ctrl_heat_corr_matrices = permute(ctrl_heat_corr_matrices,[3 2 1]);
% ctrl_neut_corr_matrices = permute(ctrl_neut_corr_matrices,[3 2 1]);
% medt_heat_corr_matrices = permute(medt_heat_corr_matrices,[3 2 1]);
% medt_neut_corr_matrices = permute(medt_neut_corr_matrices,[3 2 1]);
ctrl_premanip_neut_corr_matrices = permute(ctrl_premanip_neut_corr_matrices,[3 2 1]);
ctrl_premanip_heat_corr_matrices = permute(ctrl_premanip_heat_corr_matrices,[3 2 1]);
ctrl_posmanip_neut_corr_matrices = permute(ctrl_posmanip_neut_corr_matrices,[3 2 1]);
ctrl_posmanip_heat_corr_matrices = permute(ctrl_posmanip_heat_corr_matrices,[3 2 1]);
medt_premanip_neut_corr_matrices = permute(medt_premanip_neut_corr_matrices,[3 2 1]);
medt_premanip_heat_corr_matrices = permute(medt_premanip_heat_corr_matrices,[3 2 1]);
medt_posmanip_neut_corr_matrices = permute(medt_posmanip_neut_corr_matrices,[3 2 1]);
medt_posmanip_heat_corr_matrices = permute(medt_posmanip_heat_corr_matrices,[3 2 1]);

% average into single group corr matrix

% ctrl_neut_gp_matrix = mean(ctrl_neut_corr_matrices,3);
% ctrl_heat_gp_matrix = mean(ctrl_heat_corr_matrices,3);
% medt_neut_gp_matrix = mean(medt_neut_corr_matrices,3);
% medt_heat_gp_matrix = mean(medt_heat_corr_matrices,3);
ctrl_premanip_neut_gp_matrix = mean(ctrl_premanip_neut_corr_matrices,3);
ctrl_premanip_heat_gp_matrix = mean(ctrl_premanip_heat_corr_matrices,3);
ctrl_posmanip_neut_gp_matrix = mean(ctrl_posmanip_neut_corr_matrices,3);
ctrl_posmanip_heat_gp_matrix = mean(ctrl_posmanip_heat_corr_matrices,3);
medt_premanip_neut_gp_matrix = mean(medt_premanip_neut_corr_matrices,3);
medt_premanip_heat_gp_matrix = mean(medt_premanip_heat_corr_matrices,3);
medt_posmanip_neut_gp_matrix = mean(medt_posmanip_neut_corr_matrices,3);
medt_posmanip_heat_gp_matrix = mean(medt_posmanip_heat_corr_matrices,3);

% make diagonal zeroes

% ctrl_neut_gp_matrix = ctrl_neut_gp_matrix - diag(diag(ctrl_neut_gp_matrix));
% ctrl_heat_gp_matrix = ctrl_heat_gp_matrix - diag(diag(ctrl_heat_gp_matrix));
% medt_neut_gp_matrix = medt_neut_gp_matrix - diag(diag(medt_neut_gp_matrix));
% medt_heat_gp_matrix = medt_heat_gp_matrix - diag(diag(medt_heat_gp_matrix));
ctrl_premanip_neut_gp_matrix = ctrl_premanip_neut_gp_matrix - diag(diag(ctrl_premanip_neut_gp_matrix));
ctrl_premanip_heat_gp_matrix = ctrl_premanip_heat_gp_matrix - diag(diag(ctrl_premanip_heat_gp_matrix));
ctrl_posmanip_neut_gp_matrix = ctrl_posmanip_neut_gp_matrix - diag(diag(ctrl_posmanip_neut_gp_matrix));
ctrl_posmanip_heat_gp_matrix = ctrl_posmanip_heat_gp_matrix - diag(diag(ctrl_posmanip_heat_gp_matrix));
medt_premanip_neut_gp_matrix = medt_premanip_neut_gp_matrix - diag(diag(medt_premanip_neut_gp_matrix));
medt_premanip_heat_gp_matrix = medt_premanip_heat_gp_matrix - diag(diag(medt_premanip_heat_gp_matrix));
medt_posmanip_neut_gp_matrix = medt_posmanip_neut_gp_matrix - diag(diag(medt_posmanip_neut_gp_matrix));
medt_posmanip_heat_gp_matrix = medt_posmanip_heat_gp_matrix - diag(diag(medt_posmanip_heat_gp_matrix));

% plot visual test
figure; tiledlayout(2,4); sgtitle('Correlation Matrices');
nexttile; imagesc(ctrl_premanip_neut_gp_matrix); title('Controls, Premanip, Neut'); colormap('jet'); 
nexttile; imagesc(ctrl_premanip_heat_gp_matrix); title('Controls, Premanip, Heat'); colormap('jet');
nexttile; imagesc(ctrl_posmanip_neut_gp_matrix); title('Controls, Posmanip, Neut'); colormap('jet');
nexttile; imagesc(ctrl_posmanip_heat_gp_matrix); title('Controls, Posmanip, Heat'); colormap('jet');
nexttile; imagesc(medt_premanip_neut_gp_matrix); title('Meditatr, Premanip, Neut'); colormap('jet'); xticks(1:48); xticklabels(labels); xtickangle(45); yticks(1:48); yticklabels(labels); ytickangle(45);
nexttile; imagesc(medt_premanip_heat_gp_matrix); title('Meditatr, Premanip, Heat'); colormap('jet'); 
nexttile; imagesc(medt_posmanip_neut_gp_matrix); title('Meditatr, Posmanip, Neut'); colormap('jet'); 
nexttile; imagesc(medt_posmanip_heat_gp_matrix); title('Meditatr, Posmanip, Heat'); colormap('jet'); 

%% threshold 

% Proportional thresholding - controls, neutral, pre+post
W_thr_prop_ctrl_premanip_neut = threshold_proportional(ctrl_premanip_neut_gp_matrix, thresh_prop);
W_thr_prop_ctrl_posmanip_neut = threshold_proportional(ctrl_posmanip_neut_gp_matrix, thresh_prop);
% Absolute thresholding - controls, neutral, pre+post
W_thr_abs_ctrl_premanip_neut = threshold_absolute(ctrl_premanip_neut_gp_matrix, thresh_abs);
W_thr_abs_ctrl_posmanip_neut = threshold_absolute(ctrl_posmanip_neut_gp_matrix, thresh_abs);
% Proportional thresholding - medit, neutral, pre+post
W_thr_prop_medt_premanip_neut = threshold_proportional(medt_premanip_neut_gp_matrix, thresh_prop);
W_thr_prop_medt_posmanip_neut = threshold_proportional(medt_posmanip_neut_gp_matrix, thresh_prop);
% Absolute thresholding - medit, neutral, pre+post
W_thr_abs_medt_premanip_neut = threshold_absolute(medt_premanip_neut_gp_matrix, thresh_abs);
W_thr_abs_medt_posmanip_neut = threshold_absolute(medt_posmanip_neut_gp_matrix, thresh_abs);

% Proportional thresholding - controls, heat, pre+post
W_thr_prop_ctrl_premanip_heat = threshold_proportional(ctrl_premanip_heat_gp_matrix, thresh_prop);
W_thr_prop_ctrl_posmanip_heat = threshold_proportional(ctrl_posmanip_heat_gp_matrix, thresh_prop);
% Absolute thresholding - controls, heat, pre+post
W_thr_abs_ctrl_premanip_heat = threshold_absolute(ctrl_premanip_heat_gp_matrix, thresh_abs);
W_thr_abs_ctrl_posmanip_heat = threshold_absolute(ctrl_posmanip_heat_gp_matrix, thresh_abs);
% Proportional thresholding - medit, heat, pre+post
W_thr_prop_medt_premanip_heat = threshold_proportional(medt_premanip_heat_gp_matrix, thresh_prop);
W_thr_prop_medt_posmanip_heat = threshold_proportional(medt_posmanip_heat_gp_matrix, thresh_prop);
% Absolute thresholding - medit, heat, pre+post
W_thr_abs_medt_premanip_heat = threshold_absolute(medt_premanip_heat_gp_matrix, thresh_abs);
W_thr_abs_medt_posmanip_heat = threshold_absolute(medt_posmanip_heat_gp_matrix, thresh_abs);

figure; tiledlayout(4,4); 
str = sprintf('Thresholded Correlation Matrices, all sub, all conditions at %.2f Abs. & %.2f Prop. Correlation Thresholds',thresh_abs, thresh_prop);
sgtitle(str);
%figure; tiledlayout(2,2); sgtitle('Controls, Neutral Scans'); % controls, neutral
nexttile; imagesc(W_thr_prop_ctrl_premanip_neut); colormap('jet'); title('Ctrl Neut PreManip PropThresh');
nexttile; imagesc(W_thr_prop_ctrl_posmanip_neut); colormap('jet'); title('Ctrl Neut PosManip PropThresh');
nexttile; imagesc(W_thr_abs_ctrl_premanip_neut); colormap('jet'); title('Ctrl Neut PreManip AbsThresh');
nexttile; imagesc(W_thr_abs_ctrl_posmanip_neut); colormap('jet'); title('Ctrl Neut PosManip AbsThresh');
%figure; tiledlayout(2,2); sgtitle('Meditators, Neutral Scans'); % meditators, neutral
nexttile; imagesc(W_thr_prop_medt_premanip_neut); colormap('jet'); title('Medt Neut PreManip PropThresh');
nexttile; imagesc(W_thr_prop_medt_posmanip_neut); colormap('jet'); title('Medt Neut PosManip PropThresh');
nexttile; imagesc(W_thr_abs_medt_premanip_neut); colormap('jet'); title('Medt Neut PreManip AbsThresh');
nexttile; imagesc(W_thr_abs_medt_posmanip_neut); colormap('jet'); title('Medt Neut PosManip AbsThresh');
%figure; tiledlayout(2,2); sgtitle('Controls, Heat Scans'); % controls, heat
nexttile; imagesc(W_thr_prop_ctrl_premanip_heat); colormap('jet'); title('Ctrl Heat PreManip PropThresh');
nexttile; imagesc(W_thr_prop_ctrl_posmanip_heat); colormap('jet'); title('Ctrl Heat PosManip PropThresh');
nexttile; imagesc(W_thr_abs_ctrl_premanip_heat); colormap('jet'); title('Ctrl Heat PreManip AbsThresh');
nexttile; imagesc(W_thr_abs_ctrl_posmanip_heat); colormap('jet'); title('Ctrl Heat PosManip AbsThresh');
%figure; tiledlayout(2,2); sgtitle('Meditators, Heat Scans'); % meditators, heat
nexttile; imagesc(W_thr_prop_medt_premanip_heat); colormap('jet'); title('Medt Heat PreManip PropThresh');
nexttile; imagesc(W_thr_prop_medt_posmanip_heat); colormap('jet'); title('Medt Heat PosManip PropThresh');
nexttile; imagesc(W_thr_abs_medt_premanip_heat); colormap('jet'); title('Medt Heat PreManip AbsThresh');
nexttile; imagesc(W_thr_abs_medt_posmanip_heat); colormap('jet'); title('Medt Heat PosManip AbsThresh');

%% binarize

% Controls - Proportional - Neutral
W_bin_prop_ctrl_premanip_neut = weight_conversion(W_thr_prop_ctrl_premanip_neut, 'binarize');
W_bin_prop_ctrl_posmanip_neut = weight_conversion(W_thr_prop_ctrl_posmanip_neut, 'binarize');
% Controls - Absolute - Neutral
W_bin_abs_ctrl_premanip_neut = weight_conversion(W_thr_abs_ctrl_premanip_neut, 'binarize');
W_bin_abs_ctrl_posmanip_neut = weight_conversion(W_thr_abs_ctrl_posmanip_neut, 'binarize');
% Meditators - Proportional - Neutral
W_bin_prop_medt_premanip_neut = weight_conversion(W_thr_prop_medt_premanip_neut, 'binarize');
W_bin_prop_medt_posmanip_neut = weight_conversion(W_thr_prop_medt_posmanip_neut, 'binarize');
% Meditators - Absolute - Neutral
W_bin_abs_medt_premanip_neut = weight_conversion(W_thr_abs_medt_premanip_neut, 'binarize');
W_bin_abs_medt_posmanip_neut = weight_conversion(W_thr_abs_medt_posmanip_neut, 'binarize');
% Controls - Proportional - heat
W_bin_prop_ctrl_premanip_heat = weight_conversion(W_thr_prop_ctrl_premanip_heat, 'binarize');
W_bin_prop_ctrl_posmanip_heat = weight_conversion(W_thr_prop_ctrl_posmanip_heat, 'binarize');
% Controls - Absolute - heat
W_bin_abs_ctrl_premanip_heat = weight_conversion(W_thr_abs_ctrl_premanip_heat, 'binarize');
W_bin_abs_ctrl_posmanip_heat = weight_conversion(W_thr_abs_ctrl_posmanip_heat, 'binarize');
% Meditators - Proportional - heat
W_bin_prop_medt_premanip_heat = weight_conversion(W_thr_prop_medt_premanip_heat, 'binarize');
W_bin_prop_medt_posmanip_heat = weight_conversion(W_thr_prop_medt_posmanip_heat, 'binarize');
% Meditators - Absolute - heat
W_bin_abs_medt_premanip_heat = weight_conversion(W_thr_abs_medt_premanip_heat, 'binarize');
W_bin_abs_medt_posmanip_heat = weight_conversion(W_thr_abs_medt_posmanip_heat, 'binarize');

% figure; tiledlayout(2,2); sgtitle('Controls, Neutral Scans');
% nexttile; imagesc(W_bin_prop_ctrl_premanip_neut); colormap('jet'); title('Ctr Pre Prop');
% nexttile; imagesc(W_bin_prop_ctrl_posmanip_neut); colormap('jet'); title('Ctr Pos Prop');
% nexttile; imagesc(W_bin_abs_ctrl_premanip_neut); colormap('jet'); title('Ctr Pre Abs');
% nexttile; imagesc(W_bin_abs_ctrl_posmanip_neut); colormap('jet'); title('Ctr Pos Abs');
% figure; tiledlayout(2,2); sgtitle('Meditators, Neutral Scans');
% nexttile; imagesc(W_bin_prop_medt_premanip_neut); colormap('jet'); title('Medt Pre Prop');
% nexttile; imagesc(W_bin_prop_medt_posmanip_neut); colormap('jet'); title('Medt Pos Prop');
% nexttile; imagesc(W_bin_abs_medt_premanip_neut); colormap('jet'); title('Medt Pre Abs');
% nexttile; imagesc(W_bin_abs_medt_posmanip_neut); colormap('jet'); title('Medt Pos Abs');
% figure; tiledlayout(2,2); sgtitle('Controls, Heat Scans');
% nexttile; imagesc(W_bin_prop_ctrl_premanip_heat); colormap('jet'); title('Ctr Pre Prop');
% nexttile; imagesc(W_bin_prop_ctrl_posmanip_heat); colormap('jet'); title('Ctr Pos Prop');
% nexttile; imagesc(W_bin_abs_ctrl_premanip_heat); colormap('jet'); title('Ctr Pre Abs');
% nexttile; imagesc(W_bin_abs_ctrl_posmanip_heat); colormap('jet'); title('Ctr Pos Abs');
% figure; tiledlayout(2,2); sgtitle('Meditators, Heat Scans');
% nexttile; imagesc(W_bin_prop_medt_premanip_heat); colormap('jet'); title('Medt Pre Prop');
% nexttile; imagesc(W_bin_prop_medt_posmanip_heat); colormap('jet'); title('Medt Pos Prop');
% nexttile; imagesc(W_bin_abs_medt_premanip_heat); colormap('jet'); title('Medt Pre Abs');
% nexttile; imagesc(W_bin_abs_medt_posmanip_heat); colormap('jet'); title('Medt Pos Abs');

figure; tiledlayout(4,4); 
str = sprintf('Adjacency Matrices (Prop/Abs) (Control/Meditator) (Pre/Post) (Heat/Neutral) at %.2f Abs. & %.2f Prop. Corr. Thresholds',thresh_abs, thresh_prop);
sgtitle(str);
nexttile; imagesc(W_bin_prop_ctrl_premanip_neut); colormap('jet'); title('Ctrl Neut Pre Prop');
nexttile; imagesc(W_bin_prop_ctrl_posmanip_neut); colormap('jet'); title('Ctrl Neut Pos Prop');
nexttile; imagesc(W_bin_abs_ctrl_premanip_neut); colormap('jet'); title('Ctrl Neut Pre Abs');
nexttile; imagesc(W_bin_abs_ctrl_posmanip_neut); colormap('jet'); title('Ctrl Neut Pos Abs');
nexttile; imagesc(W_bin_prop_medt_premanip_neut); colormap('jet'); title('Medt Neut Pre Prop');
nexttile; imagesc(W_bin_prop_medt_posmanip_neut); colormap('jet'); title('Medt Neut Pos Prop');
nexttile; imagesc(W_bin_abs_medt_premanip_neut); colormap('jet'); title('Medt Neut Pre Abs');
nexttile; imagesc(W_bin_abs_medt_posmanip_neut); colormap('jet'); title('Medt Neut Pos Abs');
nexttile; imagesc(W_bin_prop_ctrl_premanip_heat); colormap('jet'); title('Ctr Heat Pre Prop');
nexttile; imagesc(W_bin_prop_ctrl_posmanip_heat); colormap('jet'); title('Ctr Heat Pos Prop');
nexttile; imagesc(W_bin_abs_ctrl_premanip_heat); colormap('jet'); title('Ctr Heat Pre Abs');
nexttile; imagesc(W_bin_abs_ctrl_posmanip_heat); colormap('jet'); title('Ctr Heat Pos Abs');
nexttile; imagesc(W_bin_prop_medt_premanip_heat); colormap('jet'); title('Medt Heat Pre Prop');
nexttile; imagesc(W_bin_prop_medt_posmanip_heat); colormap('jet'); title('Medt Heat Pos Prop');
nexttile; imagesc(W_bin_abs_medt_premanip_heat); colormap('jet'); title('Medt Heat Pre Abs');
nexttile; imagesc(W_bin_abs_medt_posmanip_heat); colormap('jet'); title('Medt Heat Pos Abs');

% cmap = gray(2);
% figure; imagesc(W_bin_abs_medt_posmanip_heat); colormap(cmap);
%% degree distribution

% Controls - Proportional - Neutral
ddist_prop_ctrl_premanip_neut = degrees_und(W_bin_prop_ctrl_premanip_neut);
ddist_prop_ctrl_posmanip_neut = degrees_und(W_bin_prop_ctrl_posmanip_neut);
% Controls - Absolute - Neutral
ddist_abs_ctrl_premanip_neut = degrees_und(W_bin_abs_ctrl_premanip_neut);
ddist_abs_ctrl_posmanip_neut = degrees_und(W_bin_abs_ctrl_posmanip_neut);
% Meditators - Proportional - Neutral
ddist_prop_medt_premanip_neut = degrees_und(W_bin_prop_medt_premanip_neut);
ddist_prop_medt_posmanip_neut = degrees_und(W_bin_prop_medt_posmanip_neut);
% Meditators - Absolute - Neutral
ddist_abs_medt_premanip_neut = degrees_und(W_bin_abs_medt_premanip_neut);
ddist_abs_medt_posmanip_neut = degrees_und(W_bin_abs_medt_posmanip_neut);
% Controls - Proportional - heat
ddist_prop_ctrl_premanip_heat = degrees_und(W_bin_prop_ctrl_premanip_heat);
ddist_prop_ctrl_posmanip_heat = degrees_und(W_bin_prop_ctrl_posmanip_heat);
% Controls - Absolute - heat
ddist_abs_ctrl_premanip_heat = degrees_und(W_bin_abs_ctrl_premanip_heat);
ddist_abs_ctrl_posmanip_heat = degrees_und(W_bin_abs_ctrl_posmanip_heat);
% Meditators - Proportional - heat
ddist_prop_medt_premanip_heat = degrees_und(W_bin_prop_medt_premanip_heat);
ddist_prop_medt_posmanip_heat = degrees_und(W_bin_prop_medt_posmanip_heat);
% Meditators - Absolute - heat
ddist_abs_medt_premanip_heat = degrees_und(W_bin_abs_medt_premanip_heat);
ddist_abs_medt_posmanip_heat = degrees_und(W_bin_abs_medt_posmanip_heat);

% % Plot Proportional Neutral - pre/post ctrl/medt
% figure; tiledlayout('flow');
% nexttile
% histogram(ddist_prop_medt_premanip_neut,'BinMethod','integers','FaceAlpha',0.9);
% hold on
% histogram(ddist_prop_medt_posmanip_neut,'BinMethod','integers','FaceAlpha',0.7);
% title('DDist Meditators Proportional Pre/Post-Manip');
% xlabel('Degree of node');
% ylabel('Number of nodes');
% legend('pre','post','Location','Northeast');
% %axis([0 100 0 250]);
% hold off
% nexttile
% histogram(ddist_prop_ctrl_premanip_neut,'BinMethod','integers','FaceAlpha',0.9);
% hold on
% histogram(ddist_prop_ctrl_posmanip_neut,'BinMethod','integers','FaceAlpha',0.7);
% title('DDist Controls Proportional Pre/Post-Manip');
% xlabel('Degree of node');
% ylabel('Number of nodes');
% legend('pre','post','Location','Northeast');
% %axis([0 100 0 250]);
% hold off

% Plot Absolute Neutral - pre/post ctrl/medt
figure
tiledlayout(2,2)
nexttile
histogram(ddist_abs_ctrl_premanip_neut,'BinMethod','integers','FaceAlpha',0.9);
hold on
histogram(ddist_abs_ctrl_posmanip_neut,'BinMethod','integers','FaceAlpha',0.7);
title('Controls, Neutral Scans');
xlabel('Degree of node');
ylabel('Number of nodes');
legend('pre','post','Location','Northeast');
axis([0 21 0 10]);
hold off
nexttile
histogram(ddist_abs_medt_premanip_neut,'BinMethod','integers','FaceAlpha',0.9);
hold on
histogram(ddist_abs_medt_posmanip_neut,'BinMethod','integers','FaceAlpha',0.7);
title('Meditators, Neutral Scans');
xlabel('Degree of node');
ylabel('Number of nodes');
legend('pre','post','Location','Northeast');
axis([0 21 0 10]);
hold off
str = sprintf('Degree Distribution - (Absolute Threshhold) at %.2f Abs. Corr. Threshold',thresh_abs);
sgtitle(str);

% Plot Absolute Heat - pre/post ctrl/medt
nexttile
histogram(ddist_abs_ctrl_premanip_heat,'BinMethod','integers','FaceAlpha',0.9);
hold on
histogram(ddist_abs_ctrl_posmanip_heat,'BinMethod','integers','FaceAlpha',0.7);
title('Controls, Heat Scans');
xlabel('Degree of node');
ylabel('Number of nodes');
legend('pre','post','Location','Northeast');
axis([0 21 0 10]);
hold off
nexttile
histogram(ddist_abs_medt_premanip_heat,'BinMethod','integers','FaceAlpha',0.9);
hold on
histogram(ddist_abs_medt_posmanip_heat,'BinMethod','integers','FaceAlpha',0.7);
title('Meditators, Heat Scans');
xlabel('Degree of node');
ylabel('Number of nodes');
legend('pre','post','Location','Northeast');
axis([0 21 0 10]);
hold off
str = sprintf('Degree Distribution - (Absolute Threshhold) at %.2f Abs. Corr. Threshold',thresh_abs);
sgtitle(str);

%% group-level entropy

% Absolute Threshold - all 8 tasks/groups
entropy_bins = 120;
% Controls, Pre, Neutral
f = figure('visible','off'); %initialize fig invisible
hist1 = histogram(ddist_abs_ctrl_premanip_neut,entropy_bins);
hist1.Normalization = 'probability';
%bincenters = (hist1.BinEdges(2:end) + hist1.BinEdges(1:end-1))/2;
binvalues = hist1.Values;
nS_abs_ctrl_premanip_neut = sum(-(binvalues(binvalues>0).*(log(binvalues(binvalues>0)))));
clear f; clear hist1;
% Controls, Post, Neutral
f = figure('visible','off'); %initialize fig invisible
hist1 = histogram(ddist_abs_ctrl_posmanip_neut,entropy_bins);
hist1.Normalization = 'probability';
%bincenters = (hist1.BinEdges(2:end) + hist1.BinEdges(1:end-1))/2;
binvalues = hist1.Values;
nS_abs_ctrl_posmanip_neut = sum(-(binvalues(binvalues>0).*(log(binvalues(binvalues>0)))));
clear f; clear hist1;
% Controls, Pre, Heat
f = figure('visible','off'); %initialize fig invisible
hist1 = histogram(ddist_abs_ctrl_premanip_heat,entropy_bins);
hist1.Normalization = 'probability';
%bincenters = (hist1.BinEdges(2:end) + hist1.BinEdges(1:end-1))/2;
binvalues = hist1.Values;
nS_abs_ctrl_premanip_heat = sum(-(binvalues(binvalues>0).*(log(binvalues(binvalues>0)))));
clear f; clear hist1;
% Controls, Post, Heat
f = figure('visible','off'); %initialize fig invisible
hist1 = histogram(ddist_abs_ctrl_posmanip_heat,entropy_bins);
hist1.Normalization = 'probability';
%bincenters = (hist1.BinEdges(2:end) + hist1.BinEdges(1:end-1))/2;
binvalues = hist1.Values;
nS_abs_ctrl_posmanip_heat = sum(-(binvalues(binvalues>0).*(log(binvalues(binvalues>0)))));
clear f; clear hist1;
% Meditators, Pre, Neutral
f = figure('visible','off'); %initialize fig invisible
hist1 = histogram(ddist_abs_medt_premanip_neut,entropy_bins);
hist1.Normalization = 'probability';
%bincenters = (hist1.BinEdges(2:end) + hist1.BinEdges(1:end-1))/2;
binvalues = hist1.Values;
nS_abs_medt_premanip_neut = sum(-(binvalues(binvalues>0).*(log(binvalues(binvalues>0)))));
clear f; clear hist1;
% Meditators, Post, Neutral
f = figure('visible','off'); %initialize fig invisible
hist1 = histogram(ddist_abs_medt_posmanip_neut,entropy_bins);
hist1.Normalization = 'probability';
%bincenters = (hist1.BinEdges(2:end) + hist1.BinEdges(1:end-1))/2;
binvalues = hist1.Values;
nS_abs_medt_posmanip_neut = sum(-(binvalues(binvalues>0).*(log(binvalues(binvalues>0)))));
clear f; clear hist1;
% Meditators, Pre, Heat
f = figure('visible','off'); %initialize fig invisible
hist1 = histogram(ddist_abs_medt_premanip_heat,entropy_bins);
hist1.Normalization = 'probability';
%bincenters = (hist1.BinEdges(2:end) + hist1.BinEdges(1:end-1))/2;
binvalues = hist1.Values;
nS_abs_medt_premanip_heat = sum(-(binvalues(binvalues>0).*(log(binvalues(binvalues>0)))));
clear f; clear hist1;
% Meditators, Post, Heat
f = figure('visible','off'); %initialize fig invisible
hist1 = histogram(ddist_abs_medt_posmanip_heat,entropy_bins);
hist1.Normalization = 'probability';
%bincenters = (hist1.BinEdges(2:end) + hist1.BinEdges(1:end-1))/2;
binvalues = hist1.Values;
nS_abs_medt_posmanip_heat = sum(-(binvalues(binvalues>0).*(log(binvalues(binvalues>0)))));
clear f; clear hist1;



%% plot shannon entropy (group level)

% x = categorical({'Controls Neutral','Controls Heat','Meditators Neutral','Meditators Heat'});
% y = [nS_abs_ctrl_premanip_neut nS_abs_ctrl_premanip_heat nS_abs_medt_premanip_neut nS_abs_medt_premanip_heat; ...
%     nS_abs_ctrl_posmanip_neut nS_abs_ctrl_posmanip_heat nS_abs_medt_posmanip_neut nS_abs_medt_posmanip_heat];
% figure; 
% b=bar(x,y); legend('pre','post','Location','north'); 
% str = sprintf('Shannon Entropy, Group-level at %.2f Abs. Corr. Threshold',thresh_abs);
% sgtitle(str); ylabel('S');
% xtips1 = b(1).XEndPoints;
% ytips1 = b(1).YEndPoints;
% labels1 = string(b(1).YData);
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
% xtips2 = b(2).XEndPoints;
% ytips2 = b(2).YEndPoints;
% labels2 = string(b(2).YData);
% text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom')

%% Individual sub matrices, ddist, entropy

% get # rois
nrois = size(ctrl_premanip_heat_corr_matrices);
nrois = nrois(1);

% Create condition matrices w/1 matrix per subject (i.e. flatten 2 into 1)
% Controls (20 sub, 160 scans)
% premanip-heat
ctrl_premanip_heat_sub = zeros(nrois,nrois,20);
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,ctrl_premanip_heat_corr_matrices(:,:,i),ctrl_premanip_heat_corr_matrices(:,:,i+1));
    sub_ctrl_premanip_heat(:,:,i) = mean(sub,3);
    sub_ctrl_premanip_heat(:,:,i) = sub_ctrl_premanip_heat(:,:,i) - diag(diag(sub_ctrl_premanip_heat(:,:,i)));
end
sub_ctrl_premanip_heat = sub_ctrl_premanip_heat(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% posmanip-heat
ctrl_posmanip_heat_sub = zeros(nrois,nrois,20);
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,ctrl_posmanip_heat_corr_matrices(:,:,i),ctrl_posmanip_heat_corr_matrices(:,:,i+1));
    sub_ctrl_posmanip_heat(:,:,i) = mean(sub,3);
    sub_ctrl_posmanip_heat(:,:,i) = sub_ctrl_posmanip_heat(:,:,i) - diag(diag(sub_ctrl_posmanip_heat(:,:,i)));
end
sub_ctrl_posmanip_heat = sub_ctrl_posmanip_heat(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% premanip-neut
ctrl_premanip_neut_sub = zeros(nrois,nrois,20);
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,ctrl_premanip_neut_corr_matrices(:,:,i),ctrl_premanip_neut_corr_matrices(:,:,i+1));
    sub_ctrl_premanip_neut(:,:,i) = mean(sub,3);
    sub_ctrl_premanip_neut(:,:,i) = sub_ctrl_premanip_neut(:,:,i) - diag(diag(sub_ctrl_premanip_neut(:,:,i)));
end
sub_ctrl_premanip_neut = sub_ctrl_premanip_neut(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% posmanip-neut
ctrl_posmanip_neut_sub = zeros(nrois,nrois,20);
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,ctrl_posmanip_neut_corr_matrices(:,:,i),ctrl_posmanip_neut_corr_matrices(:,:,i+1));
    sub_ctrl_posmanip_neut(:,:,i) = mean(sub,3);
    sub_ctrl_posmanip_neut(:,:,i) = sub_ctrl_posmanip_neut(:,:,i) - diag(diag(sub_ctrl_posmanip_neut(:,:,i)));
end
sub_ctrl_posmanip_neut = sub_ctrl_posmanip_neut(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% Meditators (20 sub, 160 scans)
% premanip-heat
medt_premanip_heat_sub = zeros(nrois,nrois,20);
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,medt_premanip_heat_corr_matrices(:,:,i),medt_premanip_heat_corr_matrices(:,:,i+1));
    sub_medt_premanip_heat(:,:,i) = mean(sub,3);
    sub_medt_premanip_heat(:,:,i) = sub_medt_premanip_heat(:,:,i) - diag(diag(sub_medt_premanip_heat(:,:,i)));
end
sub_medt_premanip_heat = sub_medt_premanip_heat(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% posmanip-heat
medt_posmanip_heat_sub = zeros(nrois,nrois,20);
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,medt_posmanip_heat_corr_matrices(:,:,i),medt_posmanip_heat_corr_matrices(:,:,i+1));
    sub_medt_posmanip_heat(:,:,i) = mean(sub,3);
    sub_medt_posmanip_heat(:,:,i) = sub_medt_posmanip_heat(:,:,i) - diag(diag(sub_medt_posmanip_heat(:,:,i)));
end
sub_medt_posmanip_heat = sub_medt_posmanip_heat(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% premanip-neut
medt_premanip_neut_sub = zeros(nrois,nrois,20);
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,medt_premanip_neut_corr_matrices(:,:,i),medt_premanip_neut_corr_matrices(:,:,i+1));
    sub_medt_premanip_neut(:,:,i) = mean(sub,3);
    sub_medt_premanip_neut(:,:,i) = sub_medt_premanip_neut(:,:,i) - diag(diag(sub_medt_premanip_neut(:,:,i)));
end
sub_medt_premanip_neut = sub_medt_premanip_neut(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);
% posmanip-neut
medt_posmanip_neut_sub = zeros(nrois,nrois,20);
for i = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]
    sub = cat(3,medt_posmanip_neut_corr_matrices(:,:,i),medt_posmanip_neut_corr_matrices(:,:,i+1));
    sub_medt_posmanip_neut(:,:,i) = mean(sub,3);
    sub_medt_posmanip_neut(:,:,i) = sub_medt_posmanip_neut(:,:,i) - diag(diag(sub_medt_posmanip_neut(:,:,i)));
end
sub_medt_posmanip_neut = sub_medt_posmanip_neut(:,:,[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]);

% get subject-level thresholded matrices, ddist, and nS

% controls, premanip, heat
sub_bin_abs_ctrl_premanip_heat = zeros(nrois,nrois,20);
sub_ddist_ctrl_premanip_heat = zeros(1,nrois);
sub_nS_ctrl_premanip_heat = zeros(1,20);
for i = 1:20
    unthresh_matrix = sub_ctrl_premanip_heat(:,:,i);
    [w_bin_abs, ddist, nS] = preproc_corr_matrix(unthresh_matrix, thresh_abs);
    sub_bin_abs_ctrl_premanip_heat(:,:,i) = w_bin_abs;
    sub_ddist_ctrl_premanip_heat(i,:) = ddist;
    sub_nS_ctrl_premanip_heat(1,i) = nS;
end
clear w_bin_abs; clear ddist; clear nS; clear unthresh_matrix; 
% controls, posmanip, heat
sub_bin_abs_ctrl_posmanip_heat = zeros(nrois,nrois,20);
sub_ddist_ctrl_posmanip_heat = zeros(1,nrois);
sub_nS_ctrl_posmanip_heat = zeros(1,20);
for i = 1:20
    unthresh_matrix = sub_ctrl_posmanip_heat(:,:,i);
    [w_bin_abs, ddist, nS] = preproc_corr_matrix(unthresh_matrix, thresh_abs);
    sub_bin_abs_ctrl_posmanip_heat(:,:,i) = w_bin_abs;
    sub_ddist_ctrl_posmanip_heat(i,:) = ddist;
    sub_nS_ctrl_posmanip_heat(1,i) = nS;
end
% controls, premanip, neut
sub_bin_abs_ctrl_premanip_neut = zeros(nrois,nrois,20);
sub_ddist_ctrl_premanip_neut = zeros(1,nrois);
sub_nS_ctrl_premanip_neut = zeros(1,20);
for i = 1:20
    unthresh_matrix = sub_ctrl_premanip_neut(:,:,i);
    [w_bin_abs, ddist, nS] = preproc_corr_matrix(unthresh_matrix, thresh_abs);
    sub_bin_abs_ctrl_premanip_neut(:,:,i) = w_bin_abs;
    sub_ddist_ctrl_premanip_neut(i,:) = ddist;
    sub_nS_ctrl_premanip_neut(1,i) = nS;
end
% controls, posmanip, neut
sub_bin_abs_ctrl_posmanip_neut = zeros(nrois,nrois,20);
sub_ddist_ctrl_posmanip_neut = zeros(1,nrois);
sub_nS_ctrl_posmanip_neut = zeros(1,20);
for i = 1:20
    unthresh_matrix = sub_ctrl_posmanip_neut(:,:,i);
    [w_bin_abs, ddist, nS] = preproc_corr_matrix(unthresh_matrix, thresh_abs);
    sub_bin_abs_ctrl_posmanip_neut(:,:,i) = w_bin_abs;
    sub_ddist_ctrl_posmanip_neut(i,:) = ddist;
    sub_nS_ctrl_posmanip_neut(1,i) = nS;
end
% meditators, premanip, heat
sub_bin_abs_medt_premanip_heat = zeros(nrois,nrois,20);
sub_ddist_medt_premanip_heat = zeros(1,nrois);
sub_nS_medt_premanip_heat = zeros(1,20);
for i = 1:20
    unthresh_matrix = sub_medt_premanip_heat(:,:,i);
    [w_bin_abs, ddist, nS] = preproc_corr_matrix(unthresh_matrix, thresh_abs);
    sub_bin_abs_medt_premanip_heat(:,:,i) = w_bin_abs;
    sub_ddist_medt_premanip_heat(i,:) = ddist;
    sub_nS_medt_premanip_heat(1,i) = nS;
end
% meditators, posmanip, heat
sub_bin_abs_medt_posmanip_heat = zeros(nrois,nrois,20);
sub_ddist_medt_posmanip_heat = zeros(1,nrois);
sub_nS_medt_posmanip_heat = zeros(1,20);
for i = 1:20
    unthresh_matrix = sub_medt_posmanip_heat(:,:,i);
    [w_bin_abs, ddist, nS] = preproc_corr_matrix(unthresh_matrix, thresh_abs);
    sub_bin_abs_medt_posmanip_heat(:,:,i) = w_bin_abs;
    sub_ddist_medt_posmanip_heat(i,:) = ddist;
    sub_nS_medt_posmanip_heat(1,i) = nS;
end
% meditators, premanip, neut
sub_bin_abs_medt_premanip_neut = zeros(nrois,nrois,20);
sub_ddist_medt_premanip_neut = zeros(1,nrois);
sub_nS_medt_premanip_neut = zeros(1,20);
for i = 1:20
    unthresh_matrix = sub_medt_premanip_neut(:,:,i);
    [w_bin_abs, ddist, nS] = preproc_corr_matrix(unthresh_matrix, thresh_abs);
    sub_bin_abs_medt_premanip_neut(:,:,i) = w_bin_abs;
    sub_ddist_medt_premanip_neut(i,:) = ddist;
    sub_nS_medt_premanip_neut(1,i) = nS;
end
% meditators, posmanip, neut
sub_bin_abs_medt_posmanip_neut = zeros(nrois,nrois,20);
sub_ddist_medt_posmanip_neut = zeros(1,nrois);
sub_nS_medt_posmanip_neut = zeros(1,20);
for i = 1:20
    unthresh_matrix = sub_medt_posmanip_neut(:,:,i);
    [w_bin_abs, ddist, nS] = preproc_corr_matrix(unthresh_matrix, thresh_abs);
    sub_bin_abs_medt_posmanip_neut(:,:,i) = w_bin_abs;
    sub_ddist_medt_posmanip_neut(i,:) = ddist;
    sub_nS_medt_posmanip_neut(1,i) = nS;
end

%% plot ddist of all indiv subjects

mx = max(max(sub_ddist_medt_posmanip_neut));

% Mean, Median, Max
% PRE
% Neutral Controls
mean_NC_pre = mean(mean(sub_ddist_ctrl_premanip_neut, 'omitnan'));
median_NC_pre = mean(median(sub_ddist_ctrl_premanip_neut,2, 'omitnan'));
max_NC_pre = max(max(sub_ddist_ctrl_premanip_neut));
% Neutral Meditators
mean_NM_pre = mean(mean(sub_ddist_medt_premanip_neut, 'omitnan'));
median_NM_pre = mean(median(sub_ddist_medt_premanip_neut,2, 'omitnan'));
max_NM_pre = max(max(sub_ddist_medt_premanip_neut));
% Heat Controls
mean_HC_pre = mean(mean(sub_ddist_ctrl_premanip_heat, 'omitnan'));
median_HC_pre = mean(median(sub_ddist_ctrl_premanip_heat,2, 'omitnan'));
max_HC_pre = max(max(sub_ddist_ctrl_premanip_heat));
% Heat Meditators
mean_HM_pre = mean(mean(sub_ddist_medt_premanip_heat, 'omitnan'));
median_HM_pre = mean(median(sub_ddist_medt_premanip_heat,2, 'omitnan'));
max_HM_pre = max(max(sub_ddist_medt_premanip_heat));
% POS
% Neutral Controls
mean_NC_pos = mean(mean(sub_ddist_ctrl_posmanip_neut, 'omitnan'));
median_NC_pos = mean(median(sub_ddist_ctrl_posmanip_neut,2, 'omitnan'));
max_NC_pos = max(max(sub_ddist_ctrl_posmanip_neut));
% Neutral Meditators
mean_NM_pos = mean(mean(sub_ddist_medt_posmanip_neut, 'omitnan'));
median_NM_pos = mean(median(sub_ddist_medt_posmanip_neut,2, 'omitnan'));
max_NM_pos = max(max(sub_ddist_medt_posmanip_neut));
% Heat Controls
mean_HC_pos = mean(mean(sub_ddist_ctrl_posmanip_heat, 'omitnan'));
median_HC_pos = mean(median(sub_ddist_ctrl_posmanip_heat,2, 'omitnan'));
max_HC_pos = max(max(sub_ddist_ctrl_posmanip_heat));
% Heat Meditators
mean_HM_pos = mean(mean(sub_ddist_medt_posmanip_heat, 'omitnan'));
median_HM_pos = mean(median(sub_ddist_medt_posmanip_heat,2, 'omitnan'));
max_HM_pos = max(max(sub_ddist_medt_posmanip_heat));

f=figure; 
tiledlayout(2,2);  

nexttile % Neutral, controls
         histogram(sub_ddist_ctrl_premanip_neut,'BinMethod','integers','FaceAlpha',0.7, 'Normalization','pdf');
hold on; histogram(sub_ddist_ctrl_posmanip_neut,'BinMethod','integers','FaceAlpha',0.7, 'Normalization','pdf');
         scatter(mean_NC_pre,.01,40,'p','b','filled')
         scatter(median_NC_pre,0,40,'^','b','filled')
         scatter(max_NC_pre,0,40,'v','b','filled')
         scatter(mean_NC_pos,.01,40,'p','r','filled')
         scatter(median_NC_pos,0,40,'^','r','filled')
         scatter(max_NC_pos,0,40,'v','r','filled')
hold off; title('Controls, Neutral','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
legend('pre','post');  xlabel('Node Degree'); ylabel('Probability'); axis([0 mx+3 0 .15]);
 
nexttile % Neutral, meditators
         histogram(sub_ddist_medt_premanip_neut,'BinMethod','integers','FaceAlpha',0.7, 'Normalization','pdf');
hold on; histogram(sub_ddist_medt_posmanip_neut,'BinMethod','integers','FaceAlpha',0.7, 'Normalization','pdf');
         scatter(mean_NM_pre,.01,40,'p','b','filled')
         scatter(median_NM_pre,0,40,'^','b','filled')
         scatter(max_NM_pre,0,40,'v','b','filled')
         scatter(mean_NM_pos,.01,40,'p','r','filled')
         scatter(median_NM_pos,0,40,'^','r','filled')
         scatter(max_NM_pos,0,40,'v','r','filled')
hold off; title('Meditators, Neutral','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
xlabel('Node Degree'); ylabel('Probability'); axis([0 mx+3 0 .15]); legend('pre','post','pre, mean','pre, median','pre, max',...
    'post, mean','post, median','post, max','Location','eastoutside');
 
nexttile % Heat, controls
         histogram(sub_ddist_ctrl_premanip_heat,'BinMethod','integers','FaceAlpha',0.7, 'Normalization','pdf');
hold on; histogram(sub_ddist_ctrl_posmanip_heat,'BinMethod','integers','FaceAlpha',0.7, 'Normalization','pdf');
         scatter(mean_HC_pre,.01,40,'p','b','filled')
         scatter(median_HC_pre,0,40,'^','b','filled')
         scatter(max_HC_pre,0,40,'v','b','filled')
         scatter(mean_HC_pos,.01,40,'p','r','filled')
         scatter(median_HC_pos,0,40,'^','r','filled')
         scatter(max_HC_pos,0,40,'v','r','filled')
hold off; title('Controls, Heat','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
legend('pre','post'); xlabel('Node Degree'); ylabel('Probability'); axis([0 mx+3 0 .15]);
 
nexttile % Heat, meditators
         histogram(sub_ddist_medt_premanip_heat,'BinMethod','integers','FaceAlpha',0.7, 'Normalization','pdf');
hold on; histogram(sub_ddist_medt_posmanip_heat,'BinMethod','integers','FaceAlpha',0.7, 'Normalization','pdf'); 
         scatter(mean_HM_pre,.01,40,'p','b','filled');
         scatter(median_HM_pre,0,40,'^','b','filled')
         scatter(max_HM_pre,0,40,'v','b','filled')
         scatter(mean_HM_pos,.01,40,'p','r','filled')
         scatter(median_HM_pos,0,40,'^','r','filled')
         scatter(max_HM_pos,0,40,'v','r','filled')
hold off; title('Meditators, Heat','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri'); 
xlabel('Node Degree'); ylabel('Probability'); axis([0 mx+3 0 .15]); legend('pre','post','pre, mean','pre, median','pre, max',...
    'post, mean','post, median','post, max','Location','eastoutside');

% str = sprintf('Network Degree Distribution, Group-Level, %.2f Threshold',thresh_abs);
% sgtitle(str);

set(gcf,'Position',[100,100,1000,700]);
print -depsc2 ddist_sum.eps  % export

clear mx mean* median* max*
%figure; histogram(sub_ddist_ctrl_premanip_neut,'BinMethod','integers','FaceAlpha',0.7, 'Normalization','pdf'); xlabel('Node Degree'); ylabel('Probability');
%% variance & kurtosis, sub level
for sub = 1:size(sub_ddist_ctrl_premanip_neut,1)
    sub_kur_ctrl_premanip_neut(1,sub) = kurtosis(sub_ddist_ctrl_premanip_neut(sub,:));
    sub_var_ctrl_premanip_neut(1,sub) = var(sub_ddist_ctrl_premanip_neut(sub,:));
end
for sub = 1:size(sub_ddist_ctrl_posmanip_neut,1)
    sub_kur_ctrl_posmanip_neut(1,sub) = kurtosis(sub_ddist_ctrl_posmanip_neut(sub,:));
    sub_var_ctrl_posmanip_neut(1,sub) = var(sub_ddist_ctrl_posmanip_neut(sub,:));
end
for sub = 1:size(sub_ddist_ctrl_premanip_heat,1)
    sub_kur_ctrl_premanip_heat(1,sub) = kurtosis(sub_ddist_ctrl_premanip_heat(sub,:));
    sub_var_ctrl_premanip_heat(1,sub) = var(sub_ddist_ctrl_premanip_heat(sub,:));
end
for sub = 1:size(sub_ddist_ctrl_posmanip_heat,1)
    sub_kur_ctrl_posmanip_heat(1,sub) = kurtosis(sub_ddist_ctrl_posmanip_heat(sub,:));
    sub_var_ctrl_posmanip_heat(1,sub) = var(sub_ddist_ctrl_posmanip_heat(sub,:));
end
for sub = 1:size(sub_ddist_medt_premanip_neut,1)
    sub_kur_medt_premanip_neut(1,sub) = kurtosis(sub_ddist_medt_premanip_neut(sub,:));
    sub_var_medt_premanip_neut(1,sub) = var(sub_ddist_medt_premanip_neut(sub,:));
end
for sub = 1:size(sub_ddist_medt_posmanip_neut,1)
    sub_kur_medt_posmanip_neut(1,sub) = kurtosis(sub_ddist_medt_posmanip_neut(sub,:));
    sub_var_medt_posmanip_neut(1,sub) = var(sub_ddist_medt_posmanip_neut(sub,:));
end
for sub = 1:size(sub_ddist_medt_premanip_heat,1)
    sub_kur_medt_premanip_heat(1,sub) = kurtosis(sub_ddist_medt_premanip_heat(sub,:));
    sub_var_medt_premanip_heat(1,sub) = var(sub_ddist_medt_premanip_heat(sub,:));
end
for sub = 1:size(sub_ddist_medt_posmanip_heat,1)
    sub_kur_medt_posmanip_heat(1,sub) = kurtosis(sub_ddist_medt_posmanip_heat(sub,:));
    sub_var_medt_posmanip_heat(1,sub) = var(sub_ddist_medt_posmanip_heat(sub,:));
end

%% plot indiv sub kur/var
% Controls Pre + Post
x_ctrl = [1:20];    % sample size controls
x_medt = [21:40];   % adjust x to sample size meditators

figure % Kurtosis, Neutral
p1 = plot(x_ctrl,sub_kur_ctrl_premanip_neut,'ko',...
    'LineWidth',1,...
    'MarkerSize',6,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor','w');
hold on
p2 = plot(x_ctrl,sub_kur_ctrl_posmanip_neut,'ko',...
    'LineWidth',1,...
    'MarkerSize',6,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor','g');
p3 = plot(x_medt,sub_kur_medt_premanip_neut,'ko',...
    'LineWidth',1,...
    'MarkerSize',6,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','w');
p4 = plot(x_medt,sub_kur_medt_posmanip_neut,'ko',...
    'LineWidth',1,...
    'MarkerSize',6,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b');
hold off;
x_fullaxis = linspace(0,40,50);         %change number of subjects
str = sprintf('Kurtosis of Degree Distribution, Pre|Post Meditators|Controls, Neutral at %.2f Abs. Corr. Threshold',thresh_abs);
title(str);
xlabel('Controls      -      (subjects)      -      Meditators');
ylabel('Kurtosis');
legend({'Ctrl Pre','Ctrl Post','Medt Pre','Medt Post'},'Location','northeast');
xticks(1:40);
% y axis height set
kurtosis_all = [sub_kur_ctrl_premanip_neut,sub_kur_ctrl_posmanip_neut,sub_kur_medt_premanip_neut,sub_kur_medt_posmanip_neut];
YMax = max(max(kurtosis_all)) * 1.1;
axis([0 40 0 YMax]);

figure % Kurtosis, Heat
p1 = plot(x_ctrl,sub_kur_ctrl_premanip_heat,'ko',...
    'LineWidth',1,...
    'MarkerSize',6,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor','w');
hold on
p2 = plot(x_ctrl,sub_kur_ctrl_posmanip_heat,'ko',...
    'LineWidth',1,...
    'MarkerSize',6,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor','g');
p3 = plot(x_medt,sub_kur_medt_premanip_heat,'ko',...
    'LineWidth',1,...
    'MarkerSize',6,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','w');
p4 = plot(x_medt,sub_kur_medt_posmanip_heat,'ko',...
    'LineWidth',1,...
    'MarkerSize',6,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b');
hold off;
x_fullaxis = linspace(0,40,50);         %change number of subjects
str = sprintf('Kurtosis of Degree Distribution, Pre|Post Meditators|Controls, Heat at %.2f Abs. Corr. Threshold',thresh_abs);
title(str);
xlabel('Controls      -      (subjects)      -      Meditators');
ylabel('Kurtosis');
legend({'Ctrl Pre','Ctrl Post','Medt Pre','Medt Post'},'Location','northeast');
xticks(1:40);
% y axis height set
kurtosis_all = [sub_kur_ctrl_premanip_heat,sub_kur_ctrl_posmanip_heat,sub_kur_medt_premanip_heat,sub_kur_medt_posmanip_heat];
YMax = max(max(kurtosis_all)) * 1.1;
axis([0 40 0 YMax]);

figure % Variance, Neutral
p1 = plot(x_ctrl,sub_var_ctrl_premanip_neut,'ko',...
    'LineWidth',1,...
    'MarkerSize',6,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor','w');
hold on
p2 = plot(x_ctrl,sub_var_ctrl_posmanip_neut,'ko',...
    'LineWidth',1,...
    'MarkerSize',6,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor','g');
p3 = plot(x_medt,sub_var_medt_premanip_neut,'ko',...
    'LineWidth',1,...
    'MarkerSize',6,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','w');
p4 = plot(x_medt,sub_var_medt_posmanip_neut,'ko',...
    'LineWidth',1,...
    'MarkerSize',6,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b');
hold off;
x_fullaxis = linspace(0,40,50);         %change number of subjects
str = sprintf('Variance of Degree Distribution, Pre|Post Meditators|Controls, Neutral at %.2f Abs. Corr. Threshold',thresh_abs);
title(str);
xlabel('Controls      -      (subjects)      -      Meditators');
ylabel('variance');
legend({'Ctrl Pre','Ctrl Post','Medt Pre','Medt Post'},'Location','northeast');
xticks(1:40);
% y axis height set
variance_all = [sub_var_ctrl_premanip_neut,sub_var_ctrl_posmanip_neut,sub_var_medt_premanip_neut,sub_var_medt_posmanip_neut];
YMax = max(max(variance_all)) * 1.1;
axis([0 40 0 YMax]);

figure % Variance, Heat
p1 = plot(x_ctrl,sub_var_ctrl_premanip_heat,'ko',...
    'LineWidth',1,...
    'MarkerSize',6,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor','w');
hold on
p2 = plot(x_ctrl,sub_var_ctrl_posmanip_heat,'ko',...
    'LineWidth',1,...
    'MarkerSize',6,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor','g');
p3 = plot(x_medt,sub_var_medt_premanip_heat,'ko',...
    'LineWidth',1,...
    'MarkerSize',6,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','w');
p4 = plot(x_medt,sub_var_medt_posmanip_heat,'ko',...
    'LineWidth',1,...
    'MarkerSize',6,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b');
hold off;
x_fullaxis = linspace(0,40,50);         %change number of subjects
str = sprintf('Variance of Degree Distribution, Pre|Post Meditators|Controls, Heat at %.2f Abs. Corr. Threshold',thresh_abs);
title(str);
xlabel('Controls      -      (subjects)      -      Meditators');
ylabel('variance');
legend({'Ctrl Pre','Ctrl Post','Medt Pre','Medt Post'},'Location','northeast');
xticks(1:40);
% y axis height set
variance_all = [sub_var_ctrl_premanip_heat,sub_var_ctrl_posmanip_heat,sub_var_medt_premanip_heat,sub_var_medt_posmanip_heat];
YMax = max(max(variance_all)) * 1.1;
axis([0 40 0 YMax]);

%% group-level kur/var
ave_kur_ctrl_pos_heat = mean(sub_kur_ctrl_posmanip_heat);
ave_kur_ctrl_pos_neut = mean(sub_kur_ctrl_posmanip_neut);
ave_kur_ctrl_pre_heat = mean(sub_kur_ctrl_premanip_heat);
ave_kur_ctrl_pre_neut = mean(sub_kur_ctrl_premanip_neut);
ave_kur_medt_pos_heat = mean(sub_kur_medt_posmanip_heat);
ave_kur_medt_pos_neut = mean(sub_kur_medt_posmanip_neut);
ave_kur_medt_pre_heat = mean(sub_kur_medt_premanip_heat);
ave_kur_medt_pre_neut = mean(sub_kur_medt_premanip_neut);
ave_var_ctrl_pos_heat = mean(sub_var_ctrl_posmanip_heat);
ave_var_ctrl_pos_neut = mean(sub_var_ctrl_posmanip_neut);
ave_var_ctrl_pre_heat = mean(sub_var_ctrl_premanip_heat);
ave_var_ctrl_pre_neut = mean(sub_var_ctrl_premanip_neut);
ave_var_medt_pos_heat = mean(sub_var_medt_posmanip_heat);
ave_var_medt_pos_neut = mean(sub_var_medt_posmanip_neut);
ave_var_medt_pre_heat = mean(sub_var_medt_premanip_heat);
ave_var_medt_pre_neut = mean(sub_var_medt_premanip_neut);

% plot kurtosis
x = categorical({'Ctrl Neut','Ctrl Heat','Medt Neut','Medt Heat'});
y = [ave_kur_ctrl_pre_neut ave_kur_ctrl_pre_heat ave_kur_medt_pre_neut ave_kur_medt_pre_heat; ...
    ave_kur_ctrl_pos_neut ave_kur_ctrl_pos_heat ave_kur_medt_pos_neut ave_kur_medt_pos_heat];
figure; b=bar(x,y); legend('pre','post','Location','north'); 
str = sprintf('Kurtosis of Network Degree Distribution, Group-level at %.2f Abs. Corr. Threshold',thresh_abs);
sgtitle(str); ylabel('Kurtosis');
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom')
% plot variance
x = categorical({'Ctrl Neut','Ctrl Heat','Medt Neut','Medt Heat'});
y = [ave_var_ctrl_pre_neut ave_var_ctrl_pre_heat ave_var_medt_pre_neut ave_var_medt_pre_heat; ...
    ave_var_ctrl_pos_neut ave_var_ctrl_pos_heat ave_var_medt_pos_neut ave_var_medt_pos_heat];
figure; b=bar(x,y); legend('pre','post','Location','north'); 
str = sprintf('Variance of Network Degree Distribution, Group-level at %.2f Abs. Corr. Threshold',thresh_abs);
sgtitle(str); ylabel('Variance');
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom')

%% plot sub-level pre vs post entropy 
x = 1:1:20;
y = [sub_nS_medt_premanip_heat ; sub_nS_medt_posmanip_heat]; %change this for 1 of 4 groups/tasks

figure; tiledlayout('flow'); 
str = sprintf('Shannon Entropy, Subject-level at %.2f Abs. Corr. Threshold',thresh_abs);
sgtitle(str)

nexttile
b=bar(x,y); xlabel('subjects'); ylabel('entropy'); legend('pre','post','Location','south');
xtips1 = b(1).XEndPoints; ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
xtips2 = b(2).XEndPoints; ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom')

%% Regress chg(VAS-int/unp) on chg(nS) & level(nS) at subject level

% chg(nS) per sub (post-pre)
chg_nS_ctrl_neut = sub_nS_ctrl_posmanip_neut - sub_nS_ctrl_premanip_neut;
chg_nS_ctrl_heat = sub_nS_ctrl_posmanip_heat - sub_nS_ctrl_premanip_heat;
chg_nS_medt_neut = sub_nS_medt_posmanip_neut - sub_nS_medt_premanip_neut;
chg_nS_medt_heat = sub_nS_medt_posmanip_heat - sub_nS_medt_premanip_heat;

% demeaned nS
dm_nS_ctrl_pre_neut = sub_nS_ctrl_premanip_neut - mean(sub_nS_ctrl_premanip_neut);
dm_nS_ctrl_pos_neut = sub_nS_ctrl_posmanip_neut - mean(sub_nS_ctrl_posmanip_neut);
dm_nS_ctrl_pre_heat = sub_nS_ctrl_premanip_heat - mean(sub_nS_ctrl_premanip_heat);
dm_nS_ctrl_pos_neut = sub_nS_ctrl_posmanip_heat - mean(sub_nS_ctrl_posmanip_heat);
dm_nS_medt_pre_neut = sub_nS_medt_premanip_neut - mean(sub_nS_medt_premanip_neut);
dm_nS_medt_pos_neut = sub_nS_medt_posmanip_neut - mean(sub_nS_medt_posmanip_neut);
dm_nS_medt_pre_heat = sub_nS_medt_premanip_heat - mean(sub_nS_medt_premanip_heat);
dm_nS_medt_pos_neut = sub_nS_medt_posmanip_heat - mean(sub_nS_medt_posmanip_heat);

%% get VASint+unp
% % load
% [id int_Ah1 unp_Ah1 int_Ah2 unp_Ah2 int_Ah3 unp_Ah3 int_Ah4 unp_Ah4 ...
%     int_Bh1 unp_Bh1 int_Bh2 unp_Bh2 int_Bh3 unp_Bh3 int_Bh4 unp_Bh4] = readvars('mfc_2018-2.csv');
% % split into ctrl/medt  pre/pos
% int_ctrl_pre = (int_Ah1([1:20],:) + int_Ah2([1:20],:) + int_Ah3([1:20],:) + int_Ah4([1:20],:)) ./ 4;
% int_ctrl_pos = (int_Bh1([1:20],:) + int_Bh2([1:20],:) + int_Bh3([1:20],:) + int_Bh4([1:20],:)) ./ 4;
% unp_ctrl_pre = (unp_Ah1([1:20],:) + unp_Ah2([1:20],:) + unp_Ah3([1:20],:) + unp_Ah4([1:20],:)) ./ 4;
% unp_ctrl_pos = (unp_Bh1([1:20],:) + unp_Bh2([1:20],:) + unp_Bh3([1:20],:) + unp_Bh4([1:20],:)) ./ 4;
% int_medt_pre = (int_Ah1([21:40],:) + int_Ah2([21:40],:) + int_Ah3([21:40],:) + int_Ah4([21:40],:)) ./ 4;
% int_medt_pos = (int_Bh1([21:40],:) + int_Bh2([21:40],:) + int_Bh3([21:40],:) + int_Bh4([21:40],:)) ./ 4;
% unp_medt_pre = (unp_Ah1([21:40],:) + unp_Ah2([21:40],:) + unp_Ah3([21:40],:) + unp_Ah4([21:40],:)) ./ 4;
% unp_medt_pos = (unp_Bh1([21:40],:) + unp_Bh2([21:40],:) + unp_Bh3([21:40],:) + unp_Bh4([21:40],:)) ./ 4;

% load
[id	v6_int_h1 v6_unpl_h1 v6_int_h2 v6_unpl_h2 v6_int_h3 v6_unpl_h3 v6_int_h4 v6_unpl_h4] = readvars('pain_ratings.xlsx');

% average h1h2 h3h4 and split into ctrl/medt  pre/pos vars
int_ctrl_pre = (v6_int_h1([1:20],:) + v6_int_h2([1:20],:)) ./ 2;
int_ctrl_pos = (v6_int_h3([1:20],:) + v6_int_h4([1:20],:)) ./ 2;
int_medt_pre = (v6_int_h1([21:40],:) + v6_int_h2([21:40],:)) ./ 2;
int_medt_pos = (v6_int_h3([21:40],:) + v6_int_h4([21:40],:)) ./ 2;
unp_ctrl_pre = (v6_unpl_h1([1:20],:) + v6_unpl_h2([1:20],:)) ./ 2;
unp_ctrl_pos = (v6_unpl_h3([1:20],:) + v6_unpl_h4([1:20],:)) ./ 2;
unp_medt_pre = (v6_unpl_h1([21:40],:) + v6_unpl_h2([21:40],:)) ./ 2;
unp_medt_pos = (v6_unpl_h3([21:40],:) + v6_unpl_h4([21:40],:)) ./ 2;

% diff post-pre
chg_int_ctrl = int_ctrl_pos - int_ctrl_pre;
chg_unp_ctrl = unp_ctrl_pos - int_ctrl_pre;
chg_int_medt = int_medt_pos - int_medt_pre;
chg_unp_medt = unp_medt_pos - int_medt_pre;
% demean/normalize?  -later: redo everything per run (nS and vas)

% plot scatter chg vas vs chg nS
figure; tiledlayout(2,2); 
str = sprintf('\DeltaS vs \DeltaPain at %.2f Abs. Corr. Threshold',thresh_abs);
sgtitle(str);

nexttile
scatter(chg_int_medt, chg_nS_medt_heat); axis([-5 5 -1 1]); hold on
scatter(chg_int_medt, chg_nS_medt_neut); axis([-5 5 -1 1]); hold off
xlabel('\DeltaVAS'); ylabel('\DeltaS');
legend('\DeltaS heat','\DeltaS neut')
title('Meditators, \DeltaS vs \DeltaVAS Intensity')

nexttile
scatter(chg_unp_medt, chg_nS_medt_heat); axis([-5 5 -1 1]); hold on
scatter(chg_unp_medt, chg_nS_medt_neut); axis([-5 5 -1 1]); hold off
xlabel('\DeltaVAS'); ylabel('\DeltaS');
legend('\DeltaS heat','\DeltaS neut')
title('Meditators, \DeltaS vs \DeltaVAS Unpleasantness')

nexttile;
scatter(chg_int_ctrl, chg_nS_ctrl_heat); axis([-5 5 -1 1]); hold on
scatter(chg_int_ctrl, chg_nS_ctrl_neut); axis([-5 5 -1 1]); hold off
xlabel('\DeltaVAS'); ylabel('\DeltaS');
legend('\DeltaS heat','\DeltaS neut')
title('Controls, \DeltaS vs \DeltaVAS Intensity')

nexttile;
scatter(chg_unp_ctrl, chg_nS_ctrl_heat); axis([-5 5 -1 1]); hold on
scatter(chg_unp_ctrl, chg_nS_ctrl_neut); axis([-5 5 -1 1]); hold off
xlabel('\DeltaVAS'); ylabel('\DeltaS');
legend('\DeltaS heat','\DeltaS neut')
title('Controls, \DeltaS vs \DeltaVAS Unpleasantness')

%% Modularity for Group-averaged Weighted Correlation Matrices by Condition
% inputs: weighted corr matrices
% asymmetric treatment of negative weights - Rubinov & Sporns (2011) Neuroimage 56:2068-79.

[M, Q] = community_louvain(ctrl_premanip_neut_gp_matrix,[],[],'negative_asym'); mod_ctrl_pre_neut = Q;
[M, Q] = community_louvain(ctrl_premanip_heat_gp_matrix,[],[],'negative_asym'); mod_ctrl_pre_heat = Q;
[M, Q] = community_louvain(ctrl_posmanip_neut_gp_matrix,[],[],'negative_asym'); mod_ctrl_pos_neut = Q;
[M, Q] = community_louvain(ctrl_posmanip_heat_gp_matrix,[],[],'negative_asym'); mod_ctrl_pos_heat = Q;
[M, Q] = community_louvain(medt_premanip_neut_gp_matrix,[],[],'negative_asym'); mod_medt_pre_neut = Q;
[M, Q] = community_louvain(medt_premanip_heat_gp_matrix,[],[],'negative_asym'); mod_medt_pre_heat = Q;
[M, Q] = community_louvain(medt_posmanip_neut_gp_matrix,[],[],'negative_asym'); mod_medt_pos_neut = Q;
[M, Q] = community_louvain(medt_posmanip_heat_gp_matrix,[],[],'negative_asym'); mod_medt_pos_heat = Q;
clear M Q

% Change (Post-Pre)
mod_chg_ctrl_neut = mod_ctrl_pos_neut - mod_ctrl_pre_neut;
mod_chg_ctrl_heat = mod_ctrl_pos_heat - mod_ctrl_pre_heat;
mod_chg_medt_neut = mod_medt_pos_neut - mod_medt_pre_neut;
mod_chg_medt_heat = mod_medt_pos_heat - mod_medt_pre_heat;

% Plot
figure
tiledlayout(1,2)
sgtitle('Modularity of Functional Network')

nexttile
x = categorical({'Controls Neutral','Controls Heat','Meditators Neutral','Meditators Heat'});
y = [mod_ctrl_pre_neut mod_ctrl_pre_heat mod_medt_pre_neut mod_medt_pre_heat; ...
    mod_ctrl_pos_neut mod_ctrl_pos_heat mod_medt_pos_neut mod_medt_pos_heat];
b=bar(x,y); legend('pre','post','Location','north'); 
%b.FaceColor = 'flat'; b.CData(3) = [0.6350 0.0780 0.1840]; b.CData(4) = [0.6350 0.0780 0.1840];
ylabel('Modularity'); title('Modularity, Pre vs Post');
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom')

nexttile
x = categorical({'Controls Neutral','Controls Heat','Meditators Neutral','Meditators Heat'});
y = [mod_chg_ctrl_neut mod_chg_ctrl_heat mod_chg_medt_neut mod_chg_medt_heat];
b=bar(x,y); ylabel('\Delta Modularity'); title('Change in Modularity, Post - Pre');
b.FaceColor = 'flat'; b.CData(3,:) = [0.6350 0.0780 0.1840]; b.CData(4,:) = [0.6350 0.0780 0.1840];
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')

%% functions (in separate .m file)

% % takes unthresholded corr matrix and returns thresholded binarized matrix,
% % its degree distribution, and its entropy
% function [w_bin_abs, ddist, nS] = preproc_corr_matrix(unthresh_matrix, thresh_abs)
%    abs_thresh_matrix = threshold_absolute(unthresh_matrix, thresh_abs);
%    w_bin_abs = weight_conversion(abs_thresh_matrix, 'binarize');
%    ddist = degrees_und(w_bin_abs);
%    
%    entropy_bins = 120;
%    f = figure('visible','off'); %initialize fig invisible
%    hist1 = histogram(ddist,entropy_bins);
%    hist1.Normalization = 'probability';
%    binvalues = hist1.Values;
%    nS = sum(-(binvalues(binvalues>0).*(log(binvalues(binvalues>0)))));
%    clear f; clear hist1;
% end