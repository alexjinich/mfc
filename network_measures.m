clc;
close all;
clearvars;
addpath '/Volumes/Seagate_Desktop_Drive/mfc/code'

%%%%%%%%%%%%%%%% ENTROPY PER sub V2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This version takes matrices per individual sub, averages over 2 subs and then
% computes ddist, S per sub to make sub & group comparisons.
% Threshold_robustness_2 and adjmatrix.m average 2 sub matrices 1st, then
% compute ddist, S on the averaged sub-condition matrix.
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
thresh = 0.40;
%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET 1 IF LOG(S), 0 NO LOG
%log_S = 0;
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
%clear ctrl* medt*

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
clear i unthresh_matrix w_bin_abs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modularity (per subject)
% inputs: weighted corr matrices
% asymmetric treatment of negative weights - Rubinov & Sporns (2011) Neuroimage 56:2068-79.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mod_ctrl_pos_heat = zeros(1,20);
mod_ctrl_pos_neut = zeros(1,20);
mod_ctrl_pre_heat = zeros(1,20);
mod_ctrl_pre_neut = zeros(1,20);
mod_medt_pos_heat = zeros(1,20);
mod_medt_pos_neut = zeros(1,20);
mod_medt_pre_heat = zeros(1,20);
mod_medt_pre_neut = zeros(1,20);

for i = 1:20
    [M, Q] = community_louvain(sub_ctrl_posmanip_heat(:,:,i),[],[],'negative_asym'); mod_ctrl_pos_heat(i) = Q;
    [M, Q] = community_louvain(sub_ctrl_posmanip_neut(:,:,i),[],[],'negative_asym'); mod_ctrl_pos_neut(i) = Q;
    [M, Q] = community_louvain(sub_ctrl_premanip_heat(:,:,i),[],[],'negative_asym'); mod_ctrl_pre_heat(i) = Q;
    [M, Q] = community_louvain(sub_ctrl_premanip_neut(:,:,i),[],[],'negative_asym'); mod_ctrl_pre_neut(i) = Q;
    [M, Q] = community_louvain(sub_medt_posmanip_heat(:,:,i),[],[],'negative_asym'); mod_medt_pos_heat(i) = Q;
    [M, Q] = community_louvain(sub_medt_posmanip_neut(:,:,i),[],[],'negative_asym'); mod_medt_pos_neut(i) = Q;
    [M, Q] = community_louvain(sub_medt_premanip_heat(:,:,i),[],[],'negative_asym'); mod_medt_pre_heat(i) = Q;
    [M, Q] = community_louvain(sub_medt_premanip_neut(:,:,i),[],[],'negative_asym'); mod_medt_pre_neut(i) = Q;
end
clear M Q

% change in modularity post-pre
mod_chg_ctrl_heat = mod_ctrl_pos_heat - mod_ctrl_pre_heat;
mod_chg_ctrl_neut = mod_ctrl_pos_neut - mod_ctrl_pre_neut;
mod_chg_medt_heat = mod_medt_pos_heat - mod_medt_pre_heat;
mod_chg_medt_neut = mod_medt_pos_neut - mod_medt_pre_neut;
 
% SE of change in modularity
SE_dmod_ctrl_heat = std(mod_chg_ctrl_heat) / sqrt(20);
SE_dmod_ctrl_neut = std(mod_chg_ctrl_neut) / sqrt(20);
SE_dmod_medt_heat = std(mod_chg_medt_heat) / sqrt(20);
SE_dmod_medt_neut = std(mod_chg_medt_neut) / sqrt(20);

% Plot change in modularity per group
x = categorical({'Controls Neutral','Controls Heat','Meditators Neutral','Meditators Heat'});
y = [mean(mod_chg_ctrl_neut) mean(mod_chg_ctrl_heat) mean(mod_chg_medt_neut) mean(mod_chg_medt_heat)];
err = [SE_dmod_ctrl_neut SE_dmod_ctrl_heat SE_dmod_medt_neut SE_dmod_medt_heat];
figure; 
b=bar(x,y); 
b.FaceColor = 'flat';
b.CData(3,:) = [0.6350 0.0780 0.1840]; b.CData(4,:) = [0.6350 0.0780 0.1840];
hold on
er = errorbar(x,y,err); er.Color = 'k'; er.LineStyle = 'none';
hold off
str = sprintf('Change in Modularity (Post - Pre), Group-level at %.2f Abs. Corr. Threshold',thresh);
sgtitle(str); ylabel('\Delta Modularity Post - Pre');
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
clear x y b xtips1 ytips1

% All sub, chg in modularity
x1=1:20; x2=21:40;
% Heat
figure
str = sprintf('Change in Modularity Per Subject at %.2f Corr. Threshold',thresh);
sgtitle(str);
tiledlayout(1,2)
nexttile
stem(x1,mod_chg_ctrl_heat, 'filled')
hold on
stem(x2,mod_chg_medt_heat, 'filled')
hold off
title('Heat Runs'); ylabel('\Delta Modularity Post - Pre'); xlabel('subject'); legend('Control','Meditator');
% Neut
nexttile
stem(x1,mod_chg_ctrl_neut, 'filled')
hold on
stem(x2,mod_chg_medt_neut, 'filled')
hold off
title('Neutral Runs'); ylabel('\Delta Modularity Post - Pre'); xlabel('subject'); legend('Control','Meditator');

clear x* y 

% RANOVA modularity
mod_ctrl_pre_neut = permute(mod_ctrl_pre_neut,[2,1]);
mod_ctrl_pos_neut = permute(mod_ctrl_pos_neut,[2,1]);
mod_ctrl_pre_heat = permute(mod_ctrl_pre_heat,[2,1]);
mod_ctrl_pos_heat = permute(mod_ctrl_pos_heat,[2,1]);
mod_medt_pre_neut = permute(mod_medt_pre_neut,[2,1]);
mod_medt_pos_neut = permute(mod_medt_pos_neut,[2,1]);
mod_medt_pre_heat = permute(mod_medt_pre_heat,[2,1]);
mod_medt_pos_heat = permute(mod_medt_pos_heat,[2,1]);

% NEUTRAL
% Construct datatable with the 4 clusters of data and name the variables.
datatable = table(mod_ctrl_pre_neut,mod_ctrl_pos_neut,mod_medt_pre_neut,mod_medt_pos_neut);
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
disp('Repeated Measures ANOVA: Modularity, Neutral, 2 by 2:');
ranovatable = ranova(rm,'WithinModel','Treatment*Group')
%disp('Coefficients:'); 
%rm.Coefficients
%disp('Covariance:'); 
%rm.Covariance
% Clear temp vars.
clear datatable WithinStructure rm ranovatable

% HEAT
% Construct datatable with the 4 clusters of data and name the variables.
datatable = table(mod_ctrl_pre_heat,mod_ctrl_pos_heat,mod_medt_pre_heat,mod_medt_pos_heat);
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
disp('Repeated Measures ANOVA: Modularity, Heat, 2 by 2:');
ranovatable = ranova(rm,'WithinModel','Treatment*Group')
%disp('Coefficients:'); 
%rm.Coefficients
%disp('Covariance:'); 
%rm.Covariance
% Clear temp vars.
clear datatable WithinStructure rm ranovatable

% T-tests per subject for sig diff in S pre vs post

disp('T-tests for Difference in Modularity Pre vs Post');

[h,p] = ttest(mod_ctrl_pre_neut, mod_ctrl_pos_neut); % Controls, Neutral
if h == 0
    h = '*not* rejected';
else
    h = '***rejected***';
end
d = sprintf('\nControls,   Neutral - Null H (no mean difference) is %s with p = %.3f', h, p); disp(d); 
 
[h,p] = ttest(mod_ctrl_pre_heat, mod_ctrl_pos_heat); % Controls, Heat
if h == 0
    h = '*not* rejected';
else
    h = '***rejected***';
end
d = sprintf('Controls,   Heat    - Null H (no mean difference) is %s with p = %.3f', h, p); disp(d); 
 
[h,p] = ttest(mod_medt_pre_neut, mod_medt_pos_neut); % Meditators, Neutral
if h == 0
    h = '*not* rejected';
else
    h = '***rejected***';
end
d = sprintf('Meditators, Neutral - Null H (no mean difference) is %s with p = %.3f', h, p); disp(d); 
 
[h,p] = ttest(mod_medt_pre_heat, mod_medt_pos_heat); % Meditators, Heat
if h == 0
    h = '*not* rejected';
else
    h = '***rejected***';
end
d = sprintf('Meditators, Heat    - Null H (no mean difference) is %s with p = %.3f', h, p); disp(d); 
 
clear h p d


% % Change in Modularity VS Change in VAS per sub
% 
% % Heat runs (nS of h3h4 - nS of h1h2)
% x_tick = [-.5:.1:.6]; y_tick = [-6:1:4];
% figure; tiledlayout(2,2); 
% str = sprintf('Change in Modularity (During Heat Runs) vs Change in VAS Pain at %.2f Correlation Threshold',thresh);
% sgtitle(str)
% nexttile
% scatter(vas_chg_int_ctrl,mod_chg_ctrl_heat,'filled'); title('Controls, VAS Int');   
% xlabel('\Delta pain int'); ylabel('\Delta Modularity'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% nexttile
% scatter(vas_chg_int_medt,mod_chg_medt_heat,'filled'); title('Meditators, VAS Int'); 
% xlabel('\Delta pain int'); ylabel('\Delta Modularity'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% nexttile
% scatter(vas_chg_unp_ctrl,mod_chg_ctrl_heat,'filled'); title('Controls, VAS Unp');   
% xlabel('\Delta pain unp'); ylabel('\Delta Modularity'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% nexttile
% scatter(vas_chg_unp_medt,mod_chg_medt_heat,'filled'); title('Meditators, VAS Unp'); 
% xlabel('\Delta pain unp'); ylabel('\Delta Modularity'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% clear x_tick y_tick
% 
% % Neut runs (nS of h3h4 - nS of h1h2)
% x_tick = [-.5:.1:.6]; y_tick = [-6:1:4];
% figure; tiledlayout(2,2); 
% str = sprintf('Change in Modularity (During neut Runs) vs Change in VAS Pain at %.2f Correlation Threshold',thresh);
% sgtitle(str)
% nexttile
% scatter(vas_chg_int_ctrl,mod_chg_ctrl_neut,'filled'); title('Controls, VAS Int');   
% xlabel('\Delta pain int'); ylabel('\Delta Modularity'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% nexttile
% scatter(vas_chg_int_medt,mod_chg_medt_neut,'filled'); title('Meditators, VAS Int'); 
% xlabel('\Delta pain int'); ylabel('\Delta Modularity'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% nexttile
% scatter(vas_chg_unp_ctrl,mod_chg_ctrl_neut,'filled'); title('Controls, VAS Unp');   
% xlabel('\Delta pain unp'); ylabel('\Delta Modularity'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% nexttile
% scatter(vas_chg_unp_medt,mod_chg_medt_neut,'filled'); title('Meditators, VAS Unp'); 
% xlabel('\Delta pain unp'); ylabel('\Delta Modularity'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% clear x_tick y_tick

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Global Efficiency (per subject)
% inputs: weighted corr matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Eglob_ctrl_pos_heat = zeros(1,20);
Eglob_ctrl_pos_neut = zeros(1,20);
Eglob_ctrl_pre_heat = zeros(1,20);
Eglob_ctrl_pre_neut = zeros(1,20);
Eglob_medt_pos_heat = zeros(1,20);
Eglob_medt_pos_neut = zeros(1,20);
Eglob_medt_pre_heat = zeros(1,20);
Eglob_medt_pre_neut = zeros(1,20);

for i = 1:20
    Eglob = efficiency_bin(sub_bin_abs_ctrl_pos_heat(:,:,i)); Eglob_ctrl_pos_heat(i) = Eglob;
    Eglob = efficiency_bin(sub_bin_abs_ctrl_pos_neut(:,:,i)); Eglob_ctrl_pos_neut(i) = Eglob;
    Eglob = efficiency_bin(sub_bin_abs_ctrl_pre_heat(:,:,i)); Eglob_ctrl_pre_heat(i) = Eglob;
    Eglob = efficiency_bin(sub_bin_abs_ctrl_pre_neut(:,:,i)); Eglob_ctrl_pre_neut(i) = Eglob;
    Eglob = efficiency_bin(sub_bin_abs_medt_pos_heat(:,:,i)); Eglob_medt_pos_heat(i) = Eglob;
    Eglob = efficiency_bin(sub_bin_abs_medt_pos_neut(:,:,i)); Eglob_medt_pos_neut(i) = Eglob;
    Eglob = efficiency_bin(sub_bin_abs_medt_pre_heat(:,:,i)); Eglob_medt_pre_heat(i) = Eglob;
    Eglob = efficiency_bin(sub_bin_abs_medt_pre_neut(:,:,i)); Eglob_medt_pre_neut(i) = Eglob;
end
clear Eglob

% change in modularity post-pre
Eglob_chg_ctrl_heat = Eglob_ctrl_pos_heat - Eglob_ctrl_pre_heat;
Eglob_chg_ctrl_neut = Eglob_ctrl_pos_neut - Eglob_ctrl_pre_neut;
Eglob_chg_medt_heat = Eglob_medt_pos_heat - Eglob_medt_pre_heat;
Eglob_chg_medt_neut = Eglob_medt_pos_neut - Eglob_medt_pre_neut;

% SE of change in global efficiency
SE_dEglob_ctrl_heat = std(Eglob_chg_ctrl_heat) / sqrt(20);
SE_dEglob_ctrl_neut = std(Eglob_chg_ctrl_neut) / sqrt(20);
SE_dEglob_medt_heat = std(Eglob_chg_medt_heat) / sqrt(20);
SE_dEglob_medt_neut = std(Eglob_chg_medt_neut) / sqrt(20);

% Plot change in global efficiency per group
x = categorical({'Controls Neutral','Controls Heat','Meditators Neutral','Meditators Heat'});
y = [mean(Eglob_chg_ctrl_neut) mean(Eglob_chg_ctrl_heat) mean(Eglob_chg_medt_neut) mean(Eglob_chg_medt_heat)];
err = [SE_dEglob_ctrl_neut SE_dEglob_ctrl_heat SE_dEglob_medt_neut SE_dEglob_medt_heat];
figure; 
b=bar(x,y); 
b.FaceColor = 'flat';
b.CData(3,:) = [0.6350 0.0780 0.1840]; b.CData(4,:) = [0.6350 0.0780 0.1840];
hold on
er = errorbar(x,y,err); er.Color = 'k'; er.LineStyle = 'none';
hold off
str = sprintf('Change in Global Efficiency (Post - Pre), Group-level at %.2f Abs. Corr. Threshold',thresh);
sgtitle(str); ylabel('\Delta Global Efficiency Post - Pre');
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
clear x y b xtips1 ytips1

% All sub, chg in global efficiency
x1=1:20; x2=21:40;
% Heat
figure
str = sprintf('Change in Global Efficiency Per Subject at %.2f Corr. Threshold',thresh);
sgtitle(str);
tiledlayout(1,2)
nexttile
stem(x1,Eglob_chg_ctrl_heat, 'filled')
hold on
stem(x2,Eglob_chg_medt_heat, 'filled')
hold off
title('Heat Runs'); ylabel('\Delta global efficiency Post - Pre'); xlabel('subject'); legend('Control','Meditator');
% Neut
nexttile
stem(x1,Eglob_chg_ctrl_neut, 'filled')
hold on
stem(x2,Eglob_chg_medt_neut, 'filled')
hold off
title('Neutral Runs'); ylabel('\Delta global efficiency Post - Pre'); xlabel('subject'); legend('Control','Meditator');

clear x* y 

% RANOVA modularity
Eglob_ctrl_pre_neut = permute(Eglob_ctrl_pre_neut,[2,1]);
Eglob_ctrl_pos_neut = permute(Eglob_ctrl_pos_neut,[2,1]);
Eglob_ctrl_pre_heat = permute(Eglob_ctrl_pre_heat,[2,1]);
Eglob_ctrl_pos_heat = permute(Eglob_ctrl_pos_heat,[2,1]);
Eglob_medt_pre_neut = permute(Eglob_medt_pre_neut,[2,1]);
Eglob_medt_pos_neut = permute(Eglob_medt_pos_neut,[2,1]);
Eglob_medt_pre_heat = permute(Eglob_medt_pre_heat,[2,1]);
Eglob_medt_pos_heat = permute(Eglob_medt_pos_heat,[2,1]);

% NEUTRAL
% Construct datatable with the 4 clusters of data and name the variables.
datatable = table(Eglob_ctrl_pre_neut,Eglob_ctrl_pos_neut,Eglob_medt_pre_neut,Eglob_medt_pos_neut);
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
disp('Repeated Measures ANOVA: Global Efficiency, Neutral, 2 by 2:');
ranovatable = ranova(rm,'WithinModel','Treatment*Group')
%disp('Coefficients:'); 
%rm.Coefficients
%disp('Covariance:'); 
%rm.Covariance
% Clear temp vars.
clear datatable WithinStructure rm ranovatable

% HEAT
% Construct datatable with the 4 clusters of data and name the variables.
datatable = table(Eglob_ctrl_pre_heat,Eglob_ctrl_pos_heat,Eglob_medt_pre_heat,Eglob_medt_pos_heat);
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
disp('Repeated Measures ANOVA: Global Efficiency, Heat, 2 by 2:');
ranovatable = ranova(rm,'WithinModel','Treatment*Group')
%disp('Coefficients:'); 
%rm.Coefficients
%disp('Covariance:'); 
%rm.Covariance
% Clear temp vars.
clear datatable WithinStructure rm ranovatable

% T-tests per subject for sig diff in S pre vs post

disp('T-tests for Difference in Global Efficiency Pre vs Post');
[h,p] = ttest(Eglob_ctrl_pre_neut, Eglob_ctrl_pos_neut); % Controls, Neutral
if h == 0
    h = '*not* rejected';
else
    h = '***rejected***';
end
d = sprintf('\nControls,   Neutral - Null H (no mean difference) is %s with p = %.3f', h, p); disp(d); 
 
[h,p] = ttest(Eglob_ctrl_pre_heat, Eglob_ctrl_pos_heat); % Controls, Heat
if h == 0
    h = '*not* rejected';
else
    h = '***rejected***';
end
d = sprintf('Controls,   Heat    - Null H (no mean difference) is %s with p = %.3f', h, p); disp(d); 
 
[h,p] = ttest(Eglob_medt_pre_neut, Eglob_medt_pos_neut); % Meditators, Neutral
if h == 0
    h = '*not* rejected';
else
    h = '***rejected***';
end
d = sprintf('Meditators, Neutral - Null H (no mean difference) is %s with p = %.3f', h, p); disp(d); 
 
[h,p] = ttest(Eglob_medt_pre_heat, Eglob_medt_pos_heat); % Meditators, Heat
if h == 0
    h = '*not* rejected';
else
    h = '***rejected***';
end
d = sprintf('Meditators, Heat    - Null H (no mean difference) is %s with p = %.3f \n\n', h, p); disp(d); 
 
clear h p d

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clustering Coefficient (per ROI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clustering coefficient (weighted undirected corr matrices)

cc_ctrl_pos_heat = zeros(nrois,20);
cc_ctrl_pos_neut = zeros(nrois,20);
cc_ctrl_pre_heat = zeros(nrois,20);
cc_ctrl_pre_neut = zeros(nrois,20);
cc_medt_pos_heat = zeros(nrois,20);
cc_medt_pos_neut = zeros(nrois,20);
cc_medt_pre_heat = zeros(nrois,20);
cc_medt_pre_neut = zeros(nrois,20);

for i = 1:20
    C = clustering_coef_bu(sub_bin_abs_ctrl_pos_heat(:,:,i));
    cc_ctrl_pos_heat(:,i) = C;
    C = clustering_coef_bu(sub_bin_abs_ctrl_pos_neut(:,:,i));
    cc_ctrl_pos_neut(:,i) = C;
    C = clustering_coef_bu(sub_bin_abs_ctrl_pre_heat(:,:,i));
    cc_ctrl_pre_heat(:,i) = C;
    C = clustering_coef_bu(sub_bin_abs_ctrl_pre_neut(:,:,i));
    cc_ctrl_pre_neut(:,i) = C;
    C = clustering_coef_bu(sub_bin_abs_medt_pos_heat(:,:,i));
    cc_medt_pos_heat(:,i) = C;
    C = clustering_coef_bu(sub_bin_abs_medt_pos_neut(:,:,i));
    cc_medt_pos_neut(:,i) = C;
    C = clustering_coef_bu(sub_bin_abs_medt_pre_heat(:,:,i));
    cc_medt_pre_heat(:,i) = C;
    C = clustering_coef_bu(sub_bin_abs_medt_pre_neut(:,:,i));
    cc_medt_pre_neut(:,i) = C;
end 
clear C

%change in cc
cc_chg_ctrl_heat = cc_ctrl_pos_heat - cc_ctrl_pre_heat;
cc_chg_ctrl_neut = cc_ctrl_pos_neut - cc_ctrl_pre_neut;
cc_chg_medt_heat = cc_medt_pos_heat - cc_medt_pre_heat;
cc_chg_medt_neut = cc_medt_pos_neut - cc_medt_pre_neut;

%mean change per roi
cc_mchg_ctrl_heat = mean(cc_chg_ctrl_heat,2);
cc_mchg_ctrl_neut = mean(cc_chg_ctrl_neut,2);
cc_mchg_medt_heat = mean(cc_chg_medt_heat,2);
cc_mchg_medt_neut = mean(cc_chg_medt_neut,2);

% RMAnova all rois
d = sprintf('Clustering Coefficient: RMAnova for each ROI:\n');
disp(d);
% group var of size = #sub
A = [0]; X = 20; C = repmat(A,X,1); % controls
A = [1]; X = 20; M = repmat(A,X,1); % meditators

% loop all rois: rmanova on each for all sub, store p-vals in array
for roi = 1:nrois
   % extract roi
   cc_ctrl_pre_heat_roi = cc_ctrl_pre_heat(roi,:);
   cc_ctrl_pre_neut_roi = cc_ctrl_pre_neut(roi,:);
   cc_ctrl_pos_heat_roi = cc_ctrl_pos_heat(roi,:);
   cc_ctrl_pos_neut_roi = cc_ctrl_pos_neut(roi,:);
   cc_medt_pre_heat_roi = cc_medt_pre_heat(roi,:);
   cc_medt_pre_neut_roi = cc_medt_pre_neut(roi,:);
   cc_medt_pos_heat_roi = cc_medt_pos_heat(roi,:);
   cc_medt_pos_neut_roi = cc_medt_pos_neut(roi,:);
   % permute
   cc_ctrl_pre_heat_roi = cc_ctrl_pre_heat_roi';
   cc_ctrl_pre_neut_roi = cc_ctrl_pre_neut_roi';
   cc_ctrl_pos_heat_roi = cc_ctrl_pos_heat_roi';
   cc_ctrl_pos_neut_roi = cc_ctrl_pos_neut_roi';
   cc_medt_pre_heat_roi = cc_medt_pre_heat_roi';
   cc_medt_pre_neut_roi = cc_medt_pre_neut_roi';
   cc_medt_pos_heat_roi = cc_medt_pos_heat_roi';
   cc_medt_pos_neut_roi = cc_medt_pos_neut_roi';
   % ranova

   % NEUTRAL-only
   
   % define data 
    dataTable = array2table([cc_ctrl_pre_neut_roi,cc_ctrl_pos_neut_roi,C;...
        cc_medt_pre_neut_roi,cc_medt_pos_neut_roi,M]);
    % convert to table format
    dataTable.Properties.VariableNames = {'Pre','Post','Group'};
    % convert group to categorical
    dataTable.Group = categorical(dataTable.Group);
    % define levels of within subject factor
    wsVariable = table([0 1]','VariableNames',{'Treatment'});
    %run RANOVA (within subjects)
    rm = fitrm(dataTable,'Pre,Post~Group','WithinDesign',wsVariable);    
    ranovatbl = ranova(rm);
    %run ANOVA (between subjects)
    anovatbl = anova(rm);
    % store p values
    clusteringcoef_pvals(:,roi) = [ranovatbl.pValue; anovatbl.pValue];
    % clear temp vars
    clear dataTable;
    clear rm;
    clear anovatbl;
    clear ranovatbl;
end

% remove cols 3,4,6 to leave only 3 p-values (treatment, interaction, group)
clusteringcoef_pvals(3,:) = [];   %remove 3rd row
clusteringcoef_pvals(3,:) = [];   %remove 4th row
clusteringcoef_pvals(4,:) = [];   %remove 6th row
clear A X C M roi;

% Flag ROIs with p<0.05 and build subset array and index vars with them

% 1) which rois are significant? build 2 arrays with them:
% clusteringcoef_sigrois = zeroes for non-sig and the p-values for sig rois
% clusteringcoef_sigrois_index = index of sig rois

% 1a) For p significant for treatment, group, or t*g (any)
clusteringcoef_sigrois_any = zeros(3, nrois);
clusteringcoef_sigrois_any_index = [];
for i = 1:nrois
    if min(clusteringcoef_pvals(:,i)) <= 0.05
        clusteringcoef_sigrois_any(:,i) = clusteringcoef_pvals(:,i);
        clusteringcoef_sigrois_any_index = [clusteringcoef_sigrois_any_index, i];
    else
        clusteringcoef_sigrois_any(:,i) = zeros(3,1);
    end
end
% display how many significant rois were found
y = size(clusteringcoef_sigrois_any_index);
y = y(2);
x = ['Clustering Coefficients RM-ANOVA: ', num2str(y), ' significant ROIs out of ', ...
    num2str(nrois), ' for any effects'];
disp(x);
clear x;
clear y;
clear i;

% 1b) For p significant for treatment-only
clusteringcoef_sigrois_treatment = zeros(1, nrois);
clusteringcoef_sigrois_treatment_index = [];
for i = 1:nrois
    if clusteringcoef_pvals(1,i) <= 0.05
        clusteringcoef_sigrois_treatment(1,i) = clusteringcoef_pvals(1,i);
        clusteringcoef_sigrois_treatment_index = [clusteringcoef_sigrois_treatment_index, i];
    else
        clusteringcoef_sigrois_treatment(1,i) = zeros(1,1);
    end
end
% display how many significant rois were found
y = size(clusteringcoef_sigrois_treatment_index);
y = y(2);
x = ['Clustering Coefficients RM-ANOVA: ', num2str(y), ' significant ROIs out of ', ...
    num2str(nrois), ' for treatment effects'];
disp(x);
clear x;
clear y;
clear i;

% 1c) For p significant for treatment-by-group-interaction only
clusteringcoef_sigrois_tbyg = zeros(1, nrois);
clusteringcoef_sigrois_tbyg_index = [];
for i = 1:nrois
    if clusteringcoef_pvals(2,i) <= 0.05
        clusteringcoef_sigrois_tbyg(1,i) = clusteringcoef_pvals(2,i);
        clusteringcoef_sigrois_tbyg_index = [clusteringcoef_sigrois_tbyg_index, i];
    else
        clusteringcoef_sigrois_tbyg(1,i) = zeros(1,1);
    end
end
% display how many significant rois were found
y = size(clusteringcoef_sigrois_tbyg_index);
y = y(2);
x = ['Clustering Coefficients RM-ANOVA: ', num2str(y), ' significant ROIs out of ', ...
    num2str(nrois), ' for treatment-by-group interaction effects'];
disp(x);
clear x;
clear y;
clear i;

% 1d) For p significant for group-only
clusteringcoef_sigrois_group = zeros(1, nrois);
clusteringcoef_sigrois_group_index = [];
for i = 1:nrois
    if clusteringcoef_pvals(3,i) <= 0.05
        clusteringcoef_sigrois_group(1,i) = clusteringcoef_pvals(3,i);
        clusteringcoef_sigrois_group_index = [clusteringcoef_sigrois_group_index, i];
    else
        clusteringcoef_sigrois_group(1,i) = zeros(1,1);
    end
end
% display how many significant rois were found
y = size(clusteringcoef_sigrois_group_index);
y = y(2);
x = ['Clustering Coefficients RM-ANOVA: ', num2str(y), ' significant ROIs out of ', ...
    num2str(nrois), ' for between group effects'];
disp(x);

clear x y i

% Clustering Coefficient sig-ROI vectors, p<0.05:
%
% (a) all effects:
%     clusteringcoef_sigrois_any
%     clusteringcoef_sigrois_any_index
%
% (b) treatment effects:
%     clusteringcoef_sigrois_treatment
%     clusteringcoef_sigrois_treatment_index
%
% (b) treatment-group interaction effects:
%     clusteringcoef_sigrois_tbyg
%     clusteringcoef_sigrois_tbyg_index
%
% (b) group effects:
%     clusteringcoef_sigrois_group
%     clusteringcoef_sigrois_group_index

% 2) Get mean change for sig rois post-pre for both groups

% take index of sig p-values and change clusteringcoef_meandiff for both
% groups such that all non-sig columns are null.

% initialize vectors for 3 sig conditions/effect-types
% treatment-only
clusteringcoef_meandiff_controls_sigtreatment = cc_mchg_ctrl_neut;
clusteringcoef_meandiff_meditators_sigtreatment = cc_mchg_medt_neut;
% interaction-only
clusteringcoef_meandiff_controls_sigtbyg = cc_mchg_ctrl_neut;
clusteringcoef_meandiff_meditators_sigtbyg = cc_mchg_medt_neut;
% group-only
clusteringcoef_meandiff_controls_siggroup = cc_mchg_ctrl_neut;
clusteringcoef_meandiff_meditators_siggroup = cc_mchg_medt_neut;

% check if roi has significance; if so copy mean difference #, else null
% (a) treatment-only
for i = 1:nrois
    if ismember(i,clusteringcoef_sigrois_treatment_index) == 0
        clusteringcoef_meandiff_controls_sigtreatment(i) = nan;
        clusteringcoef_meandiff_meditators_sigtreatment(i) = nan;
    else
        clusteringcoef_meandiff_controls_sigtreatment(i) = cc_mchg_ctrl_neut(i);
        clusteringcoef_meandiff_meditators_sigtreatment(i) = cc_mchg_medt_neut(i);
    end
end
clear i

% (b) interaction-only
for i = 1:nrois
    if ismember(i,clusteringcoef_sigrois_tbyg_index) == 0
        clusteringcoef_meandiff_controls_sigtbyg(i) = nan;
        clusteringcoef_meandiff_meditators_sigtbyg(i) = nan;
    else
        clusteringcoef_meandiff_controls_sigtbyg(i) = cc_mchg_ctrl_neut(i);
        clusteringcoef_meandiff_meditators_sigtbyg(i) = cc_mchg_medt_neut(i);
    end
end
clear i

% (c) group-only
for i = 1:nrois
    if ismember(i,clusteringcoef_sigrois_group_index) == 0
        clusteringcoef_meandiff_controls_siggroup(i) = nan;
        clusteringcoef_meandiff_meditators_siggroup(i) = nan;
    else
        clusteringcoef_meandiff_controls_siggroup(i) = cc_mchg_ctrl_neut(i);
        clusteringcoef_meandiff_meditators_siggroup(i) = cc_mchg_medt_neut(i);
    end
end
clear i

% 3) Create stem plot for just significant rois (post-pre mean chg for both groups)

%ROInames2 = ROInames(:);

[~,labels] = readvars('atlas_labels_harvardoxford_cort-maxprob-thr25-2mm.csv');
labels = labels(2:49,:);

%get # of ROIs and set x 
%sizeROI = size(ROInames);
%sizeROI = sizeROI(2);
% x = linspace(1,sizeROI,sizeROI);

% plot
figure;
%treatment-only
stem(clusteringcoef_meandiff_controls_sigtreatment, 'b<');
xticks([1:nrois]);
xticklabels(labels);
xtickangle(45);
xlim([1 nrois]);
title({['{\bf\fontsize{14} Clustering Coefficient for ROIs with p<0.05 for Any Sig. Effect, Mean Change Post-Pre (Neutral)}' ,...
    ', corr ', num2str(thresh)]; 'clustering coefficient: fraction of triangles around a node and is equivalent to the fraction of node neighbors that are neighbors of each other'},'FontWeight','Normal')
ylabel('Change in Clustering Coefficient');
hold on;
stem(clusteringcoef_meandiff_meditators_sigtreatment, 'r<', 'filled');
%interaction-only
stem(clusteringcoef_meandiff_controls_sigtbyg, 'b*');
stem(clusteringcoef_meandiff_meditators_sigtbyg, 'r*');
%group-only
stem(clusteringcoef_meandiff_controls_siggroup, 'b>');
stem(clusteringcoef_meandiff_meditators_siggroup, 'r>', 'filled');
legend({'controls - sig. treatment effect','meditators - sig. treatment effect',...
    'controls - sig. T*G effect','meditators - sig. T*G effect',...
    'controls - sig. group effect','meditators - sig. group effect'},'Location','northeast');
hold off;

clear cc_precontrols_roi cc_postcontrols_roi cc_premeditators_roi cc_postmeditators_roi clusteringcoef_pvals clusteringcoef_sigrois_any;
clear clusteringcoef_sigrois_any_index clusteringcoef_sigrois_group clusteringcoef_sigrois_group_index clusteringcoef_sigrois_tbyg clusteringcoef_sigrois_tbyg_index clusteringcoef_sigrois_treatment clusteringcoef_sigrois_treatment_index;
clear clusteringcoef_meandiff_controls clusteringcoef_meandiff_controls_siggroup clusteringcoef_meandiff_controls_sigtbyg clusteringcoef_meandiff_controls_sigtreatment clusteringcoef_meandiff_meditators clusteringcoef_meandiff_meditators_siggroup clusteringcoef_meandiff_meditators_sigtbyg clusteringcoef_meandiff_meditators_sigtreatment;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local Efficiency (binary undirected matrices w/thresh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Eloc_ctrl_pos_heat = zeros(nrois,20);
Eloc_ctrl_pos_neut = zeros(nrois,20);
Eloc_ctrl_pre_heat = zeros(nrois,20);
Eloc_ctrl_pre_neut = zeros(nrois,20);
Eloc_medt_pos_heat = zeros(nrois,20);
Eloc_medt_pos_neut = zeros(nrois,20);
Eloc_medt_pre_heat = zeros(nrois,20);
Eloc_medt_pre_neut = zeros(nrois,20);
 
for i = 1:20
    C = efficiency_bin(sub_bin_abs_ctrl_pos_heat(:,:,i),1);
    Eloc_ctrl_pos_heat(:,i) = C;
    C = efficiency_bin(sub_bin_abs_ctrl_pos_neut(:,:,i),1);
    Eloc_ctrl_pos_neut(:,i) = C;
    C = efficiency_bin(sub_bin_abs_ctrl_pre_heat(:,:,i),1);
    Eloc_ctrl_pre_heat(:,i) = C;
    C = efficiency_bin(sub_bin_abs_ctrl_pre_neut(:,:,i),1);
    Eloc_ctrl_pre_neut(:,i) = C;
    C = efficiency_bin(sub_bin_abs_medt_pos_heat(:,:,i),1);
    Eloc_medt_pos_heat(:,i) = C;
    C = efficiency_bin(sub_bin_abs_medt_pos_neut(:,:,i),1);
    Eloc_medt_pos_neut(:,i) = C;
    C = efficiency_bin(sub_bin_abs_medt_pre_heat(:,:,i),1);
    Eloc_medt_pre_heat(:,i) = C;
    C = efficiency_bin(sub_bin_abs_medt_pre_neut(:,:,i),1);
    Eloc_medt_pre_neut(:,i) = C;
end 
clear C

%change in Eloc
Eloc_chg_ctrl_heat = Eloc_ctrl_pos_heat - Eloc_ctrl_pre_heat;
Eloc_chg_ctrl_neut = Eloc_ctrl_pos_neut - Eloc_ctrl_pre_neut;
Eloc_chg_medt_heat = Eloc_medt_pos_heat - Eloc_medt_pre_heat;
Eloc_chg_medt_neut = Eloc_medt_pos_neut - Eloc_medt_pre_neut;
 
%mean change per roi
Eloc_mchg_ctrl_heat = mean(Eloc_chg_ctrl_heat,2);
Eloc_mchg_ctrl_neut = mean(Eloc_chg_ctrl_neut,2);
Eloc_mchg_medt_heat = mean(Eloc_chg_medt_heat,2);
Eloc_mchg_medt_neut = mean(Eloc_chg_medt_neut,2);

% RMAnova all rois
d = sprintf('\n\nLocal Efficiency: RMAnova for each ROI:\n');
disp(d);
% group var of size = #sub
A = [0]; X = 20; C = repmat(A,X,1); % controls
A = [1]; X = 20; M = repmat(A,X,1); % meditators
 
% loop all rois: rmanova on each for all sub, store p-vals in array
for roi = 1:nrois
   % extract roi
   Eloc_ctrl_pre_heat_roi = Eloc_ctrl_pre_heat(roi,:);
   Eloc_ctrl_pre_neut_roi = Eloc_ctrl_pre_neut(roi,:);
   Eloc_ctrl_pos_heat_roi = Eloc_ctrl_pos_heat(roi,:);
   Eloc_ctrl_pos_neut_roi = Eloc_ctrl_pos_neut(roi,:);
   Eloc_medt_pre_heat_roi = Eloc_medt_pre_heat(roi,:);
   Eloc_medt_pre_neut_roi = Eloc_medt_pre_neut(roi,:);
   Eloc_medt_pos_heat_roi = Eloc_medt_pos_heat(roi,:);
   Eloc_medt_pos_neut_roi = Eloc_medt_pos_neut(roi,:);
   % permute
   Eloc_ctrl_pre_heat_roi = Eloc_ctrl_pre_heat_roi';
   Eloc_ctrl_pre_neut_roi = Eloc_ctrl_pre_neut_roi';
   Eloc_ctrl_pos_heat_roi = Eloc_ctrl_pos_heat_roi';
   Eloc_ctrl_pos_neut_roi = Eloc_ctrl_pos_neut_roi';
   Eloc_medt_pre_heat_roi = Eloc_medt_pre_heat_roi';
   Eloc_medt_pre_neut_roi = Eloc_medt_pre_neut_roi';
   Eloc_medt_pos_heat_roi = Eloc_medt_pos_heat_roi';
   Eloc_medt_pos_neut_roi = Eloc_medt_pos_neut_roi';
   % ranova
 
   % NEUTRAL-only
   
   % define data 
    dataTable = array2table([Eloc_ctrl_pre_neut_roi,Eloc_ctrl_pos_neut_roi,C;...
        Eloc_medt_pre_neut_roi,Eloc_medt_pos_neut_roi,M]);
    % convert to table format
    dataTable.Properties.VariableNames = {'Pre','Post','Group'};
    % convert group to categorical
    dataTable.Group = categorical(dataTable.Group);
    % define levels of within subject factor
    wsVariable = table([0 1]','VariableNames',{'Treatment'});
    %run RANOVA (within subjects)
    rm = fitrm(dataTable,'Pre,Post~Group','WithinDesign',wsVariable);    
    ranovatbl = ranova(rm);
    %run ANOVA (between subjects)
    anovatbl = anova(rm);
    % store p values
    localeff_pvals(:,roi) = [ranovatbl.pValue; anovatbl.pValue];
    % clear temp vars
    clear dataTable;
    clear rm;
    clear anovatbl;
    clear ranovatbl;
end
 
% remove cols 3,4,6 to leave only 3 p-values (treatment, interaction, group)
localeff_pvals(3,:) = [];   %remove 3rd row
localeff_pvals(3,:) = [];   %remove 4th row
localeff_pvals(4,:) = [];   %remove 6th row
clear A X C M roi;
 
% Flag ROIs with p<0.05 and build subset array and index vars with them
 
% 1) which rois are significant? build 2 arrays with them:
% localeff_sigrois = zeroes for non-sig and the p-values for sig rois
% localeff_sigrois_index = index of sig rois
 
% 1a) For p significant for treatment, group, or t*g (any)
localeff_sigrois_any = zeros(3, nrois);
localeff_sigrois_any_index = [];
for i = 1:nrois
    if min(localeff_pvals(:,i)) <= 0.05
        localeff_sigrois_any(:,i) = localeff_pvals(:,i);
        localeff_sigrois_any_index = [localeff_sigrois_any_index, i];
    else
        localeff_sigrois_any(:,i) = zeros(3,1);
    end
end
% display how many significant rois were found
y = size(localeff_sigrois_any_index);
y = y(2);
x = ['Local Efficiency RM-ANOVA: ', num2str(y), ' significant ROIs out of ', ...
    num2str(nrois), ' for any effects'];
disp(x);
clear x;
clear y;
clear i;
 
% 1b) For p significant for treatment-only
localeff_sigrois_treatment = zeros(1, nrois);
localeff_sigrois_treatment_index = [];
for i = 1:nrois
    if localeff_pvals(1,i) <= 0.05
        localeff_sigrois_treatment(1,i) = localeff_pvals(1,i);
        localeff_sigrois_treatment_index = [localeff_sigrois_treatment_index, i];
    else
        localeff_sigrois_treatment(1,i) = zeros(1,1);
    end
end
% display how many significant rois were found
y = size(localeff_sigrois_treatment_index);
y = y(2);
x = ['Local Efficiency RM-ANOVA: ', num2str(y), ' significant ROIs out of ', ...
    num2str(nrois), ' for treatment effects'];
disp(x);
clear x;
clear y;
clear i;
 
% 1c) For p significant for treatment-by-group-interaction only
localeff_sigrois_tbyg = zeros(1, nrois);
localeff_sigrois_tbyg_index = [];
for i = 1:nrois
    if localeff_pvals(2,i) <= 0.05
        localeff_sigrois_tbyg(1,i) = localeff_pvals(2,i);
        localeff_sigrois_tbyg_index = [localeff_sigrois_tbyg_index, i];
    else
        localeff_sigrois_tbyg(1,i) = zeros(1,1);
    end
end
% display how many significant rois were found
y = size(localeff_sigrois_tbyg_index);
y = y(2);
x = ['Local Efficiency RM-ANOVA: ', num2str(y), ' significant ROIs out of ', ...
    num2str(nrois), ' for treatment-by-group interaction effects'];
disp(x);
clear x;
clear y;
clear i;
 
% 1d) For p significant for group-only
localeff_sigrois_group = zeros(1, nrois);
localeff_sigrois_group_index = [];
for i = 1:nrois
    if localeff_pvals(3,i) <= 0.05
        localeff_sigrois_group(1,i) = localeff_pvals(3,i);
        localeff_sigrois_group_index = [localeff_sigrois_group_index, i];
    else
        localeff_sigrois_group(1,i) = zeros(1,1);
    end
end
% display how many significant rois were found
y = size(localeff_sigrois_group_index);
y = y(2);
x = ['Local Efficiency RM-ANOVA: ', num2str(y), ' significant ROIs out of ', ...
    num2str(nrois), ' for between group effects'];
disp(x);
 
clear x y i
 
% Clustering Coefficient sig-ROI vectors, p<0.05:
%
% (a) all effects:
%     localeff_sigrois_any
%     localeff_sigrois_any_index
%
% (b) treatment effects:
%     localeff_sigrois_treatment
%     localeff_sigrois_treatment_index
%
% (b) treatment-group interaction effects:
%     localeff_sigrois_tbyg
%     localeff_sigrois_tbyg_index
%
% (b) group effects:
%     localeff_sigrois_group
%     localeff_sigrois_group_index
 
% 2) Get mean change for sig rois post-pre for both groups
 
% take index of sig p-values and change localeff_meandiff for both
% groups such that all non-sig columns are null.
 
% initialize vectors for 3 sig conditions/effect-types
% treatment-only
localeff_meandiff_controls_sigtreatment = Eloc_mchg_ctrl_neut;
localeff_meandiff_meditators_sigtreatment = Eloc_mchg_medt_neut;
% interaction-only
localeff_meandiff_controls_sigtbyg = Eloc_mchg_ctrl_neut;
localeff_meandiff_meditators_sigtbyg = Eloc_mchg_medt_neut;
% group-only
localeff_meandiff_controls_siggroup = Eloc_mchg_ctrl_neut;
localeff_meandiff_meditators_siggroup = Eloc_mchg_medt_neut;
 
% check if roi has significance; if so copy mean difference #, else null
% (a) treatment-only
for i = 1:nrois
    if ismember(i,localeff_sigrois_treatment_index) == 0
        localeff_meandiff_controls_sigtreatment(i) = nan;
        localeff_meandiff_meditators_sigtreatment(i) = nan;
    else
        localeff_meandiff_controls_sigtreatment(i) = Eloc_mchg_ctrl_neut(i);
        localeff_meandiff_meditators_sigtreatment(i) = Eloc_mchg_medt_neut(i);
    end
end
clear i
 
% (b) interaction-only
for i = 1:nrois
    if ismember(i,localeff_sigrois_tbyg_index) == 0
        localeff_meandiff_controls_sigtbyg(i) = nan;
        localeff_meandiff_meditators_sigtbyg(i) = nan;
    else
        localeff_meandiff_controls_sigtbyg(i) = Eloc_mchg_ctrl_neut(i);
        localeff_meandiff_meditators_sigtbyg(i) = Eloc_mchg_medt_neut(i);
    end
end
clear i
 
% (c) group-only
for i = 1:nrois
    if ismember(i,localeff_sigrois_group_index) == 0
        localeff_meandiff_controls_siggroup(i) = nan;
        localeff_meandiff_meditators_siggroup(i) = nan;
    else
        localeff_meandiff_controls_siggroup(i) = Eloc_mchg_ctrl_neut(i);
        localeff_meandiff_meditators_siggroup(i) = Eloc_mchg_medt_neut(i);
    end
end
clear i
 
% 3) Create stem plot for just significant rois (post-pre mean chg for both groups)
 
%ROInames2 = ROInames(:);
 
[~,labels] = readvars('atlas_labels_harvardoxford_cort-maxprob-thr25-2mm.csv');
labels = labels(2:49,:);
 
%get # of ROIs and set x 
%sizeROI = size(ROInames);
%sizeROI = sizeROI(2);
% x = linspace(1,sizeROI,sizeROI);
 
% plot
figure;
%treatment-only
stem(localeff_meandiff_controls_sigtreatment, 'b<');
xticks([1:nrois]);
xticklabels(labels);
xtickangle(45);
xlim([1 nrois]);
title({['{\bf\fontsize{14} Local Efficiency for ROIs with p<0.05 for Any Sig. Effect, Mean Change Post-Pre (Neutral)}' ,...
    ', corr ', num2str(thresh)]; 'local efficiency: global efficiency (average inverse shortest path length) computed on neighborhood of node'},'FontWeight','Normal')
ylabel('Change in Local Efficiency');
hold on;
stem(localeff_meandiff_meditators_sigtreatment, 'r<', 'filled');
%interaction-only
stem(localeff_meandiff_controls_sigtbyg, 'b*');
stem(localeff_meandiff_meditators_sigtbyg, 'r*');
%group-only
stem(localeff_meandiff_controls_siggroup, 'b>');
stem(localeff_meandiff_meditators_siggroup, 'r>', 'filled');
legend({'controls - sig. treatment effect','meditators - sig. treatment effect',...
    'controls - sig. T*G effect','meditators - sig. T*G effect',...
    'controls - sig. group effect','meditators - sig. group effect'},'Location','northeast');
hold off;
 
clear Eloc_precontrols_roi Eloc_postcontrols_roi Eloc_premeditators_roi Eloc_postmeditators_roi localeff_pvals localeff_sigrois_any;
clear localeff_sigrois_any_index localeff_sigrois_group localeff_sigrois_group_index localeff_sigrois_tbyg localeff_sigrois_tbyg_index localeff_sigrois_treatment localeff_sigrois_treatment_index;
clear localeff_meandiff_controls localeff_meandiff_controls_siggroup localeff_meandiff_controls_sigtbyg localeff_meandiff_controls_sigtreatment localeff_meandiff_meditators localeff_meandiff_meditators_siggroup localeff_meandiff_meditators_sigtbyg localeff_meandiff_meditators_sigtreatment;

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
