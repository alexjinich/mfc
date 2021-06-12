clc;
close all;
clearvars;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% S or logS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if log_S == 1
    % transform S values by log to help visualize changes in S, with convention log(0)=0.
    sub_nS_ctrl_pos_heat = log(sub_nS_ctrl_pos_heat); sub_nS_ctrl_pos_heat(~isfinite(sub_nS_ctrl_pos_heat)) = 0;
    sub_nS_ctrl_pre_heat = log(sub_nS_ctrl_pre_heat); sub_nS_ctrl_pre_heat(~isfinite(sub_nS_ctrl_pre_heat)) = 0;
    sub_nS_ctrl_pos_neut = log(sub_nS_ctrl_pos_neut); sub_nS_ctrl_pos_neut(~isfinite(sub_nS_ctrl_pos_neut)) = 0;
    sub_nS_ctrl_pre_neut = log(sub_nS_ctrl_pre_neut); sub_nS_ctrl_pre_neut(~isfinite(sub_nS_ctrl_pre_neut)) = 0;
    sub_nS_medt_pos_heat = log(sub_nS_medt_pos_heat); sub_nS_medt_pos_heat(~isfinite(sub_nS_medt_pos_heat)) = 0;
    sub_nS_medt_pre_heat = log(sub_nS_medt_pre_heat); sub_nS_medt_pre_heat(~isfinite(sub_nS_medt_pre_heat)) = 0;
    sub_nS_medt_pos_neut = log(sub_nS_medt_pos_neut); sub_nS_medt_pos_neut(~isfinite(sub_nS_medt_pos_neut)) = 0;
    sub_nS_medt_pre_neut = log(sub_nS_medt_pre_neut); sub_nS_medt_pre_neut(~isfinite(sub_nS_medt_pre_neut)) = 0;
end
    
% change in entropy post-pre
nS_chg_ctrl_heat = sub_nS_ctrl_pos_heat - sub_nS_ctrl_pre_heat;
nS_chg_ctrl_neut = sub_nS_ctrl_pos_neut - sub_nS_ctrl_pre_neut;
nS_chg_medt_heat = sub_nS_medt_pos_heat - sub_nS_medt_pre_heat;
nS_chg_medt_neut = sub_nS_medt_pos_neut - sub_nS_medt_pre_neut;

% SE of change in S
SE_dS_ctrl_heat = std(nS_chg_ctrl_heat) / sqrt(20);
SE_dS_ctrl_neut = std(nS_chg_ctrl_neut) / sqrt(20);
SE_dS_medt_heat = std(nS_chg_medt_heat) / sqrt(20);
SE_dS_medt_neut = std(nS_chg_medt_neut) / sqrt(20);

% Plot change in S per group
x = categorical({'Controls Neutral','Controls Heat','Meditators Neutral','Meditators Heat'});
x = reordercats(x,{'Controls Neutral' 'Meditators Neutral' 'Controls Heat' 'Meditators Heat'});
y = [mean(nS_chg_ctrl_neut) mean(nS_chg_ctrl_heat) mean(nS_chg_medt_neut) mean(nS_chg_medt_heat)];
err = [SE_dS_ctrl_neut SE_dS_ctrl_heat SE_dS_medt_neut SE_dS_medt_heat];
figure; 
b=bar(x,y); 
b.FaceColor = 'flat';
b.CData(3,:) = [0.6350 0.0780 0.1840]; b.CData(4,:) = [0.6350 0.0780 0.1840];
hold on
er = errorbar(x,y,err); er.Color = 'k'; er.LineStyle = 'none';
hold off
%ylim([-.1 .4])
%str = sprintf('Change in Shannon Entropy (Post - Pre), Group-level at %.2f Abs. Corr. Threshold',thresh);
%sgtitle(str); 
ylabel('\Delta S','Rotation',0,'FontUnits','points','FontWeight','normal','FontSize',14)
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
clear x y b xtips1 ytips1

print -depsc2 chgS_groups.eps  % export


%% Plot S (pre|post) per individual
x1 = 1:20;
x2 = 21:40;
x = 1:40;
y = 4 * ones(1,40);

% find max and min, make min 0 if >0
mx1 = cat(2,sub_nS_ctrl_pre_heat,sub_nS_ctrl_pos_heat,sub_nS_medt_pre_heat,sub_nS_medt_pos_heat); 
mx1a = max(mx1); % find max1 value
mx2 = cat(2,sub_nS_ctrl_pre_neut,sub_nS_ctrl_pos_neut,sub_nS_medt_pre_neut,sub_nS_medt_pos_neut); 
mx2a = max(mx2); % find max2 value
mx = max(mx1a,mx2a); clear mx1* mx2*
mn1 = cat(2,sub_nS_ctrl_pre_heat,sub_nS_ctrl_pos_heat,sub_nS_medt_pre_heat,sub_nS_medt_pos_heat); 
mn1a = min(mn1); % find min1 value
mn2 = cat(2,sub_nS_ctrl_pre_neut,sub_nS_ctrl_pos_neut,sub_nS_medt_pre_neut,sub_nS_medt_pos_neut); 
mn2a = min(mn2); % find min2 value
mn = min(mn1a,mn2a); clear mn1* mn2*
% if mn > 0
%     mn = 0;
% end

% S_ctrl_pre_heat = mean(sub_nS_ctrl_pre_heat,2);
% S_ctrl_pos_heat = mean(sub_nS_ctrl_pos_heat,2);
% S_medt_pre_heat = mean(sub_nS_medt_pre_heat,2);
% S_medt_pos_heat = mean(sub_nS_medt_pos_heat,2);
%  
% S_ctrl_pre_neut = mean(sub_nS_ctrl_pre_neut,2);
% S_ctrl_pos_neut = mean(sub_nS_ctrl_pos_neut,2);
% S_medt_pre_neut = mean(sub_nS_medt_pre_neut,2);
% S_medt_pos_neut = mean(sub_nS_medt_pos_neut,2);

% Heat 
figure
%str = sprintf('Shannon Entropy Per Subject at %.2f Abs. Corr. Threshold',thresh);
%sgtitle(str);
tiledlayout(1,2)
nexttile
s1 = scatter(x1,sub_nS_ctrl_pre_heat,'v','b');
hold on
s2 = scatter(x1,sub_nS_ctrl_pos_heat,'^','b','filled');
s3 = scatter(x2,sub_nS_medt_pre_heat,'v','r');
s4 = scatter(x2,sub_nS_medt_pos_heat,'^','r','filled');
stem(x,y,':','MarkerEdgeColor','w')  % draw vertical ref lines
stem(x,-y,':','MarkerEdgeColor','w')
ylim([mn-.3 mx+.1])
hold off
ylabel('S','Rotation',0,'FontSize',12); xlabel('Subjects'); legend('Control Pre','Control Post','Meditator Pre','Meditator Post','Location','south','NumColumns',2);
title('Heat Runs','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');

% Neutral 
nexttile
s1 = scatter(x1,sub_nS_ctrl_pre_neut,'v','b');
hold on
s2 = scatter(x1,sub_nS_ctrl_pos_neut,'^','b','filled');
s3 = scatter(x2,sub_nS_medt_pre_neut,'v','r');
s4 = scatter(x2,sub_nS_medt_pos_neut,'^','r','filled');
stem(x,y,':','MarkerEdgeColor','w')
stem(x,-y,':','MarkerEdgeColor','w')
ylim([mn-.3 mx+.1])
hold off
ylabel('S','Rotation',0,'FontSize',12); xlabel('Subjects'); legend('Control Pre','Control Post','Meditator Pre','Meditator Post','Location','south','NumColumns',2);
title('Neutral Runs','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');

clear mx mn

set(gcf,'Position',[100,100,1300,400]);
print -depsc2 S_per_sub.eps  % export

%% All sub, chg in logS
x = 1:40;
x1 = 1:20;
x2 = 21:40;
% Heat
figure
%str = sprintf('Change in Shannon Entropy Per Subject at %.2f Corr. Threshold',thresh);
%sgtitle(str);
tiledlayout(1,2)
nexttile
stem(x1,nS_chg_ctrl_heat,'filled')
hold on
stem(x2,nS_chg_medt_heat,'filled')
hold off
title('Heat Runs','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
ylabel('\Delta S','Rotation',0,'FontSize',12); xlabel('Subjects','FontSize',12); legend('Control','Meditator','Location','south');
ylim([-.8 1]); xlim([0 41]);

% Neut
nexttile
stem(x1,nS_chg_ctrl_neut,'filled')
hold on
stem(x2,nS_chg_medt_neut,'filled')
hold off
title('Neutral Runs','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Calibri');
ylabel('\Delta S','Rotation',0,'FontSize',12); xlabel('Subjects','FontSize',12); legend('Control','Meditator','Location','south');
ylim([-.8 1]); xlim([0 41]);

clear x* y 
set(gcf,'Position',[100,100,1300,500]);

print -depsc2 chgS_per_sub.eps  % export


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
%% VAS Pain VS S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Scatter 1: abs S vs abs VAS pain per heat sub
% 
% % intensity
% str = sprintf('Entropy Level vs VAS Pain Intensity at %.2f Correlation Threshold',thresh);
% figure; tiledlayout(2,2); sgtitle(str)
% nexttile
% scatter(vas_int_h12_medt,sub_nS_medt_pre_heat,'*'); hold on
%  title('Meditators, Pre-Manip Heat'); 
% xlabel('pain intensity'); ylabel('S'); %axis([1.5 3 0 9])
% nexttile
% scatter(vas_int_h34_medt,sub_nS_medt_pos_heat,'*'); hold on
%  title('Meditators, Post-Manip Heat'); 
% xlabel('pain intensity'); ylabel('S'); %axis([1.5 3 0 9])
% nexttile
% scatter(vas_int_h12_ctrl,sub_nS_ctrl_pre_heat,'*'); hold on
%  title('Controls, Pre-Manip Heat'); 
% xlabel('pain intensity'); ylabel('S'); %axis([1.5 3 0 9])
% nexttile
% scatter(vas_int_h34_ctrl,sub_nS_ctrl_pos_heat,'*'); hold on
%  title('Controls, Post-Manip Heat'); 
% xlabel('pain intensity'); ylabel('S'); %axis([1.5 3 0 9])
% 
% % unpleasantness
% figure; tiledlayout(2,2); 
% str = sprintf('Entropy Level vs VAS Pain Unpleasantness at %.2f Correlation Threshold',thresh);
% sgtitle(str)
% nexttile
% scatter(vas_unp_h12_medt,sub_nS_medt_pre_heat,'*'); hold on
%  title('Meditators, Pre-Manip Heat'); 
% xlabel('pain unpleasantness'); ylabel('S'); %axis([1.6 3 0 10])
% nexttile
% scatter(vas_unp_h34_medt,sub_nS_medt_pos_heat,'*'); hold on
%  title('Meditators, Post-Manip Heat'); 
% xlabel('pain unpleasantness'); ylabel('S'); %axis([1.6 3 0 10])
% nexttile
% scatter(vas_unp_h12_ctrl,sub_nS_ctrl_pre_heat,'*'); hold on
%  title('Controls, Pre-Manip Heat'); 
% xlabel('pain unpleasantness'); ylabel('S'); %axis([1.6 3 0 10])
% nexttile
% scatter(vas_unp_h34_ctrl,sub_nS_ctrl_pos_heat,'*'); hold on
%  title('Controls, Post-Manip Heat'); 
% xlabel('pain unpleasantness'); ylabel('S'); %axis([1.6 3 0 10])
% 
% % regress S vs VAS

% Scatter 2: change in pain vs change in entropy

% Change in Entropy of Heat subs (nS of h3h4 - nS of h1h2)
x_tick = [-.5:.1:.6]; y_tick = [-6:1:4];
figure; tiledlayout(2,2); 
str = sprintf('Change in Entropy (Heat Runs) vs Change in VAS Pain at %.2f Correlation Threshold',thresh);
sgtitle(str)
nexttile
scatter(vas_chg_int_ctrl,nS_chg_ctrl_heat,'filled'); title('Controls, VAS Int');   
xlabel('\Delta pain int'); ylabel('\Delta S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
nexttile
scatter(vas_chg_int_medt,nS_chg_medt_heat,'filled'); title('Meditators, VAS Int'); 
xlabel('\Delta pain int'); ylabel('\Delta S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
nexttile
scatter(vas_chg_unp_ctrl,nS_chg_ctrl_heat,'filled'); title('Controls, VAS Unp');   
xlabel('\Delta pain unp'); ylabel('\Delta S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
nexttile
scatter(vas_chg_unp_medt,nS_chg_medt_heat,'filled'); title('Meditators, VAS Unp'); 
xlabel('\Delta pain unp'); ylabel('\Delta S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
clear x_tick y_tick

% Change in Entropy of Neutral subs (nS of n3n4 - nS of n1n2)
x_tick = [-.5:.1:.6]; y_tick = [-6:1:4];
figure; tiledlayout(2,2); 
str = sprintf('Change in Entropy (Neut Runs) vs Change in VAS Pain at %.2f Correlation Threshold',thresh);
sgtitle(str)
nexttile
scatter(vas_chg_int_ctrl,nS_chg_ctrl_neut,'filled'); title('Controls, VAS Int');   
xlabel('\Delta pain int'); ylabel('\Delta S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
nexttile
scatter(vas_chg_int_medt,nS_chg_medt_neut,'filled'); title('Meditators, VAS Int'); 
xlabel('\Delta pain int'); ylabel('\Delta S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
nexttile
scatter(vas_chg_unp_ctrl,nS_chg_ctrl_neut,'filled'); title('Controls, VAS Unp');   
xlabel('\Delta pain unp'); ylabel('\Delta S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
nexttile
scatter(vas_chg_unp_medt,nS_chg_medt_neut,'filled'); title('Meditators, VAS Unp'); 
xlabel('\Delta pain unp'); ylabel('\Delta S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
clear x_tick y_tick

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

%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Optional removal of 3 meditator-outliers based on qual
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

vas_chg_int_medt_outl = cat(1,vas_chg_int_medt(1:3), vas_chg_int_medt(6:10), vas_chg_int_medt(12:20));
vas_chg_unp_medt_outl = cat(1,vas_chg_unp_medt(1:3), vas_chg_unp_medt(6:10), vas_chg_unp_medt(12:20));
nS_chg_medt_neut_outl = cat(2,nS_chg_medt_neut(1:3), nS_chg_medt_neut(6:10), nS_chg_medt_neut(12:20));
nS_chg_medt_heat_outl = cat(2,nS_chg_medt_heat(1:3), nS_chg_medt_heat(6:10), nS_chg_medt_heat(12:20));

% Scatter 2: change in pain vs change in entropy without outliers

% Change in Entropy of Heat subs (nS of h3h4 - nS of h1h2)
x_tick = [-.5:.1:.6]; y_tick = [-6:1:4];
figure; tiledlayout(2,2); 
str = sprintf('Change in Entropy (Heat Runs) vs Change in VAS Pain at %.2f Correlation Threshold (no outliers)',thresh);
sgtitle(str)
nexttile
scatter(vas_chg_int_ctrl,nS_chg_ctrl_heat,'filled'); title('Controls, VAS Int');   
xlabel('\Delta pain int'); ylabel('\Delta S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
nexttile
scatter(vas_chg_int_medt_outl,nS_chg_medt_heat_outl,'filled'); title('Meditators, VAS Int'); 
xlabel('\Delta pain int'); ylabel('\Delta S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
nexttile
scatter(vas_chg_unp_ctrl,nS_chg_ctrl_heat,'filled'); title('Controls, VAS Unp');   
xlabel('\Delta pain unp'); ylabel('\Delta S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
nexttile
scatter(vas_chg_unp_medt_outl,nS_chg_medt_heat_outl,'filled'); title('Meditators, VAS Unp'); 
xlabel('\Delta pain unp'); ylabel('\Delta S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
clear x_tick y_tick

% Change in Entropy of Neutral subs (nS of n3n4 - nS of n1n2)
x_tick = [-.5:.1:.6]; y_tick = [-6:1:4];
figure; tiledlayout(2,2); 
str = sprintf('Change in Entropy (Neut Runs) vs Change in VAS Pain at %.2f Correlation Threshold (no outliers)',thresh);
sgtitle(str)
nexttile
scatter(vas_chg_int_ctrl,nS_chg_ctrl_neut,'filled'); title('Controls, VAS Int');   
xlabel('\Delta pain int'); ylabel('\Delta S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
nexttile
scatter(vas_chg_int_medt_outl,nS_chg_medt_neut_outl,'filled'); title('Meditators, VAS Int'); 
xlabel('\Delta pain int'); ylabel('\Delta S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
nexttile
scatter(vas_chg_unp_ctrl,nS_chg_ctrl_neut,'filled'); title('Controls, VAS Unp');   
xlabel('\Delta pain unp'); ylabel('\Delta S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
nexttile
scatter(vas_chg_unp_medt_outl,nS_chg_medt_neut_outl,'filled'); title('Meditators, VAS Unp'); 
xlabel('\Delta pain unp'); ylabel('\Delta S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
clear x_tick y_tick

% Heat subs - regress chgS vs chg VAS
%fprintf('\nChange in Entropy during Heat vs Change in VAS Pain during Heat\n');
%fprintf('\nVAS INT, Controls\n'); 
dS_heat_dInt_ctrl = fitlm(vas_chg_int_ctrl,nS_chg_ctrl_heat);
%fprintf('\nVAS INT, Meditators\n'); 
dS_heat_dInt_medt_outl = fitlm(vas_chg_int_medt_outl,nS_chg_medt_heat_outl);
%fprintf('\nVAS UNP, Controls\n'); 
dS_heat_dUnp_ctrl = fitlm(vas_chg_unp_ctrl,nS_chg_ctrl_heat);
%fprintf('\nVAS UNP, Meditators\n'); 
dS_heat_dUnp_medt_outl = fitlm(vas_chg_unp_medt_outl,nS_chg_medt_heat_outl);

% plot regressions
figure; tiledlayout(2,2); 
str = sprintf('Change in Entropy (Heat Runs) vs Change in VAS Pain at %.2f Correlation Threshold (no outliers)',thresh);
sgtitle(str)

nexttile
plot(dS_heat_dInt_ctrl); title({'Controls, VAS Int';''}); 
xlabel('Change in pain intensity'); ylabel('Change in S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_heat_dInt_ctrl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_heat_dInt_ctrl.Rsquared.Ordinary)];
annotation('textbox',[.13 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

nexttile
plot(dS_heat_dInt_medt_outl); title({'Meditators, VAS Int';''}); 
xlabel('Change in pain intensity'); ylabel('Change in S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_heat_dInt_medt_outl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_heat_dInt_medt_outl.Rsquared.Ordinary)];
annotation('textbox',[.59 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

nexttile
plot(dS_heat_dUnp_ctrl); title({'Controls, VAS Unp';''});   
xlabel('Change in pain unpleasantness'); ylabel('Change in S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_heat_dUnp_ctrl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_heat_dUnp_ctrl.Rsquared.Ordinary)];
annotation('textbox',[.13 0.44 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

nexttile
plot(dS_heat_dUnp_medt_outl); title({'Meditators, VAS Unp';''}); 
xlabel('Change in pain unpleasantness'); ylabel('Change in S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_heat_dUnp_medt_outl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_heat_dUnp_medt_outl.Rsquared.Ordinary)];
annotation('textbox',[.59 0.44 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

% Neut subs - regress chg S vs chg VAS
%fprintf('\nChange in Entropy during Neut vs Change in VAS Pain during Heat\n');
%fprintf('\nVAS INT, Controls\n'); 
dS_neut_dInt_ctrl = fitlm(vas_chg_int_ctrl,nS_chg_ctrl_neut);
%fprintf('\nVAS INT, Meditators\n'); 
dS_neut_dInt_medt_outl = fitlm(vas_chg_int_medt_outl,nS_chg_medt_neut_outl);
%fprintf('\nVAS UNP, Controls\n'); 
dS_neut_dUnp_ctrl = fitlm(vas_chg_unp_ctrl,nS_chg_ctrl_neut);
%fprintf('\nVAS UNP, Meditators\n'); 
dS_neut_dUnp_medt_outl = fitlm(vas_chg_unp_medt_outl,nS_chg_medt_neut_outl);

% plot regressions
figure; tiledlayout(2,2); 
str = sprintf('Change in Entropy (Neutral Runs) vs Change in VAS Pain at %.2f Correlation Threshold (no outliers)',thresh);
sgtitle(str)

nexttile
plot(dS_neut_dInt_ctrl); title({'Controls, VAS Int';''}); 
xlabel('Change in pain intensity'); ylabel('Change in S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_neut_dInt_ctrl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_neut_dInt_ctrl.Rsquared.Ordinary)];
annotation('textbox',[.13 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

nexttile
plot(dS_neut_dInt_medt_outl); title({'Meditators, VAS Int';''}); 
xlabel('Change in pain intensity'); ylabel('Change in S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_neut_dInt_medt_outl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_neut_dInt_medt_outl.Rsquared.Ordinary)];
annotation('textbox',[.59 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

nexttile
plot(dS_neut_dUnp_ctrl); title({'Controls, VAS Unp';''});   
xlabel('Change in pain unpleasantness'); ylabel('Change in S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_neut_dUnp_ctrl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_neut_dUnp_ctrl.Rsquared.Ordinary)];
annotation('textbox',[.13 0.44 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

nexttile
plot(dS_neut_dUnp_medt_outl); title({'Meditators, VAS Unp';''}); 
xlabel('Change in pain unpleasantness'); ylabel('Change in S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_neut_dUnp_medt_outl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_neut_dUnp_medt_outl.Rsquared.Ordinary)];
annotation('textbox',[.59 0.44 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  


%%
% % scatter 3: change in pain vs abs entropy(pre & post)
% 
% % Abs entropy in Heat subs
% %x_tick = [-.5:.1:.6]; y_tick = [-6:1:4];
% figure; tiledlayout(2,2); 
% str = sprintf('Entropy (During Heat Runs) vs Change in VAS Pain at %.2f Correlation Threshold',thresh);
% sgtitle(str)
% nexttile
% scatter(vas_chg_int_ctrl,sub_nS_ctrl_pre_heat,'filled'); title('Controls, VAS Int'); hold on
% scatter(vas_chg_int_ctrl,sub_nS_ctrl_pos_heat,'filled'); legend('S pre','S pos');    hold off
% xlabel('\Delta pain int'); ylabel('ln(S)'); %axis([1.6 3 -4 4]); %xticks(x_tick); yticks(y_tick);
% nexttile
% scatter(vas_chg_int_medt,sub_nS_medt_pre_heat,'filled'); title('Meditators, VAS Int'); hold on
% scatter(vas_chg_int_medt,sub_nS_medt_pos_heat,'filled'); legend('S pre','S pos');    hold off
% xlabel('\Delta pain int'); ylabel('ln(S)'); %axis([1.6 3 -4 4]); %xticks(x_tick); yticks(y_tick);
% nexttile
% scatter(vas_chg_unp_ctrl,sub_nS_ctrl_pre_heat,'filled'); title('Controls, VAS unp'); hold on
% scatter(vas_chg_unp_ctrl,sub_nS_ctrl_pos_heat,'filled'); legend('S pre','S pos');    hold off
% xlabel('\Delta pain unp'); ylabel('ln(S)'); %axis([1.6 3 -4 4]); %xticks(x_tick); yticks(y_tick);
% nexttile
% scatter(vas_chg_unp_medt,sub_nS_medt_pre_heat,'filled'); title('Meditators, VAS unp'); hold on
% scatter(vas_chg_unp_medt,sub_nS_medt_pos_heat,'filled'); legend('S pre','S pos');    hold off
% xlabel('\Delta pain unp'); ylabel('ln(S)'); %axis([1.6 3 -4 4]); %xticks(x_tick); yticks(y_tick);
% % clear x_tick y_tick
% 
% % Abs Entropy in Neut subs
% % x_tick = [-.5:.1:.6]; y_tick = [-6:1:4];
% figure; tiledlayout(2,2); 
% str = sprintf('Entropy (During Neut Runs) vs Change in VAS Pain at %.2f Correlation Threshold',thresh);
% sgtitle(str)
% nexttile
% scatter(vas_chg_int_ctrl,sub_nS_ctrl_pre_neut,'filled'); title('Controls, VAS Int'); hold on
% scatter(vas_chg_int_ctrl,sub_nS_ctrl_pos_heat,'filled'); legend('S pre','S pos');    hold off
% xlabel('\Delta pain int'); ylabel('S'); %axis([1.6 3 -4 4]); %xticks(x_tick); yticks(y_tick);
% nexttile
% scatter(vas_chg_int_medt,sub_nS_medt_pre_neut,'filled'); title('Meditators, VAS Int'); hold on
% scatter(vas_chg_int_medt,sub_nS_medt_pos_heat,'filled'); legend('S pre','S pos');    hold off
% xlabel('\Delta pain int'); ylabel('S'); %axis([1.6 3 -4 4]); %xticks(x_tick); yticks(y_tick);
% nexttile
% scatter(vas_chg_unp_ctrl,sub_nS_ctrl_pre_neut,'filled'); title('Controls, VAS unp'); hold on
% scatter(vas_chg_unp_ctrl,sub_nS_ctrl_pos_heat,'filled'); legend('S pre','S pos');    hold off
% xlabel('\Delta pain unp'); ylabel('S'); %axis([1.6 3 -4 4]); %xticks(x_tick); yticks(y_tick);
% nexttile
% scatter(vas_chg_unp_medt,sub_nS_medt_pre_neut,'filled'); title('Meditators, VAS unp'); hold on
% scatter(vas_chg_unp_medt,sub_nS_medt_pos_heat,'filled'); legend('S pre','S pos');    hold off
% xlabel('\Delta pain unp'); ylabel('S'); %axis([1.6 3 -4 4]); %xticks(x_tick); yticks(y_tick);
% % clear x_tick y_tick
% 
% % regress S vs chg VAS
% % heat subs - regress S (pre) vs chg VAS
% %fprintf('\nEntropy during Heat Pre vs Change in VAS Pain during Heat\n');
% %fprintf('\nVAS INT, Controls\n'); 
% S_heat_dInt_ctrl = fitlm(vas_chg_int_ctrl,sub_nS_ctrl_pre_heat);
% %fprintf('\nVAS INT, Meditators\n'); 
% S_heat_dInt_medt = fitlm(vas_chg_int_medt,sub_nS_medt_pre_heat);
% %fprintf('\nVAS UNP, Controls\n'); 
% S_heat_dUnp_ctrl = fitlm(vas_chg_unp_ctrl,sub_nS_ctrl_pre_heat);
% %fprintf('\nVAS UNP, Meditators\n'); 
% S_heat_dUnp_medt = fitlm(vas_chg_unp_medt,sub_nS_medt_pre_heat);
% 
% % plot regressions
% figure; tiledlayout(2,2); 
% str = sprintf('Entropy (Heat Pre Runs) vs Change in VAS Pain at %.2f Correlation Threshold',thresh);
% sgtitle(str)
% 
% nexttile
% plot(S_heat_dInt_ctrl); title({'Controls, VAS Int';''}); 
% xlabel('Change in pain intensity'); ylabel('ln S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% str=['p = ',sprintf('%.2f',S_heat_dInt_ctrl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',S_heat_dInt_ctrl.Rsquared.Ordinary)];
% annotation('textbox',[.13 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  
% 
% nexttile
% plot(S_heat_dInt_medt); title({'Meditators, VAS Int';''}); 
% xlabel('Change in pain intensity'); ylabel('ln S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% str=['p = ',sprintf('%.2f',S_heat_dInt_medt.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',S_heat_dInt_medt.Rsquared.Ordinary)];
% annotation('textbox',[.59 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  
% 
% nexttile
% plot(S_heat_dUnp_ctrl); title({'Controls, VAS Unp';''});   
% xlabel('Change in pain unpleasantness'); ylabel('ln S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% str=['p = ',sprintf('%.2f',S_heat_dUnp_ctrl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',S_heat_dUnp_ctrl.Rsquared.Ordinary)];
% annotation('textbox',[.13 0.44 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  
% 
% nexttile
% plot(S_heat_dUnp_medt); title({'Meditators, VAS Unp';''}); 
% xlabel('Change in pain unpleasantness'); ylabel('ln S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% str=['p = ',sprintf('%.2f',S_heat_dUnp_medt.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',S_heat_dUnp_medt.Rsquared.Ordinary)];
% annotation('textbox',[.59 0.44 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  
% 
% % Neut subs - regress S vs chg VAS
% %fprintf('\nEntropy during Neut Pre vs Change in VAS Pain during Heat\n');
% %fprintf('\nVAS INT, Controls\n'); 
% S_neut_dInt_ctrl = fitlm(vas_chg_int_ctrl,sub_nS_ctrl_pre_neut);     % regress
% %fprintf('\nVAS INT, Meditators\n'); 
% S_neut_dInt_medt = fitlm(vas_chg_int_medt,sub_nS_medt_pre_neut);
% %fprintf('\nVAS UNP, Controls\n'); 
% S_neut_dUnp_ctrl = fitlm(vas_chg_unp_ctrl,sub_nS_ctrl_pre_neut);
% %fprintf('\nVAS UNP, Meditators\n'); 
% S_neut_dUnp_medt = fitlm(vas_chg_unp_medt,sub_nS_medt_pre_neut);
% 
% % plot regressions
% figure; tiledlayout(2,2); 
% str = sprintf('Entropy (Neut Pre Runs) vs Change in VAS Pain at %.2f Correlation Threshold',thresh);
% sgtitle(str)
% 
% nexttile
% plot(S_neut_dInt_ctrl); title({'Controls, VAS Int';''}); 
% xlabel('Change in pain intensity'); ylabel('ln S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% str=['p = ',sprintf('%.2f',S_neut_dInt_ctrl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',S_neut_dInt_ctrl.Rsquared.Ordinary)];
% annotation('textbox',[.13 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  
% 
% nexttile
% plot(S_neut_dInt_medt); title({'Meditators, VAS Int';''}); 
% xlabel('Change in pain intensity'); ylabel('ln S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% str=['p = ',sprintf('%.2f',S_neut_dInt_medt.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',S_neut_dInt_medt.Rsquared.Ordinary)];
% annotation('textbox',[.59 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  
% 
% nexttile
% plot(S_neut_dUnp_ctrl); title({'Controls, VAS Unp';''});   
% xlabel('Change in pain unpleasantness'); ylabel('ln S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% str=['p = ',sprintf('%.2f',S_neut_dUnp_ctrl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',S_neut_dUnp_ctrl.Rsquared.Ordinary)];
% annotation('textbox',[.13 0.44 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  
% 
% nexttile
% plot(S_neut_dUnp_medt); title({'Meditators, VAS Unp';''}); 
% xlabel('Change in pain unpleasantness'); ylabel('ln S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% str=['p = ',sprintf('%.2f',S_neut_dUnp_medt.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',S_neut_dUnp_medt.Rsquared.Ordinary)];
% annotation('textbox',[.59 0.44 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% S VS NPS (Pain Rating)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load
[subject group nps_h1 nps_h2 nps_h3 nps_h4 nps_pre nps_post nps_diff] = readvars('signature_responses-2.csv');

%Permute
nps_diff = permute(nps_diff,[2,1]);
nps_h1 = permute(nps_h1,[2,1]);
nps_h2 = permute(nps_h2,[2,1]);
nps_h3 = permute(nps_h3,[2,1]);
nps_h4 = permute(nps_h4,[2,1]);
nps_pre = permute(nps_pre,[2,1]);
nps_pos = permute(nps_post,[2,1]); clear nps_post;

%Ctrl & Medt Split
nps_diff_ctrl = nps_diff(1,1:20);
nps_diff_medt = nps_diff(1,21:40);
nps_pre_ctrl = nps_pre(1,1:20);
nps_pre_medt = nps_pre(1,21:40);
nps_pos_ctrl = nps_pos(1,1:20);
nps_pos_medt = nps_pos(1,21:40);
nps_h1_ctrl = nps_h1(1,1:20);
nps_h1_medt = nps_h1(1,21:40);
nps_h2_ctrl = nps_h2(1,1:20);
nps_h2_medt = nps_h2(1,21:40);
nps_h3_ctrl = nps_h3(1,1:20);
nps_h3_medt = nps_h3(1,21:40);
nps_h4_ctrl = nps_h4(1,1:20);
nps_h4_medt = nps_h4(1,21:40);

% scatter 2A, nps_diff vs dlogS fitlm

% Change in Entropy of Heat Runs (nS of h3h4 - nS of h1h2)
x_tick = [-.5:.1:.6]; y_tick = [-6:1:4];
figure; tiledlayout(2,2); 
str = sprintf('Change in Entropy vs Change in NPS Pain at %.2f Correlation Threshold',thresh);
sgtitle(str)
nexttile
scatter(nps_diff_ctrl,nS_chg_ctrl_heat,'filled'); title('Heat S, Controls');   
xlabel('\Delta NPS'); ylabel('\Delta log(S)'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
nexttile
scatter(nps_diff_medt,nS_chg_medt_heat,'filled'); title('Heat S, Meditators'); 
xlabel('\Delta NPS'); ylabel('\Delta log(S)'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% Change in Entropy of Neutral Runs (nS of n3n4 - nS of n1n2)
nexttile
scatter(nps_diff_ctrl,nS_chg_ctrl_neut,'filled'); title('Neut S, Controls');   
xlabel('\Delta NPS'); ylabel('\Delta log(S)'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
nexttile
scatter(nps_diff_medt,nS_chg_medt_neut,'filled'); title('Neut S, Meditators'); 
xlabel('\Delta NPS'); ylabel('\Delta log(S)'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
clear x_tick y_tick

% Heat Runs - regress chgS vs chg NPS

%fprintf('\nChange in Entropy during Heat vs Change in VAS Pain during Heat\n');
%fprintf('\nVAS INT, Controls\n'); 
dS_heat_dNPS_ctrl = fitlm(nps_diff_ctrl,nS_chg_ctrl_heat);
%fprintf('\nVAS INT, Meditators\n'); 
dS_heat_dNPS_medt = fitlm(nps_diff_medt,nS_chg_medt_heat);

% Neut Runs - regress chg S vs chg NPS

%fprintf('\nChange in Entropy during Neut vs Change in VAS Pain during Heat\n');
%fprintf('\nVAS INT, Controls\n'); 
dS_neut_dNPS_ctrl = fitlm(nps_diff_ctrl,nS_chg_ctrl_neut);
%fprintf('\nVAS INT, Meditators\n'); 
dS_neut_dNPS_medt = fitlm(nps_diff_medt,nS_chg_medt_neut);

% plot regressions
figure; tiledlayout(2,2); 
str = sprintf('Change in Entropy vs Change in NPS Pain at %.2f Correlation Threshold',thresh);
sgtitle(str)

nexttile
plot(dS_heat_dNPS_ctrl); title({'Controls, Heat S'}); 
xlabel('Change in pain'); ylabel('Change in Heat ln(S)'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_heat_dNPS_ctrl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_heat_dNPS_ctrl.Rsquared.Ordinary)];
annotation('textbox',[.13 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

nexttile
plot(dS_heat_dNPS_medt); title({'Meditators, Heat S'}); 
xlabel('Change in pain'); ylabel('Change in Heat ln(S)'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_heat_dNPS_medt.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_heat_dNPS_medt.Rsquared.Ordinary)];
annotation('textbox',[.59 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

nexttile
plot(dS_neut_dNPS_ctrl); title({'Controls, Neut S'}); 
xlabel('Change in pain'); ylabel('Change in Neutral ln(S)'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_neut_dNPS_ctrl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_neut_dNPS_ctrl.Rsquared.Ordinary)];
annotation('textbox',[.13 0.44 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

nexttile
plot(dS_neut_dNPS_medt); title({'Meditators, Neut S'}); 
xlabel('Change in pain'); ylabel('Change in Neutral ln(S)'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_neut_dNPS_medt.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_neut_dNPS_medt.Rsquared.Ordinary)];
annotation('textbox',[.59 0.44 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

% scatter 2B, nps_pos vs dlogS fitlm

% Scatters: Change in Entropy of Heat Runs vs NPS (post)
x_tick = [-.5:.1:.6]; y_tick = [-6:1:4];
figure; tiledlayout(2,2); 
str = sprintf('Change in Entropy vs NPS Pain (Post) at %.2f Correlation Threshold',thresh);
sgtitle(str)
nexttile
scatter(nps_pos_ctrl,nS_chg_ctrl_heat,'filled'); title('Heat dS, Controls');   
xlabel('NPS (post)'); ylabel('\Delta log(S)'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
nexttile
scatter(nps_pos_medt,nS_chg_medt_heat,'filled'); title('Heat dS, Meditators'); 
xlabel('NPS (post)'); ylabel('\Delta log(S)'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% Change in Entropy of Neutral Runs vs NPS (post)
nexttile
scatter(nps_pos_ctrl,nS_chg_ctrl_neut,'filled'); title('Neut dS, Controls');   
xlabel('NPS (post)'); ylabel('\Delta log(S)'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
nexttile
scatter(nps_pos_medt,nS_chg_medt_neut,'filled'); title('Neut dS, Meditators'); 
xlabel('NPS (post)'); ylabel('\Delta log(S)'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
clear x_tick y_tick

% Regress
dS_heat_NPS_pos_ctrl = fitlm(nps_pos_ctrl,nS_chg_ctrl_heat);
dS_heat_NPS_pos_medt = fitlm(nps_pos_medt,nS_chg_medt_heat);
dS_neut_NPS_pos_ctrl = fitlm(nps_pos_ctrl,nS_chg_ctrl_neut);
dS_neut_NPS_pos_medt = fitlm(nps_pos_medt,nS_chg_medt_neut);

% plot regressions
figure; tiledlayout(2,2); 
str = sprintf('Change in Entropy vs NPS Pain (Post) at %.2f Correlation Threshold',thresh);
sgtitle(str)

nexttile
plot(dS_heat_NPS_pos_ctrl); title({'Controls, Heat S'}); 
xlabel('NPS post'); ylabel('Change in Heat ln(S)'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_heat_NPS_pos_ctrl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_heat_NPS_pos_ctrl.Rsquared.Ordinary)];
annotation('textbox',[.13 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  
nexttile
plot(dS_heat_NPS_pos_medt); title({'Meditators, Heat S'}); 
xlabel('NPS post'); ylabel('Change in Heat ln(S)'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_heat_NPS_pos_medt.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_heat_NPS_pos_medt.Rsquared.Ordinary)];
annotation('textbox',[.59 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  
nexttile
plot(dS_neut_NPS_pos_ctrl); title({'Controls, Neut S'}); 
xlabel('NPS post'); ylabel('Change in Neutral ln(S)'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_neut_NPS_pos_ctrl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_neut_NPS_pos_ctrl.Rsquared.Ordinary)];
annotation('textbox',[.13 0.44 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  
nexttile
plot(dS_neut_NPS_pos_medt); title({'Meditators, Neut S'}); 
xlabel('NPS post'); ylabel('Change in Neutral ln(S)'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
str=['p = ',sprintf('%.3f',dS_neut_NPS_pos_medt.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_neut_NPS_pos_medt.Rsquared.Ordinary)];
annotation('textbox',[.59 0.44 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  

% scatter 2C: NPS vs logS fitlm

% Entropy vs NPS 
%x_tick = [-.5:.1:.6]; y_tick = [-6:1:4];
figure; tiledlayout(2,2); 
str = sprintf('Entropy vs NPS Pain at %.2f Correlation Threshold',thresh);
sgtitle(str)

nexttile
scatter(nps_pos_ctrl,sub_nS_ctrl_pos_heat,'filled'); hold on
scatter(nps_pre_ctrl,sub_nS_ctrl_pre_heat,'filled'); hold off
title('Heat S, Controls'); legend('post','pre');
xlabel('NPS'); ylabel('log(S) Heat'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
nexttile
scatter(nps_pos_medt,sub_nS_medt_pos_heat,'filled'); hold on
scatter(nps_pre_medt,sub_nS_medt_pre_heat,'filled'); hold off
title('Heat S, Meditators'); legend('post','pre');
xlabel('NPS'); ylabel('log(S) Heat'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
nexttile
scatter(nps_pos_ctrl,sub_nS_ctrl_pos_neut,'filled'); hold on
scatter(nps_pre_ctrl,sub_nS_ctrl_pre_neut,'filled'); hold off
title('Neutral S, Controls'); legend('post','pre');
xlabel('NPS'); ylabel('log(S) Neut'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
nexttile
scatter(nps_pos_medt,sub_nS_medt_pos_neut,'filled'); hold on
scatter(nps_pre_medt,sub_nS_medt_pre_neut,'filled'); hold off
title('Neutral S, Meditators'); legend('post','pre');
xlabel('NPS'); ylabel('log(S) Neut'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
%
%heatscatter(nps_pos_medt',sub_nS_medt_pos_heat,' ','test1.png')
%heatscatter(nps_pre_medt',sub_nS_medt_pre_heat,' ','test2.png')

%% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NPS vs S - minus outliers
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% nps_diff_medt_outl = cat(2,nps_diff_medt(1:3), nps_diff_medt(6:10), nps_diff_medt(12:20));
% %nS_chg_medt_neut_outl = cat(2,nS_chg_medt_neut(1:3), nS_chg_medt_neut(6:10), nS_chg_medt_neut(12:20));
% %nS_chg_medt_heat_outl = cat(2,nS_chg_medt_heat(1:3), nS_chg_medt_heat(6:10), nS_chg_medt_heat(12:20));
% nps_pos_medt_outl = cat(2,nps_pos_medt(1:3), nps_pos_medt(6:10), nps_pos_medt(12:20));
% 
% % Scatters: Change in Entropy of Heat Runs vs NPS (post)
% x_tick = [-.5:.1:.6]; y_tick = [-6:1:4];
% figure; tiledlayout(2,2); 
% str = sprintf('Change in Entropy vs NPS Pain (Post) at %.2f Correlation Threshold (no outliers)',thresh);
% sgtitle(str)
% nexttile
% scatter(nps_pos_ctrl,nS_chg_ctrl_heat,'filled'); title('Heat dS, Controls');   
% xlabel('NPS (post)'); ylabel('\Delta S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% nexttile
% scatter(nps_pos_medt_outl,nS_chg_medt_heat_outl,'filled'); title('Heat dS, Meditators'); 
% xlabel('NPS (post)'); ylabel('\Delta S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% % Change in Entropy of Neutral Runs vs NPS (post)
% nexttile
% scatter(nps_pos_ctrl,nS_chg_ctrl_neut,'filled'); title('Neutral dS, Controls');   
% xlabel('NPS (post)'); ylabel('\Delta S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% nexttile
% scatter(nps_pos_medt_outl,nS_chg_medt_neut_outl,'filled'); title('Neutral dS, Meditators'); 
% xlabel('NPS (post)'); ylabel('\Delta S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% clear x_tick y_tick
% 
% % Regress
% dS_heat_NPS_pos_ctrl = fitlm(nps_pos_ctrl,nS_chg_ctrl_heat);
% dS_heat_NPS_pos_medt_outl = fitlm(nps_pos_medt_outl,nS_chg_medt_heat_outl);
% dS_neut_NPS_pos_ctrl = fitlm(nps_pos_ctrl,nS_chg_ctrl_neut);
% dS_neut_NPS_pos_medt_outl = fitlm(nps_pos_medt_outl,nS_chg_medt_neut_outl);
% 
% % plot regressions
% figure; tiledlayout(2,2); 
% str = sprintf('Change in Entropy vs NPS Pain (Post) at %.2f Correlation Threshold (no outliers)',thresh);
% sgtitle(str)
% 
% nexttile
% plot(dS_heat_NPS_pos_ctrl); title({'Controls, Heat S'}); 
% xlabel('NPS post'); ylabel('Change in Heat S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% str=['p = ',sprintf('%.3f',dS_heat_NPS_pos_ctrl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_heat_NPS_pos_ctrl.Rsquared.Ordinary)];
% annotation('textbox',[.13 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  
% nexttile
% plot(dS_heat_NPS_pos_medt_outl); title({'Meditators, Heat S'}); 
% xlabel('NPS post'); ylabel('Change in Heat S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% str=['p = ',sprintf('%.3f',dS_heat_NPS_pos_medt_outl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_heat_NPS_pos_medt_outl.Rsquared.Ordinary)];
% annotation('textbox',[.59 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  
% nexttile
% plot(dS_neut_NPS_pos_ctrl); title({'Controls, Neut S'}); 
% xlabel('NPS post'); ylabel('Change in Neutral S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% str=['p = ',sprintf('%.3f',dS_neut_NPS_pos_ctrl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_neut_NPS_pos_ctrl.Rsquared.Ordinary)];
% annotation('textbox',[.13 0.44 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  
% nexttile
% plot(dS_neut_NPS_pos_medt_outl); title({'Meditators, Neut S'}); 
% xlabel('NPS post'); ylabel('Change in Neutral S'); %axis([-.5 .6 -4 4]); xticks(x_tick); yticks(y_tick);
% str=['p = ',sprintf('%.3f',dS_neut_NPS_pos_medt_outl.Coefficients.pValue(2)),', R^2 = ',sprintf('%.2f',dS_neut_NPS_pos_medt_outl.Rsquared.Ordinary)];
% annotation('textbox',[.59 0.44 0 0],'string',str,'FitBoxToText','on','EdgeColor','white')  
% 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Trait Mindfulness
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% mfc_spss = readtable('spss.csv');
% 
% % trait mindfulness table vars: (access as mfc_spss.var_name)
% % FMI:  v1_fmi v6_fmi 
% % MAAS: v1_maas v6_maas 
% % FFMQ: v1_ffmq_obs v1_ffmq_desc v1_ffmq_act_aware v1_ffmq_nonjudge v1_ffmq_nonreact
% %       v6_ffmq_obs v6_ffmq_desc v6_ffmq_act_aware v6_ffmq_nonjudge v6_ffmq_nonreact
% 
% % Split (ctrl/medt)
% TM.v6_fmi_ctrl = mfc_spss.v6_fmi(1:20,1);
% TM.v6_fmi_medt = mfc_spss.v6_fmi(21:40,1);
% TM.v1_fmi_ctrl = mfc_spss.v1_fmi(1:20,1);
% TM.v1_fmi_medt = mfc_spss.v1_fmi(21:40,1);
% TM.v6_maas_ctrl = mfc_spss.v6_maas(1:20,1);
% TM.v6_maas_medt = mfc_spss.v6_maas(21:40,1);
% TM.v1_maas_ctrl = mfc_spss.v1_maas(1:20,1);
% TM.v1_maas_medt = mfc_spss.v1_maas(21:40,1);
% TM.v6_ffmq_obs_ctrl = mfc_spss.v6_ffmq_obs(1:20,1);
% TM.v6_ffmq_obs_medt = mfc_spss.v6_ffmq_obs(21:40,1);
% TM.v1_ffmq_obs_ctrl = mfc_spss.v1_ffmq_obs(1:20,1);
% TM.v1_ffmq_obs_medt = mfc_spss.v1_ffmq_obs(21:40,1);
% TM.v6_ffmq_desc_ctrl = mfc_spss.v6_ffmq_desc(1:20,1);
% TM.v6_ffmq_desc_medt = mfc_spss.v6_ffmq_desc(21:40,1);
% TM.v1_ffmq_desc_ctrl = mfc_spss.v1_ffmq_desc(1:20,1);
% TM.v1_ffmq_desc_medt = mfc_spss.v1_ffmq_desc(21:40,1);
% TM.v6_ffmq_act_aware_ctrl = mfc_spss.v6_ffmq_act_aware(1:20,1);
% TM.v6_ffmq_act_aware_medt = mfc_spss.v6_ffmq_act_aware(21:40,1);
% TM.v1_ffmq_act_aware_ctrl = mfc_spss.v1_ffmq_act_aware(1:20,1);
% TM.v1_ffmq_act_aware_medt = mfc_spss.v1_ffmq_act_aware(21:40,1);
% TM.v6_ffmq_nonjudge_ctrl = mfc_spss.v6_ffmq_nonjudge(1:20,1);
% TM.v6_ffmq_nonjudge_medt = mfc_spss.v6_ffmq_nonjudge(21:40,1);
% TM.v1_ffmq_nonjudge_ctrl = mfc_spss.v1_ffmq_nonjudge(1:20,1);
% TM.v1_ffmq_nonjudge_medt = mfc_spss.v1_ffmq_nonjudge(21:40,1);
% TM.v6_ffmq_nonreact_ctrl = mfc_spss.v6_ffmq_nonreact(1:20,1);
% TM.v6_ffmq_nonreact_medt = mfc_spss.v6_ffmq_nonreact(21:40,1);
% TM.v1_ffmq_nonreact_ctrl = mfc_spss.v1_ffmq_nonreact(1:20,1);
% TM.v1_ffmq_nonreact_medt = mfc_spss.v1_ffmq_nonreact(21:40,1);
% 
% % Trait Mindfulness vs Entropy (Does S covary w/trait mindfulness?)
% 
% % FMI
% 
%     figure
%     tiledlayout(2,2)
% 
%     % Trait mindfulness (pre)|(post)|ave(pre,post) *VS* rs entropy (pre)|(post)|ave(pre,post)
% 
%     % TM (post) vs S (post)
%     nexttile
%     scatter(TM.v6_fmi_ctrl,sub_nS_ctrl_pos_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_fmi_medt,sub_nS_medt_pos_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (FMI) Post VS Log S (Post, Neutral)');
%     xlabel('FMI'); ylabel('log(s)'); legend('controls','meditators');
% 
%     % TM (pre) vs chg S
%     
%     nexttile
%     scatter(TM.v1_fmi_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v1_fmi_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (FMI) Pre VS Change in Log S (Neutral)');
%     xlabel('FMI'); ylabel('log(s)'); legend('controls','meditators');
%     
%     % Trait mindfulness (pre)|(post)|ave(pre,post)|chg(post-pre) *VS* change in entropy (post-pre)
% 
%     % TM (post) vs chg S
%     nexttile
%     scatter(TM.v6_fmi_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_fmi_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (FMI) Post VS Change in S (Neutral');
%     xlabel('FMI'); ylabel('\Delta log(s)'); legend('controls','meditators');
% 
%     % Chg TM (post-pre) vs chg S
%     nexttile
%     scatter(TM.v6_fmi_ctrl-TM.v1_fmi_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_fmi_medt-TM.v1_fmi_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Change in Trait Mindfulness (FMI) VS Change in S (Neutral)');
%     xlabel('\Delta FMI'); ylabel('\Delta log(s)'); legend('controls','meditators');
% 
% % maas
%  
%     figure
%     tiledlayout(2,2)
%  
%     % Trait mindfulness (pre)|(post)|ave(pre,post) *VS* rs entropy (pre)|(post)|ave(pre,post)
%  
%     % TM (post) vs S (post)
%     nexttile
%     scatter(TM.v6_maas_ctrl,sub_nS_ctrl_pos_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_maas_medt,sub_nS_medt_pos_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (maas) Post VS Log S (Post, Neutral)');
%     xlabel('maas'); ylabel('log(s)'); legend('controls','meditators');
%  
%     % TM (pre) vs chg S
%     
%     nexttile
%     scatter(TM.v1_maas_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v1_maas_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (maas) Pre VS Change in Log S (Neutral)');
%     xlabel('FMI'); ylabel('log(s)'); legend('controls','meditators');
%     
%     % Trait mindfulness (pre)|(post)|ave(pre,post)|chg(post-pre) *VS* change in entropy (post-pre)
%  
%     % TM (post) vs chg S
%     nexttile
%     scatter(TM.v6_maas_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_maas_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (maas) Post VS Change in S (Neutral');
%     xlabel('maas'); ylabel('\Delta log(s)'); legend('controls','meditators');
%  
%     % Chg TM (post-pre) vs chg S
%     nexttile
%     scatter(TM.v6_maas_ctrl-TM.v1_maas_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_maas_medt-TM.v1_maas_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Change in Trait Mindfulness (maas) VS Change in S (Neutral)');
%     xlabel('\Delta maas'); ylabel('\Delta log(s)'); legend('controls','meditators');
% 
% % ffmq_obs
%  
%     figure
%     tiledlayout(2,2)
%  
%     % Trait mindfulness (pre)|(post)|ave(pre,post) *VS* rs entropy (pre)|(post)|ave(pre,post)
%  
%     % TM (post) vs S (post)
%     nexttile
%     scatter(TM.v6_ffmq_obs_ctrl,sub_nS_ctrl_pos_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_ffmq_obs_medt,sub_nS_medt_pos_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (ffmq_obs) Post VS Log S (Post, Neutral)');
%     xlabel('ffmq_obs'); ylabel('log(s)'); legend('controls','meditators');
%  
%     % TM (pre) vs chg S
%     
%     nexttile
%     scatter(TM.v1_ffmq_obs_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v1_ffmq_obs_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (ffmq_obs) Pre VS Change in Log S (Neutral)');
%     xlabel('FMI'); ylabel('log(s)'); legend('controls','meditators');
%     
%     % Trait mindfulness (pre)|(post)|ave(pre,post)|chg(post-pre) *VS* change in entropy (post-pre)
%  
%     % TM (post) vs chg S
%     nexttile
%     scatter(TM.v6_ffmq_obs_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_ffmq_obs_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (ffmq_obs) Post VS Change in S (Neutral');
%     xlabel('ffmq_obs'); ylabel('\Delta log(s)'); legend('controls','meditators');
%  
%     % Chg TM (post-pre) vs chg S
%     nexttile
%     scatter(TM.v6_ffmq_obs_ctrl-TM.v1_ffmq_obs_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_ffmq_obs_medt-TM.v1_ffmq_obs_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Change in Trait Mindfulness (ffmq_obs) VS Change in S (Neutral)');
%     xlabel('\Delta ffmq_obs'); ylabel('\Delta log(s)'); legend('controls','meditators');
% 
% % ffmq_desc
%  
%     figure
%     tiledlayout(2,2)
%  
%     % Trait mindfulness (pre)|(post)|ave(pre,post) *VS* rs entropy (pre)|(post)|ave(pre,post)
%  
%     % TM (post) vs S (post)
%     nexttile
%     scatter(TM.v6_ffmq_desc_ctrl,sub_nS_ctrl_pos_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_ffmq_desc_medt,sub_nS_medt_pos_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (ffmq_desc) Post VS Log S (Post, Neutral)');
%     xlabel('ffmq_desc'); ylabel('log(s)'); legend('controls','meditators');
%  
%     % TM (pre) vs chg S
%     
%     nexttile
%     scatter(TM.v1_ffmq_desc_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v1_ffmq_desc_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (ffmq_desc) Pre VS Change in Log S (Neutral)');
%     xlabel('FMI'); ylabel('log(s)'); legend('controls','meditators');
%     
%     % Trait mindfulness (pre)|(post)|ave(pre,post)|chg(post-pre) *VS* change in entropy (post-pre)
%  
%     % TM (post) vs chg S
%     nexttile
%     scatter(TM.v6_ffmq_desc_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_ffmq_desc_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (ffmq_desc) Post VS Change in S (Neutral');
%     xlabel('ffmq_desc'); ylabel('\Delta log(s)'); legend('controls','meditators');
%  
%     % Chg TM (post-pre) vs chg S
%     nexttile
%     scatter(TM.v6_ffmq_desc_ctrl-TM.v1_ffmq_desc_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_ffmq_desc_medt-TM.v1_ffmq_desc_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Change in Trait Mindfulness (ffmq_desc) VS Change in S (Neutral)');
%     xlabel('\Delta ffmq_desc'); ylabel('\Delta log(s)'); legend('controls','meditators');
% 
% % ffmq_act_aware
%  
%     figure
%     tiledlayout(2,2)
%  
%     % Trait mindfulness (pre)|(post)|ave(pre,post) *VS* rs entropy (pre)|(post)|ave(pre,post)
%  
%     % TM (post) vs S (post)
%     nexttile
%     scatter(TM.v6_ffmq_act_aware_ctrl,sub_nS_ctrl_pos_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_ffmq_act_aware_medt,sub_nS_medt_pos_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (ffmq_act_aware) Post VS Log S (Post, Neutral)');
%     xlabel('ffmq_act_aware'); ylabel('log(s)'); legend('controls','meditators');
%  
%     % TM (pre) vs chg S
%     
%     nexttile
%     scatter(TM.v1_ffmq_act_aware_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v1_ffmq_act_aware_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (ffmq_act_aware) Pre VS Change in Log S (Neutral)');
%     xlabel('FMI'); ylabel('log(s)'); legend('controls','meditators');
%     
%     % Trait mindfulness (pre)|(post)|ave(pre,post)|chg(post-pre) *VS* change in entropy (post-pre)
%  
%     % TM (post) vs chg S
%     nexttile
%     scatter(TM.v6_ffmq_act_aware_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_ffmq_act_aware_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (ffmq_act_aware) Post VS Change in S (Neutral');
%     xlabel('ffmq_act_aware'); ylabel('\Delta log(s)'); legend('controls','meditators');
%  
%     % Chg TM (post-pre) vs chg S
%     nexttile
%     scatter(TM.v6_ffmq_act_aware_ctrl-TM.v1_ffmq_act_aware_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_ffmq_act_aware_medt-TM.v1_ffmq_act_aware_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Change in Trait Mindfulness (ffmq_act_aware) VS Change in S (Neutral)');
%     xlabel('\Delta ffmq_act_aware'); ylabel('\Delta log(s)'); legend('controls','meditators');
% 
% % ffmq_nonjudge
%  
%     figure
%     tiledlayout(2,2)
%  
%     % Trait mindfulness (pre)|(post)|ave(pre,post) *VS* rs entropy (pre)|(post)|ave(pre,post)
%  
%     % TM (post) vs S (post)
%     nexttile
%     scatter(TM.v6_ffmq_nonjudge_ctrl,sub_nS_ctrl_pos_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_ffmq_nonjudge_medt,sub_nS_medt_pos_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (ffmq_nonjudge) Post VS Log S (Post, Neutral)');
%     xlabel('ffmq_nonjudge'); ylabel('log(s)'); legend('controls','meditators');
%  
%     % TM (pre) vs chg S
%     
%     nexttile
%     scatter(TM.v1_ffmq_nonjudge_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v1_ffmq_nonjudge_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (ffmq_nonjudge) Pre VS Change in Log S (Neutral)');
%     xlabel('FMI'); ylabel('log(s)'); legend('controls','meditators');
%     
%     % Trait mindfulness (pre)|(post)|ave(pre,post)|chg(post-pre) *VS* change in entropy (post-pre)
%  
%     % TM (post) vs chg S
%     nexttile
%     scatter(TM.v6_ffmq_nonjudge_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_ffmq_nonjudge_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (ffmq_nonjudge) Post VS Change in S (Neutral');
%     xlabel('ffmq_nonjudge'); ylabel('\Delta log(s)'); legend('controls','meditators');
%  
%     % Chg TM (post-pre) vs chg S
%     nexttile
%     scatter(TM.v6_ffmq_nonjudge_ctrl-TM.v1_ffmq_nonjudge_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_ffmq_nonjudge_medt-TM.v1_ffmq_nonjudge_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Change in Trait Mindfulness (ffmq_nonjudge) VS Change in S (Neutral)');
%     xlabel('\Delta ffmq_nonjudge'); ylabel('\Delta log(s)'); legend('controls','meditators');
% 
% % ffmq_nonreact
%  
%     figure
%     tiledlayout(2,2)
%  
%     % Trait mindfulness (pre)|(post)|ave(pre,post) *VS* rs entropy (pre)|(post)|ave(pre,post)
%  
%     % TM (post) vs S (post)
%     nexttile
%     scatter(TM.v6_ffmq_nonreact_ctrl,sub_nS_ctrl_pos_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_ffmq_nonreact_medt,sub_nS_medt_pos_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (ffmq_nonreact) Post VS Log S (Post, Neutral)');
%     xlabel('ffmq_nonreact'); ylabel('log(s)'); legend('controls','meditators');
%  
%     % TM (pre) vs chg S
%     
%     nexttile
%     scatter(TM.v1_ffmq_nonreact_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v1_ffmq_nonreact_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (ffmq_nonreact) Pre VS Change in Log S (Neutral)');
%     xlabel('FMI'); ylabel('log(s)'); legend('controls','meditators');
%     
%     % Trait mindfulness (pre)|(post)|ave(pre,post)|chg(post-pre) *VS* change in entropy (post-pre)
%  
%     % TM (post) vs chg S
%     nexttile
%     scatter(TM.v6_ffmq_nonreact_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_ffmq_nonreact_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Trait Mindfulness (ffmq_nonreact) Post VS Change in S (Neutral');
%     xlabel('ffmq_nonreact'); ylabel('\Delta log(s)'); legend('controls','meditators');
%  
%     % Chg TM (post-pre) vs chg S
%     nexttile
%     scatter(TM.v6_ffmq_nonreact_ctrl-TM.v1_ffmq_nonreact_ctrl,nS_chg_ctrl_neut','filled') % controls, neutral
%     hold on
%     scatter(TM.v6_ffmq_nonreact_medt-TM.v1_ffmq_nonreact_medt,nS_chg_medt_neut','filled') % meditatr, neutral
%     hold off
%     title('Change in Trait Mindfulness (ffmq_nonreact) VS Change in S (Neutral)');
%     xlabel('\Delta ffmq_nonreact'); ylabel('\Delta log(s)'); legend('controls','meditators');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Modularity (per subject)
% % inputs: weighted corr matrices
% % asymmetric treatment of negative weights - Rubinov & Sporns (2011) Neuroimage 56:2068-79.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% mod_ctrl_pos_heat = zeros(1,20);
% mod_ctrl_pos_neut = zeros(1,20);
% mod_ctrl_pre_heat = zeros(1,20);
% mod_ctrl_pre_neut = zeros(1,20);
% mod_medt_pos_heat = zeros(1,20);
% mod_medt_pos_neut = zeros(1,20);
% mod_medt_pre_heat = zeros(1,20);
% mod_medt_pre_neut = zeros(1,20);
% 
% for i = 1:20
%     [M, Q] = community_louvain(sub_ctrl_posmanip_heat(:,:,i),[],[],'negative_asym'); mod_ctrl_pos_heat(i) = Q;
%     [M, Q] = community_louvain(sub_ctrl_posmanip_neut(:,:,i),[],[],'negative_asym'); mod_ctrl_pos_neut(i) = Q;
%     [M, Q] = community_louvain(sub_ctrl_premanip_heat(:,:,i),[],[],'negative_asym'); mod_ctrl_pre_heat(i) = Q;
%     [M, Q] = community_louvain(sub_ctrl_premanip_neut(:,:,i),[],[],'negative_asym'); mod_ctrl_pre_neut(i) = Q;
%     [M, Q] = community_louvain(sub_medt_posmanip_heat(:,:,i),[],[],'negative_asym'); mod_medt_pos_heat(i) = Q;
%     [M, Q] = community_louvain(sub_medt_posmanip_neut(:,:,i),[],[],'negative_asym'); mod_medt_pos_neut(i) = Q;
%     [M, Q] = community_louvain(sub_medt_premanip_heat(:,:,i),[],[],'negative_asym'); mod_medt_pre_heat(i) = Q;
%     [M, Q] = community_louvain(sub_medt_premanip_neut(:,:,i),[],[],'negative_asym'); mod_medt_pre_neut(i) = Q;
% end
% clear M Q
% 
% % change in modularity post-pre
% mod_chg_ctrl_heat = mod_ctrl_pos_heat - mod_ctrl_pre_heat;
% mod_chg_ctrl_neut = mod_ctrl_pos_neut - mod_ctrl_pre_neut;
% mod_chg_medt_heat = mod_medt_pos_heat - mod_medt_pre_heat;
% mod_chg_medt_neut = mod_medt_pos_neut - mod_medt_pre_neut;
%  
% % SE of change in modularity
% SE_dmod_ctrl_heat = std(mod_chg_ctrl_heat) / sqrt(20);
% SE_dmod_ctrl_neut = std(mod_chg_ctrl_neut) / sqrt(20);
% SE_dmod_medt_heat = std(mod_chg_medt_heat) / sqrt(20);
% SE_dmod_medt_neut = std(mod_chg_medt_neut) / sqrt(20);
% 
% % Plot change in modularity per group
% x = categorical({'Controls Neutral','Controls Heat','Meditators Neutral','Meditators Heat'});
% y = [mean(mod_chg_ctrl_neut) mean(mod_chg_ctrl_heat) mean(mod_chg_medt_neut) mean(mod_chg_medt_heat)];
% err = [SE_dmod_ctrl_neut SE_dmod_ctrl_heat SE_dmod_medt_neut SE_dmod_medt_heat];
% figure; 
% b=bar(x,y); 
% b.FaceColor = 'flat';
% b.CData(3,:) = [0.6350 0.0780 0.1840]; b.CData(4,:) = [0.6350 0.0780 0.1840];
% hold on
% er = errorbar(x,y,err); er.Color = 'k'; er.LineStyle = 'none';
% hold off
% str = sprintf('Change in Modularity (Post - Pre), Group-level at %.2f Abs. Corr. Threshold',thresh);
% sgtitle(str); ylabel('\Delta Modularity Post - Pre');
% xtips1 = b(1).XEndPoints;
% ytips1 = b(1).YEndPoints;
% labels1 = string(b(1).YData);
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
% clear x y b xtips1 ytips1
% 
% % All sub, chg in modularity
% x1=1:20; x2=21:40;
% % Heat
% figure
% str = sprintf('Change in Modularity Per Subject at %.2f Corr. Threshold',thresh);
% sgtitle(str);
% tiledlayout(1,2)
% nexttile
% stem(x1,mod_chg_ctrl_heat)
% hold on
% stem(x2,mod_chg_medt_heat)
% hold off
% title('Heat Runs'); ylabel('\Delta Modularity Post - Pre'); xlabel('subject'); legend('Control','Meditator');
% % Neut
% nexttile
% stem(x1,mod_chg_ctrl_neut)
% hold on
% stem(x2,mod_chg_medt_neut)
% hold off
% title('Neutral Runs'); ylabel('\Delta Modularity Post - Pre'); xlabel('subject'); legend('Control','Meditator');
% 
% clear x* y 
% 
% % RANOVA modularity
% mod_ctrl_pre_neut = permute(mod_ctrl_pre_neut,[2,1]);
% mod_ctrl_pos_neut = permute(mod_ctrl_pos_neut,[2,1]);
% mod_ctrl_pre_heat = permute(mod_ctrl_pre_heat,[2,1]);
% mod_ctrl_pos_heat = permute(mod_ctrl_pos_heat,[2,1]);
% mod_medt_pre_neut = permute(mod_medt_pre_neut,[2,1]);
% mod_medt_pos_neut = permute(mod_medt_pos_neut,[2,1]);
% mod_medt_pre_heat = permute(mod_medt_pre_heat,[2,1]);
% mod_medt_pos_heat = permute(mod_medt_pos_heat,[2,1]);
% 
% % NEUTRAL
% % Construct datatable with the 4 clusters of data and name the variables.
% datatable = table(mod_ctrl_pre_neut,mod_ctrl_pos_neut,mod_medt_pre_neut,mod_medt_pos_neut);
% datatable.Properties.VariableNames = {'pre_controls','pre_meditators','post_controls','post_meditators'};
% % Since we have more than one repeated-measures factor, we set up a table
% % to indicate the levels on each factor for each of the different variables.
% WithinStructure = table([1 1 2 2]',[1 2 1 2]','VariableNames',{'Treatment','Group'});
% % The 4 rows of the WithinStructure table correspond to 4 different cols
% % 'pre_A','pre_B','post_A','post_B' in data table. Each 'pre_A','pre_B','post_A','post_B' column is coded as 1/2 on the Treatment factor
% % and as 1/2 on the Group factor. Now pass the WithinStructure table to fitrm so that it knows how the different
% % columns correspond to the different levels of the repeated-measures factors.
% rm = fitrm(datatable, 'pre_controls,pre_meditators,post_controls,post_meditators~1','WithinDesign',WithinStructure);
% % Finally, specify the repeated-measures factors again when you call ranova & display results:
% disp('Repeated Measures ANOVA: Modularity, Neutral, 2 by 2:');
% ranovatable = ranova(rm,'WithinModel','Treatment*Group')
% %disp('Coefficients:'); 
% %rm.Coefficients
% %disp('Covariance:'); 
% %rm.Covariance
% % Clear temp vars.
% clear datatable WithinStructure rm ranovatable
% 
% % HEAT
% % Construct datatable with the 4 clusters of data and name the variables.
% datatable = table(mod_ctrl_pre_heat,mod_ctrl_pos_heat,mod_medt_pre_heat,mod_medt_pos_heat);
% datatable.Properties.VariableNames = {'pre_controls','pre_meditators','post_controls','post_meditators'};
% % Since we have more than one repeated-measures factor, we set up a table
% % to indicate the levels on each factor for each of the different variables.
% WithinStructure = table([1 1 2 2]',[1 2 1 2]','VariableNames',{'Treatment','Group'});
% % The 4 rows of the WithinStructure table correspond to 4 different cols
% % 'pre_A','pre_B','post_A','post_B' in data table. Each 'pre_A','pre_B','post_A','post_B' column is coded as 1/2 on the Treatment factor
% % and as 1/2 on the Group factor. Now pass the WithinStructure table to fitrm so that it knows how the different
% % columns correspond to the different levels of the repeated-measures factors.
% rm = fitrm(datatable, 'pre_controls,pre_meditators,post_controls,post_meditators~1','WithinDesign',WithinStructure);
% % Finally, specify the repeated-measures factors again when you call ranova & display results:
% disp('Repeated Measures ANOVA: Modularity, Heat, 2 by 2:');
% ranovatable = ranova(rm,'WithinModel','Treatment*Group')
% %disp('Coefficients:'); 
% %rm.Coefficients
% %disp('Covariance:'); 
% %rm.Covariance
% % Clear temp vars.
% clear datatable WithinStructure rm ranovatable
% 
% % T-tests per subject for sig diff in S pre vs post
% 
% disp('T-tests for Difference in Modularity Pre vs Post');
% 
% [h,p] = ttest(mod_ctrl_pre_neut, mod_ctrl_pos_neut); % Controls, Neutral
% if h == 0
%     h = '*not* rejected';
% else
%     h = '***rejected***';
% end
% d = sprintf('\nControls,   Neutral - Null H (no mean difference) is %s with p = %.3f', h, p); disp(d); 
%  
% [h,p] = ttest(mod_ctrl_pre_heat, mod_ctrl_pos_heat); % Controls, Heat
% if h == 0
%     h = '*not* rejected';
% else
%     h = '***rejected***';
% end
% d = sprintf('Controls,   Heat    - Null H (no mean difference) is %s with p = %.3f', h, p); disp(d); 
%  
% [h,p] = ttest(mod_medt_pre_neut, mod_medt_pos_neut); % Meditators, Neutral
% if h == 0
%     h = '*not* rejected';
% else
%     h = '***rejected***';
% end
% d = sprintf('Meditators, Neutral - Null H (no mean difference) is %s with p = %.3f', h, p); disp(d); 
%  
% [h,p] = ttest(mod_medt_pre_heat, mod_medt_pos_heat); % Meditators, Heat
% if h == 0
%     h = '*not* rejected';
% else
%     h = '***rejected***';
% end
% d = sprintf('Meditators, Heat    - Null H (no mean difference) is %s with p = %.3f', h, p); disp(d); 
%  
% clear h p d
% 
% 
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