
clc;
close all;
clearvars;

cd /Volumes/Seagate_Desktop_Drive/mfc/code/nilearn_matrices_0321
addpath '/Volumes/Seagate_Desktop_Drive/mfc/code'

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

%% Levels (Group) & SEM

% Intensity
   VAS_int_ctrl_pre = mean(vas_int_h12_ctrl);         % pre, controls
VAS_int_ctrl_pre_SE = std(vas_int_h12_ctrl)/sqrt(20);
   VAS_int_ctrl_pos = mean(vas_int_h34_ctrl);         % pos, controls
VAS_int_ctrl_pos_SE = std(vas_int_h34_ctrl)/sqrt(20);
   VAS_int_medt_pre = mean(vas_int_h12_medt);         % pre, meditators
VAS_int_medt_pre_SE = std(vas_int_h12_medt)/sqrt(20);
   VAS_int_medt_pos = mean(vas_int_h34_medt);         % pos, meditators
VAS_int_medt_pos_SE = std(vas_int_h34_medt)/sqrt(20);

% Unpleasantness
   VAS_unp_ctrl_pre = mean(vas_unp_h12_ctrl);         % pre, controls
VAS_unp_ctrl_pre_SE = std(vas_unp_h12_ctrl)/sqrt(20);
   VAS_unp_ctrl_pos = mean(vas_unp_h34_ctrl);         % pos, controls
VAS_unp_ctrl_pos_SE = std(vas_unp_h34_ctrl)/sqrt(20);
   VAS_unp_medt_pre = mean(vas_unp_h12_medt);         % pre, meditators
VAS_unp_medt_pre_SE = std(vas_unp_h12_medt)/sqrt(20);
   VAS_unp_medt_pos = mean(vas_unp_h34_medt);         % pos, meditators
VAS_unp_medt_pos_SE = std(vas_unp_h34_medt)/sqrt(20);

%% Plot 
x = categorical({'ctrl pre','medt pre','ctrl pos','medt pos'});
x = reordercats(x,{'ctrl pre','medt pre','ctrl pos','medt pos'});

y1 = [VAS_int_ctrl_pre VAS_int_medt_pre VAS_int_ctrl_pos VAS_int_medt_pos]; % intensity
err1 = [VAS_int_ctrl_pre_SE VAS_int_medt_pre_SE VAS_int_ctrl_pos_SE VAS_int_medt_pos_SE];

y2 = [VAS_unp_ctrl_pre VAS_unp_medt_pre VAS_unp_ctrl_pos VAS_unp_medt_pos]; % unpleasantness
err2 = [VAS_unp_ctrl_pre_SE VAS_unp_medt_pre_SE VAS_unp_ctrl_pos_SE VAS_unp_medt_pos_SE]; 

figure
tiledlayout(1,2)

nexttile % intensity
b = bar(x,y1);
b.FaceColor = 'flat';
b.CData(2,:) = [0.6350 0.0780 0.1840]; b.CData(4,:) = [0.6350 0.0780 0.1840];
hold on
er = errorbar(x,y1,err1); er.Color = 'k'; er.LineStyle = 'none';
hold off
ylim([0 10])
ylabel('VAS Pain Intensity Rating','Rotation',90,'FontUnits','points','FontWeight','normal','FontSize',14)

nexttile % unpleasantness
b = bar(x,y2);
b.FaceColor = 'flat';
b.CData(2,:) = [0.6350 0.0780 0.1840]; b.CData(4,:) = [0.6350 0.0780 0.1840];
hold on
er = errorbar(x,y2,err2); er.Color = 'k'; er.LineStyle = 'none';
hold off
ylim([0 10])
ylabel('VAS Pain Unp[leasantness Rating','Rotation',90,'FontUnits','points','FontWeight','normal','FontSize',14)

set(gcf,'Position',[500,500,500,450]); %xl = xlabel('d'); xl.FontSize = 20;

%% Change Calc

% Pain Intensity Reduction (Meditators)
Pntg_Int_Medt = (VAS_int_medt_pos - VAS_int_medt_pre) / VAS_int_medt_pre

% Pain Unpleasantness Reduction (Meditators)
Pntg_Unp_Medt = (VAS_unp_medt_pos - VAS_unp_medt_pre) / VAS_unp_medt_pre

%% 

