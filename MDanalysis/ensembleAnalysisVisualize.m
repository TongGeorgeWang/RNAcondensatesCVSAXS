
%% Visualize results of ensemble analysis on polymers in CG-RNA simulation 
%   Wrote a separate script for this so all the parameters do not need to get recomputed every run 
%
%   First, load results of ensembleAnalysis then run this script
%   If >1 dataset is desired to be plotted, you can load a dataset, generate Rg_mean, rename the variable or save in separate .mat file,
%       and then run it again for second dataset, then alter the plotting section to plot both in one figure 
%
%   GW - July 2024
%

close all
framesToPlot = (1:50:10000)'; % rebinned increments of 20-50 frames makes for a less cluttered graph 

colors = {[46 80 122]./255, [146 0 0]./255, [197 192 0]./255}; %rU blue, rA red, rC deep yellow
lw = 2; %plot linewidth 


%% Calculate parameters to plot (comment out if already calculated)
% Rg_mean_rU30alone = mean(Rg_ind,2); Rg_err_rU30alone = std(Rg_ind,0,2)./numel(Rg_ind(1,:));
% K2_mean_rU30alone = mean(K2_ind,2); K2_err_rU30alone = std(K2_ind,0,2)./numel(K2_ind(1,:));
% FA_mean_rU30alone = mean(FA_ind,2); FA_err_rU30alone = std(FA_ind,0,2)./numel(FA_ind(1,:));
% Ree_mean_rU30alone = mean(Ree_ind,2); Ree_err_rU30alone = std(Ree_ind,0,2)./numel(Ree_ind(1,:));
% T2_mean_rU30alone = mean(T2_ind,2); T2_err_rU30alone = std(T2_ind,0,2)./numel(T2_ind(1,:));
% b_mean_rU30alone = mean(b_ind,2); b_err_rU30alone = std(b_ind,0,2)./numel(b_ind(1,:));
% lOCF_mean_rU30alone = mean(l_OCF_ind,2); lOCF_mean_err_rU30alone = std(l_OCF_ind,0,2)./numel(l_OCF_ind(1,:));


%% Make plots across entire trajectory 

time = framesToPlot/10000 ; % in usec

figure; hold all
set(gcf,'color','white')
set(gca,'LineWidth',1.5,'FontSize',22)
grid on; box on
xlabel('Simulation time ($\mu$s)','Interpreter','latex')
ylabel('$R_{g} (\AA)$','Interpreter','latex')

%errorbar(time, Rg_mean_rU30alone(framesToPlot), Rg_err_rU30alone(framesToPlot), 'Color',colors{1},'LineWidth',lw)
%errorbar(time, Rg_mean_rA30alone(framesToPlot), Rg_err_rA30alone(framesToPlot), 'Color',colors{2},'LineWidth',lw)

errorbar(time, Rg_mean_rU30(framesToPlot), Rg_err_rU30(framesToPlot), 'Color',colors{1},'LineWidth',lw)
errorbar(time, Rg_mean_rA30(framesToPlot), Rg_err_rA30(framesToPlot), 'Color',colors{2},'LineWidth',lw)
errorbar(time, Rg_mean_rC30(framesToPlot), Rg_err_rC30(framesToPlot), 'Color',colors{3},'LineWidth',lw)
legend('rU30','rA30','rC30','Location','southeast')
%legend('rU30','rA30','Location','southeast')


figure; hold all
set(gcf,'color','white')
set(gca,'LineWidth',1.5,'FontSize',22)
grid on; box on
xlabel('Simulation time ($\mu$s)','Interpreter','latex')
ylabel('$\kappa^{2}$','Interpreter','latex')

%errorbar(time, K2_mean_rU30alone(framesToPlot), K2_err_rU30alone(framesToPlot), 'Color',colors{1},'LineWidth',lw)
%errorbar(time, K2_mean_rA30alone(framesToPlot), K2_err_rA30alone(framesToPlot), 'Color',colors{2},'LineWidth',lw)

errorbar(time, K2_mean_rU30(framesToPlot), K2_err_rU30(framesToPlot), 'Color',colors{1},'LineWidth',lw)
errorbar(time, K2_mean_rA30(framesToPlot), K2_err_rA30(framesToPlot), 'Color',colors{2},'LineWidth',lw)
errorbar(time, K2_mean_rC30(framesToPlot), K2_err_rC30(framesToPlot), 'Color',colors{3},'LineWidth',lw)
legend('rU30','rA30','rC30','Location','southeast')
%legend('rU30','rA30','Location','southeast')


figure; hold all
set(gcf,'color','white')
set(gca,'LineWidth',1.5,'FontSize',22)
grid on; box on
xlabel('Simulation time ($\mu$s)','Interpreter','latex')
ylabel('$R_{EE} (\AA)$','Interpreter','latex')

%errorbar(time, Ree_mean_rU30alone(framesToPlot), Ree_err_rU30alone(framesToPlot), 'Color',colors{1},'LineWidth',lw)
%errorbar(time, Ree_mean_rA30alone(framesToPlot), Ree_err_rA30alone(framesToPlot), 'Color',colors{2},'LineWidth',lw)

errorbar(time, Ree_mean_rU30(framesToPlot), Ree_err_rU30(framesToPlot), 'Color',colors{1},'LineWidth',lw)
errorbar(time, Ree_mean_rA30(framesToPlot), Ree_err_rA30(framesToPlot), 'Color',colors{2},'LineWidth',lw)
errorbar(time, Ree_mean_rC30(framesToPlot), Ree_err_rC30(framesToPlot), 'Color',colors{3},'LineWidth',lw)
legend('rU30','rA30','rC30','Location','southeast')
%legend('rU30','rA30','Location','southeast')



%% Plot early vs late chain conformational properties - RNA + peptide

earlyFrames = 50:1:150;
midFrames = 2000:1:2100;
lateFrames = 9900:1:10000;

% Rg
Rg_early_rU30 = Rg_mean_rU30(earlyFrames); % Used std as errors  
Rg_early_rA30 = Rg_mean_rA30(earlyFrames);
Rg_early_rC30 = Rg_mean_rC30(earlyFrames); 

Rg_mid_rU30 = Rg_mean_rU30(midFrames);
Rg_mid_rA30 = Rg_mean_rA30(midFrames);
Rg_mid_rC30 = Rg_mean_rC30(midFrames);

Rg_late_rU30 = Rg_mean_rU30(lateFrames);
Rg_late_rA30 = Rg_mean_rA30(lateFrames);
Rg_late_rC30 = Rg_mean_rC30(lateFrames);
 
% Ree
Ree_early_rU30 = Ree_mean_rU30(earlyFrames);
Ree_early_rA30 = Ree_mean_rA30(earlyFrames);
Ree_early_rC30 = Ree_mean_rC30(earlyFrames);

Ree_mid_rU30 = Ree_mean_rU30(midFrames);
Ree_mid_rA30 = Ree_mean_rA30(midFrames);
Ree_mid_rC30 = Ree_mean_rC30(midFrames);

Ree_late_rU30 = Ree_mean_rU30(lateFrames);
Ree_late_rA30 = Ree_mean_rA30(lateFrames);
Ree_late_rC30 = Ree_mean_rC30(lateFrames);

% K2
K2_early_rU30 = K2_mean_rU30(earlyFrames);
K2_early_rA30 = K2_mean_rA30(earlyFrames);
K2_early_rC30 = K2_mean_rC30(earlyFrames);

K2_mid_rU30 = K2_mean_rU30(midFrames);
K2_mid_rA30 = K2_mean_rA30(midFrames);
K2_mid_rC30 = K2_mean_rC30(midFrames);

K2_late_rU30 = K2_mean_rU30(lateFrames);
K2_late_rA30 = K2_mean_rA30(lateFrames);
K2_late_rC30 = K2_mean_rC30(lateFrames);


%% 

%% Box plot 
%Rgs
RgsBox = [Rg_early_rU30;Rg_late_rU30;Rg_early_rA30;Rg_late_rA30;Rg_early_rC30;Rg_late_rC30];
indices = [zeros(length(Rg_early_rU30), 1); ones(length(Rg_late_rU30), 1);...
            2*ones(length(Rg_early_rA30), 1); 3*ones(length(Rg_late_rA30), 1);...
            4*ones(length(Rg_early_rC30), 1); 5*ones(length(Rg_late_rC30), 1);];

figure; hold all
grid on; box on
set(gcf,'color','w')
ylabel('$R_{g} (\AA)$','Interpreter','latex')

colors = {[197 192 0]./255, [197 192 0]./255, [146 0 0]./255,...
          [146 0 0]./255, [46 80 122]./255, [46 80 122]./255};

bp = boxplot(RgsBox,indices,...
    'Labels',{'RU30 early','RU30 late','RA30 early','RA30 late','RC30 early','RC30 late'},...
    'Symbol','.','Whisker',1,'OutlierSize',8);

set(bp,{'linew'},{3})

box on
set(gcf,'color','w')
set(gca,'LineWidth',3)
set(gca,'FontSize',25)
ylim([15.9 25.1])

h = findobj(gca,'Tag','Box');
a = get(get(gca,'children'),'children');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors{j},'FaceAlpha',.7);
    a(j).MarkerEdgeColor = colors{j};
end

% T-tests
[h_RgrU,p_RgrU] = ttest2(Rg_early_rU30,Rg_late_rU30,'Vartype','unequal','Alpha',0.05);
[h_RgrA,p_RgrA] = ttest2(Rg_early_rA30,Rg_late_rA30,'Vartype','unequal','Alpha',0.05);
[h_RgrC,p_RgrC] = ttest2(Rg_early_rC30,Rg_late_rC30,'Vartype','unequal','Alpha',0.05);





%% Plot early vs late chain conformational properties - RNA alone (control)

earlyFrames = 50:1:150;
lateFrames = 9900:1:10000;

Rg_early_rU30alone = Rg_mean_rU30alone(earlyFrames); %Rg_err_rC30 = std(Rg_ind,0,2)./numel(Rg_ind(1,:));
Rg_early_rA30alone = Rg_mean_rA30alone(earlyFrames);

Rg_late_rU30alone = Rg_mean_rU30alone(lateFrames);
Rg_late_rA30alone = Rg_mean_rA30alone(lateFrames);

%% Box plot 
%Rgs
RgsBox = [Rg_early_rU30alone;Rg_late_rU30alone;Rg_early_rA30alone;Rg_late_rA30alone];
indices = [zeros(length(Rg_early_rU30alone), 1); ones(length(Rg_late_rU30alone), 1);...
            2*ones(length(Rg_early_rA30alone), 1); 3*ones(length(Rg_late_rA30alone), 1)];

figure; hold all
grid on; box on
set(gcf,'color','w')
ylabel('$R_{g} (\AA)$','Interpreter','latex')

colors = {[146 0 0]./255,...
          [146 0 0]./255, [46 80 122]./255, [46 80 122]./255};

bp = boxplot(RgsBox,indices,...
    'Labels',{'RU30 early','RU30 late','RA30 early','RA30 late'},...
    'Symbol','.','Whisker',1,'OutlierSize',8);

set(bp,{'linew'},{3})

box on
set(gcf,'color','w')
set(gca,'LineWidth',3)
set(gca,'FontSize',25)
ylim([15.9 25.1])

h = findobj(gca,'Tag','Box');
a = get(get(gca,'children'),'children');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors{j},'FaceAlpha',.7);
    a(j).MarkerEdgeColor = colors{j};
end

% T-tests
[h_RgrU,p_RgrU] = ttest2(Rg_early_rU30alone,Rg_late_rU30alone,'Vartype','unequal','Alpha',0.05);
[h_RgrA,p_RgrA] = ttest2(Rg_early_rA30alone,Rg_late_rA30alone,'Vartype','unequal','Alpha',0.05);


