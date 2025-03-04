
%close all

%% Plot similarities 
%   Import J and D values (Jaccard Indices and Dice-Sorensen Coefficients) resulting from running clusterGeometry.m 


%% Time plot
frames = 1:1:1000;

figure; hold all
plot(frames,J_rU30_5000_1_6000,'.-','Color',[46 80 122]./255)   
plot(frames,J_rA30_5000_1_6000,'.-','Color',[146 0 0]./255)
plot(frames,J_rC30_5000_1_6000,'.-','Color',[197 192 0]./255)



%% Box plot 
%JACCARD
Jaccards = [J_rU30_5000_1_6000';J_rA30_5000_1_6000';J_rC30_5000_1_6000'];
indices = [zeros(length(J_rU30_5000_1_6000), 1); ones(length(J_rA30_5000_1_6000), 1); ones(length(J_rC30_5000_1_6000), 1).*2];

figure; hold all
grid on; box on
set(gcf,'color','w')
ylabel('Jaccard index')

colors = { [197 192 0]./255, [146 0 0]./255, [46 80 122]./255};

b = boxplot(Jaccards,indices,...
    'Labels',{'RU30','RA30','RC30'},...
    'Symbol','.','Whisker',1,'OutlierSize',8);

set(b,{'linew'},{3})

box on
set(gcf,'color','w')
set(gca,'LineWidth',3)
set(gca,'FontSize',25)

h = findobj(gca,'Tag','Box');
a = get(get(gca,'children'),'children');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors{j},'FaceAlpha',.7);
    a(j).MarkerEdgeColor = colors{j};
end

% T-tests
[h_JrUrA,p_JrUrA] = ttest2(J_rU30_5000_1_6000,J_rA30_5000_1_6000,'Vartype','unequal','Alpha',0.05)
[h_JrCrA,p_JrCrA] = ttest2(J_rC30_5000_1_6000,J_rA30_5000_1_6000,'Vartype','unequal','Alpha',0.05)


%% DICE
Dices = [D_rU30_5000_1_6000';D_rA30_5000_1_6000';D_rC30_5000_1_6000'];
indices = [zeros(length(D_rU30_5000_1_6000), 1); ones(length(D_rA30_5000_1_6000), 1); ones(length(D_rC30_5000_1_6000), 1).*2];

figure; hold all
grid on; box on
set(gcf,'color','w')
ylabel('Dice-Sorensen coefficient')

colors = { [197 192 0]./255, [146 0 0]./255, [46 80 122]./255};

b = boxplot(Dices,indices,...
    'Labels',{'RU30','RA30','RC30'},...
    'Symbol','.','Whisker',1,'OutlierSize',8);

set(b,{'linew'},{3})

box on
set(gcf,'color','w')
set(gca,'LineWidth',3)
set(gca,'FontSize',25)

h = findobj(gca,'Tag','Box');
a = get(get(gca,'children'),'children');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors{j},'FaceAlpha',.7);
    a(j).MarkerEdgeColor = colors{j};
end

% T-tests
[h_DrUrA,p_DrUrA] = ttest2(D_rU30_5000_1_6000,D_rA30_5000_1_6000,'Vartype','unequal','Alpha',0.05)
[h_DrCrA,p_DrCrA] = ttest2(J_rC30_5000_1_6000,J_rA30_5000_1_6000,'Vartype','unequal','Alpha',0.05)

