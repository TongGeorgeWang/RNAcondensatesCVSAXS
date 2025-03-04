
%% Quantify condensate (cluster) growth and constituent exchange across MD simulation frames
%
%   Designed for ssRNA-peptide Mpipi simulations, focusing on RNA condensation
%   Methods: 
%       Cluster growth: 
%           only consider clusters with >=0 constituents 
%           tabulate the number of RNAs in each cluster, across frames of simulation        
%           do the same for Rg         
%   
%       Constituent exchange: 
%           use the cluster graph from getRNA_v5 to determine what nodes are connected in each cluster
%           do the same thing for the next frame (i+1)
%           take the difference in adjacency matrices 
%
%   GW - October 2024
%

close all 

%% First, load the coordsAndClusters.mat and structuralParameters.mat workspaces of the molecule you wish to analyze

frames = 1:100:10000; % use this if you only want to view a subset of the total frames 

% Primary color and arrival color
%plotColor = [46/255, 80/255, 122/255]; %rU blue
plotColor = [146/255, 0/255, 0/255]; %rA red 
%plotColor = [197/255, 192/255, 0/255]; %rC yellow

% Lighter version of primary color for departure color
%plotColor2 = [134/255 167/255 219/255]; %rU light blue
plotColor2 = [214/255 129/255 129/255]; %rA light red 
%plotColor2 = [255/255 245/255 179/255]; %rC30 light yellow 

%colors = {[46 80 122]./255, [146 0 0]./255, [197 192 0]./255}; %rU blue, rA red, rC deep yellow



%% Quantify cluster growth

%nFrames = numel(clusterIndices(:,1));
nFrames = numel(frames);
numelClusters = cellfun(@numel, clusterIndices); % how many constituents each cluster has
time = frames/10000 ; % in usec


% Boxplot/swarmplot approach 
boxData = [];
boxDataRg = [];
labels = [];
labelsRg = [];
for i = 1:nFrames
    nConstit = numelClusters(frames(i),:);
    nConstit = nConstit(~ismember(nConstit,[0])); % remove clusters w/ 0 constituents (just matrix filler, unphysical)
    clusterRgs = Rg(frames(i),:);
    clusterRgs = clusterRgs(~ismember(clusterRgs,[0]));

    boxData = [boxData; nConstit'];
    labels = [labels; i*ones(length(nConstit), 1)];
    boxDataRg = [boxDataRg; clusterRgs'];
    labelsRg = [labelsRg; i*ones(length(clusterRgs), 1)];

    barY(i) = max(nConstit);
    nClusters(i) = numel(nConstit); 
end
figure; hold all
set(gcf,'color','white')
set(gca,'LineWidth',1.5,'FontSize',12)
grid on; box on
ylabel('# of RNA in cluster','FontSize',20)
xlabel('Simulation progress (%)','FontSize',20)
% set(gca,'XTick',[0:0.01:1])

% xticks(0:0.01:1);
% ax.XAxis.MinorTick = 'on';
% ax.XAxis.MinorTickValues = 0:0.01:1;

%labels = labels./100;
boxplot(boxData, labels,'Colors','k','MedianStyle','line')
bx2 = findobj(gca,'Tag','boxplot');
set(bx2.Children,'LineWidth',1.5)

bar(unique(labels),barY','LineWidth',1.5,'FaceColor',plotColor,'FaceAlpha',0.2,'EdgeAlpha',0.4)
swarmchart(labels,boxData, 20,plotColor,'filled','MarkerEdgeColor','k')

%xticks(0:0.01:1);
ylim([0 80]);  xlim([0 nFrames])

figure; hold all
plot(time,nClusters,'o-','Color',plotColor,'LineWidth',1.5)
set(gcf,'color','white')
set(gca,'LineWidth',2,'FontSize',12)
grid on; box on
ylabel('# of RNA in cluster','FontSize',20)
xlabel('Simulation time ($\mu$s)','FontSize',20,'Interpreter','latex')



%% Quantify cluster exchange 

% Use node distance as more robust way to tell if a certain node is still present or not between frames
for k = 1:numel(adjacencyMatrix)
    distanceMatrix{k} = distances(G{k});
    distanceMatrix{k}(distanceMatrix{k} == Inf) = 0; % If distance=Inf, two nodes are not in the same cluster 
    distanceMatrix{k}(distanceMatrix{k} > 0) = 1;
end

% Count the number of departures and arrivals
%   this part may be confusing, and that's because it is very difficult to avoid overcounting here; had to take an ABSOLUTE ROUNDABOUT approach
%   I've tried to comment what the fuck I'm doing as best as possible but it probably won't make sense unfortunately 
nTotal = numel(distanceMatrix);
for j = 1:nTotal-1
    differenceMatrix{j} = distanceMatrix{j+1} - distanceMatrix{j};

    % The general method: only consider a unique arrival/departure as unique pairs of RNA1,RNA2 in the pairwise difference matrices
    [a1,a2] = find(differenceMatrix{j}==1); 
    arrivalPairs = [a1,a2];
    if ~isempty(arrivalPairs)
        % Initial removal of repeated rows (ie [1 19] and [19 1])
        arrivalPairs = unique(sort(arrivalPairs')','rows','stable');
        % Detect RNA indices that are repeated:
        uniqueVals = unique(arrivalPairs);
        valCount = hist( arrivalPairs , uniqueVals );
        valCount2 = valCount(:,1) + valCount(:,2);
        % Identify repeat arrivers:
        repeatArriversi = find(valCount2>5); % May need to increase this threshold 
        repeatArrivers = uniqueVals(repeatArriversi);
        % Remove repeat arrivers in arrivalPairs matrix (set to NaN):
        for p = 1:numel(repeatArrivers)
            arrivalPairs(arrivalPairs==repeatArrivers(p))=NaN;
        end
        % Count the rows with no NaNs:
        nRestOfArrivals = sum(sum(~isnan(arrivalPairs),2)==2);
        % nArrivals is the number of repeat (non-unique) arrivers + the rest of the non-repeat arrivals
        nArrivals = numel(repeatArrivers) + nRestOfArrivals;
    else
        nArrivals = 0;
    end

    % Repeat the exact same method as above, only now for departures:
    [d1,d2] = find(differenceMatrix{j}==-1);
    departurePairs = [d1,d2];
    if ~isempty(departurePairs)
        departurePairs = unique(sort(departurePairs')','rows','stable');
        uniqueVals = unique(departurePairs);
        valCount = hist( departurePairs , uniqueVals );
        valCount2 = valCount(:,1) + valCount(:,2);
        repeatDepartersi = find(valCount2>3);
        repeatDeparters = uniqueVals(repeatDepartersi);
        for p = 1:numel(repeatDeparters)
            departurePairs(departurePairs==repeatDeparters(p))=NaN;
        end
        nRestOfDepartures = sum(sum(~isnan(departurePairs),2)==2);
        nDepartures = numel(repeatDeparters) + nRestOfDepartures;
    else
        nDepartures = 0;
    end

    nArrivalsTot(j) = nArrivals; nDeparturesTot(j) = nDepartures;

end

% round out these arrays for rebinning 
nArrivalsTot(nTotal) = 0; nDeparturesTot(nTotal) = 0;


% Count the number of departures and arrivals - OLD method that does not avoid overcounting
% for j = 1:numel(distanceMatrix)-1
%     differenceMatrix{j} = distanceMatrix{j+1} - distanceMatrix{j}; 
%     nArrivals = 0; nDepartures = 0;
%     for m = 1:numel(differenceMatrix{j}(:,1))
%         columnToAnalyze = differenceMatrix{j}(m,:);
%         if sum(columnToAnalyze) > 0
%             nArrivals = nArrivals + 1;
%         end
%         if sum(columnToAnalyze) < 0
%             nDepartures = nDepartures + 1;
%         end
%     end
%     nArrivalsTot(j) = nArrivals; nDeparturesTot(j) = nDepartures;
%     %nArrivalsTot(j) = sum(differenceMatrix{j}(:) == 1) /2; % this will not consider uniqueness
%     %nDeparturesTot(j) = sum(differenceMatrix{j}(:) == -1) /2; % this will not consider uniqueness
% end


% Rebin data according to number of frames initially declared
rebinInterval = 10000/numel(time);
nArrivalsReshape = reshape(nArrivalsTot, rebinInterval, nTotal/rebinInterval)';
nArrivalsPlot = sum(nArrivalsReshape, 2);
nArrivalsPlot(nArrivalsPlot==0) = NaN; 
nDeparturesReshape = reshape(nDeparturesTot, rebinInterval, nTotal/rebinInterval)';
nDeparturesPlot = sum(nDeparturesReshape, 2);
nDeparturesPlot(nDeparturesPlot==0) = NaN;

figure; hold all
set(gcf,'color','white')
set(gca,'LineWidth',2,'FontSize',12)
grid on; box on
ylim([0 100])
xlim([0 1])
ylabel('Number of events','FontSize',20)
xlabel('Simulation time ($\mu$s)','FontSize',20,'Interpreter','latex')

%plot(time,nArrivalsPlot,'o-','Color',[29/232 131/232 72/232],'LineWidth',2)
stem(time,nArrivalsPlot,'LineWidth',1,'Color',plotColor,'LineWidth',1.5)
%plot(time,nDeparturesPlot,'ro-','LineWidth',2)
stem(time,nDeparturesPlot,'LineWidth',1,'Color',plotColor2,'LineWidth',1.5)

legend('Arrivals','Departures','FontSize',20) 

disp(['Total number of arrivals: ',num2str(sum(nArrivalsTot))])
disp(['Total number of departures: ',num2str(sum(nDeparturesTot))])


%% Visualize a difference matrix and graphs on either side (this is largely for debugging)
%  toVisualize = triu(differenceMatrix{12},1);
%  figure
%  heatmap(toVisualize)
%  colormap('jet');
%  clim([-1 1])
% % 
% % 
% % figure
% % heatmap(adjacencyMatrix{2})
%  figure
%  h = plot(G{12},'-db','LineWidth',1,'MarkerSize',5);
% % 
% % figure
% % heatmap(adjacencyMatrix{3})
%  figure
%  h = plot(G{13},'-db','LineWidth',1,'MarkerSize',5);
% 


