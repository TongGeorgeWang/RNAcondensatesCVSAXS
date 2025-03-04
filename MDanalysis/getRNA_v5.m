function getRNA_v5(dataBasename,analSubfolder)
%% Identify RNA molecules and their clusters in a given simulation frame 
%
%   Stripped coarse grain PDB file format:  
%       ATOM      1   41         1     784.238  81.441  11.110  
%
%   In order, the information for these columns is:
%       Atom type ('ATOM'); Serial# (1+i*1); Species (eg N, O, or for CG: UNK for unknown, or LAMMPS atom ID number);
%       residue seq ID; X-coord; Y-coord; Z-coord
%
%   These data columns are the minimum required to represent the data and
%   be used by programs that read pdb files.
%
%
%   This script works by isolating individual chains, then IDing clusters,
%   then performing periodic boundary unification, then isolating RNA
%   chains. 
%
%   Things that are saved for each processed frame: 
%       - graph for cluster identification, along with spectral clusters 
%       - non PBC corrected atomic coordinates
%       - graph for PBC correction, along with spectral clusters 
%       - PBD corrected atomic coordinates
%
%   
%   GW June 2024
%       v1: shitty version for getting things sorted 
%       v2: updated method of PBC correction (made different versions in case something got screwed up)
%       v3: added parallel computing functionality to run faster (involved changing the way variables are indexed)
%       v4: clustering and PBC correction done in Ovito instead of here; start by reading that output
%         - time testing: 10000 frames took 7-8 hours to run
%       v5: culled out the things we no longer want to do (after deciding on the things that we want to include in the final paper) - October 2024        
%         - time testing: 10000 frames took 0.5 hours to run 

%clear; close all
%% input parameters 

frames = 1:1:10000; % what MD frames you want to run this script over
%frames = 1:10:3000;

RNAnumber = [41 42 43 44]; % atom ID# of RNA (this is derived from the Mpipi atom ID#)
bondLength = 5; % RNA bond length used in the simulation (Anstroms)
%distanceThreshold = bondLength*3; % what distance to consider a cutoff threshold in graph formation (may need to refine this)
distanceThreshold = bondLength*1.5;

dataBasename = 'inputFiles/rU30alone_250uM_150Na_Ovito/rU30alone_.';
analSubfolder = 'rU30PR30_250uM_150Na'; % Subfolder in 'Analysis' where the results will be saved. 


%% Preallocate cells and matrices that are amenable to preallocation
nFrames = numel(frames);
indRNAcoords_allFrames = cell([nFrames 1]); clustersToSave = cell([nFrames 1]); 
G = cell([nFrames 1]); r_all = cell([nFrames 1]); pr_all = cell([nFrames 1]); gr_all = cell([nFrames 1]); I_all = cell([nFrames 1]);
P_RNA_clusters_Coords = cell([nFrames 1]); clusterIndices = cell([nFrames 1]); adjacencyMatrix = cell([nFrames 1]);
gyr = cell([nFrames 1]); Rg = zeros([nFrames 1]); b = zeros([nFrames 1]); K2 = zeros([nFrames 1]); FA = zeros([nFrames 1]); VecMax = cell([nFrames 1]); Lmean = zeros([nFrames 1]); L = cell([nFrames 1]);


%%
for iter = 1:nFrames
    %% Get RNA coordinates
    filename = [dataBasename,num2str(frames(iter))];

    raw = readmatrix(filename,'FileType','text');

    % Isolate RNA and remove amino acids
    raw(~ismember(raw(:,2),RNAnumber),:) = [];

    % Split master matrix into a cell array P where P{i} contains the {x,y,z} coords of the i'th polymer in the frame
    [~,~,Pindices] = unique(raw(:,3));
    P_RNA = accumarray(Pindices,1:size(raw,1),[],@(r){raw(r,4:6)});
    indRNAcoords_allFrames{iter} = P_RNA; % individual RNA coordinates for each frame, for saving

    %% Determine which polymers are in the same cluster, and which are connected, using graph theory

    % Compute difference in space between all polymer pairs
    distances = zeros([numel(P_RNA) numel(P_RNA)]); nContacts = zeros([numel(P_RNA) numel(P_RNA)]);
    numReps = numel(P_RNA); % # of chains within the frame

    for i = 1:numReps
        for j = 1:numReps
            if i >= j % (*) only need to compute for one half of symmetric matrix; populate other half with mirror symmetry in later line
                P_diff = pdist2(P_RNA{i},P_RNA{j});
                %P_diff(P_diff==0)=inf; minDistance(i,j) = min(min(P_diff)); % check minimum distances, for debugging purposes
                if any(any(P_diff < distanceThreshold)) % if any monomers are close together between a pair of polymers,
                    distances(i,j) = 1; % consider them to be in the same cluster
                    nContacts(i,j) = sum(P_diff < distanceThreshold, 'all'); % record the number of contacts present; divide by 2 due to matrix symmetry
                    
                else
                    distances(i,j) = 0; % else, not in the same cluster
                end
            else
                % piss off to next iteration (*)
            end
        end
    end
    adjacencyMatrix{iter} = ceil((distances + distances.')/2); % (*) populate the other half of distances, which is a symmetric binary matrix

    % Spectrally cluster graph:
    C = SpecClust_v2(adjacencyMatrix{iter}, numReps);
    clusters = [1:numel(C); C']';
    clusters = sortrows(clusters,2);
    clustersToSave{iter} = clusters; % For saving in final lines
    nClusters = max(clusters(:,2))+1;

    G{iter} = graph(adjacencyMatrix{iter},'omitselfloops','upper');


    %% Plot graph of clusters (for visualization)
    % figure('Name','Clustered Graph'); hold all
    % h = plot(G{iter},'-db','LineWidth',1,'MarkerSize',5);
    % set(gcf,'color','w')
    % set(gca,'FontSize',20)
    % set(gca,'LineWidth',2)
    % box on
    % grid off
    % h.EdgeColor = [0.5 0.5 0.5];
    % h.Marker = 'O';
    % 
    % %Color graph nodes by cluster
    % colors = colormap(jet);
    % colorSpacing = floor(numel(colors(:,1)) / (nClusters+1) * 0.9);
    % for j = 1:nClusters
    %     clusterIndicess = clusters(clusters(:,2)==j-1);
    %     highlight(h,clusterIndicess,'NodeColor',colors(j*colorSpacing,:))
    % end


    %% Plot coordinates (for visualization)
    % figure
    % hold all
    % xlabel('X'); ylabel('Y'); zlabel('Z')

    for n = 1:nClusters

        clusterIndices{iter,n} = clusters(clusters(:,2)==n-1);
        P_RNA_clusters = P_RNA(clusterIndices{iter,n});

        thisClusterCoords = cell2mat(P_RNA_clusters); % unfurl RNA coordinates        
        P_RNA_clusters_Coords{iter,n} = thisClusterCoords;  
        
        %plot3(thisClusterCoords(:,1),thisClusterCoords(:,2),thisClusterCoords(:,3),'-o')

        %[r{iter,n},pr{iter,n},Dmax{iter,n}, gr{iter,n}, gyr{iter,n},Rg(iter,n),b(iter,n),K2(iter,n),FA(iter,n),VecMax{iter,n},Lmean(iter,n),L{iter,n}] = lengthCorrelations(P_RNA_clusters); % perform on coordinates within single clusters
        [gyr{iter,n},Rg(iter,n),b(iter,n),K2(iter,n),FA(iter,n),VecMax{iter,n},Lmean(iter,n),L{iter,n}] = lengthCorrelations(P_RNA_clusters); % perform on coordinates within single clusters
                                                            
    end

end


%% Save output results en masse (all frames in a single, celled variable within a .mat file, within 'Analysis' subfolder)
%   This is so downstream visualization and further analysis may be modularized to avoid one long chain of misfortune

if ~exist('Analysis','dir')
    mkdir('Analysis') 
end
if ~exist(['Analysis/',analSubfolder],'dir')
    mkdir(['Analysis/',analSubfolder]) 
end

% P(r) - pair distance distribution functions | G(r) - radial distribution functions
%save(['Analysis/',analSubfolder,'/structuralParameters'],'r','pr','Dmax','gr','gyr','Rg','b','K2','FA','VecMax','Lmean','L')
save(['Analysis/',analSubfolder,'/structuralParameters'],'gyr','Rg','b','K2','FA','VecMax','Lmean','L')

% G - Graphs for clustering, with accompanying spectral clusters | G2 - Graphs for subclustering (PBC correction), with accompanying spectral clusters | Uncorrected coordinates | Corrected coordinates
save(['Analysis/',analSubfolder,'/coordsAndClusters'],'G','clustersToSave','indRNAcoords_allFrames','clusterIndices','adjacencyMatrix','P_RNA_clusters_Coords','frames')


end