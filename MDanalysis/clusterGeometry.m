
%% Quantify geometric properties of condensates (clusters) across MD simulation frames
%
%   Designed for ssRNA-peptide Mpipi simulations, focusing on RNA condensation
%   Methods:
%       Convex hull projection: a way to visualize shape changes/deformations
%
%       2D mesh projections / cross sections: a way to simplify and visualize mesh structures
%
%   Some functions from this package are used: 
%       Matt J (2021). Analyze N-dimensional Convex Polyhedra 
%       (https://www.mathworks.com/matlabcentral/fileexchange/30892-analyze-n-dimensional-convex-polyhedra)
%       MATLAB Central File Exchange.
%
%   GW - October 2024
%

close all

%%
%color = [46 80 122]./255; %rU blue
color = [146 0 0]./255; %rA red
%color = [197 192 0]./255; %rC deep yellow 

% For similarity calculations:
%frames = 5000:1:6000;

% For visualizations: 
%frames = 1000:10:1100;
%frames = 5000:10:5100;
frames = 9000:10:9100;

convexSum = 0; % if 1: sum over all frames and plot overall convex hull | 0 for plotting all frames together
% I've found 1 doesn't work very well 

%% Unravel RNA coordinates
% thisFrame = indRNAcoords_allFrames{frame};
% RNA = [];
% for i = 1:numel(thisFrame)
%     RNA = [RNA; thisFrame{i}];
% end

nFrames = numel(frames);

figure; hold all; set(gcf,'color','white')
set(gca,'XColor', 'none','YColor','none','ZColor','none'); set(gca, 'color', 'none')
set(gcf, 'Position',  [100, 100, 400, 400])

%% Get RNA coordinates of largest cluster in frame

if convexSum == 0


    %% Compute difference index 

 % Commented out lines here are shit that didn't work 
    % for N = 1:1:(nFrames-1)
    %     D(N) = sum(K{N+1} - K{N});
    % end

    % k1 = K{1}(:,1); k2 = K{1}(:,2); k3 = K{1}(:,3);
    % figure
    % T=delaunay(k1,k2);
    % trisurf(T,k1,k2,k3);
    % xy = [k1,k2];
    % a = xy(T(:,2),:)-xy(T(:,1),:);
    % b = xy(T(:,3),:)-xy(T(:,1),:);
    % V = ((a(:,1).*b(:,2)-a(:,2).*b(:,1))' * sum(k3(T),2))/6;
    % jaccardIndex = height(intersect(K{1},K{2},'rows')) / height( union(K{1},K{2},'rows') );

    %figure; hold all
    i = frames(1);

    thisFrame = P_RNA_clusters_Coords(i,:);
    [max_size, clusteridx] = max(cellfun('size', thisFrame, 1)); % determine largest cluster
    RNA = P_RNA_clusters_Coords{i,clusteridx};
    RNA_cm = [mean(RNA(:,1)),mean(RNA(:,2)),mean(RNA(:,3))];
    RNA_centred = RNA - RNA_cm;
    x = RNA_centred(:,1); y = RNA_centred(:,2); z = RNA_centred(:,3);
    V{1} = [x,y,z];
    [C{1},vol(1)] = convhull(x,y,z);


    for N = 2:1:nFrames
        i = frames(N);

        thisFrame = P_RNA_clusters_Coords(i,:);
        [max_size, clusteridx] = max(cellfun('size', thisFrame, 1)); % determine largest cluster
        RNA = P_RNA_clusters_Coords{i,clusteridx};


        %% Equalize centres of mass by shifting to origin
        RNA_cm = [mean(RNA(:,1)),mean(RNA(:,2)),mean(RNA(:,3))];
        RNA_centred = RNA - RNA_cm;


        %% Compute convex hull of RNA condensate coordinates
        x = RNA_centred(:,1); y = RNA_centred(:,2); z = RNA_centred(:,3);
        V{N} = [x,y,z];
        [C{N}, vol(N)] = convhull(x,y,z);

        trisurf(C{N},x,y,z,'FaceColor',color,'FaceAlpha',0.05,'EdgeAlpha',0.2)

        %% Compute 3D Jaccard index as a metric of similarity 
        IntHull = intersectionHull('vert',V{N-1},'vert',V{N});
        xInt = IntHull.vert(:,1); yInt = IntHull.vert(:,2); zInt = IntHull.vert(:,3);
        UnionHull = unionHull('vert',V{N-1},'vert',V{N});
        xUn = UnionHull.vert(:,1); yUn = UnionHull.vert(:,2); zUn = UnionHull.vert(:,3);
        [CInt, VInt] = convhull(xInt,yInt,zInt);
        [CUn, VUn] = convhull(xUn,yUn,zUn);
        J(N-1) = VInt./VUn; % Jaccard Index
        D(N-1) = 2*VInt./(vol(N-1)+vol(N)); % Dice-Sorensen Coefficient
    end
    
    %% Visualize the intersection and union in an example 
    IntHull = intersectionHull('vert',V{1},'vert',V{2});
    xInt = IntHull.vert(:,1); yInt = IntHull.vert(:,2); zInt = IntHull.vert(:,3);
    UnionHull = unionHull('vert',V{1},'vert',V{2});
    xUn = UnionHull.vert(:,1); yUn = UnionHull.vert(:,2); zUn = UnionHull.vert(:,3);

    [CInt, VInt] = convhull(xInt,yInt,zInt);
    plot3(xInt,yInt,zInt,'ro')
    trisurf(CInt,xInt,yInt,zInt,'FaceColor','r','FaceAlpha',0.05,'EdgeAlpha',0.2)

    [CUn, VUn] = convhull(xUn,yUn,zUn);
    plot3(xUn,yUn,zUn,'bo')
    trisurf(CUn,xUn,yUn,zUn,'FaceColor','b','FaceAlpha',0.05,'EdgeAlpha',0.2)

    J = VInt./VUn;
    
     
    

%% Visualize dissimilarities by plotting convex hulls of ~10 frames spaced ~10 frames apart (for clarity) 
    for N = 1:1:nFrames
        
        i = frames(N);

        thisFrame = P_RNA_clusters_Coords(i,:);
        [max_size, clusteridx] = max(cellfun('size', thisFrame, 1)); % determine largest cluster
        RNA = P_RNA_clusters_Coords{i,clusteridx}; 


        %% Equalize centres of mass by shifting to origin
        RNA_cm = [mean(RNA(:,1)),mean(RNA(:,2)),mean(RNA(:,3))];
        RNA_centred = RNA - RNA_cm;


        %% Compute convex hull of RNA condensate coordinates
        x = RNA_centred(:,1); y = RNA_centred(:,2); z = RNA_centred(:,3);
        K{N} = convhull(x,y,z);


        %% Plot 2D convex hull projections - simply view plot from each axis head on

        
        trisurf(K{N},x,y,z,'FaceColor',color,'FaceAlpha',0.05,'EdgeAlpha',0.2)
        axis equal
        xlim([-180 180]); ylim([-180 180]); zlim([-180 180])
        xlabel('X');ylabel('Y');zlabel('Z');
        view(45,45) % isometric

        trisurf(K{N},x,y,z,'FaceColor',color,'FaceAlpha',0.05,'EdgeAlpha',0.2)
        xlim([-180 180]); ylim([-180 180]); zlim([-180 180])
        xlabel('X');ylabel('Y');zlabel('Z');
        view(90,0) % yz projection

        subplot(2,2,3); hold a ll
        trisurf(K{N},x,y,z,'FaceColor',color,'FaceAlpha',0.05,'EdgeAlpha',0.2)
        xlim([-180 180]); ylim([-180 180]); zlim([-180 180])
        xlabel('X');ylabel('Y');zlabel('Z');
        view(0,90) % xy projection

        subplot(2,2,4); hold all
        trisurf(K{N},x,y,z,'FaceColor',color,'FaceAlpha',0.05,'EdgeAlpha',0.2)
        xlim([-180 180]); ylim([-180 180]); zlim([-180 180])
        xlabel('X');ylabel('Y');zlabel('Z');
        view(0,0) % xz projection



    end


%%
elseif convexSum == 1

    RNA_centred_tot = [0 0 0];
    for i = frames

        thisFrame = P_RNA_clusters_Coords(i,:);
        [max_size, clusteridx] = max(cellfun('size', thisFrame, 1)); % determine largest cluster
        RNA = P_RNA_clusters_Coords{i,clusteridx};


        %% Equalize centres of mass by shifting to origin
        RNA_cm = [mean(RNA(:,1)),mean(RNA(:,2)),mean(RNA(:,3))];
        RNA_centred = RNA - RNA_cm;

        RNA_centred_tot = [RNA_centred_tot;RNA_centred];

    end

    %% Compute convex hull of RNA condensate coordinates
    x = RNA_centred_tot(:,1); y = RNA_centred_tot(:,2); z = RNA_centred_tot(:,3);
    K = convhull(x,y,z);

    %% Plot 2D convex hull projections - simply view plot from each axis head on

    trisurf(K,x,y,z,'FaceColor',color,'FaceAlpha',0.05,'EdgeAlpha',0.2)
    axis equal
    xlim([-180 180]); ylim([-180 180]); zlim([-180 180])
    xlabel('X');ylabel('Y');zlabel('Z');
    view(45,45) % isometric

    trisurf(K,x,y,z,'FaceColor',color,'FaceAlpha',0.05,'EdgeAlpha',0.2)
    xlim([-180 180]); ylim([-180 180]); zlim([-180 180])
    xlabel('X');ylabel('Y');zlabel('Z');
    view(90,0) % yz projection

    subplot(2,2,3); hold all
    trisurf(K,x,y,z,'FaceColor',color,'FaceAlpha',0.05,'EdgeAlpha',0.2)
    xlim([-180 180]); ylim([-180 180]); zlim([-180 180])
    xlabel('X');ylabel('Y');zlabel('Z');
    view(0,90) % xy projection

    subplot(2,2,4); hold all
    trisurf(K,x,y,z,'FaceColor',color,'FaceAlpha',0.05,'EdgeAlpha',0.2)
    xlim([-180 180]); ylim([-180 180]); zlim([-180 180])
    xlabel('X');ylabel('Y');zlabel('Z');
    view(0,0) % xz projection

end
