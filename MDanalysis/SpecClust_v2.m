function C = SpecClust_v2(AdjacencyMatrix, nTotal)

warning('off','all') 

%% Perform spectral clustering on a graph, whose information is contained in an adjacency matrix 
%   AdjacencyMatrix - Input adjacency matrix (in this case of binary macromolecular pairings, needs to be square
%   nKmeans - Number of K-means clusters to compute
%   nClustersMax - if nKmeans=0, program will determine optimal nKmeans, using nClustersMax as a search maxima
%
%   References: 
%       - Ng, Jordan, Weiss: On spectral clustering and an algorighm, NIPS 2001
%       - von Luxburg: A tutorial on spectral clustering, Statistics and Computing 2007
%       some code implementations adapted from Ingo Burk, 'areslp', Github, on Academic License
% 
%   v2 (Aug 2024): used alternate way of estimating the number of clusters, that doesn't involve an iterative sweep
%

%% Compute Laplacian matrix and normalize
degrees = sum(AdjacencyMatrix, 2);
D = sparse(1:size(AdjacencyMatrix, 1), 1:size(AdjacencyMatrix,2), degrees);
L = D - AdjacencyMatrix; 
D = spdiags(1./(degrees.^0.5), 0, size(D, 1), size(D, 2)); 
L = D * L * D;
%degrees(degrees == 0) = eps; 
[V, eigenvalues] = eigs(L);

%% Compute nKmeans smallest eigenvectors and perform K-means

diff  = eps;
evs = eigs(L,nTotal,diff);
nClusters = numel(find(evs<1E-10 == 1)); % # of K-means clusters is estimated by computing evs (eigenvalues) of L and counting # of evs=0 (or functionally 0, very small)
  
[V, eigenvalues] = eigs(L, nClusters, diff);
%eigenvalues = diag(eigenvalues);
V = bsxfun(@rdivide, V, sqrt(sum(V.^2, 2))); %normalize eigenvectors

% C: n-by-1 matrix containing the cluster number for each data point, minus 1 for the consistant with the true label
C = kmeans(V, nClusters) - 1;
%C = sparse(1:size(D, 1), C, 1); % now convert C to a n-by-k matrix containing the k indicator vectors as columns




end
