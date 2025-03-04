
%% Ensemble analysis on polymers in CG-RNA simulation 
%   Computes the prescribed quantities for each individual RNA strand in the input cell array
%   After computing the quantities, assigns RNAs based on their cluster indices 
%
%   GW - July 2024
%
%   Running this on 10000 frames took ~500-800seconds.

%clear; close all
%load Analysis\dev\coordsAndClusters.mat
frames = 1:1:10000;                                                                  
analSubfolder = 'rU30PR30_250uM_150Na'; %name of analysis subfolder where results will be saved


%% Preallocate matrices and cells 
nFrames = numel(frames); 
nRNAs = numel(indRNAcoords_allFrames{1,1}); % Note: if for some reason # of RNA changes between frames, this preallocation will not work 
gyr_ind = cell([nFrames,nRNAs]); Rg_ind = zeros([nFrames,nRNAs]); b_ind = zeros([nFrames,nRNAs]); K2_ind = zeros([nFrames,nRNAs]);  
OCF_ind = cell([nFrames,nRNAs]); l_OCF_ind = zeros([nFrames,nRNAs]); l_OCF_err_ind = zeros([nFrames,nRNAs]); T1_ind = zeros([nFrames,nRNAs]); T2_ind = zeros([nFrames,nRNAs]);
Ree_ind = zeros([nFrames,nRNAs]);
gyr_ind = cell([nFrames 1]); Rg_ind = zeros([nFrames 1]); b_ind = zeros([nFrames 1]); K2_ind = zeros([nFrames 1]); FA_ind = zeros([nFrames 1]); VecMax_ind = cell([nFrames 1]); Lmean_ind = zeros([nFrames 1]); L_ind = cell([nFrames 1]);


%% Compute metrics  
for iter = 1:nFrames
    for nRNA = 1:numel(indRNAcoords_allFrames{frames(iter),1})

        indRNAcoords = indRNAcoords_allFrames{frames(iter),1}{nRNA}; % get coordinates of current RNA

        %% Rgs of individual chains
        [gyr_ind{iter,nRNA},Rg_ind(iter,nRNA),b_ind(iter,nRNA),K2_ind(iter,nRNA),FA_ind(iter,nRNA),VecMax_ind{iter,nRNA},Lmean_ind(iter,nRNA),L_ind{iter,nRNA}] = gyrationTensor(indRNAcoords);

        %% OCF and tortuosity
        [OCF_ind{iter,nRNA},l_OCF_ind(iter,nRNA),l_OCF_err_ind(iter,nRNA), T1_ind(iter,nRNA),T2_ind(iter,nRNA)] = tortuosityAndOCF(indRNAcoords);

        %% End-end distances
        [Ree_ind(iter,nRNA)] = EEdist(indRNAcoords);

        
    end
end



%% Save results when satisfied 
save(['Analysis/',analSubfolder,'/ensembleAnalysisResults'],'gyr_ind','Rg_ind','b_ind','K2_ind','FA_ind','VecMax_ind','Lmean_ind','L_ind','OCF_ind','l_OCF_ind','l_OCF_err_ind','T1_ind','T2_ind','Ree_ind')





