
function [OCF,l_OCF,l_OCF_err, T1,T2] = tortuosityAndOCF(coords)
%% Compute tortuosity of PAR sequence 
%
%   OCF: 
%   The code for computing OCF was written by Alex Plumridge 
%   Also added computation of correlation length 
%
%
%   Tortuosity method 1: sample coordinates of P atoms (PA & PB) and O atoms in between the two
%   bases (O1D) as 'points'
%       Take pairwise lengths of adjacent points
%       C = total arc length = total length of pairwise vectors
%       L = end-to-end distance = length of 1st & last coordinates
%
%       Tortuosity (T) = C/L
%
%   Tortuosity method 2: sample coordinates of P atoms (PA & PB) and O atoms in between the two
%   bases (O1D)
%       Take pairwise lengths of adjacent points
%       SA = sum of all angles between pairwise vectors
%       L = end-to-end distance 
%       
%       Tortuosity (T) = SA/L
%
%
%   Derived from analysis in Nat Comm PAR paper, but the input is a coordinate set from 
%       CG models, so no atom sampling needs to be done
%
%   GW - 2022
%


%% Set up vectors

S_1 = [coords; 0 0 0];
S_2 = [0 0 0; coords];
v = S_2 - S_1;
v = v(2:end-1,:); % truncate first and last points

% Normalize 
vl = sqrt(sum(v.^2,2));
vn = v./repmat(vl,[1,size(v,2),1]);
nv = zeros(size(vn,1)-1,size(vn,3));


%% Compute OCF
% Compute P-O displacements (b value in correlation length computation)
bVec = sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);
bMean = mean(bVec);

% Compute average of cos(theta) vs separation, for each chain
for j = 1:(size(vn,1)-1)
    nv(j,:) = squeeze(mean((vn(1:(end-j),1,:).*vn((j+1):end,1,:) + ...
        vn(1:(end-j),2,:).*vn((j+1):end,2,:) + ...
        vn(1:(end-j),3,:).*vn((j+1):end,3,:)),1)); 
end

OCF = nv(:);


%% Compute correlation length

OCF_ensembleMean = mean(OCF,1); % arithmetic mean
OCF_ensembleVar = var(OCF,0,1);
b = mean(bMean); % arithmetic mean
b_err = var(bMean);
l_OCF = b.*sum(OCF_ensembleMean);
l_OCF_err = l_OCF.*sqrt( (b_err./b).^2 + (mean(OCF_ensembleVar)./l_OCF).^2 ); % propagate error multiplicatively 


%% Compute tortuosity index (T) via method 1

magnitudes = sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);
C = sum(magnitudes);
vEE = coords(end,:)-coords(1,:); % end-to-end vector
L = sqrt(vEE(1).^2 + vEE(2).^2 +vEE(3).^2);

T1 = C./L ;


%% Compute tortuosity index (T) via method 2

v1 = v(1:end-1,:);
v2 = v(2:end,:);
angles = acos(dot(v1,v2,2)./...
    (sqrt(v1(:,1).^2 + v1(:,2).^2 +v1(:,3).^2) .* sqrt(v2(:,1).^2 + v2(:,2).^2 +v2(:,3).^2) )); % Based on definition of dot product
%angles = atan2(norm(cross(v1,v2)),dot(v1,v2)); % More robust calculation at small angles
anglesDegrees = angles .* (180./pi);
SA = sum(anglesDegrees);

T2 = SA./L ;


