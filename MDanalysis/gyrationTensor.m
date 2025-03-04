
function [gyr, Rg,b,K2, FA,VecMax,Lmean,L,Evecs,Evals] = gyrationTensor(coords) 
%% Compute gyration tensor and associated shape parameters
%
%   GW - July 2024
%
%   See for example https://www.physik.uni-leipzig.de/~janke/Paper/jcp138_054904_2013.pdf
%

%% Shift coordinates such that the origin is at the centre of mass 
r_cm = mean(coords,1);
r = coords - repmat(r_cm,[numel(coords(:,1)) 1]);
N = numel(coords(:,1)); % # of particles


%% Compute gyration tensor, component-wise 
Gxx = (1/N).*sum(r(:,1).*r(:,1));
Gxy = (1/N).*sum(r(:,1).*r(:,2));
Gxz = (1/N).*sum(r(:,1).*r(:,3));

Gyx = (1/N).*sum(r(:,2).*r(:,1));
Gyy = (1/N).*sum(r(:,2).*r(:,2));
Gyz = (1/N).*sum(r(:,2).*r(:,3));

Gzx = (1/N).*sum(r(:,3).*r(:,1));
Gzy = (1/N).*sum(r(:,3).*r(:,2));
Gzz = (1/N).*sum(r(:,3).*r(:,3));

gyr = [Gxx,Gxy,Gxz ; Gyx,Gyy,Gyz ; Gzx,Gzy,Gzz]; 


%% Shape parameters derived from eigendecomposition of G
[V, L] = eig(gyr);
LorderPreserve = diag(L); % save a snapshot of eigenvalues before sorting for correct delineation of max eigenvector
L = diag(L); 
L = sort(L,'descend');


% Eigenvalue-derived physical metrics: 
Rg = sqrt(L(1)+L(2)+L(3)); %radius of gyration 
b = L(1) - (1/2)*(L(2)+L(3)); %asphericity
%K2 = (3/2) * ((L(1)^2 + L(2)^2 + L(3)^2) / (L(1) + L(2) + L(3)).^2) - (1/2); %shape anisotropy
K2 = 1 - 3*((L(1)*L(2) + L(2)*L(3) + L(3)*L(1)) ./ (L(1) + L(2) + L(3)).^2); %shape anisotropy (written in more intuitive form)

Lmean = (1/3).*(L(1)+L(2)+L(3)); %mean mass distribution 
FA = sqrt((3/2).* ((L(1)-Lmean).^2 + (L(2)-Lmean).^2 + (L(3)-Lmean).^2) ./ (L(1).^2+L(2).^2+L(3).^2) ); %fractional anisotropy (borrowed from diffusion tensor imaging)
VecMax = V(:, LorderPreserve == max(LorderPreserve)); %direction of principle eigenvector

% Eigenvectors and eigenvalues for outputting
Evecs = V;
Evals = LorderPreserve; 

