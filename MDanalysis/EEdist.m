
%% Compute end to end distance 
%
%   GW - 2024
%

function [Ree] = EEdist(coords)
    
    P1 = coords(1,:);
    Plast = coords(end,:);
    
    Ree = pdist([P1;Plast],'euclidean');

end



