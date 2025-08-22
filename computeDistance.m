%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function computeDistance will return the geodesic/chordal distance 
% between two points on the Riemann sphere
%
% INPUT:
%      z: Complex number (real + j imag) 
% OUTPUT:
%      r: 3D coordinate of point in Riemann sphere
%
% Copyrights @ 2025
% 
% Authors: Venkatraman Renganathan 
%          Cranfield University, United Kingdom.
%
% Email: v.renganathan@cranfield.ac.uk
%
% Date last updated: 07 August, 2025.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function distance = computeDistance(r1, r2)

% Set distance choice: 1-chordal, 2-geodesic
distanceChoice = 1;

if(distanceChoice == 1)
    % chordal distance
    distance = norm(r1 - r2);
else
    % geodesic distance
    sphereRadius = [0;0;0.5];
    distance = 0.5 * real(acos(4*dot(r1-sphereRadius, r2-sphereRadius)));
    %geoDistance = 0.5 * real(acos(dot(r1, r2) ./ (vecnorm(r1) .* vecnorm(r2))));
end

end