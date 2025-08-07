%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function computeGeodesicDistance will return the geodesic distance 
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
function geoDistance = computeGeodesicDistance(r1, r2)

geoDistance = 0.5 * real(acos(dot(r1, r2) ./ (vecnorm(r1) .* vecnorm(r2))));

end