%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function inverseStereographicProjection will return the coordinates of
% points on the Riemann sphere given a complex number
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
function r = inverseStereographicProjection(z)

r = [real(z) ./ (1 + abs(z).^2); ...
     imag(z) ./ (1 + abs(z).^2); ...
     abs(z).^2 ./ (1 + abs(z).^2)];

end