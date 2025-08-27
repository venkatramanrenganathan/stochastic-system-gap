%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function projectToRiemannSphere will return the coordinates of
% points on the Riemann sphere given a complex number
%
% INPUT:
%      % z: N x 1 complex number vector; 
% OUTPUT:
%      r: 3 x N matrix of 3D points
%
% Copyrights @ 2025
% 
% Authors: Venkatraman Renganathan 
%          Cranfield University, United Kingdom.
%
% Email: v.renganathan@cranfield.ac.uk
%
% Date last updated: 27 August, 2025.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = projectToRiemannSphere(z)
    x = real(z(:)).'; 
    y = imag(z(:)).'; 
    r2 = x.^2 + y.^2;
    X = x ./ (1 + r2);
    Y = y ./ (1 + r2);
    Z = r2 ./ (1 + r2);
    R = [X; Y; Z]; 
end