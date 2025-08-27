%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function generateRandomDirectionData will return random direction
% matrices for perturbing the nominal A, B, C state space data
%
% INPUT:
%      n,m,l,p: dimensions of state, input, output, theta
%      scale  : scaling factor to ensure perturbation is small
%      A0,B0,C0: nominal state space data
%      Aj,Bj,Cj: random directions data
%      theta   : theta parameter
%
% OUTPUT:
%      Aj, Bj, Cj: random direction matrices
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
function [Aj,Bj,Cj] = generateRandomDirectionData(n,m,l,p,scale)
    Aj = zeros(n,n,p);
    Bj = zeros(n,m,p);
    Cj = zeros(l,n,p);
    for r = 1:p
        Aj(:,:,r) = scale*(randn(n)*0.5);
        Bj(:,:,r) = scale*(randn(n,m)*0.5);
        Cj(:,:,r) = scale*(randn(l,n)*0.5);
    end
end