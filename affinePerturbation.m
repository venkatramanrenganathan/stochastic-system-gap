%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function affinePerturbation will return stable perturbed state space data
% given the nominal state space data, random directions data and theta
%
% INPUT:
%      A0,B0,C0: nominal state space data
%      Aj,Bj,Cj: random directions data
%      theta   : theta parameter
%
% OUTPUT:
%      A,B,C: perturbed state space data
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
function [A,B,C] = affinePerturbation(A0,B0,C0, Aj,Bj,Cj, theta)
    % Store the nominal values
    A = A0; 
    B = B0; 
    C = C0;
    % Add the perturbation to the nominal
    for r = 1:numel(theta)
        A = A + theta(r)*Aj(:,:,r);
        B = B + theta(r)*Bj(:,:,r);
        C = C + theta(r)*Cj(:,:,r);
    end
    % Check & ensure that the matrix is Hurwitz
    maxRe = max(real(eig(A)));
    if(maxRe >= -0.2)
        % Shift any right-half-plane eigenvalues to ensure stability.
        A = A - (maxRe + 0.3)*eye(size(A));
    end
end