function [M2, P2_P1, T2_T1, rho2_rho1, P02_P01] = normal_shock(M1, gamma)
% NORMAL_SHOCK Compute flow properties across a normal shock.
%
%   [M2, P2_P1, T2_T1, rho2_rho1, P02_P01] = normal_shock(M1, gamma)
%
%   Inputs:
%       M1    - upstream Mach number (must be >= 1)
%       gamma - ratio of specific heats (default 1.4)
%
%   Outputs:
%       M2       - downstream Mach number
%       P2_P1    - static pressure ratio across shock
%       T2_T1    - static temperature ratio across shock
%       rho2_rho1 - density ratio across shock
%       P02_P01  - stagnation pressure ratio (< 1, entropy increase)

    if nargin < 2
        gamma = 1.4;
    end

    g = gamma;

    % Normal shock relations (Anderson, Modern Compressible Flow, Ch. 3)
    M2 = sqrt((1 + (g-1)/2 .* M1.^2) ./ (g .* M1.^2 - (g-1)/2));

    P2_P1 = 1 + 2*g/(g+1) .* (M1.^2 - 1);

    T2_T1 = P2_P1 .* (2 + (g-1) .* M1.^2) ./ ((g+1) .* M1.^2);

    rho2_rho1 = (g+1) .* M1.^2 ./ (2 + (g-1) .* M1.^2);

    % Stagnation pressure ratio (Rayleigh pitot formula)
    P02_P01 = ((((g+1).*M1.^2) ./ ((g-1).*M1.^2 + 2)).^(g/(g-1))) ...
              .* ((2*g.*M1.^2 - (g-1)) ./ (g+1)).^(-1/(g-1));
end
