function [T_ratio, P_ratio, rho_ratio, A_ratio] = isentropic_flow(M, gamma)
% ISENTROPIC_FLOW Compute isentropic flow ratios for a given Mach number.
%
%   [T_ratio, P_ratio, rho_ratio, A_ratio] = isentropic_flow(M, gamma)
%
%   Inputs:
%       M     - Mach number (scalar or array)
%       gamma - ratio of specific heats (default 1.4)
%
%   Outputs:
%       T_ratio   - T/T0   (static-to-stagnation temperature ratio)
%       P_ratio   - P/P0   (static-to-stagnation pressure ratio)
%       rho_ratio - rho/rho0 (static-to-stagnation density ratio)
%       A_ratio   - A/A*   (area ratio relative to sonic throat)

    if nargin < 2
        gamma = 1.4;
    end

    % Isentropic relations (Anderson, Modern Compressible Flow, Ch. 3)
    factor = 1 + (gamma - 1) / 2 .* M.^2;

    T_ratio   = 1 ./ factor;
    P_ratio   = factor .^ (-gamma / (gamma - 1));
    rho_ratio = factor .^ (-1 / (gamma - 1));

    % Area-Mach relation: A/A* = (1/M) * [(2/(gamma+1)) * (1 + (gamma-1)/2 * M^2)]^((gamma+1)/(2*(gamma-1)))
    exponent = (gamma + 1) / (2 * (gamma - 1));
    A_ratio  = (1 ./ M) .* ((2 / (gamma + 1)) .* factor) .^ exponent;
end
