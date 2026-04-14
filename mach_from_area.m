function M = mach_from_area(A_ratio, gamma, branch)
% MACH_FROM_AREA Solve the area-Mach relation for Mach number.
%
%   M = mach_from_area(A_ratio, gamma, branch)
%
%   Inputs:
%       A_ratio - A/A* (must be >= 1)
%       gamma   - ratio of specific heats (default 1.4)
%       branch  - 'sub' for subsonic root, 'sup' for supersonic root
%
%   Output:
%       M - Mach number (scalar or array, same size as A_ratio)
%
%   Uses fzero with appropriate initial guesses for each branch.

    if nargin < 2 || isempty(gamma)
        gamma = 1.4;
    end
    if nargin < 3 || isempty(branch)
        branch = 'sup';
    end

    M = zeros(size(A_ratio));

    for i = 1:numel(A_ratio)
        AR = A_ratio(i);

        if AR < 1 - 1e-10
            error('A/A* = %.4f is less than 1. Not physical.', AR);
        end

        if abs(AR - 1) < 1e-10
            M(i) = 1.0;
            continue;
        end

        % Define residual: f(M) = A/A*(M) - target A/A*
        exponent = (gamma + 1) / (2 * (gamma - 1));
        f = @(m) (1/m) * ((2/(gamma+1)) * (1 + (gamma-1)/2 * m^2))^exponent - AR;

        if strcmpi(branch, 'sub')
            % Subsonic root: M in (0, 1)
            M(i) = fzero(f, [1e-6, 1 - 1e-10]);
        else
            % Supersonic root: M in (1, ~large)
            % Upper bound estimate: for large M, A/A* ~ M, so M_max ~ 2*AR
            M_upper = max(2, 2 * AR);
            M(i) = fzero(f, [1 + 1e-10, M_upper]);
        end
    end
end
