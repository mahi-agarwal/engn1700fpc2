function [x, r, A] = nozzle_geometry(geom_type, params, r_t, N)
% NOZZLE_GEOMETRY Generate nozzle wall radius r(x) and area A(x).
%
%   [x, r, A] = nozzle_geometry(geom_type, params, r_t, N)
%
%   Inputs:
%       geom_type - 'conical' or 'bell'
%       params    - geometry parameters:
%                   conical: [theta_deg, L]  (half-angle in degrees, length)
%                   bell:    [a, b, L]       (curvature, initial slope, length)
%       r_t       - throat radius (default 1)
%       N         - number of x-stations (default 200)
%
%   Outputs:
%       x - axial coordinate along diverging section [0, L]
%       r - wall radius r(x)
%       A - cross-sectional area A(x) = pi * r(x)^2
%
%   Note: Only the diverging section is modeled (x=0 is the throat).

    if nargin < 3 || isempty(r_t)
        r_t = 1;
    end
    if nargin < 4 || isempty(N)
        N = 200;
    end

    switch lower(geom_type)
        case 'conical'
            theta_deg = params(1);
            L = params(2);
            theta = deg2rad(theta_deg);
            x = linspace(0, L, N);
            r = r_t + x .* tan(theta);

        case 'bell'
            r_e         = params(1);                    % exit radius
            theta_exit  = deg2rad(params(2));           % exit half-angle (deg input)
            L           = params(3);
            a = (L * tan(theta_exit) - (r_e - r_t)) / L^2;
            b = tan(theta_exit) - 2 * a * L;
            x = linspace(0, L, N);
            r = r_t + a .* x.^2 + b .* x;
    

        otherwise
            error('Unknown geometry type: %s. Use ''conical'' or ''bell''.', geom_type);
    end

    A = pi .* r.^2;
end

