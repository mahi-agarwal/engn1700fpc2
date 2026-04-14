function sol = solve_nozzle(x, A, P_b, P_0, T_0, gamma, R)
% SOLVE_NOZZLE Solve quasi-1D steady-state flow through a C-D nozzle.
%
%   sol = solve_nozzle(x, A, P_b, P_0, T_0, gamma, R)
%
%   Inputs:
%       x     - axial coordinate vector (diverging section, x=0 at throat)
%       A     - area distribution A(x), same size as x
%       P_b   - back pressure (Pa)
%       P_0   - stagnation pressure (Pa)
%       T_0   - stagnation temperature (K)
%       gamma - ratio of specific heats (default 1.4)
%       R     - specific gas constant (default 287 J/kg/K)
%
%   Output: struct sol with fields:
%       M, P, T, rho, V  - flow properties along x
%       regime            - 'subsonic', 'shock_in_nozzle', 'overexpanded',
%                           'perfectly_expanded', 'underexpanded'
%       shock_x           - shock location (NaN if no shock inside nozzle)
%       shock_idx         - index of shock location

    if nargin < 6 || isempty(gamma), gamma = 1.4; end
    if nargin < 7 || isempty(R), R = 287; end

    N = length(x);
    A_t = A(1);          % throat area (x=0)
    A_e = A(end);        % exit area
    A_star = A_t;        % sonic throat
    AR = A ./ A_star;    % area ratio distribution
    AR_e = A_e / A_star; % exit area ratio

    Pb_P0 = P_b / P_0;

    % --- Compute critical pressure ratios to classify the flow ---

    % Case 1: Fully subsonic (no supersonic flow in diverging section)
    M_sub_exit = mach_from_area(AR_e, gamma, 'sub');
    [~, P_sub_exit] = isentropic_flow(M_sub_exit, gamma);
    % P_sub_exit is P_e/P_0 for fully subsonic case

    % Case 2: Normal shock at exit plane
    M_sup_exit = mach_from_area(AR_e, gamma, 'sup');
    [~, P_sup_exit_isen] = isentropic_flow(M_sup_exit, gamma);
    [M_post_shock, P2_P1, ~, ~, ~] = normal_shock(M_sup_exit, gamma);
    P_shock_at_exit = P_sup_exit_isen * P2_P1;
    % P_shock_at_exit is P_e/P_0 if normal shock sits right at exit

    % Case 3: Fully supersonic isentropic (design condition)
    % P_sup_exit_isen is P_e/P_0 for fully supersonic isentropic

    % --- Classify regime based on back pressure ---

    if Pb_P0 >= P_sub_exit
        % Back pressure too high: fully subsonic in diverging section
        % (or subsonic with some deceleration)
        sol = solve_subsonic(x, AR, P_0, T_0, gamma, R);
        sol.regime = 'subsonic';
        sol.shock_x = NaN;
        sol.shock_idx = NaN;

    elseif Pb_P0 >= P_shock_at_exit
        % Normal shock inside the diverging section
        sol = solve_with_shock(x, AR, Pb_P0, P_0, T_0, gamma, R);
        sol.regime = 'shock_in_nozzle';

    elseif Pb_P0 >= P_sup_exit_isen
        % Overexpanded: supersonic throughout, but P_e < P_b outside
        % (oblique shocks form outside the nozzle)
        sol = solve_supersonic(x, AR, P_0, T_0, gamma, R);
        sol.regime = 'overexpanded';
        sol.shock_x = NaN;
        sol.shock_idx = NaN;

    elseif abs(Pb_P0 - P_sup_exit_isen) < 1e-6
        % Perfectly expanded: design condition
        sol = solve_supersonic(x, AR, P_0, T_0, gamma, R);
        sol.regime = 'perfectly_expanded';
        sol.shock_x = NaN;
        sol.shock_idx = NaN;

    else
        % Underexpanded: P_e > P_b, expansion fans outside
        sol = solve_supersonic(x, AR, P_0, T_0, gamma, R);
        sol.regime = 'underexpanded';
        sol.shock_x = NaN;
        sol.shock_idx = NaN;
    end

    sol.x = x;
    sol.A = A;
end


function sol = solve_subsonic(x, AR, P_0, T_0, gamma, R)
% Fully subsonic flow in diverging section
    N = length(x);
    M = zeros(1, N);
    for i = 1:N
        M(i) = mach_from_area(AR(i), gamma, 'sub');
    end
    [T_ratio, P_ratio, rho_ratio, ~] = isentropic_flow(M, gamma);
    sol.M   = M;
    sol.P   = P_ratio .* P_0;
    sol.T   = T_ratio .* T_0;
    sol.rho = sol.P ./ (R .* sol.T);
    sol.V   = M .* sqrt(gamma .* R .* sol.T);
end


function sol = solve_supersonic(x, AR, P_0, T_0, gamma, R)
% Fully supersonic isentropic flow in diverging section
    N = length(x);
    M = zeros(1, N);
    M(1) = 1.0; % sonic at throat
    for i = 2:N
        M(i) = mach_from_area(AR(i), gamma, 'sup');
    end
    [T_ratio, P_ratio, rho_ratio, ~] = isentropic_flow(M, gamma);
    sol.M   = M;
    sol.P   = P_ratio .* P_0;
    sol.T   = T_ratio .* T_0;
    sol.rho = sol.P ./ (R .* sol.T);
    sol.V   = M .* sqrt(gamma .* R .* sol.T);
end


function sol = solve_with_shock(x, AR, Pb_P0, P_0, T_0, gamma, R)
% Normal shock inside the diverging section.
% Strategy: search for the shock location where the post-shock subsonic
% flow, when decelerated isentropically to the exit, matches P_b.

    N = length(x);
    AR_e = AR(end);

    % Search for shock location: try each station in the diverging section
    % At station i, supersonic Mach is M1. After shock, M2 and P02/P01.
    % Post-shock flow decelerates subsonically to exit with new A/A*.
    % Find where exit pressure matches P_b.

    % Build supersonic Mach distribution first
    M_sup = zeros(1, N);
    M_sup(1) = 1.0;
    for i = 2:N
        M_sup(i) = mach_from_area(AR(i), gamma, 'sup');
    end

    % Compute exit pressure for shock at each station
    P_exit_ratio = zeros(1, N);
    for i = 2:N-1
        M1 = M_sup(i);
        [M2, ~, ~, ~, P02_P01] = normal_shock(M1, gamma);

        % New A* after shock (A*_new = A_shock / A_ratio(M2))
        [~, ~, ~, A_ratio_M2] = isentropic_flow(M2, gamma);
        A_star_new = AR(i) * (pi * 1^2) / A_ratio_M2; % in terms of A_t
        % Actually, we need A/A*_new at the exit
        AR_exit_new = AR(end) * (pi * 1^2) / A_star_new;

        % Subsonic Mach at exit with new A*
        if AR_exit_new >= 1
            M_exit_sub = mach_from_area(AR_exit_new, gamma, 'sub');
            [~, P_exit_isen] = isentropic_flow(M_exit_sub, gamma);
            % P_exit / P_0 = (P_exit/P02) * (P02/P01) * (P01/P_0)
            % P01 = P_0 (isentropic upstream of shock)
            P_exit_ratio(i) = P_exit_isen * P02_P01;
        else
            P_exit_ratio(i) = NaN;
        end
    end

    % Find shock location by interpolation: where P_exit_ratio = Pb_P0
    valid = ~isnan(P_exit_ratio(2:end-1));
    idx_range = 2:N-1;
    idx_valid = idx_range(valid);

    if isempty(idx_valid)
        % Fallback: no valid shock location found, treat as supersonic
        sol = solve_supersonic(x, AR, P_0, T_0, gamma, R);
        sol.shock_x = NaN;
        sol.shock_idx = NaN;
        return;
    end

    % Interpolate to find exact shock index
    P_valid = P_exit_ratio(idx_valid);
    % P_exit_ratio increases as shock moves toward throat (higher M1 = stronger shock)
    % Find where P_valid crosses Pb_P0
    shock_idx = idx_valid(1); % default
    for k = 1:length(P_valid)-1
        if (P_valid(k) - Pb_P0) * (P_valid(k+1) - Pb_P0) <= 0
            % Linear interpolation for better accuracy
            frac = (Pb_P0 - P_valid(k)) / (P_valid(k+1) - P_valid(k));
            shock_idx = idx_valid(k) + round(frac);
            break;
        end
    end
    shock_idx = max(2, min(shock_idx, N-1));

    % --- Build flow field ---
    M1_shock = M_sup(shock_idx);
    [M2_shock, ~, ~, ~, P02_P01_shock] = normal_shock(M1_shock, gamma);

    % Upstream of shock: supersonic isentropic
    M = zeros(1, N);
    M(1) = 1.0;
    for i = 2:shock_idx
        M(i) = M_sup(i);
    end

    % New A* after shock
    [~, ~, ~, A_ratio_M2] = isentropic_flow(M2_shock, gamma);
    A_star_new_ratio = AR(shock_idx) / A_ratio_M2; % A*_new / A_t

    % Downstream of shock: subsonic isentropic with new A*
    M(shock_idx) = M2_shock; % post-shock Mach (overwrites pre-shock)
    for i = shock_idx+1:N
        AR_new = AR(i) / A_star_new_ratio;
        if AR_new >= 1
            M(i) = mach_from_area(AR_new, gamma, 'sub');
        else
            M(i) = mach_from_area(max(AR_new, 1.0), gamma, 'sub');
        end
    end

    % Compute thermodynamic properties
    % Upstream of shock: use P_0
    % Downstream of shock: use P_02 = P_0 * P02_P01
    P_02 = P_0 * P02_P01_shock;

    sol.M = M;
    sol.P = zeros(1, N);
    sol.T = zeros(1, N);

    % Upstream of shock (isentropic from P_0, T_0)
    [T_r_up, P_r_up, ~, ~] = isentropic_flow(M(1:shock_idx-1), gamma);
    sol.P(1:shock_idx-1) = P_r_up .* P_0;
    sol.T(1:shock_idx-1) = T_r_up .* T_0;

    % Downstream of shock (isentropic from P_02, T_0 -- T_0 unchanged across shock)
    [T_r_dn, P_r_dn, ~, ~] = isentropic_flow(M(shock_idx:N), gamma);
    sol.P(shock_idx:N) = P_r_dn .* P_02;
    sol.T(shock_idx:N) = T_r_dn .* T_0;

    sol.rho = sol.P ./ (R .* sol.T);
    sol.V   = M .* sqrt(gamma .* R .* sol.T);
    sol.shock_x   = x(shock_idx);
    sol.shock_idx = shock_idx;
end
