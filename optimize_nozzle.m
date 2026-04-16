function results = optimize_nozzle(geom_type, param1_vec, param2_vec, ...
                                    P_0_cruise, T_0_cruise, P_inf_cruise, ...
                                    P_0_takeoff, T_0_takeoff, P_inf_takeoff, ...
                                    gamma, R, r_t, Ve_limit, vane_CF_penalty)
% OPTIMIZE_NOZZLE Grid search over nozzle parameter space.
%
%   Models a variable-area SST nozzle:
%   - CRUISE: full nozzle geometry, core exhaust (hot, high P_0)
%   - TAKEOFF: adapted (smaller) exit area, mixed bypass+core exhaust
%     (cooler T_0, lower P_0). The noise constraint V_e <= Ve_limit
%     applies here.
%
%   The optimization finds the nozzle shape that maximizes cruise C_F
%   while keeping the adapted takeoff V_e within the FAA noise limit.

    N1 = length(param1_vec);
    N2 = length(param2_vec);

    CF_cruise       = NaN(N1, N2);
    CF_takeoff      = NaN(N1, N2);
    Ve_takeoff      = NaN(N1, N2);
    Pe_Pinf_cruise  = NaN(N1, N2);
    Pe_Pinf_takeoff = NaN(N1, N2);
    Me_cruise       = NaN(N1, N2);
    Me_takeoff      = NaN(N1, N2);
    AR_cruise       = NaN(N1, N2);
    AR_takeoff_opt  = NaN(N1, N2);
    regime_cruise   = strings(N1, N2);
    regime_takeoff  = strings(N1, N2);

    L_bell = 5 * r_t;
    A_t = pi * r_t^2;

    % --- Takeoff analysis: find noise-limited exit Mach ---
    NPR_takeoff = P_0_takeoff / P_inf_takeoff;

    % Perfect expansion Mach at takeoff NPR
    Me_perfect_to = sqrt(2/(gamma-1) * (NPR_takeoff^((gamma-1)/gamma) - 1));
    [~,~,~, AR_perfect_to] = isentropic_flow(Me_perfect_to, gamma);
    Te_perfect = T_0_takeoff / (1 + (gamma-1)/2 * Me_perfect_to^2);
    Ve_perfect_to = Me_perfect_to * sqrt(gamma * R * Te_perfect);

    fprintf('Takeoff: NPR=%.2f, T_0=%.0f K\n', NPR_takeoff, T_0_takeoff);
    fprintf('  Perfect expansion: Me=%.3f, A/A*=%.3f, Ve=%.1f m/s\n', ...
            Me_perfect_to, AR_perfect_to, Ve_perfect_to);

    % Find max Mach that keeps Ve <= Ve_limit
    % Ve(Me) = Me * sqrt(gamma*R*T_0/(1+(g-1)/2*Me^2))
    f_Ve = @(Me) Me * sqrt(gamma*R*T_0_takeoff / (1 + (gamma-1)/2 * Me^2)) - Ve_limit;

    % Check if even M=0.01 is under the limit (it should be)
    Ve_at_low_M = 0.01 * sqrt(gamma*R*T_0_takeoff);
    if Ve_at_low_M > Ve_limit
        error('T_0_takeoff too high: even M=0.01 exceeds Ve limit.');
    end

    % Check if perfectly expanded is within limit
    if Ve_perfect_to <= Ve_limit
        Me_max_noise = Me_perfect_to;
        AR_max_noise = AR_perfect_to;
        fprintf('  Perfect expansion is within noise limit.\n');
        noise_binding = false;
    else
        % Find the maximum Mach number that satisfies the noise constraint
        % Ve(Me) increases monotonically for reasonable Me, peaks around Me~2-3
        % Search in [0.01, Me_perfect_to]
        Me_max_noise = fzero(f_Ve, [0.01, Me_perfect_to]);
        [~,~,~, AR_max_noise] = isentropic_flow(Me_max_noise, gamma);
        fprintf('  Noise-limited: Me_max=%.3f, A/A*=%.3f, Ve=%.1f m/s\n', ...
                Me_max_noise, AR_max_noise, Ve_limit);
        noise_binding = true;
    end

    fprintf('Optimizing %s nozzle: %d x %d = %d evaluations\n', ...
            geom_type, N1, N2, N1*N2);

    for i = 1:N1
        for j = 1:N2
            try
                % Build geometry
                if strcmpi(geom_type, 'conical')
                    params = [param1_vec(i), param2_vec(j)];
                else
                    params = [param1_vec(i), param2_vec(j), L_bell];
                end

                [x, r, A] = nozzle_geometry(geom_type, params, r_t);
                A_e = A(end);
                A_t_local = A(1);

                % Validity checks
                if A_e <= A_t_local, continue; end
                if any(r <= 0), continue; end
                if any(diff(A) < -1e-10 * A_t_local), continue; end
                if A_e / A_t_local > 10, continue; end

                AR_cruise(i,j) = A_e / A_t_local;

                % ============================================
                % CRUISE: full geometry, core exhaust
                % ============================================
                sol_cruise = solve_nozzle(x, A, P_inf_cruise, P_0_cruise, ...
                                          T_0_cruise, gamma, R);
                qoi_cruise = compute_qois(sol_cruise, P_inf_cruise, ...
                                          P_0_cruise, gamma, R);

                % Apply divergence loss correction (lambda factor)
                % Conical nozzles lose thrust because exit flow is not
                % axial — the velocity has a radial component.
                % lambda = (1 + cos(theta))/2 for a conical nozzle
                % Bell nozzles have nearly axial exit flow: lambda ~ 0.99
                % (Anderson, Modern Compressible Flow, Ch. 10; Sutton, Rocket
                %  Propulsion Elements, Ch. 3)
                if strcmpi(geom_type, 'conical')
                    theta_rad = atan((r(end) - r(1)) / (x(end) - x(1)));
                    lambda = (1 + cos(theta_rad)) / 2;
                else
                    % Bell nozzle: exit angle IS the parameter
                    theta_exit = deg2rad(param2_vec(j));
                    lambda = 0.99;
                end

                CF_cruise(i,j)      = qoi_cruise.C_F * lambda;
                Pe_Pinf_cruise(i,j) = qoi_cruise.Pe_Pinf;
                Me_cruise(i,j)      = qoi_cruise.M_e;
                regime_cruise(i,j)  = string(qoi_cruise.regime);

                % ============================================
                % TAKEOFF: adapted nozzle, mixed exhaust
                % ============================================
                % Variable geometry closes the exit to a smaller area ratio.
                % Use the noise-limited area ratio (or perfect expansion if
                % that's already within the limit).
                % The nozzle can only close DOWN from the cruise geometry,
                % so if cruise AR < takeoff target, use cruise AR.

                if noise_binding
                    AR_to_target = min(AR_max_noise, A_e/A_t_local);
                else
                    AR_to_target = min(AR_perfect_to, A_e/A_t_local);
                end
                AR_takeoff_opt(i,j) = AR_to_target;

                % Compute takeoff V_e and C_F directly from isentropic
                % relations at the adapted exit Mach number.
                % This avoids building a nozzle geometry just for takeoff.
                if AR_to_target <= 1.0
                    % Convergent-only: choked throat, exit is sonic
                    Me_to = 1.0;
                else
                    % Find exit Mach (subsonic or supersonic depends on NPR)
                    % For the adapted nozzle, we WANT subsonic exit
                    % if the back pressure can support it
                    Me_to_sub = mach_from_area(AR_to_target, gamma, 'sub');
                    Me_to_sup = mach_from_area(AR_to_target, gamma, 'sup');
                    [~, P_sub] = isentropic_flow(Me_to_sub, gamma);
                    [~, P_sup] = isentropic_flow(Me_to_sup, gamma);
                    Pb_P0_to = P_inf_takeoff / P_0_takeoff;

                    if Pb_P0_to >= P_sub
                        Me_to = Me_to_sub;
                    elseif Pb_P0_to <= P_sup
                        Me_to = Me_to_sup;
                    else
                        % Shock in nozzle — use subsonic exit
                        Me_to = Me_to_sub;
                    end
                end

                Te_to = T_0_takeoff / (1 + (gamma-1)/2 * Me_to^2);
                Ve_to = Me_to * sqrt(gamma * R * Te_to);
                [~, Pe_P0_to] = isentropic_flow(Me_to, gamma);
                Pe_to = Pe_P0_to * P_0_takeoff;
                rho_to = Pe_to / (R * Te_to);

                % C_F at takeoff
                m_dot_to = rho_to * Ve_to * A_t_local * AR_to_target;
                thrust_to = m_dot_to * Ve_to + (Pe_to - P_inf_takeoff) * A_t_local * AR_to_target;
                CF_to = thrust_to / (P_0_takeoff * A_t_local);

                % Apply mixing vane thrust penalty at takeoff
                CF_takeoff(i,j)      = CF_to * (1 - vane_CF_penalty);
                Ve_takeoff(i,j)      = Ve_to;
                Pe_Pinf_takeoff(i,j) = Pe_to / P_inf_takeoff;
                Me_takeoff(i,j)      = Me_to;
                regime_takeoff(i,j)  = "adapted";

            catch ME
                continue;
            end
        end

        if mod(i, 10) == 0
            fprintf('  %d / %d rows complete\n', i, N1);
        end
    end

    % --- Find optimal: max CF_cruise where Ve_takeoff <= Ve_limit ---
    feasible = (Ve_takeoff - Ve_limit) < 2.0 & ~isnan(CF_cruise); % 2 m/s tolerance
    if any(feasible(:))
        CF_feasible = CF_cruise;
        CF_feasible(~feasible) = -Inf;
        [~, lin_idx] = max(CF_feasible(:));
        [opt_i, opt_j] = ind2sub([N1, N2], lin_idx);
    else
        warning('No feasible point (Ve <= %.0f m/s). Returning best CF.', Ve_limit);
        [~, lin_idx] = max(CF_cruise(:));
        [opt_i, opt_j] = ind2sub([N1, N2], lin_idx);
    end

    % Pack results
    results.param1_vec      = param1_vec;
    results.param2_vec      = param2_vec;
    results.CF_cruise       = CF_cruise;
    results.CF_takeoff      = CF_takeoff;
    results.Ve_takeoff      = Ve_takeoff;
    results.Pe_Pinf_cruise  = Pe_Pinf_cruise;
    results.Pe_Pinf_takeoff = Pe_Pinf_takeoff;
    results.Me_cruise       = Me_cruise;
    results.Me_takeoff      = Me_takeoff;
    results.AR_cruise       = AR_cruise;
    results.AR_takeoff_opt  = AR_takeoff_opt;
    results.regime_cruise   = regime_cruise;
    results.regime_takeoff  = regime_takeoff;
    results.feasible        = feasible;
    results.opt_idx         = [opt_i, opt_j];
    results.opt_CF          = CF_cruise(opt_i, opt_j);
    results.opt_Ve          = Ve_takeoff(opt_i, opt_j);
    results.opt_CF_takeoff  = CF_takeoff(opt_i, opt_j);
    results.opt_params      = [param1_vec(opt_i), param2_vec(opt_j)];
    results.geom_type       = geom_type;
    results.Ve_perfect_takeoff = Ve_perfect_to;
    results.AR_perfect_takeoff = AR_perfect_to;
    results.noise_binding   = noise_binding;

    fprintf('Optimal %s: param1=%.3f, param2=%.3f\n', ...
            geom_type, results.opt_params(1), results.opt_params(2));
    fprintf('  Cruise:  CF=%.4f, Me=%.3f, AR=%.3f\n', ...
            results.opt_CF, Me_cruise(opt_i,opt_j), AR_cruise(opt_i,opt_j));
    fprintf('  Takeoff: CF=%.4f, Ve=%.1f m/s, AR_adapted=%.3f\n', ...
            results.opt_CF_takeoff, results.opt_Ve, AR_takeoff_opt(opt_i,opt_j));
end
