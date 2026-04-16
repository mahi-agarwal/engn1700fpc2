function plot_results(results_conical, results_bell, ...
                      P_0_cruise, T_0_cruise, P_inf_cruise, ...
                      P_0_takeoff, T_0_takeoff, P_inf_takeoff, ...
                      gamma, R, r_t, Ve_limit)
% PLOT_RESULTS Generate all presentation figures for FPC2 Part 3.
%
%   Figures produced:
%     1. Cruise C_F contour map (per geometry) with regime overlay
%     2. Cruise exit Mach contour map (per geometry)
%     3. P(x), T(x), rho(x) profiles at cruise for optimal designs
%     4. P(x), T(x), rho(x) profiles at takeoff (full geometry w/ shock)
%     5. Nozzle geometry comparison (both optimal shapes)
%     6. Conical vs Bell comparison bar chart

    geom_names = {'Conical', 'Bell'};
    all_results = {results_conical, results_bell};
    param_labels = {{'Half-angle \theta (deg)', 'Length L / r_t'}, ...
                    {'Curvature a', 'Initial slope b'}};

    % =====================================================================
    % Figures 1-2: Parameter space contour maps for each geometry
    % =====================================================================
    for g = 1:2
        res = all_results{g};
        [P1, P2] = meshgrid(res.param1_vec, res.param2_vec);

        % --- (a) C_F at cruise with regime boundaries ---
        fig = figure('Name', sprintf('%s: C_F at Cruise', geom_names{g}), ...
               'Color', 'w', 'Position', [50+600*(g-1) 500 560 450]);

        % Identify regimes numerically for overlay
        regime_num = NaN(size(res.regime_cruise));
        regime_num(res.regime_cruise == "underexpanded")    = 1;
        regime_num(res.regime_cruise == "perfectly_expanded") = 2;
        regime_num(res.regime_cruise == "overexpanded")      = 3;
        regime_num(res.regime_cruise == "shock_in_nozzle")   = 4;
        regime_num(res.regime_cruise == "subsonic")          = 5;

        contourf(P1', P2', res.CF_cruise, 20, 'LineStyle', 'none');
        hold on;
        % Overlay regime boundaries as contour lines
        if any(~isnan(regime_num(:)))
            contour(P1', P2', regime_num, 'LineColor', [0.3 0.3 0.3], ...
                    'LineWidth', 1, 'LineStyle', ':');
        end
        % Mark optimal
        plot(res.opt_params(1), res.opt_params(2), 'rp', ...
             'MarkerSize', 18, 'MarkerFaceColor', 'r', 'LineWidth', 2);
        hold off;
        cb = colorbar; ylabel(cb, 'C_F', 'FontName', 'Arial', 'FontSize', 11);
        colormap(fig, parula);
        xlabel(param_labels{g}{1}, 'FontName', 'Arial', 'FontSize', 12);
        ylabel(param_labels{g}{2}, 'FontName', 'Arial', 'FontSize', 12);
        title(sprintf('%s Nozzle: Cruise C_F (with divergence correction)', ...
              geom_names{g}), 'FontName', 'Arial', 'FontSize', 13);
        set(gca, 'FontName', 'Arial');

        % --- (b) Exit Mach at cruise ---
        fig2 = figure('Name', sprintf('%s: M_e at Cruise', geom_names{g}), ...
               'Color', 'w', 'Position', [50+600*(g-1) 50 560 450]);

        contourf(P1', P2', res.Me_cruise, 20, 'LineStyle', 'none');
        hold on;
        % Mark sonic line (Me = 1)
        if any(res.Me_cruise(:) < 1) && any(res.Me_cruise(:) > 1)
            contour(P1', P2', res.Me_cruise, [1 1], ...
                    'LineColor', 'w', 'LineWidth', 2, 'LineStyle', '--');
        end
        plot(res.opt_params(1), res.opt_params(2), 'rp', ...
             'MarkerSize', 18, 'MarkerFaceColor', 'r', 'LineWidth', 2);
        hold off;
        cb = colorbar; ylabel(cb, 'M_e', 'FontName', 'Arial', 'FontSize', 11);
        colormap(fig2, turbo);
        xlabel(param_labels{g}{1}, 'FontName', 'Arial', 'FontSize', 12);
        ylabel(param_labels{g}{2}, 'FontName', 'Arial', 'FontSize', 12);
        title(sprintf('%s Nozzle: Exit Mach Number at Cruise', geom_names{g}), ...
              'FontName', 'Arial', 'FontSize', 13);
        set(gca, 'FontName', 'Arial');

        % --- (c) Area ratio map ---
        fig3 = figure('Name', sprintf('%s: Area Ratio', geom_names{g}), ...
               'Color', 'w', 'Position', [50+600*(g-1) 50 560 450]);

        contourf(P1', P2', res.AR_cruise, 20, 'LineStyle', 'none');
        hold on;
        plot(res.opt_params(1), res.opt_params(2), 'rp', ...
             'MarkerSize', 18, 'MarkerFaceColor', 'r', 'LineWidth', 2);
        hold off;
        cb = colorbar; ylabel(cb, 'A_e / A_t', 'FontName', 'Arial', 'FontSize', 11);
        colormap(fig3, parula);
        xlabel(param_labels{g}{1}, 'FontName', 'Arial', 'FontSize', 12);
        ylabel(param_labels{g}{2}, 'FontName', 'Arial', 'FontSize', 12);
        title(sprintf('%s Nozzle: Exit Area Ratio A_e/A_t', geom_names{g}), ...
              'FontName', 'Arial', 'FontSize', 13);
        set(gca, 'FontName', 'Arial');
    end

    % =====================================================================
    % Figures 3-4: P(x), T(x), rho(x) profiles for optimal geometries
    % =====================================================================
    line_colors = {[0 0.447 0.741], [0.850 0.325 0.098]}; % blue, orange

    for g = 1:2
        res = all_results{g};
        opt = res.opt_params;

        if strcmpi(res.geom_type, 'conical')
            params = [opt(1), opt(2)];
        else
            params = [opt(1), opt(2), 5 * r_t];
        end

        [x, r, A] = nozzle_geometry(res.geom_type, params, r_t);

        % --- CRUISE profiles ---
        sol_c = solve_nozzle(x, A, P_inf_cruise, P_0_cruise, T_0_cruise, gamma, R);

        figure('Name', sprintf('%s Optimal: Cruise Flow', geom_names{g}), ...
               'Color', 'w', 'Position', [100 100 750 650]);

        subplot(3,1,1);
        plot(x/r_t, sol_c.P / P_0_cruise, '-', 'Color', line_colors{g}, 'LineWidth', 2);
        ylabel('P / P_0', 'FontName', 'Arial', 'FontSize', 12);
        title(sprintf('%s Optimal at Cruise  (C_F = %.4f, M_e = %.2f, %s)', ...
              geom_names{g}, res.opt_CF, res.Me_cruise(res.opt_idx(1), res.opt_idx(2)), ...
              sol_c.regime), 'FontName', 'Arial', 'FontSize', 13);
        grid on; set(gca, 'FontName', 'Arial');

        subplot(3,1,2);
        plot(x/r_t, sol_c.T, '-', 'Color', [0.8 0.1 0.1], 'LineWidth', 2);
        ylabel('T (K)', 'FontName', 'Arial', 'FontSize', 12);
        yline(1573, 'k--', 'Inconel 718 limit (1573 K)', ...
              'LineWidth', 1.5, 'FontName', 'Arial', 'LabelHorizontalAlignment', 'left');
        grid on; set(gca, 'FontName', 'Arial');

        subplot(3,1,3);
        plot(x/r_t, sol_c.rho, '-', 'Color', [0.2 0.6 0.2], 'LineWidth', 2);
        xlabel('x / r_t', 'FontName', 'Arial', 'FontSize', 12);
        ylabel('\rho (kg/m^3)', 'FontName', 'Arial', 'FontSize', 12);
        grid on; set(gca, 'FontName', 'Arial');

        % --- TAKEOFF profiles (full cruise geometry at takeoff conditions) ---
        % This shows the shock inside the nozzle when the variable geometry
        % is NOT adapted — illustrating WHY you need variable area.
        sol_t = solve_nozzle(x, A, P_inf_takeoff, P_0_takeoff, T_0_takeoff, gamma, R);

        figure('Name', sprintf('%s Optimal: Takeoff (Unadapted)', geom_names{g}), ...
               'Color', 'w', 'Position', [100 100 750 650]);

        subplot(3,1,1);
        plot(x/r_t, sol_t.P / P_0_takeoff, '-', 'Color', line_colors{g}, 'LineWidth', 2);
        hold on;
        if ~isnan(sol_t.shock_x)
            xline(sol_t.shock_x / r_t, 'r--', 'Normal Shock', ...
                  'LineWidth', 2, 'FontName', 'Arial', 'FontSize', 10, ...
                  'LabelHorizontalAlignment', 'left');
        end
        hold off;
        ylabel('P / P_0', 'FontName', 'Arial', 'FontSize', 12);
        title(sprintf('%s at Takeoff (unadapted, %s) — shows need for variable area', ...
              geom_names{g}, sol_t.regime), 'FontName', 'Arial', 'FontSize', 13);
        grid on; set(gca, 'FontName', 'Arial');

        subplot(3,1,2);
        plot(x/r_t, sol_t.T, '-', 'Color', [0.8 0.1 0.1], 'LineWidth', 2);
        hold on;
        if ~isnan(sol_t.shock_x)
            xline(sol_t.shock_x / r_t, 'r--', 'LineWidth', 2);
        end
        hold off;
        ylabel('T (K)', 'FontName', 'Arial', 'FontSize', 12);
        yline(1573, 'k--', 'Inconel 718 limit', 'LineWidth', 1.5, 'FontName', 'Arial');
        grid on; set(gca, 'FontName', 'Arial');

        subplot(3,1,3);
        plot(x/r_t, sol_t.rho, '-', 'Color', [0.2 0.6 0.2], 'LineWidth', 2);
        hold on;
        if ~isnan(sol_t.shock_x)
            xline(sol_t.shock_x / r_t, 'r--', 'LineWidth', 2);
        end
        hold off;
        xlabel('x / r_t', 'FontName', 'Arial', 'FontSize', 12);
        ylabel('\rho (kg/m^3)', 'FontName', 'Arial', 'FontSize', 12);
        grid on; set(gca, 'FontName', 'Arial');
    end

    % =====================================================================
    % Figure 5: Nozzle geometry comparison
    % =====================================================================
    figure('Name', 'Nozzle Geometry Comparison', 'Color', 'w', ...
           'Position', [100 100 700 400]);
    for g = 1:2
        res = all_results{g};
        opt = res.opt_params;
        if strcmpi(res.geom_type, 'conical')
            params = [opt(1), opt(2)];
        else
            params = [opt(1), opt(2), 5 * r_t];
        end
        [x_div, r_div, ~] = nozzle_geometry(res.geom_type, params, r_t);
        
        % Add a cosine converging section for visualization
        r_inlet = 2.0 * r_t;       % inlet radius
        L_conv  = 1.5 * r_t;       % converging length
        N_conv  = 50;
        x_conv  = linspace(-L_conv, 0, N_conv);
        r_conv  = (r_inlet + r_t)/2 + (r_inlet - r_t)/2 * cos(pi * (x_conv + L_conv)/L_conv);
        
        % Concatenate (drop the duplicate throat point)
        x_full = [x_conv, x_div(2:end)];
        r_full = [r_conv, r_div(2:end)];
        
        plot(x_full/r_t, r_full/r_t, '-', 'Color', line_colors{g}, 'LineWidth', 2.5, ...
             'DisplayName', sprintf('%s  (C_F=%.4f)', geom_names{g}, res.opt_CF));
        hold on;
        plot(x_full/r_t, -r_full/r_t, '-', 'Color', line_colors{g}, 'LineWidth', 2.5, ...
             'HandleVisibility', 'off');
    end
    % Throat marker
    plot(0, r_t, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', ...
         'HandleVisibility', 'off');
    plot(0, -r_t, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', ...
         'HandleVisibility', 'off');
    xline(0, 'k:', 'Throat', 'LineWidth', 1, 'FontName', 'Arial', ...
          'LabelVerticalAlignment', 'top');
    hold off;
    xlabel('x / r_t', 'FontName', 'Arial', 'FontSize', 12);
    ylabel('r / r_t', 'FontName', 'Arial', 'FontSize', 12);
    title('Optimal Nozzle Geometries (Diverging Section)', ...
          'FontName', 'Arial', 'FontSize', 14);
    legend('Location', 'southeast', 'FontName', 'Arial', 'FontSize', 11);
    grid on; axis equal;
    set(gca, 'FontName', 'Arial');

    % =====================================================================
    % Figure 6: Comparison bar chart
    % =====================================================================
    figure('Name', 'Performance Comparison', 'Color', 'w', ...
           'Position', [100 100 650 500]);

    metrics = {'C_F (cruise)', 'C_F (takeoff)', 'M_e (cruise)', 'A_e/A_t (cruise)'};
    vals = [results_conical.opt_CF,       results_bell.opt_CF; ...
            results_conical.opt_CF_takeoff, results_bell.opt_CF_takeoff; ...
            results_conical.Me_cruise(results_conical.opt_idx(1), results_conical.opt_idx(2)), ...
            results_bell.Me_cruise(results_bell.opt_idx(1), results_bell.opt_idx(2)); ...
            results_conical.AR_cruise(results_conical.opt_idx(1), results_conical.opt_idx(2)), ...
            results_bell.AR_cruise(results_bell.opt_idx(1), results_bell.opt_idx(2))];

    b = bar(vals);
    b(1).FaceColor = line_colors{1};
    b(2).FaceColor = line_colors{2};
    set(gca, 'XTickLabel', metrics, 'FontName', 'Arial', 'FontSize', 11);
    ylabel('Value', 'FontName', 'Arial', 'FontSize', 12);
    title('Conical vs Bell: Optimal Design Comparison', ...
          'FontName', 'Arial', 'FontSize', 14);
    legend('Conical', 'Bell', 'Location', 'northwest', 'FontName', 'Arial');
    grid on;

    % =====================================================================
    % Print summary
    % =====================================================================
    fprintf('\n=== Results Summary ===\n');
    for g = 1:2
        res = all_results{g};
        fprintf('\n%s Nozzle Optimal:\n', geom_names{g});
        fprintf('  Parameters: [%.4f, %.4f]\n', res.opt_params);
        fprintf('  Cruise:   C_F = %.4f, M_e = %.2f, A_e/A_t = %.2f\n', ...
                res.opt_CF, ...
                res.Me_cruise(res.opt_idx(1), res.opt_idx(2)), ...
                res.AR_cruise(res.opt_idx(1), res.opt_idx(2)));
        fprintf('  Takeoff:  C_F = %.4f, V_e = %.1f m/s (limit: %.0f m/s)\n', ...
                res.opt_CF_takeoff, res.opt_Ve, Ve_limit);
        fprintf('  Feasible: %s\n', mat2str(res.opt_Ve < Ve_limit + 2));
    end
end
