%% FPC2 Part 3: Supersonic Civilian Transport Nozzle Optimization
%  Team Nozzle Be There — Diego Delgado, Faz Zaidi, Mahi Agarwal
%  Application: Supersonic Civilian Transport (Boom Overture class)
%  Design Mach range: 1.5 - 2.5
%
%  This script runs the full optimization and produces all figures.
%  Run this single file to reproduce all results.

clear; clc; close all;

%% ===== Gas Properties =====
gamma = 1.4;         % ratio of specific heats (air)
R     = 287;         % specific gas constant for air (J/kg/K)

%% ===== Operating Conditions =====
% The SST uses a medium-bypass turbofan (BPR ~ 3, like Boom Symphony).
% At cruise, the nozzle handles primarily core exhaust (hot, high pressure).
% At takeoff, core and bypass streams mix before the nozzle, producing
% a cooler, lower-pressure mixed exhaust. This is the key to noise compliance.

% --- CRUISE: 18 km altitude, core-dominated exhaust ---
P_0_cruise   = 400000;    % stagnation pressure (Pa)
T_0_cruise   = 1000;      % stagnation temperature (K) — core exhaust
P_inf_cruise = 7505;      % ambient at 18 km (Pa)

% --- TAKEOFF: sea level, mixed bypass+core exhaust (BPR ~ 3) ---
% T_0_mixed = (T_core + BPR * T_bypass) / (1 + BPR)
%           = (1000 + 3 * 300) / 4 = 475 K
% P_0 is also lower due to mixing and reduced throttle setting
P_0_takeoff   = 250000;   % stagnation pressure (Pa) — mixed flow
T_0_takeoff   = 475;      % stagnation temperature (K) — mixed exhaust
P_inf_takeoff = 101325;   % ambient at sea level (Pa)

fprintf('=== Operating Conditions ===\n');
fprintf('Cruise:  P_0=%.0f kPa, T_0=%.0f K, P_inf=%.1f kPa, NPR=%.1f\n', ...
        P_0_cruise/1e3, T_0_cruise, P_inf_cruise/1e3, P_0_cruise/P_inf_cruise);
fprintf('Takeoff: P_0=%.0f kPa, T_0=%.0f K, P_inf=%.1f kPa, NPR=%.1f\n', ...
        P_0_takeoff/1e3, T_0_takeoff, P_inf_takeoff/1e3, P_0_takeoff/P_inf_takeoff);

%% ===== Noise Constraint =====
% FAA Stage 5 / ICAO Chapter 14 compliance requires Ve <= 400 m/s
% for a bare nozzle at the sideline certification point.
% (NASA HSR program, 14 CFR Part 36, FAA NPRM Docket FAA-2020-0316)
%
% NASA's variable mixing nozzle uses deployable vanes/slots at takeoff
% that enhance jet mixing with ambient air, reducing effective acoustic
% velocity. NASA TM-2005-213894 reports 3-5 EPNdB reduction.
% Since acoustic power ~ Ve^8 (Lighthill), a dB reduction raises the
% allowable physical Ve:
%   delta_dB = 80 * log10(Ve_actual / Ve_bare)
%   => Ve_actual = Ve_bare * 10^(delta_dB / 80)
%
% With 4 dB mixing benefit (mid-range of NASA 3-5 dB):
%   Ve_with_vanes = 400 * 10^(4/80) = 400 * 1.122 = 449 m/s
%
% Vane deployment costs ~1% thrust penalty at takeoff (NASA Glenn data).

Ve_bare_limit = 400;       % m/s — bare nozzle FAA limit
mixing_dB     = 4;         % dB noise reduction from vanes (NASA: 3-5)
vane_CF_penalty = 0.01;    % 1% thrust loss when vanes deployed

% Effective Ve limit with mixing vanes deployed
Ve_limit = Ve_bare_limit * 10^(mixing_dB / 80);
fprintf('Noise: bare limit = %d m/s, with %.0f dB mixing vanes = %.0f m/s\n', ...
        Ve_bare_limit, mixing_dB, Ve_limit);
fprintf('Vane thrust penalty at takeoff: %.1f%%\n', vane_CF_penalty * 100);

%% ===== Nozzle Parameters =====
r_t = 1; % normalized throat radius

% --- Conical nozzle: r(x) = r_t + x*tan(theta) ---
% Sweep: half-angle theta (deg) and length L (multiples of r_t)
theta_vec = linspace(3, 25, 40);   % degrees
L_vec     = linspace(1.5, 8, 40);  % r_t units

% --- Bell (parabolic) nozzle: r(x) = r_t + a*x^2 + b*x ---
% Sweep: curvature a and initial slope b, fixed L = 5*r_t
% Range chosen so exit radius r_e = 1 + 25a + 5b stays in [~1.1, ~2.5]
% (area ratios ~1.2 to ~6, comparable to conical range)
a_vec = linspace(-0.01, 0.03, 40);   % curvature (negative = bell shape that turns inward)
b_vec = linspace(0.05, 0.40, 40);    % initial expansion slope

%% ===== Validation =====
fprintf('\n--- Validation checks ---\n');

% Check 1: Isentropic at M=2, gamma=1.4
[T_r, P_r, rho_r, A_r] = isentropic_flow(2.0, 1.4);
fprintf('M=2.0: T/T0=%.4f (expect 0.5556), P/P0=%.4f (expect 0.1278), A/A*=%.4f (expect 1.6875)\n', ...
        T_r, P_r, A_r);

% Check 2: Normal shock at M=2, gamma=1.4
[M2, P21, T21, ~, P0201] = normal_shock(2.0, 1.4);
fprintf('Normal shock M1=2.0: M2=%.4f (expect 0.5774), P2/P1=%.4f (expect 4.5000), P02/P01=%.4f (expect 0.7209)\n', ...
        M2, P21, P0201);

% Check 3: Area-Mach inversion
M_test = mach_from_area(1.6875, 1.4, 'sup');
fprintf('A/A*=1.6875 -> M=%.4f (expect 2.0000)\n', M_test);

M_test_sub = mach_from_area(1.6875, 1.4, 'sub');
fprintf('A/A*=1.6875 (subsonic) -> M=%.4f (expect 0.3722)\n', M_test_sub);

%% ===== Run Optimization: Conical Nozzle =====
fprintf('\n========================================\n');
fprintf('CONICAL NOZZLE OPTIMIZATION\n');
fprintf('========================================\n');
results_conical = optimize_nozzle('conical', theta_vec, L_vec, ...
    P_0_cruise, T_0_cruise, P_inf_cruise, ...
    P_0_takeoff, T_0_takeoff, P_inf_takeoff, ...
    gamma, R, r_t, Ve_limit, vane_CF_penalty);

%% ===== Run Optimization: Bell Nozzle =====
fprintf('\n========================================\n');
fprintf('BELL NOZZLE OPTIMIZATION\n');
fprintf('========================================\n');
results_bell = optimize_nozzle('bell', a_vec, b_vec, ...
    P_0_cruise, T_0_cruise, P_inf_cruise, ...
    P_0_takeoff, T_0_takeoff, P_inf_takeoff, ...
    gamma, R, r_t, Ve_limit, vane_CF_penalty);

%% ===== Generate All Figures =====
fprintf('\n========================================\n');
fprintf('GENERATING FIGURES\n');
fprintf('========================================\n');
plot_results(results_conical, results_bell, ...
             P_0_cruise, T_0_cruise, P_inf_cruise, ...
             P_0_takeoff, T_0_takeoff, P_inf_takeoff, ...
             gamma, R, r_t, Ve_limit);

%% ===== Thermal Check =====
fprintf('\n--- Thermal Analysis ---\n');
fprintf('Stagnation temperature T_0 (cruise) = %.0f K\n', T_0_cruise);
fprintf('Inconel 718 melting point  = 1573 K\n');
fprintf('T_0 < T_melt: %s\n', mat2str(T_0_cruise < 1573));
fprintf('Since flow expands isentropically, T(x) <= T_0 everywhere.\n');
fprintf('The nozzle housing will NOT reach melting temperature.\n');

%% ===== Print Summary =====
fprintf('\n========================================\n');
fprintf('FINAL SUMMARY\n');
fprintf('========================================\n');
fprintf('Application: Supersonic Civilian Transport\n');
fprintf('Noise constraint: V_e <= %.0f m/s at takeoff (FAA Stage 5 / Ch.14 + %.0f dB vanes)\n', Ve_limit, mixing_dB);
fprintf('\nConical optimal: theta=%.1f deg, L=%.2f r_t\n', ...
        results_conical.opt_params);
fprintf('  C_F(cruise) = %.4f, V_e(takeoff) = %.1f m/s\n', ...
        results_conical.opt_CF, results_conical.opt_Ve);
fprintf('\nBell optimal: a=%.4f, b=%.4f\n', results_bell.opt_params);
fprintf('  C_F(cruise) = %.4f, V_e(takeoff) = %.1f m/s\n', ...
        results_bell.opt_CF, results_bell.opt_Ve);

if results_bell.opt_CF > results_conical.opt_CF
    fprintf('\n>> Bell nozzle outperforms conical by %.2f%% in C_F at cruise.\n', ...
            (results_bell.opt_CF - results_conical.opt_CF) / results_conical.opt_CF * 100);
else
    fprintf('\n>> Conical nozzle outperforms bell by %.2f%% in C_F at cruise.\n', ...
            (results_conical.opt_CF - results_bell.opt_CF) / results_bell.opt_CF * 100);
end
fprintf('========================================\n');
