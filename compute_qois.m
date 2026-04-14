function qoi = compute_qois(sol, P_inf, P_0, gamma, R)
% COMPUTE_QOIS Extract quantities of interest from a solved nozzle flow.
%
%   qoi = compute_qois(sol, P_inf, P_0, gamma, R)
%
%   Inputs:
%       sol   - struct from solve_nozzle
%       P_inf - ambient pressure (Pa)
%       P_0   - stagnation pressure (Pa)
%       gamma - ratio of specific heats
%       R     - specific gas constant (J/kg/K)
%
%   Output: struct qoi with fields:
%       C_F       - thrust coefficient
%       V_e       - exit velocity (m/s)
%       Pe_Pinf   - exit pressure ratio P_e / P_inf
%       regime    - flow regime string
%       M_e       - exit Mach number
%       T_e       - exit temperature (K)

    if nargin < 4, gamma = 1.4; end
    if nargin < 5, R = 287; end

    A_t = sol.A(1);
    A_e = sol.A(end);
    P_e = sol.P(end);
    V_e = sol.V(end);
    M_e = sol.M(end);
    T_e = sol.T(end);
    rho_e = sol.rho(end);

    % Thrust coefficient: C_F = (m_dot * V_e + (P_e - P_inf) * A_e) / (P_0 * A_t)
    % Mass flow from choked throat: m_dot = P_0 * A_t * sqrt(gamma/(R*T_0)) * (2/(gamma+1))^((gamma+1)/(2*(gamma-1)))
    % But we can also compute it directly:
    %   m_dot = rho_e * V_e * A_e  (continuity)
    % Use momentum thrust + pressure thrust form:
    m_dot = rho_e * V_e * A_e;
    thrust = m_dot * V_e + (P_e - P_inf) * A_e;
    C_F = thrust / (P_0 * A_t);

    % Exit pressure ratio
    Pe_Pinf = P_e / P_inf;

    % Pack output
    qoi.C_F     = C_F;
    qoi.V_e     = V_e;
    qoi.Pe_Pinf = Pe_Pinf;
    qoi.regime  = sol.regime;
    qoi.M_e     = M_e;
    qoi.T_e     = T_e;
end
