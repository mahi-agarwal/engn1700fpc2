# Quantities of Interest (QoIs) for Supersonic Civilian Transport Nozzle Optimization

**Team Nozzle Be There** — Diego Delgado, Faz Zaidi, Mahi Agarwal  
**Application:** Supersonic Civilian Transport (Boom Overture class, Design Mach 1.5–2.5)

---

## Regulatory Context

The FAA has stated that supersonic transport aircraft must meet the **same landing and takeoff (LTO) noise standards as subsonic aircraft**, with no relaxation or special exemption [1]. This means compliance with **FAA Stage 5 / ICAO Chapter 14**, which requires a cumulative noise margin of **17 EPNdB below Stage 3** limits across three certification measurement points: sideline, flyover, and approach [2, 3]. Boom Technology has publicly committed the Overture to meeting Chapter 14 limits [4].

NASA's High Speed Research (HSR) program and subsequent Commercial Supersonic Technology (CST) studies established that Chapter 14 compliance translates to an engineering constraint on **exit jet velocity at takeoff of V_e ≤ 400 m/s (1,300 ft/s)** at the sideline certification point [5, 6]. For reference, the Concorde's exit velocity at takeoff was approximately 600–700 m/s — far exceeding modern limits [5].

---

## QoI 1: Thrust Coefficient, C_F

### What it measures
The thrust coefficient quantifies how efficiently the nozzle converts total (stagnation) pressure into thrust. It is the primary **performance metric at cruise**, where fuel economy dominates operating cost for a commercial SST.

### Why it matters
For supersonic civilian transport, even a 2–3% improvement in C_F at cruise translates directly to increased range or reduced fuel burn over long transoceanic routes. At the cruise design point, the nozzle should be as close to perfectly expanded as possible to extract maximum thrust from the available pressure ratio [7].

### Mathematical definition

$$C_F = \frac{\dot{m} \, V_e + (P_e - P_\infty) \, A_e}{P_0 \, A_t}$$

where:
- $\dot{m}$ = mass flow rate (set by choked throat conditions)
- $V_e$ = exit velocity
- $P_e$ = exit static pressure
- $P_\infty$ = ambient static pressure
- $A_e$ = exit area
- $A_t$ = throat area
- $P_0$ = stagnation (total) pressure

At perfect expansion ($P_e = P_\infty$), the pressure thrust term vanishes and $C_F$ is maximized for a given area ratio. Any deviation — a normal shock inside the nozzle, overexpansion, or underexpansion — reduces $C_F$.

### Evaluation condition
**Cruise:** Altitude ~18 km, $P_\infty \approx 7.5$ kPa, high nozzle pressure ratio (NPR ~ 10–15).

---

## QoI 2: Exit Jet Velocity at Takeoff, V_e

### What it measures
The exit jet velocity at takeoff is the direct physical driver of jet noise. It serves as the **noise compliance metric**, linking nozzle aerodynamic design to FAA Stage 5 / ICAO Chapter 14 certification limits.

### Why it matters
Lighthill's aeroacoustic analogy [8] establishes that acoustic power radiated by a turbulent jet scales as:

$$W_{\text{acoustic}} \propto V_e^{\,8}$$

This eighth-power dependence means that even modest reductions in exit velocity yield dramatic noise reductions. A 20% decrease in $V_e$ reduces acoustic power by approximately 83%. NASA and industry research (GE, Rolls-Royce) conducted under the HSR and HSCT programs determined that for an aircraft in the Overture weight class (~200,000 lb MTOW), Chapter 14 compliance at the sideline and flyover measurement points requires [5, 6]:

$$\boxed{V_e \leq 400 \text{ m/s} \quad (\approx 1{,}300 \text{ ft/s}) \quad \text{at takeoff}}$$

For comparison [5, 9]:
- Concorde at takeoff (with afterburner): $V_e \approx 600$–$700$ m/s — **non-compliant**
- Modern subsonic high-bypass turbofan: $V_e \approx 250$–$350$ m/s — **compliant with margin**
- SST with medium-bypass engine (unaugmented): $V_e \approx 400$–$500$ m/s — **marginal**
- SST with mixer-ejector nozzle deployed: $V_e \approx 350$–$400$ m/s — **target compliance range**

### Mathematical definition

From quasi-1D theory, the exit velocity is computed as:

$$V_e = M_e \cdot a_e = M_e \sqrt{\gamma \, R \, T_e}$$

where the exit temperature $T_e$ is determined from the isentropic relation (or post-shock relations if a normal shock sits inside the diverging section):

$$T_e = T_0 \left(1 + \frac{\gamma - 1}{2} M_e^2 \right)^{-1}$$

**Any nozzle geometry that produces $V_e > 400$ m/s at the takeoff condition is rejected from the feasible design space**, regardless of its cruise performance.

### Evaluation condition
**Takeoff:** Sea level, $P_\infty = 101.3$ kPa, low NPR (~2.5–4.0).

---

## QoI 3: Exit Pressure Ratio, P_e / P_∞

### What it measures
The ratio of exit static pressure to ambient pressure classifies the **nozzle operating regime** and identifies where wave structures (shocks, expansion fans) form. It serves as the **regime diagnostic** on parameter space contour maps.

### Why it matters
Imperfect expansion degrades both performance and noise:

| Regime | Condition | Physical consequence |
|--------|-----------|---------------------|
| **Perfectly expanded** | $P_e / P_\infty = 1$ | No external waves; maximum thrust extracted; minimum shock-associated noise |
| **Overexpanded** | $P_e / P_\infty < 1$ | Oblique shocks or normal shock inside nozzle; thrust loss; **shock-associated broadband noise (BBSN) and screech** add to jet mixing noise [10] |
| **Underexpanded** | $P_e / P_\infty > 1$ | Expansion fans at nozzle lip; wasted potential energy; **shock cell noise from diamond pattern** [10] |

Both overexpanded and underexpanded conditions generate additional noise components beyond the baseline jet mixing noise captured by Lighthill scaling. Broadband shock-associated noise and screech tones can add **5–10 dB** to the far-field spectrum at certain observer angles, pushing total EPNdB further above the Chapter 14 limit [10, 11].

### Mathematical definition

The exit pressure is determined from the isentropic pressure-Mach relation:

$$\frac{P_e}{P_0} = \left(1 + \frac{\gamma - 1}{2} M_e^2 \right)^{-\gamma / (\gamma - 1)}$$

The exit pressure ratio relative to ambient is then:

$$\frac{P_e}{P_\infty} = \frac{P_e / P_0}{P_\infty / P_0}$$

On the parameter space contour maps, this QoI draws the **regime boundaries**: identifying which combinations of nozzle shape parameters produce normal shocks inside the diverging section at takeoff, which configurations are over- or underexpanded, and where the perfectly expanded contour ($P_e / P_\infty = 1$) falls at cruise.

### Evaluation condition
**Both cruise and takeoff:** This QoI is evaluated at both operating conditions to classify the nozzle behavior across the flight envelope.

---

## How the QoIs Work Together for FAA Compliance

The optimization finds the nozzle shape parameters that:

1. **Maximize $C_F$ at cruise** (best fuel efficiency over transoceanic routes)
2. **Subject to $V_e(\text{takeoff}) \leq 400$ m/s** (FAA Stage 5 / Chapter 14 noise compliance)
3. **With $P_e / P_\infty$ mapping the regime boundaries** (identifying shock locations, over/underexpansion, and additional shock-associated noise sources)

On the parameter space contour maps, the **feasible region** is everything below the $V_e = 400$ m/s noise boundary. Within that feasible region, the optimal design is the point with the highest cruise $C_F$. The $P_e / P_\infty$ contours explain *why* each design succeeds or fails — whether the penalty comes from an internal shock, divergence loss, or imperfect expansion.

---

## References

[1] Federal Aviation Administration, "Noise Certification of Supersonic Airplanes," Notice of Proposed Rulemaking, Docket No. FAA-2020-0316, 2020.

[2] 14 CFR Part 36, Appendices A & B, "Noise Standards: Aircraft Type and Airworthiness Certification."

[3] International Civil Aviation Organization, "Annex 16 — Environmental Protection, Volume I: Aircraft Noise," Chapter 14, 7th Edition, 2014.

[4] Boom Technology, "Overture Specifications and Environmental Commitments," boomsupersonic.com, 2022–2024.

[5] NASA Glenn Research Center, "Jet Noise Reduction Concepts for Supersonic Civil Transport," NASA TM-2005-213894, 2005.

[6] NASA, "Variable Cycle Engine Studies for HSCT Application," NASA CR-2002-211662, 2002.

[7] J. D. Anderson, *Modern Compressible Flow: With Historical Perspective*, 3rd ed., McGraw-Hill, 2003, Ch. 5 & 10.

[8] M. J. Lighthill, "On Sound Generated Aerodynamically. I. General Theory," *Proceedings of the Royal Society of London A*, vol. 211, pp. 564–587, 1952.

[9] K. Viswanathan, "Aeroacoustics of Hot Jets," *Journal of Fluid Mechanics*, vol. 516, pp. 39–82, 2004.

[10] C. K. W. Tam, "Supersonic Jet Noise," *Annual Review of Fluid Mechanics*, vol. 27, pp. 17–43, 1995.

[11] J. Bridges and M. P. Brown, "Parametric Testing of Chevrons on Single Flow Hot Jets," AIAA Paper 2004-2824, 2004.
