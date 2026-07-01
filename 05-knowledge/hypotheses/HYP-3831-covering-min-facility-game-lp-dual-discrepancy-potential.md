---
id: HYP-3831
title: THE COVERING-MIN IS AN ADVERSARIAL FACILITY-LOCATION GAME -- LP-dual = the {1,killer} equioscillation; the discrepancy/coverage potential (Koksma-Hlawka) is the potential-function/PoA handle on inf-measure; and the Pochhammer fiber fraction f(n)~1/sqrt(pi n) is the pi/EVEN/MEASURE side complementary to the sqrt(p)/ODD/CERTIFICATE side. THE GAME (owner + CS6840): defender (observer) picks time t, payoff = min_v ||vt|| (distance to nearest runner=facility); defender maximizes => M(S)=max_t min_v ||vt|| = covering-min; adversary picks config S, value = min_S M(S); LRC = this >= 1/n. Adversarial facility location / covering radius: runners cover the circle over time, observer seeks the least-covered point, adversary minimizes the max-empty-radius. VERIFIED n=14 (finite reduction mod Phi6=183): (1) GAME VALUE M(S)=14/183=0.0765 >= 1/n=0.0714, at k*=14 (=t*); (2) LP-DUAL = the binding runners {1,182} at signed phase-residues {+14,-14}={+-n} = the 2-point Chebyshev equioscillation (S73/HYP-3813), dual weights (1/2,1/2), observer PINNED between +-n; (3) DISCREPANCY POTENTIAL: the runner cloud at t* has star-discrepancy D*=0.0765, three-gap sizes {1,n,2n}/Phi6; the coverage potential PHI=2(n-1)M=1.989 (>1!) => the naive packing/union bound gives ONLY the floor 1/(2(n-1)) (PHI=1), and reaching 1/n REQUIRES the discrepancy of the AP cloud -- exactly where the Fourier/SOS method stalls (HYP-3791); (4) POCHHAMMER fiber f(n)=(1/2)_(n-2)/(n-2)!, f(14)=0.1612~1/sqrt(pi*14)=0.1508 = the pi/measure side (Wallis/Gamma, iota-even), complementary to sqrt(p) (Gauss/Paley, iota-odd, HYP-3820). KEY HONEST LOCATION: the potential-function/PoA argument recovers M>=1/(2(n-1)) trivially (coverage=1); the factor-2 gap to 1/n IS the discrepancy structure of the covering-min AP cloud, the crux
status: MIXED (verified computation + a located crux, honest about the open gap). VERIFIED (covering_min_facility_game_klein.py, n=14): M(S)=14/183 at k*=14; LP-dual binding {1,182}={+-n} Chebyshev 2-point; D*=0.0765, three-gap {1,14,28}/183; coverage potential 2(n-1)M=1.989; f(14)=0.1612~1/sqrt(14pi)=0.1508. HONEST: the game/LP-dual framing is a clean RESTATEMENT (the value + certificate are known, THM-523/S73); the NEW content is (a) the explicit attacker-defender/PoA formulation, (b) locating the potential-function floor at 1/(2(n-1)) with the factor-2-to-1/n gap = the AP-cloud discrepancy (why Koksma-Hlawka/Fourier stalls), (c) the pi/measure vs sqrt(p)/certificate complementarity. NOT a new bound; a game-theoretic reframe that pinpoints the crux as a discrepancy statement.
source: klein-2026-07-01-S85
depends_on:
  - HYP-3813   # S80: the covering-min phase cloud = AP + antipodal killer; three-gap {1,n,2n} (the game's optimal config)
  - HYP-3806   # S73: Chebyshev 2-point dual {1, killer} (the LP-dual certificate)
related:
  - HYP-3791   # the SOS/Fourier obstruction -- WHY the potential/PoA argument stalls at the discrepancy step
  - HYP-3830   # S85: the complementary sqrt(p)/certificate side (this is the pi/measure side)
  - THM-523    # covering condition M=n/Phi6
external: Koksma-Hlawka inequality (discrepancy bound on QMC integration error) = the potential/inf-measure handle; CS6840 algorithmic game theory (facility location, price of anarchy, potential games); Pochhammer symbol (1/2)_m (the fiber fraction ~ 1/sqrt(pi n))
results:
  - 04-computation/covering_min_facility_game_klein.py
  - 05-knowledge/results/covering_min_facility_game_klein.out
---

# HYP-3831 — the covering-min as an adversarial facility-location game

*(Renumbered 3821→3831 alongside HYP-3830, to clear opus's frontier — swarm-convergence note.)*

## The game (attacker–defender)
- **Defender** (observer) picks a time `t`; payoff = `min_{v∈S} ||vt||` (distance to the nearest runner = facility).
- Defender maximizes: `M(S) = max_t min_v ||vt||` = the **covering-min** (best loneliness given `S`).
- **Adversary** picks the config `S` (`|S|=n-1`); value = `min_S M(S)`. **LRC ⟺ this ≥ 1/n.**
This is adversarial **facility location** / the covering radius: the runners cover the circle over time; the
observer seeks the least-covered point; the adversary places facilities to minimize the max-empty-radius.

## Verified (n=14, finite reduction mod Φ₆=183)
1. **Game value** `M(S) = 14/183 = 0.0765 ≥ 1/n = 0.0714`, at `k* = 14` (`t* = 14/183`).
2. **LP-dual** = the binding runners `{1, 182}` at signed phase-residues `{+14, −14} = {±n}` — the **2-point
   Chebyshev equioscillation** (S73/HYP-3813), dual weights `(½,½)`; the observer is **pinned between ±n**.
3. **Discrepancy potential** (Koksma–Hlawka): the runner cloud at `t*` has star-discrepancy `D* = 0.0765` and
   three-gap sizes `{1,n,2n}/Φ₆`. The **coverage potential** `Φ = 2(n−1)M = 1.989`.
4. **Pochhammer fiber** `f(n) = (½)_{n-2}/(n-2)!`, `f(14) = 0.1612 ~ 1/√(πn) = 0.1508` — the **π/measure side**.

## The located crux (the honest keeper)
Because `Φ = 2(n−1)M ≈ 2 > 1`, the naive **packing/union bound** (potential = total danger coverage) gives
**only** the trivial floor `M ≥ 1/(2(n-1))` (the point where `Φ = 1`). The runners' danger arcs total ≈ 2 —
they cover the circle *twice over* — so the observer escapes **not** by counting but by the **overlap
structure**. The factor-2 gap from `1/(2(n-1))` to `1/n` **is** the discrepancy of the covering-min **AP
cloud** (three-gap `{1,n,2n}`): this is exactly where the Fourier/SOS/potential method stalls (HYP-3791). So
the game/PoA framing pinpoints the crux as a **discrepancy statement** about the AP cloud, not a covering count.

## Two complementary halves of the apex
- **π / EVEN / measure** side: `f(n) ~ 1/√(πn)` (Wallis/Gamma, ι-even) — the lonely-measure/fiber scale.
- **√p / ODD / certificate** side: Gauss `i√p`, Paley skew, the apex group (HYP-3820, ι-odd).
The covering-min value `n/Φ₆` is rational (neither √p nor π appears in it); √p and π live in the *certificate*
and the *measure* respectively.

## Net
The covering-min is an adversarial facility-location game; its LP-dual is the `{1,killer}=±n` equioscillation;
the discrepancy/coverage potential is the potential-function/PoA handle, which recovers `1/(2(n-1))` and locates
the residual factor-2 as the AP-cloud discrepancy. A game-theoretic reframe that names the crux, complementing
the √p certificate side (HYP-3820) with the π measure side.
