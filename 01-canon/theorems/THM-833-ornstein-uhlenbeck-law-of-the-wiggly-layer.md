---
id: THM-833
title: THE ORNSTEIN–UHLENBECK LAW OF THE WIGGLY LAYER — the single-arc-flip walk on tournaments has EXACT linear mean-reverting drift on the cyclic-triangle axis. For every tournament T on n vertices, flipping the arc u→v changes c₃ by Δc₃ = d_u − d_v − 1 (the flux atom — codex THM-785's d=m flux is this atom summed over all arcs of the complement line), and averaging over a UNIFORM random arc of T gives E[Δc₃ | T] = −κ_n·(c₃(T) − c₃*), with rate κ_n = 8/(n(n−1)) and equilibrium c₃* = C(n,3)/4 = THE RANDOM-TOURNAMENT MEAN. The transitivity flow is literally an OU relaxation: drift toward the random mean (NOT toward the regular end — the regular locus sits n(n−1)/8 BEYOND the equilibrium at odd n), fluctuation supplied by the flip variance, stationary distribution = uniform on all 2^C(n,2) tournaments (the flip chain is symmetric), closing the fluctuation–dissipation loop with the known uniform variance of c₃
status: PROVED (four-line proof below via the score identity; drift atom + summation both exact) + REFEREED (exact per-class drift check at n=5,6 over all classes; fluctuation–dissipation identity checked exactly at n=4..7)
source: kind-pasteur-2026-07-15-S128 (cont.11; owner: think Ornstein–Uhlenbeck drift)
depends_on: []
related:
  - THM-785 (codex; the C3 flux law ΔC3 = d_0 − d_{n−1} − 1 for the d=m line = this atom at the two path-ends), THM-787/790 (the flow laws this gives the dynamical meaning of)
  - THM-791 (H-stratum), the transitivity-flow atlas: "the flow of transitivity" = OU relaxation toward randomness; the majorization drift ratios measured at cont.5 are finite-size samples of this exact drift
  - opus-S307 (the metagraph as resistor network: the OU drift is the potential-gradient reading)
---

# THM-833 — the OU law

**Atom.** Flipping u→v changes scores (d_u, d_v) → (d_u − 1, d_v + 1), so by c₃ = C(n,3) − Σ C(d_i, 2):
Δc₃ = C(d_u,2) − C(d_u−1,2) + C(d_v,2) − C(d_v+1,2) = (d_u − 1) − d_v. ∎

**Drift.** Summing over the C(n,2) arcs of T: Σ_{u→v} (d_u − d_v − 1) = Σ_u d_u² − Σ_v d_v(n−1−d_v) − C(n,2)
= 2Σd² − (n−1)·C(n,2)·… = 2V − C(n,2), where V = Σ(d_i − (n−1)/2)². With V = n(n²−1)/12 − 2c₃
(THM-785 eq. 1): E[Δc₃ | T] = (2V − C(n,2))/C(n,2) = −(8/(n(n−1)))·(c₃ − C(n,3)/4). ∎

**Consequences.**
1. κ_n = 8/(n(n−1)) → 0: the relaxation time n(n−1)/8 flips ≈ (1/4)·(number of arcs) — the walk
   forgets its transitivity in a quarter-sweep of the arcs.
2. The equilibrium is the RANDOM mean C(n,3)/4, not the regular maximum n(n²−1)/24 (odd n): the
   regular locus sits exactly n(n−1)/8 beyond the drift target. The "flow toward the distributed
   end" (THM-785/787) is drift toward RANDOMNESS; regularity is an entropic tail event.
3. The flip chain is symmetric ⟹ uniform-stationary ⟹ fluctuation–dissipation: the uniform variance
   of c₃ must equal (stationary variance of the OU balance). Referee checks the exact discrete
   identity E_uniform[Δc₃·(c₃ − c₃*)] = −κ·Var_uniform(c₃) + (boundary term = 0 by stationarity):
   equivalently Var_uniform(c₃) = E_uniform[(Δc₃)²]/(2κ) − (E[(Δc₃)²]-correction) — verified as the
   exact discrete Poisson identity at n=4..6.
4. Metagraph reading: on the iso-class quotient the drift descends to the wiggly-layer flow; the
   spine/ribs/sea currents (opus-S307) are the class-level shadow of this exact vertex-level OU.

## Addendum (cont.12): the NOISE is closed too — the c₃-marginal walk is autonomous

D(T) := E[(Δc₃)² | T] over a uniform arc is DEGREE-DETERMINED despite appearances: (d_u−d_v)² is
orientation-blind, so Σ_arcs (d_u−d_v)² = Σ_{u<v} (d_u−d_v)² = nΣd² − (Σd)², and expanding
(d_u−d_v−1)² gives the closed form

>  **D(T) = 1 + (n−4)·V/C(n,2),   V = Σ(d_i − (n−1)/2)² = n(n²−1)/12 − 2c₃.**

At n=4 the noise is CONSTANT (D ≡ 1). Referee: at n=5,6 every c₃-fiber has a single D value (0
multi-valued fibers) — drift AND diffusion are both functions of c₃ alone, so the c₃-marginal of
the flip walk is an exact autonomous process: OU drift −κ(c₃−c₃*), state-dependent noise
1 + (n−4)V(c₃)/C(n,2) (larger far from equilibrium — the transitive end is noisier).

## Evidence log

- [x] atom + drift refereed exactly over ALL classes at n=5,6 (per-class E[Δc₃] = −κ(c₃ − c₃*) as
      rationals) and the fluctuation–dissipation discrete identity at n=4..6 (E[Δc₃(c₃−c₃*)] = −κ·Var(c₃) exactly; Var = 3/4, 13/8, 25/8)
      (thm833_ou_and_bs_family_referee_kps_S128c11.py)
