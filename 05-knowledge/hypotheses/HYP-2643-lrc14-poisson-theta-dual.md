---
id: HYP-2643
title: LRC(14) Poisson/theta dual of the seven-sector signed correction — the convergent dual IS the finite x-cell (Weyl) evaluation; no absolutely convergent dual-lattice series exists; the signed series converges FAST for wide E and slowly only at the AP
status: PARTIALLY-TRUE / structural; the convergent dual is identified and the wide/AP split is verified exactly; does NOT close the cap by itself
source: kind-pasteur-2026-06-19-S12 (workflow angle: poisson-summation / theta-transform)
depends_on:
  - THM-538   # support-6 floor
  - HYP-2606  # the signed relation-lattice Fourier identity
  - THM-534
related:
  - MISTAKE-078  # absolute envelope diverges
  - HYP-2610     # stranger-contraction (wide-spread reduction)
  - HYP-2640     # signed/coset quotient is the ruler
  - OPEN-Q-108
---

# HYP-2641 — Poisson/theta dual of the LRC(14) signed correction

## The angle (from the endgame brief)
Apply Poisson summation to convert the absolutely-divergent relation-lattice sum
`corr(E) = sum_{0!=n in Lambda(E)} K(n)` into a sum over the dual lattice, hoping the
dual function is concentrated/smooth and the dual series converges, giving a usable bound
for WIDE E.

## What Poisson actually gives (VERIFIED exactly, fractions.Fraction)
1. **The Poisson dual = the finite x-cell (Weyl) evaluation.** The orbit-character-sum
   `integral_0^1 e^{2pi i (n.e) x} dx` equals `1_{n in Lambda(E)}` exactly (because `n.e`
   is an integer), so the relation-lattice theta sum is *literally* the orbit integral:
   `sum_{n in Lambda(E)} K(n) = integral_0^1 [1_cover(x) - M7(k)] dx`,
   an exact finite sum over the `O(7*sum e_i)` breakpoint cells. The dual "function" is the
   per-cell cover-indicator deviation, bounded by 1 and finitely supported -> convergent.
   **This dual is exactly the existing engine.** Poisson does not produce a *new* analytic
   series; it certifies that the engine's finite x-integral IS the convergent representation.
2. **`M7(k)` is the iid coupon limit:** `M7(k) = P(k-1 iid uniform sector labels cover {1..6})`
   (verified vs Monte Carlo: 0.02448/0.05770/0.10491 for k=8/9/10).
3. **Koksma–Hlawka / discrepancy is the WRONG dual.** The Erdős–Turán bound routes through
   `sum_{Lambda} 1/r(n)`, which diverges harmonically (AP: 86->436->992; squares: 18->73->169)
   because the support-6 floor (THM-538) does NOT apply to the absolute discrepancy object.
4. **The absolute envelope diverges for ALL E, including wide/Sidon.** `sum_{supp>=6} prod 0.6973/|n_j|`
   grows without bound even for squares (0.39->2.9->8.2->15.7 at L=1..4). Width does not rescue
   the absolute bound (confirms MISTAKE-078).

## The genuine payoff: signed series converges FAST for wide E
The **signed** box-cumulative sum `sum_{Lambda cap [-H,H]^{k-1}} K(n)`:
- **AP (k=8):** H=1,2,3 -> 0.0244, 0.0643, 0.0978 — crawls harmonically toward 0.3027 (SLOW).
- **squares (Sidon, k=8):** H=1,2,3 -> -0.0020, -0.0020, -0.0010 — already converged to the
  true value 0.000536 by H=1 (FAST).
Relations-per-shell: AP 20/962/8354/38518 vs squares 4/152/1388/6522. Sparse lattice (wide E)
=> few relations => fast signed convergence => small, easily-bounded correction. This is the
exact Poisson trade-off the brief asked for: it works precisely where it is needed (wide), and
is slow exactly at the AP, which is handled by the finite check.

## The ruler is R6 (support-6 relation DENSITY), not shortest length
`lambda_inf(shortest support-6 relation) = 1` for essentially every E (the `0,1` prefix and
3-term coincidences supply tiny relations), so it does not discriminate. But the **count** of
support-6 relations in a box tracks `corr` monotonically (box radius 2):
| E | R6(box2) | corr |
|---|---:|---:|
| AP | 982 | 0.302731 |
| near-AP | 924 | 0.183854 |
| odd-AP | 546 | 0.213476 |
| primes-ish | 414 | 0.009649 |
| squares | 156 | 0.000536 |
Signed-vs-absolute cancellation ratio `sum_T|P_T-iid|/|corr|` improves with width: AP 8.4,
near-AP 13, primes 102, squares 1379 — the signed structure is most effective for wide sets.

## Verdict on the angle
Poisson/theta gives a **real convergent representation** of the signed correction — but it is
the finite x-cell evaluation, not a new dual-lattice series, and no absolutely convergent
dual series exists (the divergence is intrinsic to the discontinuous 1/7-sector indicators).
The usable wide bound must retain the SIGNED support-6 structure. Net effect: this
**corroborates the existing wide-spread plan** (THM-538 + stranger-contraction HYP-2610) and
rules out a Poisson short-cut around it; the tight margin still lives in the finite check.

Script: `04-computation/lrc14_poisson_theta_dual_kps_s12.py`;
output: `05-knowledge/results/lrc14_poisson_theta_dual_kps_s12.out`.
-> THM-538, HYP-2606, HYP-2610, HYP-2640, MISTAKE-078, OPEN-Q-108.
