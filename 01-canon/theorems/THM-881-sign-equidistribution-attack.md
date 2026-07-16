---
id: THM-881
title: THE SIGN-EQUIDISTRIBUTION ATTACK — ownership, periodicity, and the sup-norm descent: (P1) every R_s endpoint is a section boundary j/(7e) of a runner (owner), so all differences lie in (1/(7·lcm E))Z; (P2) hence w ↦ Q_s(w) is PERIODIC mod P = 7·lcm(E) and sup over ALL integer w is a finite exact computation — "Q_s ≤ C·diam" is DECIDABLE per cluster; full-period sups computed: clean-w sup/diam ∈ [10.4, 15.3] on the family [0..5,t], t = 6..50 (constant ≤ 16 over ALL w); (P3) DFT descent on Z_P: Q_s = 2π² Σ_n k̂_P(n)|S(n)|² with S(0) = Σε = 0 and C_P = Σ_{n≠0}|k̂_P(n)| → 1/6, so Q_s ≤ (π²/3)(1+o(1))·max_{n≠0}|S(n)|² — SHARP O(diam) IS REDUCED TO THE SUP-NORM OF ONE SIGN EXPONENTIAL SUM (measured max|S|²/M ∈ [3.7, 5.9] on the bank)
status: P1, P2, P3 PROVED (ownership: endpoints are section-boundary rationals; periodicity: denominators divide 7·lcm; descent: DFT + S(0)=0 + computed C_P); the per-instance sharp bound CLOSED by finite check on the scanned family; the UNIFORM closure = the sup-norm lemma max_{n≠0}|S(n)|² = O(M) (HYP-6994, measured O(1)·M)
source: klein-2026-07-16-S314 (cont.2); executes the backlog attack "sign-equidistribution via the 7-section difference structure" (THM-880 handoff)
depends_on: [THM-880 (the bilinear form), THM-729 (frame)]
verification: 04-computation/sign_equidistribution_owner_descent_klein_S314.py -> 05-knowledge/results/sign_equidistribution_owner_descent_klein_S314.out (6/6)
---

# THM-881 — ownership, periodicity, sup-norm descent

**(P1) Ownership.** R_s is a Boolean combination of section conditions ⌊7{ex}⌋ = c, so every
endpoint is a boundary j/(7e) of a unique-generic owner e ∈ E; every endpoint difference has
denominator dividing P = 7·lcm(E). Block decomposition by owner pairs computed; the
within-owner blocks live on the arithmetic progressions (1/(7e))Z.

**(P2) Periodicity ⟹ decidability.** {wΔ} depends only on w mod P, so Q_s(w) is P-periodic
(verified exactly). Therefore sup_{w ∈ Z} Q_s is a FINITE computation:
> "Q_s ≤ C·diam for cluster E" is a DECIDABLE statement, per cluster, over ALL integer w.
Full-period sups on [0..5,t]: clean-w (coprime to P) sup/diam = 11.9, 10.6, 10.4, 15.3 at
t = 6, 12, 25, 50 — **sharp O(diam) holds with constant ≤ 16 on every scanned instance,
uniformly in w** (resonant w inflates the sup only at small t: w = 105 at t=6, 55 at t=12;
from t = 25 the sup is already attained at coprime w — the resonant/clean distinction fades).

**(P3) The sup-norm descent (the reduction).** DFT of K(t) = {t}(1−{t}) on Z_P:
Q_s = 2π² Σ_{n mod P} k̂_P(n)|S(n)|², where S(n) = Σ_k ε_k e(n·wp_k·(P)/P) is the endpoint
sign sum on the finite group. S(0) = Σε_k = 0 kills the mean mode, and
C_P := Σ_{n≠0}|k̂_P(n)| = 0.16665…→ 1/6. Hence
> **Q_s ≤ 2π² C_P · max_{n≠0}|S(n)|² ≤ (π²/3)(1+o(1)) · max_{n≠0}|S(n)|².**
The sharp problem is now the sup-norm of ONE exponential sum with ±1 coefficients over the
M endpoint positions. Parseval forces max|S|² ≥ M; measured max|S|²/M ∈ [3.7, 5.9] on the
bank — near the Parseval floor.

**The residual, named (HYP-6994):** max_{n≠0}|S(n)|² = O(M) (or O(M log P)) for the
7-section endpoint sign patterns. This one sup-norm lemma ⟹ Q_s = O(M) = O(diam)
unconditionally ⟹ |S| = O(√diam) ⟹ the density row closes (THM-729). Also decidable
per instance by P2 (S(n) is P-periodic data). Per-owner refinement: each S_e is supported
on an AP (1/(7e))Z with M_e = O(1)-small counts on the bank (max|S_e|²/M_e ≤ 3.63) — the
descent can be run owner-by-owner if the joint sum resists.
