---
id: HYP-9425
title: "After the three-pillar re-typing of the Collatz strategy atlas: (1) the strong form of Mazur's theorem is the renewal model's statement — the density of {tau(n) <= C ln n} tends to 1/2 at C* = 3/ln(4/3), to 1 above it and to 0 below it, with the no-parameter Gaussian profile (sd of the normalised deviation -> 1 from below); (2) the fast-set exponent is the independent-3-adic E(C) (E(3.307) = 0.738, E(2.5) = 0.584) and the two-block dyadic exponent at 2^30-2^32 lies in [0.65, 0.80] at C = 3.307; (3) Mazur's mechanism is sign-blind: Theorem 1.1 and the companion hold for 3x-1 with the seeds (2^(2j+1)+1)/3, so the three 3x-1 cycles have positive-density log-time basins; (4) the dyadic oscillation of the floor set is the carry geometry of the admissible word shapes: the independent-3-adic word model reproduces the 36 quarter-dyadic floor counts (pre-registered, PASSED: log-correlation 0.992, sign agreement 24/27, total ratio 1.098); (5) a proof that never uses a pillar has a conclusion valid on the hostile varying it, the rows whose proofs use all three are counts or identities, and no proof using all three on one orbit with an unbounded pointwise conclusion exists in the atlas"
status: >
  OPEN (opus-2026-10-27-S64). (1) CONJECTURE (Lagarias-Weiss renewal model, ordinary convention); evidence: all n < 2^28,
  block sd of z = (tau - C* ln n)/(16.42 sqrt(ln n)) rising 0.949 -> 0.963 over k = 20..28, mean -0.038 -> -0.034,
  F_X(523/50) = 0.54 flat from 2^16 to 2^28 (pre-registered P1, P2 PASSED). (2) HEURISTIC model; measured two-block
  exponents 0.61-0.78 at 2^28 (mean 0.703) against E(3.307) = 0.738; the pre-registered monotone forms (P3, P4) FAILED on
  a period-two oscillation in the dyadic index (reported as failed). (3) CONJECTURE; the seed lemma is PROVED and
  FINITE-EXACT (q <= 6), the harmonic-mass identity is sign-symmetric (digest Theorem A), the transfer is not checked in
  the formal development. (4) VERIFIED: pre-registered in 9c744f56e2 after the quarter-dyadic probe (factors up to 2.85 inside one block),
  run in 918e97bfaf, PASSED (correlation 0.992, sign agreement 0.889, total ratio 1.098). (5) SYNTHESIS of the atlas re-typing (pillars as proof inputs, walls as conclusion controls); falsified by exhibiting
  a pointwise unbounded Collatz theorem whose proof uses the integer, the map and the sign on one orbit.
  Independently audited 2026-10-27 (SOUND WITH CORRECTIONS, applied; MISTAKE-685). Proposition F of the note
  (density zero below 3.307 ln n) is PROVED and is the unconditional shadow of (1)-(2).
session: opus-2026-10-27-S64
source: 05-knowledge/results/collatz_three_pillars_mazur_external_20261027.md (sections 2-3, 5)
related:
  - 05-knowledge/results/mazur_positive_density_20260928.md (Theorems A-C; the seed-1 test)
  - 05-knowledge/results/mod2_collapse_collatz_strategy_audit_20261027.md (kind-pasteur-2026-10-27-S6; the concurrent audit)
  - 05-knowledge/reference/TRUSTED-ADVANCES-2026-10-27.md (the ledger, Collatz rows added)
  - 01-canon/theorems/THM-4517-no-root-uniform-positive-density-of-collatz-basins.md (trunk basins; 3x-1 basin densities)
  - 01-canon/theorems/THM-4813-collatz-numerator-coordinate-forcing-cocycle-gauss-sum-intertwiner-and-sharp-spectral-exception-budgets.md (positive carries)
  - 01-canon/theorems/THM-4479-strategy-cube-distance-to-provability-sharp.md (the toggle space of the strategy cube)
  - 05-knowledge/results/collatz_coalescence_20261004_proof_strategies_and_bold_predictions.md (the atlas S1-S13)
scripts:
  - 04-computation/experiments/collatz_fast_convergence_density_20261027/fast_convergence_density_20261027.py (+ .out; P1-P5 pre-registered in commit 4232963739)
  - 04-computation/experiments/collatz_fast_convergence_density_20261027/sawtooth_probe_20261027.py (+ .out; post hoc)
  - 04-computation/experiments/collatz_fast_convergence_density_20261027/fast_set_exponent_model_20261027.py (+ .out)
---

# HYP-9425 — the fast set at the typical-word constant, the sign-blind mechanism, and the pillar typing

Notation: `tau(n)` is the ordinary total stopping time (`n -> n/2`, `n -> 3n+1`, every application counted);
`C* = 3/ln(4/3) = 10.42818`; `sigma = sqrt(2 (1 - C* ln 2)^2/(2 ln 2 - ln 3)) = 16.422`;
`F_X(C) = #{1 <= n < X : tau(n) <= C ln n}/(X - 1)`; `G_k(C)` the same density on `[2^(k-1), 2^k)`.

1. **Strong form of Mazur's theorem (CONJECTURE).** `lim_X F_X(C) = 0` for `C < C*`, `= 1/2` at `C = C*` and `= 1`
   for `C > C*`; more precisely `G_k(C) -> Phi((C - C*) sqrt(ln N_k)/sigma)` in the renewal scale, and the block standard
   deviation of `z(n) = (tau(n) - C* ln n)/(sigma sqrt(ln n))` tends to `1` from below; at the scales reached the
   block distribution of `z` is right-skewed (skewness `+0.6`, audit), so the Gaussian profile is a limit statement
   and the shape is still owed. Mazur's Theorem 1.1
   (positive lower density at `523/50`, T2) is the weakest form of the middle clause; its only improvable constant
   is `c`. **Pre-registered for the next scale (`2^32`, 2-bit sieve):** block sd in `[0.96, 1.00]` and mean in
   `[-0.040, -0.020]` at `k = 30, 31, 32`; `F_(2^32)(523/50)` in `[0.530, 0.550]`; the median of `tau/ln n` over
   `n < 2^32` in `[9.90, 10.05]`. What would refute (1): a block sd that stops rising below `0.97`, or a density at
   `523/50` that drifts outside `[0.50, 0.58]` with a trend.
2. **The fast-set exponent (HEURISTIC model).** `#{n < X : tau(n) <= C ln n} = X^(E(C) + o(1))` with
   `E(C) = max { kappa h(delta/kappa) - delta ln 3 : delta = (kappa ln 2 - 1)/ln 3 >= 0, kappa + delta <= C }`,
   `E(C*) = 1`. **Pre-registered:** at `2^30` and `2^32` the two-block dyadic exponents at `C = 3.307` lie in
   `[0.65, 0.80]` and at `C = 2.5` in `[0.50, 0.65]` (model `0.738`, `0.584`), and the 4-adic block exponents lie
   within `0.08` of the model at both constants. The pre-registered P3/P4 forms of the first run (monotone block
   densities; block ratios below `2`) are withdrawn as predictions; the oscillation is (4).
3. **Sign-blindness (CONJECTURE).** Mazur's Theorem 1.1 and his companion theorem hold for the `3x-1` map on the
   positive integers, with the size inequality reversed (the carries change sign: `1/x = 3^d 2^(-A) prod(1 - 1/(3x_j))`,
   so the source charge gains the factor `prod(1 - 1/(3x_j))^(-1) = 1 + O(d/16^(b_0))` on the selected histories)
   and with the seeds `r_j = (2^(2j+1) + 1)/3` (PROVED here: `3 r_j - 1 = 2^(2j+1)`,
   `v_3(2^(2j+1) + 1) = 1 + v_3(2j+1)`, the `r_j` permute `Z/3^q`; FINITE-EXACT `q <= 6`). Consequence: each of the
   three `3x-1` cycles has a positive-lower-density basin reached in logarithmic time. Test: replay the digest's
   P1-P3 for `3x-1`; then list the sign-dependent lemmas of the formal development (expected: the seed lemma and
   every use of the source-charge direction). Refuted by a sign-dependent step that does not transfer.
4. **The dyadic sawtooth is word geometry, not 3-adic conspiracy (pre-registered here before the run; PASSED:
   correlation `0.992`, sign agreement `24/27`, total ratio `1.098`; VERIFIED).** The floor set `{tau(n) < 3.307 ln n}` has quarter-dyadic densities varying by factors up
   to `2.85` inside one block (`24.1, 19.8, 8.5, 18.8` per `10^4` in `[2^27, 2^28)`). An integer of the floor set
   with T-word of length `k` and `d` odd steps is `n = (2^k - c_w)/3^d`, `c_w = sum_(odd steps j) 3^(d_j') 2^j`
   (`d_j'` the number of odd steps after `j`), so its position inside a block is fixed by the word's carry `c_w/2^k`.
   The independent-3-adic model (HYP (2)) predicts the floor count of a sub-block `B` as
   `sum_(k,d) C(k,d) 3^(-d) P_(k,d)(n in B, k + d < 3.307 ln n)`, `P` over uniformly random words of that shape ending
   with an even step. Prediction: over the 36 quarter-dyadic cells of `k = 20..28`, the Pearson correlation between
   `log(model)` and `log(observed)` exceeds `0.8`, and the model reproduces the sign of every within-block change
   between consecutive quarters in at least `70%` of the 27 comparisons. Refuted by a correlation below `0.5` or a
   sign agreement at chance level (`50%`); a pass means the oscillation is the carry geometry of the few admissible
   `(k, d)` shapes, a fail means the integrality events are correlated at this resolution (a 3-adic effect, which
   would be the interesting outcome). (The earlier draft of this item proposed a lattice-proximity test that is
   vacuous at this scale because the lattice `{k - d log_2 3}` is dense; withdrawn before any run.)
5. **Pillar completeness (SYNTHESIS, falsifiable).** Pillars are proof inputs and the atlas's walls are conclusion
   controls: a proof that never uses the integral coordinate proves a conclusion valid on the `Z_2`/rational
   hostile, one that never uses (or uses symmetrically) positivity proves a sheet-blind conclusion (Mazur), one that
   never uses the operation proves a conclusion valid for `5x+1` or under choice; the rows whose proofs use all three
   are counts or identities (S4, THM-4633, Mazur, THM-4813). Refuted by a Collatz theorem with a pointwise unbounded
   conclusion whose proof uses all three on one orbit; the candidate named by the atlas is the orbit-coupled local-time statement (THM-4506).

Mechanism notes: (1)-(2) are the Terras model with the 3-adic integrality taken as independent, which is what Tao's
Proposition 1.14 makes rigorous on average and what Mazur's (2.3) quantifies; (3) is the observation that every
layer of Mazur's argument except the seed lemma and the direction of the source charge is sign-symmetric; (4) is
the fast-side form of the "expense is Diophantine" deadline lattice, now confirmed at quarter-dyadic resolution; (5) is the two-wall verdict of the coalescence atlas read through
the owner's three pillars.
