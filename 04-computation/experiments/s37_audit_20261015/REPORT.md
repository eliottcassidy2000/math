# S37 independent audit report (2026-10-15)

Independent adversarial referee (subagent), audit of the S37 pre-audit checkpoint `378a9b6390`
(THM-4658, HYP-9282, HYP-9283, `collatz_octave_memory_bias_energy_budgets_20261015.md`,
`special_primes_curling_genus_zero_langlands_20261015.md`). The referee could not write files outside its scripts,
so this report was saved by the session from the referee's returned text. All referee code is in this folder,
one `.out` per `.py`, independent of the session's code.

Scope notes:
- `f74c2037e1` (clock_duality_oct15) independently found E1 and repaired statement 4(a) and added a caveat on
  statement 5 while the audit ran.
- `93d1ee3885` (S37 census to 2^37) was outside the audited file list and was not audited.

Method: two independent parity schemes (sequential descent into a byte array; 8-step jump table on parallel dyadic
blocks, using `tau(n) = s + tau(T^s n) mod 2`, valid because the terminal cycle has even length); popcount prefix
tables for exact `H(n)`; every run under 2.5 minutes.

## Errors

**E1. THM-4658 4(a) was circular as typed.** The step "every class sum shares S's main term with Haar weight" is
equivalent to the conclusion. The exact identity `H(y) = -[S(y) + 3S(floor(2y/3)) + 3S(floor(y/3))]` (0 failures for
y < 3 2^22) gives: if `S(N) = c N^rho + o(N^sigma)` then `H(N) = -3c(2/3)^rho D(rho) N^rho + o(N^sigma)`, so
`D(rho) = 0` iff `H = o(N^sigma)`. A single complex term with Im rho != 0 is impossible for real S, and D has no real
zeros (D(1/2) = 2.1154, D(1) = 2). HEAD (clock_duality_oct15) adds the hypothesis; the title still overclaimed
("every power-law main term (density or excess)").

**E2. HYP-9282's independence assumption for the 19 octave means is false.** The quoted 0.06 is the point ACF of g
at one octave, not the correlation of adjacent octave means. Adjacent-octave-mean correlations: H +0.36 (sample),
+0.28 (implied by its point ACF), +0.32 (session JSON); median over 363 channels +0.63 (quartiles +0.49/+0.72); D_j
+0.67..+0.87. With an AR(1) correction (pooled rho = 0.63, effective n = 4.3) only 7 channels have |t| > 5 and 23
have |t| > 3.5 (claimed 75). Only 242 of the 363 channels are free (sibling triples sum to zero). What survives is
sign persistence: 27 channels keep one sign in all 19 octaves, against 0.04, 0.45, 2.2, 9.3 expected under zero-mean
AR(1) nulls with rho = 0.3, 0.5, 0.63, 0.75 - beating short-memory nulls but not slowly varying alternatives.

**E3. "Drift at most about +-2.5% per octave, exponent 1/2 +- 0.04" is false for 5 of the 14 table rows**, all
drifting toward zero: H_(2,3) +9.1% (t +5.0), H_(2,6) -8.2% (-3.0), H_(1,1) +7.6% (+2.2), H_(4,51) +4.8% (+3.6),
H_(3,16) +4.6% (+3.0). Among the 20 channels with |t| > 8, 8 drift by more than 2.5% per octave and 12 drift toward
zero. A 9% drift corresponds to exponent about 0.37. In the octave [2^33, 2^34) (beyond the session's mod-243 data):
pre-registered H_(5,155) flips sign (+0.186); H_(3,18) = -1.04 and H_(3,9) = +0.95, about half their means;
H_(5,59) = +0.16 against +0.92.

**E4. The note's claim of "a deterministic floor ~0.26" under every octave energy is false**: e_X < 0.26 at
k = 10, 11, 12, 16, 28, 30, 33 (e = 0.0954 at k = 33); the octave mean part is 0.0002 at k = 29. A constant main
term contributes 3b^2 only on average across octaves.

## Imprecisions

- I1. THM title "O(X^(3/2))" holds only if E_4X = O(X^2); the body's `|R_X| <= 24 sqrt(X E_4X) + 48X` is correct
  (observed |R|/bound <= 0.007; in fact sum eps^2 <= 29X).
- I2. "Transport agrees to 4 decimals on every octave 2^10..2^30" is false at k = 10, 11, 12, 13, 16 (k = 10: 0.1947
  vs 0.1942); agreement to 1e-6 starts at k = 17.
- I3. kappa = c/sqrt(e delta) is an uncentred cosine, not a correlation.
- I4. Delta = H(4N) - H(N) by definition, so the empirical content of "cancelled memory" is ACF(2 octaves) ~ 0
  (slope 0.0063); the new content is that Delta is a level-2 expression. "Random affine iteration" for
  g -> g/2 + Delta-hat is misleading (Delta-hat ~ -g/2 plus a fresh term).
- I5. Statement 5's literal claims are correct; "balanced / critical / unitarity line" readings are heuristic. Exact
  level-6 norms of P*_rho: 1.398, 1.362, 1.352 at rho = 1/2, 3/4, 1 (about 1.34 on mean-zero functions). The valid
  reading is the variance recursion for independent white-sign branches. Mercator zeros sit at Re s ~ 0.998.
- I6. "Mercator modes cancel exactly from every wavelet channel": proved only for H; with floors the identity is exact
  for every y. For other channels it is numerical (26.5-cycle amplitude <= 8e-5 of that in S, levels 1-3).
- I7. HYP-9282's model (main term plus mean-zero fluctuation of the same order sqrt(N)) violates 4(b)'s o(N^sigma)
  hypothesis; a logarithmic-mean version of 4(b) would be needed and is not stated.
- I8. HYP-9283(A) surrogate probabilities are fragile: P(dyadic E_X > 2X^2 before 2^64) = 32% (reproduced), 1.0% with
  the bias removed (5.9% before 2^200), 8%/58% with innovation sd x0.8/x1.2, 24-34% for AR orders 2-12, 40% fitted on
  u >= 20. 77 quarter-octave samples for 9 parameters; kurtosis SE ~0.95; "all-phase" figures include three windows
  starting inside the census.
- I9. HYP-9283(B) understates its margin: certified budgets reach E_X <= 227.4 X^2 (constant), (4c, c) with c = 36.9,
  and slope b = 144 with a = 2.
- I10. Credit: the wavelet note's one-octave energy rule (3) is the direct analogue of statement 3; its section 2
  already reported cross-term correlations from -0.949 to +0.936 with no Q-constant closure; the amplitude note's
  sections 3 and 6.1 already proposed logarithmic budgets (S37's contributions are the interval certification and the
  Krasikov-Lagarias calibration); T1/T2 are coarse-mean (terminal-controller) certificates, not refinement-energy
  results; HYP-9283(C) is not part of THM-4658.
- I11. Tao's Fourier decay is O_A(n^-A) for every A (stronger than "polynomial"); the "would need exponential decay"
  claim is not derived.
- I12. Companion note: the first odd ancestors of 1 under 3x-1 are (2^m+1)/3 for every odd m (m = 1, 9, 15, 21, ...
  included), a strict superset of the Wagstaff numbers; W_n's initial odd run has length n + v2(n) (n <= 80), not n
  (the landing point n 3^n - 1 after n steps is correct); "S-duality ... with a fixed coupling" is wrong (S maps tau to
  -1/tau; S^2 is charge conjugation).
- I13. Minor: HYP-9282 title omits "rank-one modified"; "the descent bound is attained" at k = 22 compares 1.7688 with
  1.8715.

## Verdicts

| Claim | Verdict |
|---|---|
| S1 class-sum recursion (5.3e8 checks) | CONFIRMED |
| S2 two-octave identity, one-octave relations, identity (2) for 120 channels | CONFIRMED |
| S3 transport, R_X bound, descent, chain bound | CONFIRMED (title fix I1) |
| S4(a) as stated | REFUTED (E1); true with H = o(N^sigma) (fixed in HEAD) |
| S4(b) excess identity and fixed point | CONFIRMED |
| S5 branch norms, V*V = I/3, phi | CONFIRMED (readings overstated, I5) |
| S(2^34) = 24953018, H(2^34) = 26869 | CONFIRMED (third independent census) |
| Octave table and grid class sums (1537 points) | CONFIRMED |
| "Transport to 4 decimals" | CONFIRMED WITH FIXES (I2) |
| HYP-9282 numbers (75/363, 20/363, table rows, D_j) | CONFIRMED (exact reproduction) |
| HYP-9282 significance and independence | REFUTED (E2) |
| HYP-9282 drift and exponent claim | REFUTED (E3) |
| HYP-9282 Mercator cancellation, fixed-point reading | CONFIRMED WITH FIXES (I6, I7) |
| HYP-9283(A) statistics / surrogate probabilities | CONFIRMED / CONFIRMED WITH FIXES (I8) |
| HYP-9283(B) intervals, tail derivation, (15,4) translation | CONFIRMED (bit-identical rerun; mpmath agrees to 1e-29) |
| HYP-9283(C) implication; Krasikov-Lagarias citation | CONFIRMED |
| Note section 4 "deterministic floor" | REFUTED (E4) |
| Non-generous primes; orbit-count corollary (7746 bases; 2^40488 at p = 40487) | CONFIRMED |
| Wagstaff / Woodall / Cullen | CONFIRMED WITH FIXES (I12) |
| Leyland, supersingular/Paley, Stern primes, curling records, Gijswijt's first 4 | CONFIRMED |
| Langlands / S-duality citations | CONFIRMED WITH FIXES (I12) |
| Pre-registered tests 2^34..2^40; checkpoint 93d1ee3885 | NOT CHECKED |

## Scripts in this folder

`a1_thm4658_identities.py` (statements 1-5, S4(a) identity), `a2_census34.py` (independent census to 2^34),
`a3_compare_census.py`, `a4_bias_recompute.py`, `a4b_bias_sign_persistence.py`, `a5_ar_surrogate_review.py`,
`a6_cert_rerun.json/.out`, `a6_logbudget_tail_check.py`, `a7_special_primes.py`, `a7b_orbit40487.py`,
`a8_curling_check.py`, `a9_transfer_operator_norms.py`.

Sources cited by the referee: OEIS A094004; Krasikov-Lagarias, Acta Arith. 109 (2003) 237-258; arXiv math/0205002.

## Session response (applied 2026-10-15)

All errors and imprecisions were applied to THM-4658, HYP-9282, HYP-9283 and both notes; the ledger entry is in
`01-canon/MISTAKES.md`. Independently of the audit, the session ran the pre-registered out-of-sample test (census to
2^37, then 2^40) and a matched halving null; HYP-9282 failed its pre-registered test and is REFUTED, and the
HYP-9283(A) surrogate failed its out-of-sample check (energies far below the predictive distribution).
