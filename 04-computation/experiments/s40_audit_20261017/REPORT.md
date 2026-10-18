# Independent audit of session S40 (opus-2026-10-17-S40)

2026-10-17. Independent referee subagent; its report is transcribed here (the subagent may not write report files). State
audited: worktree math-wt-explicit-modal-38 at HEAD 8933902975 (17:58 MDT).

Files under audit, unchanged since these commits:
- THM-4666 and HYP-9287: 309af1b454.
- Main note: 96f27780ca.
- ROOT note: fb6923b86f.
- Pre-registration file: 9db25e9068.

All audit code is in this folder and was written without the session's code. The two exceptions are a sandboxed rerun of
the session's search and the use of clock_duality's field class for the rotation-field check. No tracked file was
modified.

## Overall verdict

- **THM-4666: correct, with citation fixes.**
  - Every statement and every FINITE-EXACT number reproduces with independent code. Class numbers were computed by two
    methods (a reduced-form count and an erfc series for `L(1, χ)`, agreeing to 5.7e-12). MAS(P_31) = 285 was
    recomputed by a different DP with a different symmetry reduction.
  - Kalmynin proves only an upper bound for the counting function of L+, and the density statement for primes is his
    Thm 1, not Baker-Montgomery's.
  - Repository prior art is not credited.
- **HYP-9287: refuted as stated.**
  - The pre-registration is clean: committed and pushed before the run, both JSONs reproduce byte-for-byte, and all 50
    witness orders recompute.
  - The search was too weak. Explicit rotation orders (rank `x` by `frac(αx)`) beat the natural order at 167, 263 and
    383, three of the nine P1 primes reported unbeaten. With certified orders 4 of 10 P1 primes are beaten, and the
    pre-registration's own rule then refutes the threshold form.
  - The second pre-registration fails too: 503 and 647 are beaten.
- **Main note: computations correct, text needs fixes.**
  - Weinberger's bound is "at most two" (Kani).
  - EKN cover `e = 4` and `e = 7`.
  - Gauss's cubic and quartic criteria were misquoted and are false as written.
  - The Collatz row overclaimed.
  - The 1848 tower has sixteen self-inverse classes only for `k >= 1`.
  - The idoneal status ignored the repository's adopted THM-4566.
  - Sections 3 and 6 restated clock_duality's escape_renewal section 7 without credit.
- **ROOT note: numbers exact; two prose errors.** Anti-correlation holds at 17 of 19 octaves, not 16, and `2^12` is not
  a cancellation failure. "All 34 rotations lie in K" is true, but the session's script does not test it.

## Verdict table (summary)

| Item | Verdict |
|---|---|
| THM-4666 (1) excess, shifts, multipliers | CONFIRMED: 6079 fundamental D, 1500 random D in [-10^6, -2·10^4), 9006 primes <= 2·10^5 and 400 random primes in (2·10^5, 10^6]; 30,388 multiplier/shift orders; zero failures |
| (2) half window | CONFIRMED (9405 primes) |
| (3) split primes = arcs | CONFIRMED (42,900 pairs); the h-series needs p > 3 |
| (4) h = 1 list, Rabinowitsch, lucky numbers | CONFIRMED |
| (5) insertion = -2S(t); symmetry; local optimality iff L+ | CONFIRMED (18,712 moves at 10 primes <= 71; 16 primes 7..131; symmetry at 9005 primes) |
| (5) Fekete identity and Mellin positivity | CONFIRMED (polynomial identity for p < 3000; Mellin at three s values, four primes, differences <= 2.3e-12) |
| (5) Kalmynin "counting function" | WRONG (citation): Mathematika 69 (2023), Thm 2 is only `<< π(x)/(log log x)^c`, `c ≈ 0.0368`; whether L+ is infinite is unknown |
| (5) Baker-Montgomery density 0 | CONFIRMED with fix: their theorem is for almost all fundamental discriminants; the prime statement is Kalmynin's Thm 1 |
| (6) eigenvalues, `‖S‖ = sqrt p`, trace duality, asymptotic | CONFIRMED (eigenvalues to 7e-14 for n <= 64; asymptotic to 7e-13 at n = 10^5; optimum/ceiling 0.40-0.50 for p <= 43) |
| (6) comparison with Alon | CONFIRMED with wording fix: Alon's Cor. 3.2 is also explicit (`p^(3/2) ceil(log2 p)`), 3.4-4.1x weaker; the classical bound is Alon-Spencer section 9.1 |
| (7) exact MAS 2, 14, 35, 107, 161, 285 | CONFIRMED (DP with vertex 0 last for p <= 23; brute force for p <= 11; p = 31 with 0 last and -1 penultimate, 2^29 states, 108 s) |
| (7) L+ counts and bins | CONFIRMED |
| (7) MAS(P_43) = 543 | CONFIRMED as cited (arXiv:2607.26690 Thm 7.1, exact DP modulo the automorphism group; unrefereed); the excess-183 witness verifies, so non-optimality at 43 needs no citation |
| Joshi remark typed CONDITIONAL | can be upgraded: Joshi's Thm 2 (via Brink-Moree-Osburn, Thm 3.1) covers q ≡ 7 (mod 8) |
| Pre-registration integrity | CONFIRMED: ee91f70219 committed 16:42:31, pushed 16:42:43, run output written 16:47:37, results committed 16:48:38; byte-identical reruns |
| Witness orders and P1/P2 verdicts as executed | CONFIRMED (50 of 50) |
| L(1,χ) values | CONFIRMED except p = 3: `L(1, χ_(-3)) = π/(3 sqrt 3) = 0.605`, not 1.81 (missing w = 6) |
| HYP-9287 threshold claims | **WRONG**: rotation-order certificates below |
| "never optimal for h = 1, p > 7" | CONFIRMED; PROVED (finite): 11, 19 exact; 43, 67, 163 witnesses of excess 183, 393, 1577 |
| Note section 1 literature table | CONFIRMED WITH FIXES (Weinberger: at most two; the Collatz row; everything else verified) |
| Note section 5 exponent ladder | CONFIRMED (independent composition at 3·10^5 and 10^6; rows e = 1, 2, 3, 4, 5, 7, 9 match; group structures C3², C5², C7², C9², C2³×C4²) |
| Idoneal class-number histogram | CONFIRMED |
| Note section 2: ring class field of 1848 | CONFIRMED (all 20,639 primes `x² + 1848y²` below 10^7 match the five Legendre symbols) |
| Note section 0.5: cubic and quartic criteria | **WRONG**: missing `p ≡ 1 (mod 3)` and `p ≡ 1 (mod 4)` (counterexamples p = 5, 7) |
| Note section 3: dictionary, Cox 3.22(iii), T(n), SC(n) | CONFIRMED (Burnside recount) |
| Note sections 6-7 | CONFIRMED WITH FIXES (tower k >= 1; μ14 occurrence of 463 omitted) |
| Idoneal status vs THM-4566 | CONFLICT (fixed: both layers now stated) |
| ROOT section 1 numbers | CONFIRMED (k = 8..24, integer numerators, signs, extremes; d(n) spot-checked at 2008 values) |
| ROOT section 2 fold | CONFIRMED (zero failures; inner products exact) |
| ROOT section 2 prose | WRONG in two places (17 of 19; the 2^12 remark) |
| ROOT section 3 calibration | CONFIRMED WITH FIXES (range 4 <= X; the step from the current energy D_X to E_X) |
| ROOT section 4 rotations in K | claim CONFIRMED (34 rotations, zero w-coordinates, |U| = 1, six roots of unity); not tested by the session's script |

### Certificates against HYP-9287

Orders list `x ∈ {0..p-1}` sorted by `(x·num) mod den`; every excess was recomputed in pure Python.

| p | L | natural `p h` | α | rotation excess | best known |
|---|---|---|---|---|---|
| 167 | 2.674 | 1837 | 7599/19912 | 1861 | 1863 |
| 263 | 2.518 | 3419 | 57115/72812 | 4079 | 4089 |
| 383 | 2.729 | 6511 | 38609/72644 | 7519 | 7529 |
| 439 | 2.249 | 6585 | 27601/171580 | 9313 | 9313 |
| 503 | 2.942 | 10563 | 9521/197064 | 11657 | - |
| 647 | 2.841 | 14881 | 10901/258330 | 17661 | - |

For `p <= 500` the natural order is beaten by no known order only at 3, 7, 23 (exact optima) and 47, 71, 191, 239, 311,
359, 431, 479 (all with `L >= 1.97`), after the full rotation family, insertion and window-DP polishing, and simulated
annealing with 480 restarts at p = 47 and 71. The session's later scan reports 887 and 911 beaten as well; these were not
verified by the audit.

## Fixes applied (2026-10-17, S40)

- **THM-4666.**
  - The Kalmynin and Baker-Montgomery citations; Alon-Spencer section 9.1 and Alon's Cor. 3.2.
  - `p > 3` in statement 3.
  - The preprint status of MAS(P_43), with the witness.
  - Joshi upgraded to CITED via Brink-Moree-Osburn.
  - The prime463 gloss.
  - Repository prior art credited.
  - The HYP-9287 refutation recorded.
  - Also the new statements 1b (arithmetic orders: inverse order 3h proved by a marked-coefficient argument, pair-block
    orders, residue sub-tournament) and 5b (randomised moves = Fekete positivity), which postdate the audit.
- **HYP-9287.** Rewritten: REFUTED in the threshold form, with the certificate table, the revised conjecture and the
  calibration rule for "unbeaten" claims.
- **Pre-registration file.** `L(1, χ_(-3))` corrected; audit addendum for the first test; the second-test outcome (Q1
  FAILED, Q2 PASSED).
- **Main note.**
  - Both idoneal layers stated: classical (Weinberger, Kani "at most two") and adopted (THM-4566, modulo #003).
  - Gauss's criteria corrected; the Collatz row; EKN coverage for e = 4, 6, 7, 8; homocyclic; the tower for `k >= 1`.
  - THM-1935 and THM-4658 paraphrases; the μ14 row.
  - Credits to escape_renewal section 7, trees5 and THM-871's scope repair.
  - Euler's criterion for 463 (new).
- **ROOT note.** 17 of 19; the 2^12 remark; the range `4 <= X`; the D_X to E_X step; the rotation-field check
  attribution; root-of-unity blocks; the certificate floor (new).

## MISTAKE candidates (logged as MISTAKE-626)

1. HYP-9287's threshold form refuted by rotation orders. The cause was a light, uncalibrated insertion search that never
   tries structured orders.
2. Kalmynin's upper bound reported as the counting function; the density for primes attributed to Baker-Montgomery.
3. Uncredited repository prior art.
4. The idoneal status stated against THM-4566 without comment, and Weinberger's bound misreported.
5. Gauss's cubic and quartic criteria misquoted.

## Prior art found in the repository

- **Class number as a Paley half-row imbalance:** ogg_triangular_chowla_20260927.md section 2 (opus S17, audited
  p < 2000) and special_primes_curling_genus_zero_langlands_20261015.md section 4 (opus S37).
- **`R - N = Σ (a/p) a = -p h(-p)`:** collatz_mod6_20260917_floor_sums_odd_functions.md S5 (mac-mini).
- **Both numbers at p = 463:** prime463_conductor_clock_20261017e.md P4.
- **Small Paley rankings:**
  - MFAS(P_7) = 7 with the natural order optimal (tournament_heptagon_flipdepth_kps.out);
  - MFAS(Paley_11) <= 20 (tournament_covering_code_bounds_opus_20260701.out);
  - Legendre partial sums as interval-tournament alignment (qr_character_eigenvector.out, HYP-488).
- **Idoneal and exponent-2:** THM-4566 and the oai3 reader exp2_analysis.out (the e = 2 row); golden_tail_idoneal
  section 6 (Kani; "do not regress the inherited status").
- **clock_duality's escape_renewal section 7:** the orbit law `(|S|+|Fix|)/2`, the sixteen rotational classes against
  `Cl(-7392)`, and the sedenion contrast.
- **Merged = (A000568 + A002785)/2:** trees5_tournaments4_converse_20261005.md.
- **Euler's lucky numbers and class number one:** HYP-2226, HYP-3730.
- **Not found anywhere earlier:** the insertion/walk dictionary, L+ and the Fekete link to tournaments, the exponent
  lists for e = 3, 5, 7, 9, and the idoneal-order class-number histogram.

## Audit scripts (this folder)

- `audit_common.py`; `audit_thm4666_identities.py`; `audit_walk_insertion.py`; `audit_mellin.py`; `audit_mas_exact.py`
  (with the `_p31` outputs); `audit_spectral_bound.py`.
- `audit_witnesses.py`; `rerun/` (byte-identical reruns); `audit_search.py` (with its outputs).
- `audit_rotation_witnesses.py` (with `allp500/`, `q1_second/`); `verify_witness_purepython.py`; `verify_rotation_alpha.py`;
  `summarize_rotation_vs_session.py`.
- `audit_exponent_ladder.py`; `audit_ring_class_1848.py`; `audit_tournament_counts.py`.
- `audit_root_current.py`; `audit_rotations_field.py`.
