# Pre-registration (mac-mini-2026-11-11-race)

Written and committed before the runs named below produced any output. Notation: `D_n = P(H = 1 mod 4) - P(H = 3 mod 4)`
over uniformly random labelled tournaments of order `n` (THM-5127: exact for `n <= 11`); `b_R`, `beta_R`, `rho_m` as in
THM-5127 (C); `c_k` = number of directed `k`-cycles.

In-sample facts used to form the predictions (all exact, computed today or in THM-5127):

- `D_n 4^n = 1706, 2722, 1972, 2541, 1988` at `n = 7, ..., 11`; after removing the self-converse tournaments
  (all `= 1 mod 4` at even order, THM-791) and the non-self-converse tournaments with a non-trivial automorphism,
  the rigid remainder is `1417.5, 1181.2, 1965.7, 1653.3` at `n = 7, 8, 9, 10`;
- labelled rank law of `beta_R` at `m = 9, 10` against the random alternating law: differences `0.011`, `0.0008`;
- `rho_m = 1.244, 1.156, 1.041` at `m = 8, 9, 10`; non-zero-bias frequency at corank 1: `1, 0.815, 0.636, 0.527` at
  `m = 3, 5, 7, 9`;
- `E (-1)^(c_3) = 2^(-floor((n-1)/2))` (proved today for all `n`); at `n = 6` and `n = 8` the full race is unchanged
  when `c_5` is dropped: `E (-1)^(c_3 + c_5) = E (-1)^(c_3) = 1/4` at `n = 6`,
  `E (-1)^(c_3 + c_5 + c_7) = E (-1)^(c_3 + c_7) = 1361/32768` at `n = 8`; `E (-1)^(c_5) = 11/256` and
  `E (-1)^(c_3 + c_5) = 67/256` at both `n = 7` and `n = 8`.

## Run A: order 12 by Monte Carlo with exact averaging over the last vertex (`biasmc.c`)

`N = 14,400,000` uniformly random labelled tournaments `R` of order 11 (six processes, seeds 101..106, 2,400,000
samples each); for each, the exact Gauss sum over the `2^11` out-sets. The estimator of `D_12` is the mean of
`eps_R b_R`; its standard error is reported from the sample.

- **R1 (generic lead).** `D_12 > 0` and `D_12 4^12` lies in `[1500, 2500]`. Verdict rule: CONFIRMED if the point
  estimate is in the interval and more than 3 standard errors above 0; FAILED if it is outside the interval by more
  than 2 standard errors; otherwise UNDECIDED.
- **R2 (rank law at order 11; the falsifier named in HYP-9659 (B)).** The fraction of `R` with `rank(beta_R) = 10` is
  within `0.006` of `0.83898`, and the fraction with rank at most 6 lies in `[0.0012, 0.0030]` (random alternating
  matrix with an odd number of edges: `0.83898 / 0.15965 / 0.00137`).
- **R3 (deletion-pair correlation).** `rho_11 = 2^11 E[b_R^2]` lies in `[1.015, 1.045]`.
- **R4.** The frequency of non-zero bias among forms of rank 10 (corank 1) lies in `[0.500, 0.512]`.
- **Control (must hold by THM-5127 (C)).** No odd `mu`, no odd `a_i`, Gauss sum `0` or `+-2^(11 - rank/2)` with the
  radical criterion, and an odd number of edges of `beta_R`, in every sample.

## Run B: parity races by cycle length at order 10 (`cycrace.c`)

The run over all 9,733,056 classes of order 10 was started at 08:02 MST, before this file was committed; the program
prints only when it ends, and no output existed when this text was written (the six output files were empty).

- **R5 (the pentagon is invisible at even order).** `E (-1)^(c_3 + c_7 + c_9) = E (-1)^(c_3 + c_5 + c_7 + c_9) = D_10
  = 650471/268435456` exactly.
- **R6 (pair law for short cycles).** `E (-1)^(c_5) = 23/2048` and `E (-1)^(c_3 + c_5) = 263/2048` at `n = 10` (their
  values at `n = 9`).
- **Control.** `E (-1)^(c_3) = 1/16`; `c_3` equals the score formula on every class; the full race reproduces `D_10`.

## Refutation

Any of R1-R6 outside its statement is reported as FAILED (R1 by its verdict rule). A failure of a control refutes
the program or the theorem it checks.

---

# Second pre-registration (same session): R7, R8

Written and committed after the outcome of Run B was read (R5 and R6 FAILED: the pentagon is visible at `n = 10`,
and `E[{5}]`, `E[{3,5}]` change from `n = 9` to `n = 10`) and before Run C produced any output. Run A was still
running when this was written.

In-sample (exact, `out/cycle_tables.out`): `E[{3,5}] = E (-1)^(c_3 + c_5)` equals `13/32, 1/4, 67/256, 67/256,
263/2048, 1/16` at `n = 5, ..., 10`. Two patterns survive the failure of R5-R6: at `n = 6` and `n = 10`
(`n = 2 mod 4`) the value is `2^(-(n-2)/2)`, the triangle value; and from `n = 7` to `n = 8` (`3 -> 0 mod 4`) it does
not change.

## Run C: plain Monte Carlo for the triangle and pentagon parities (`c35mc.c`)

`N = 40,000,000` uniformly random labelled tournaments at each order `n = 11, ..., 18` (seed `1000 + n`).

- **R7 (orders `2 mod 4`).** `E_14[{3,5}] = 2^(-6)` and `E_18[{3,5}] = 2^(-8)`, each within 3 standard errors of
  the estimate (`s.e. = 1.6 * 10^(-4)`). The alternative suggested by the neighbouring orders is about twice as
  large.
- **R8 (pair law `3 -> 0 mod 4`).** `E_12[{3,5}] = E_11[{3,5}]` and `E_16[{3,5}] = E_15[{3,5}]`, each difference within
  3 standard errors (`2.2 * 10^(-4)`).
- **Control.** `E_n[{3}] = 2^(-floor((n-1)/2))` within 3 standard errors at every order (Theorem T).

Exploratory, not registered: the values at `n = 13` and `n = 17`, and `E_n[{5}]`.

---

# Outcomes (appended after the runs; the text above is unchanged)

| | prediction | outcome | verdict |
|---|---|---|---|
| R1 | `D_12 > 0`, `D_12 4^12` in `[1500, 2500]` | `D_12 = (5.892 +- 0.59) * 10^(-5)`, ten standard errors above 0; `D_12 4^12 = 988.5 +- 98.5`, 5.2 standard errors below the interval | **FAILED** (the sign holds, the plateau does not) |
| R2 | rank 10 within `0.006` of `0.83898`; rank at most 6 in `[0.0012, 0.0030]` | `0.838685 +- 0.000097`; `0.001398` | **CONFIRMED** |
| R3 | `rho_11` in `[1.015, 1.045]` | `1.0172` | **CONFIRMED** |
| R4 | corank-1 bias frequency in `[0.500, 0.512]` | `0.50399 +- 0.00014` | **CONFIRMED** |
| Run A control | no failure | none in 14,400,000 samples | holds |
| R5 | `E_10[{3,7,9}] = D_10 = 650471/2^28` | `533081/2^28` | **FAILED** |
| R6 | `E_10[{5}] = 23/2048`, `E_10[{3,5}] = 263/2048` | `1/256`, `1/16` | **FAILED** |
| Run B control | triangle column, score formula, full race | all hold | holds |
| R7 | `E_14[{3,5}] = 2^(-6)`, `E_18[{3,5}] = 2^(-8)` within 3 s.e. | `0.015803` (`+1.13` s.e.), `0.003824` (`-0.52` s.e.) | **CONFIRMED** |
| R8 | `E_12 = E_11`, `E_16 = E_15` for `{3,5}` within 3 s.e. | differences `+0.73`, `+0.60` s.e. | **CONFIRMED** |
| Run C control | `E_n[{3}]` within 3 s.e. of Theorem T at every order | largest deviation `+2.7` s.e. (`n = 13`) | holds |

Files: `out/mc11_summary.out` (Run A), `out/cycle_tables.out` (Run B), `out/c35_summary.out` (Run C).

Not registered, done afterwards:

- An estimator check for Run A (the same program at order 10, where the answer is known exactly):
  `D_11 4^11 = 1981.5 +- 66.8` against `1987.8` (`out/mc10_control_summary.out`).
- A sharper test of R7 and R8 with the exact vertex-cube sum of `c_3 + c_5` (`k35cube.c`, 4,000,000 bases per order
  `m = 11, ..., 40`; relative precision `3 * 10^(-4)`): agreement at every order (`out/k35_summary.out`).
- R7 and R8 are now theorems (THM-5132 (5); README section 6): the values are `2^(-(n-2)/2)` for every
  `n = 2 (mod 4)`, and `E_n[{3,5}] = E_(n+1)[{3,5}]` for every `n = 3 (mod 4)`. The exact values at the sampled orders
  are in `out/c35_summary.out`; all 16 estimates of Run C are within 2.3 standard errors of them.

## Erratum and audit notes (appended after the independent audit; the registered text above is unchanged)

- **Erratum (audit m22).** The first block says that Run B "was started at 08:02 MST". File times show that the six
  parts of Run B were launched at 08:06:07; the run launched at 08:02:02 was the self-converse census of order 10
  (in-sample). Run B wrote its outputs at 08:14:35-08:16:20, after the commit of the first block (08:09:37), so the
  statement that no output existed when the block was written stands.
- **Integrity (audit section C).** Both blocks were committed and visible to other agents before the outputs of their
  runs existed in the repository; the registered text is byte-identical in all later versions; `biasmc.c`,
  `cycrace.c`, `c35mc.c` were not changed after registration; every verdict follows its registered rule. Not
  verifiable: that no earlier unrecorded run was looked at.
- **R1.** An independent estimator written by the auditor (9,000,000 bases of order 11, other program and seeds)
  gives `D_12 4^12 = 1039 +- 125`; pooled with Run A, `1008 +- 77`. FAILED stands.
- **R2 (audit m17).** CONFIRMED by its rule, but the tolerance `0.006` was sixty standard errors wide; the observed
  fractions are `-3.1` and `+2.8` standard errors from the reference law.
