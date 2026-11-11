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
