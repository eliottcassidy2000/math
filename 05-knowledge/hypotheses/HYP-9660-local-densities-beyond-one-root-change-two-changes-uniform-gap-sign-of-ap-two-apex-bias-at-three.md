---
id: HYP-9660
title: "Local densities beyond THM-5128: (A) by-value limit of the two-root-change family, L(0) = 1 + (p-1)(sigma p - 1)/p^3 with sigma = 2 w_1 = 1.60718 instead of the letter-count 3/2; (B) a uniform spectral gap for the twisted digit operators, no rate above 0.93204 (the rate at 7) at any odd modulus (KNOWN for moduli up to 120, S82; scanned here to 521); (C) the sign of the first-order constants A_p; (D) the two-apex grammar G(2,0) shows no local bias at 5, 7, 11, 13 (factors within 0.0004 of 1 at n = 16) and little at 2 (within 0.0032); at 3 all three classes are biased, by up to 0.7 % at n = 16, and only the class-0 part decays for n = 12..16; (E) the coprime states have density 8/pi^2"
status: >
  OPEN (opus-2026-11-10-S91, 2026-11-11). Independent audit 2026-11-11 (seven auditors; record
  05-knowledge/results/s91_f1_local_densities_20261110/audit/AUDIT_REPORT.md): PASS WITH CORRECTIONS, all applied;
  MISTAKE-761. For this file the audit compared the numbers with the outputs and with prior art. It did NOT check the
  helper derivation of (A), and above m = 101 it recomputed the scan of (B) only at the primes up to 127, at 281
  and at 1039.
  (A) is a derivation by a helper of this session whose exact ingredients (cascade identities X1-X4 of
  census_k2/cascade_exact.py) are FINITE-EXACT, but whose limit is NOT confirmed by the census: on [2^23, 2^24) the
  measured ratios exceed both candidate limits at all eight primes tested (3, 5, 7, 11, 13, 17, 23, 31). The excess
  over the conjectured limit is 0.23 to 0.44 at 3, 5, 11, 17, the primes with a balanced divisible state at
  h <= 16, where THM-5128 (D) gives a positive 1/ln X term for one change, and 0.02 to 0.06 at 7, 13, 23, 31; at 13
  and 31 the one-change term does not explain it (it is negative there; at 7 its sign is undetermined). On the
  later windows up to 2^26 the ratios stay above both candidates, and at 7 and 23 they still rise. The derivation
  is unaudited; typed CONJECTURE.
  (B) CONJECTURE. KNOWN for every a/q with q <= 120 (S82 note section 7, twisted_scan_q120.out: rate at most 0.9320,
  attained only at +-1/7 mod 1/2; next 11/93, 1/31, 1/15, 1/63), where the failure of the worst-prime-power rule is
  already visible (93, 15). NEW here: the extension to every odd m <= 521 and to m = 1023, 1039, 1143 (VERIFIED
  numerically, floating point eigenvalues; receipts rates_exact/rates_by_modulus.txt and, for 1039,
  bridges/b3_f1_data.out): no rate above 0.932039, attained exactly at the multiples of 7; further counterexamples
  (45, 51, 65, 75, 85, 381, and others above 101 listed in (B)) to the rule "a composite modulus has the rate of
  its worst prime power" (REFUTED); and the statement as a conjecture. THM-5128 gives the gap at each fixed modulus
  only: uniformity in the modulus is OPEN.
  (C), (E) are questions with NUMERICAL data; in (C) the sign of A_7 is undetermined. (D) is FINITE-EXACT data for
  all 4^n words, n <= 16, with no theory. The census limit for two root changes is the subject of HYP-9506 (A)
  (kind-pasteur; it conjectures the letter-count limit), OPEN; THM-5128 settles one root change only, where the
  by-value limit differs from the letter-count limit in the split.
session: opus-2026-11-10-S91
date: 2026-11-11
source: 05-knowledge/results/s91_f1_local_densities_20261110/README.md
scripts:
  - 04-computation/experiments/s91_f1_local_densities_20261110/census_k2/ (census_k1.py, census_k2.py, cascade_exact.py, dict_k1.py and their .out)
  - 04-computation/experiments/s91_f1_local_densities_20261110/rates_exact/rates_scan.py (+ scan_003_101.npz, scan_103_521.npz, scan_103_521.out, scan_1023.npz, scan_1143.npz), twisted_spectrum.py, bridges/b3_f1_data.py (m = 1039)
  - 04-computation/experiments/s91_f1_local_densities_20261110/rates_exact/rates_by_modulus.txt (text table, the rate of each of the 262 moduli), bridges/b3_f1_data.out (m = 1039)
  - 04-computation/experiments/s91_f1_local_densities_20261110/ap_constants.py, divisible_states.py (+ .out each)
  - 04-computation/experiments/s91_f1_local_densities_20261110/g20_local/full_hist.py, full_hist_report.py (+ .out)
  - 04-computation/experiments/chanpak_mass_exploratory_20261106/twisted_scan.py, twisted_scan_q120.out (S82, the rates for q <= 120)
  - 05-knowledge/results/s91_f1_local_densities_20261110/audit/recompute-5128/ (a5_Ap.out, a5b_Ap_extended_*.out, a4_divisible.out, a9b_rate_power.out)
  - 05-knowledge/results/s91_f1_local_densities_20261110/audit/numbers-vs-outputs/ (check_g20_mod3.out, check_rates_npz.out, verify_misc.out, check_rate_281.out)
  - 05-knowledge/results/s91_f1_local_densities_20261110/audit/proof-5128-AB/ (check_B1.out, check_B24.out)
  - 05-knowledge/results/s91_f1_local_densities_20261110/audit/verify-corrections/HYP-9660/ (verify_hyp9660.py and .out, the corrected numbers and the elementary steps of (B); verify_exceptions_above_101.py and .out, six moduli above 101 from the full operators)
related:
  - 01-canon/theorems/THM-5128-* (the K = 1 theory), THM-4916-*, 05-knowledge/hypotheses/HYP-9506-* (letter-count model, census OPEN)
  - 01-canon/theorems/THM-4878-* (grammars G(k, w)), THM-4885-*, THM-4886-* ((A1), (A3), a modulus-uniform bound for the count on full blocks)
  - 05-knowledge/results/one_apex_join_representation_mass_pi2_over_12_20261106.md (S82 note; section 7 has the rates for q <= 120)
  - 05-knowledge/hypotheses/HYP-9661-* ((D) quotes (B) as the uniform congruence gap that is missing)
---

# HYP-9660 - local densities beyond one root change

Notation of THM-5128: one-apex states `(x(h), y(h))`, digit chain on `(h mod m, h mod 2)`, twisted operators
`T_k`, rates `max_k rho(T_k)/2`, word-pair sheet shares `w_0 = 0.19640777`, `w_1 = 0.80359223`. The rate of a
modulus `m` is the maximum over all `k != 0 (mod m)`, prime to `m` or not. Receipts are under
`04-computation/experiments/s91_f1_local_densities_20261110/` unless a path says otherwise; `audit/...` is
`05-knowledge/results/s91_f1_local_densities_20261110/audit/...`.

## (A) Two root changes, by value (CONJECTURE; helper derivation, not audited, not resolved numerically)

For the family with two root changes (three one-apex blocks; the `K = 2` words of the census of HYP-9506), ordered
by value, the limiting class factors at an odd prime `p` are conjectured to be

    L(0) = 1 + (p - 1)(sigma p - 1)/p^3,     L(c) = 1 - (sigma p - 1)/p^3  (c != 0),     sigma = 2 w_1 = 1.60718...

The letter-count model has the same shape with `sigma = 3/2` (THM-4916, HYP-9506 (A)). The mechanism is a cascade
(identities X1-X4 of cascade_exact.py, exact on 61,504 words): a `p`-divisible first block either stays confined
at the second join or resets to `(D_2, U_2) = (-1, 0)`, according to one binary feature of the first join whose
letter-count share is `3/4`; the conjecture is that by value this share becomes `w_1`, as in THM-5128 (D3) for
one change.

Evidence and its limit. The census of `[2^23, 2^24)` reproduces HYP-9506 (B)'s numbers and gives ratios
`L(0)/L(c != 0)` of `1.857, 1.614, 1.240, 1.374, 1.171, 1.534, 1.124, 1.102` at `p = 3, 5, 7, 11, 13, 17, 23, 31`;
the conjectured limits are `1.495, 1.298, 1.216, 1.140, 1.119, 1.092, 1.068, 1.051` (letter-count: `1.447, 1.274,
1.199, 1.130, 1.110, 1.085, 1.064, 1.047`); receipt census_k2/census_k2_23_24.out. The measured ratios exceed both
candidate limits at all eight primes. The excess over the conjectured limit is largest at `3, 5, 11, 17`
(`0.36, 0.32, 0.23, 0.44`), the primes with a balanced divisible state at `h <= 16` (`x = y = 3, 5, 33, 17`;
THM-5128 (F)), and smallest at `7, 13, 23, 31` (`0.02, 0.05, 0.06, 0.05`; 13 has its balanced state only at
`h = 32`). At `3, 5, 11, 17` this is the direction of the `1/ln X` term of THM-5128 (D2), whose constants `A_p` are
positive there for one change. At `13` and `31` the one-change term does not account for it: `A_p < 0` there, and
the excess here is still positive. The first-order term for two changes is not derived; the same census shows a
second finite-size effect at every prime: the seed confinement at the second root change without confinement at
the first has `1.4` to `6.7` times the share `(p - 1)/p^3` that both candidate limits assign to it (columns
`P(S2,-C1)` and `(p-1)/p^3` of the cascade table of the same output). At this size the excess hides the difference
between the two candidates (`0.05` at 3, `0.02` at 5 and 7, at most `0.01` from 11 on). The two later windows of
the same census do not resolve it: on `[2^25, 2^26)` the ratios are `1.842, 1.594, 1.257, 1.361, 1.169, 1.501,
1.130, 1.100` (census_k2/census_k2_25_26.out; `[2^24, 2^25)` lies in between), still above both candidates at
every prime; from `[2^23, 2^24)` on they fall at `3, 5, 11, 13, 17, 31` and rise at `7` and `23`, that is, away
from both candidates there. The census also measures the binary feature itself, at the first join and at the
second (`kappa` and `lambda` of census_k2.py: letter-count share `3/4` each, and the derivation takes the by-value
limit `w_1 = 0.8036` for both): `0.8041` and `0.7792` on `[2^23, 2^24)`, `0.8040` and `0.7813` on `[2^25, 2^26)`.
The first agrees with the conjectured value to `0.0005`; the second lies between `3/4` and `w_1` and rises by about
`0.001` per octave (`0.7535` on `[2^12, 2^13)`): at this size the census supports one of the two shares of the
mechanism and leaves the other open.

Decisive test: derive the first-order term for two changes (the same hyperbola argument, with the constant of the
first block), then pre-register the class factors of a window above `2^27`.

## (B) A uniform gap (CONJECTURE; KNOWN for q <= 120, S82; VERIFIED numerically in the stated range)

For every odd `m > 1` and every `k` not divisible by `m`, `rho(T_k)/2 <= 0.93204`, the rate of `m = 7`.

KNOWN (S82 note `05-knowledge/results/one_apex_join_representation_mass_pi2_over_12_20261106.md`, section 7;
receipt `04-computation/experiments/chanpak_mass_exploratory_20261106/twisted_scan_q120.out`): over all `a/q` with
`q <= 120` and `theta = a/q` not in `(1/2)Z` the rate is at most `0.9320`, attained only at `theta = +-1/7
(mod 1/2)`; the next worst are `11/93` (`0.9137`), `1/31` (`0.9002`), `1/15` (`0.8922`) and `1/63` (`0.8895`).
For odd `m` the character `k` of `T_k` is `a = -k/2 (mod m)` there, a reduced fraction being a `k` prime to `m`;
so the characters prime to the modulus have the same rates there for every odd modulus up to 119, and the failure
of the worst-prime-power rule (below) is already visible in that list: 93 above 31, 15 above 3 and 5. Since S82
lists reduced fractions, `1/63` has `0.8895` there, while the modulus 63, whose rate is taken over all `k != 0`,
has the rate `0.932039` of its divisor 7 (audit/numbers-vs-outputs/verify_misc.out).

NEW here: the extension of the scan to every odd `m <= 521` and to `1023, 1039, 1143`; further exceptions to that
rule: `45, 51, 65, 75, 85` (inside S82's range, not in its printed list), `381` and the others above 101 listed
below; the statement as a conjecture.
Scanned completely: the 260 odd moduli `3 <= m <= 521` (rates_exact/scan_003_101.npz and, for `103 <= m <= 201`,
scan_103_521.npz: dense eigenvalues of every operator; for `203 <= m <= 521`, scan_103_521.npz: power-iteration
estimates of every operator and dense eigenvalues of those whose estimate is within `0.012` of the largest for the
modulus or at least `0.835`; 65,310 operators and 8,342 dense computations in `[103, 521]`, the largest difference
between estimate and dense value being `0.002`), and `m = 1039` (`0.7636`, bridges/b3_f1_data.out). Receipt for
every rate quoted here: rates_exact/rates_by_modulus.txt (262 moduli: every odd `m <= 521`, and 1023, 1143). For
`m >= 203` the maximality of the dense value rests on the accuracy of the estimates, calibrated on the dense set
(at `205, 219, 381, 395, 445` a recomputation with dense eigenvalues at every `k` gives the same maxima, see below).
The operators of the scan are the two `m x m` blocks of `T_k` under the complement symmetry (elementary: with
`w = zeta^(k/2)`, `(T_k u)(eta, beta) = u(2 eta, 0) + w^(eta + beta) u(2 eta + 1, 1)` commutes with
`(C u)(eta, beta) = w^(beta - eta) u(-1 - eta, 1 - beta)`, `C^2 = zeta^k`, and on the eigenspace `C = +- w`, where
`u(eta, 1) = +- w^(-eta) u(-1 - eta, 0)`, it acts on `f = u(., 0)` as `f(eta) -> f(2 eta) +- w^(-(eta + 1))
f(-2 eta - 2)`; `k` and `m - k` have conjugate spectra, so `1 <= k <= (m - 1)/2` is scanned, with both signs). The
rates agree with the audit's eigenvalues of the full `2m x 2m` operators at every modulus recomputed there: the
primes up to 127 and `9, 15, 21, 25, 27, 51, 63, 65, 85, 93, 281` (audit/proof-5128-AB/check_B1.out,
audit/numbers-vs-outputs/verify_misc.out, check_rate_281.out); at 1039 the audit's power iteration over all `k`
gives the same `0.7636`, at `k = 64` (audit/recompute-5128/a9b_rate_power.out). The reduction to the two blocks,
the divisor property below and the dictionary with S82 above were checked on the full operators at moduli up to
93 when the corrections were verified (audit/verify-corrections/HYP-9660/verify_hyp9660.out).

The maximum `0.932039` is attained exactly at the multiples of 7, every one of which has that rate; the largest
rate at a modulus prime to 7 is `0.91374` (`m = 93, 279, 465`), then `0.9002` (`31, 155, 341, 403`) and `0.9001`
(`381`). The rate of a modulus is at least that of each divisor `d > 1`: for `k = (m/d) k'` the functions of
`(eta mod d, beta)` form an invariant subspace of `T_k`, on which it acts as the operator `T_k'` of the modulus
`d`. So every multiple of 7 has at least the rate of 7 (exactly that rate in the scan), and (B) at a modulus
contains (B) at each of its divisors.
The rate of a composite modulus is NOT the maximum over its prime-power divisors: `15, 45, 75` have `0.8922`
(divisors `0.8660, 0.7418, 0.8153`), `51` has `0.8690`, `65` `0.8481`, `85` `0.8503`, `93` `0.9137` (`31` has
`0.9002`), `381` `0.9001` (`127` has `0.8754`), while the multiples of 7 stay at `0.932039`. (Up to 101 the
exceptions are exactly `15, 45, 51, 65, 75, 85, 93`, audit/numbers-vs-outputs/check_rates_npz.out; 45 and 75
inherit the rate of 15. Above 101 the table has 26 more: `171, 205, 219, 381, 395, 445`, with the rates `0.8776,
0.8180, 0.8786, 0.9001, 0.7715, 0.7826`, are above every proper divisor (these six were recomputed from the full
operators with dense eigenvalues at every `k` when the corrections were verified,
audit/verify-corrections/HYP-9660/verify_exceptions_above_101.out), and 20 others, among them `1023` and `1143`,
have the rate of a proper divisor (list in verify_hyp9660.out of the same directory).) So the gap is a statement
about all characters, not about primes. The slowest prime moduli of the scan are `7, 31, 127` (`0.932, 0.900,
0.875`), then `3` (`sqrt(3)/2 = 0.866`), `11` (`0.859`), `257` (`0.855`), `73` (`0.852`) and `241` (`0.851`).
`7, 31, 127, 73` divide `2^j - 1` with `j = 3, 5, 7, 9`, but a small order of 2 is not by itself a criterion: `5`
divides `2^4 - 1` and is the fastest modulus of the scan (`0.742`), and `11` and `257`, of orders 10 and 16, are
slower than `73`.

What (B) would give. THM-5128 (B1) proves `rho(T_k) < 2` for each fixed `(m, k)`, and on this THM-5128 (C), (D)
determine the residue-class side at each fixed odd prime and squarefree odd modulus (THM-5128 (E) does so at every
power of 2, by unitarity instead of (B1)), with a power saving uniform in the region but not in the modulus: the
exponent `theta_p` depends on the rate at `p`, and the constants `C_m` of THM-5128 (B2) come from Gelfand's
formula and are not explicit. Uniformity in the modulus is OPEN. (B) would make the exponent of the remainder of
THM-5128 (C2) uniform; a remainder uniform in the modulus needs in addition a bound for the constants,
`||T_k^r|| <= C m^A (2 rho)^r` with `C`, `A` and `rho < 1` independent of `m` and `k` (not scanned; the audit
measured cell constants `0.5` to `0.6` at `m = 3, 5, 7, 9, 15` for `n <= 18`, audit/proof-5128-AB/check_B24.out).
(B) with this norm form is what the other documents of the session quote as "uniformity in the modulus,
HYP-9660 (B)": a sieve or circle-method treatment of the residue classes needs it for moduli up to a power of
`X`, and it is the congruence input of the thin-orbit density theorems quoted in HYP-9661 (D), where the gap is
uniform in the modulus. For prime-power moduli the by-value class law is numerical only (THM-5128 (D3)). A
modulus-uniform bound for the count on full blocks exists (THM-4886 (A1), S83: for `n >= 7` and `r/q` in lowest
terms with `|2 theta - r/q| <= 1/q^2`, `2 <= q <= 2^((n-3)/3)`, the sum `S_n(theta)` of `e(theta H)` over the
one-apex words of length `n` satisfies `|S_n(theta)| <= 3.52 * 2^(n-1) q^(-1/2)`; (A3) covers minor arcs); what
(B) would add is a rate per digit on cells, uniform in `m`.

Two further moduli were computed because they looked able to refute (B): `1023 = 2^10 - 1 = 3 * 11 * 31` has
`0.913738` (the rate of 93) and `1143 = 9 * 127` has `0.900128` (the rate of 381); the method is that of
`m >= 203`. Able to refute: any other modulus above 521.

## (C) The first-order constants (QUESTION; NUMERICAL)

`A_p` (THM-5128 (D2)) is positive at `p = 3, 5, 11, 17, 43` and negative at `p = 13, 31` (partial sums over the
states below `2^26`, ap_constants.out; the audit's continuation, to `2^36` at 3 and 5, `2^38` at 11, 13, 17 and
`2^40` at 31, keeps these signs, audit/recompute-5128/a5b_Ap_extended_*.out, and its partial sums over the states
below `2^28` are negative at `19, 29` and positive at `23, 257`, audit/recompute-5128/a5_Ap.out). At `p = 7` the
constant is not determined by the states below `2^26`: the partial sums are `0.0058, 0.0066` there and still
falling; continued over the 7-divisible states to `2^40` they are `0.0006, -0.0002`
(audit/recompute-5128/a5b_Ap_extended_p7_40.out). The sign of `A_7` is undetermined.

Is `A_p != 0` for every `p`? Is `A_p > 0` exactly when the `p`-divisible states are early enough to beat the
deficit `-1/p^2` of the states before them? (A balanced first divisible state does not decide it: at 13 the first
one is `h = 32`, `x = y = 65`, the next is `h = 721`, and `A_13 < 0`.) Does `sum_p A_p` converge?

## (D) The two-apex grammar (FINITE-EXACT data, no theory)

Over all `4^n` words of `G(2,0)` (THM-4878), `n <= 16`, the class factors of `H` are within `0.0004` of 1 at
`p = 5, 7, 11, 13` at `n = 16`, the odd classes modulo 8 and 16 within `0.0032` (modulo 8 the largest deviation is
`0.0094, 0.0034, 0.0020, 0.0041, 0.0032` for `n = 12, ..., 16`: small, not monotone). At `p = 3` the class factors
(classes 0, 1, 2) are `(1.0102, 0.9973, 0.9925)`, `(1.0068, 1.0009, 0.9923)`, `(1.0062, 1.0020, 0.9918)`,
`(1.0052, 1.0030, 0.9918)`, `(1.0048, 1.0025, 0.9927)` for `n = 12, ..., 16`: the excess of class 0 halves, the
deficit of class 2, the class `-1` (`0.73 %` to `0.82 %`), does not move in this range, and for `n >= 13` it is the
largest deviation. Modulo 9 at `n = 16` the factors are `1.0051, 1.0046, 1.0047` at the classes 0, 3, 6 and
`0.9925, 0.9926, 0.9929` at the classes 2, 5, 8. Receipts: g20_local/full_hist_report.out (class 0 and the extreme
factors of the other classes); the factors class by class are from the audit
(audit/numbers-vs-outputs/check_g20_mod3.out), the histograms not being kept in the repository.

Question: do all three factors at 3 tend to 1, and at what rate? By THM-5128 (A) a seam of rank `r` with
equidistributed states would have a bias `p^(-r)` at one class only (class factors `1 + p^(1-r) - p^(-r)` there and
`1 - p^(-r)` at the other classes), smaller for larger `r`. That `G(2,0)`, or the grammars `G(k, w)` in general,
behave this way is conjectural: it is not shown that they are such seams. The three factors at 3 are not of that
shape, no two of them being equal (the classes 1 and 2 differ by `0.009` to `0.011` for `n = 13, ..., 16`). So
whether the bias at 3 is a mixing rate or a limit is not decided by these data.

## (E) Coprime states (QUESTION; NUMERICAL)

The proportion of `h < 2^n` with `gcd(x(h), y(h)) = 1` is `0.8262, 0.7707, 0.8024, 0.8166, 0.8118` at
`n = 10, 14, 18, 22, 24` (divisible_states.out) and `0.8104` at `n = 26` (audit/recompute-5128/a4_divisible.out).
The product over the odd primes of the densities `1 - 1/p^2` of the states not divisible by `p` (THM-5128 (F)) is
`8/pi^2 = 0.8106`. THM-5128 (F) proves the density only for the `h` with `g(h)` coprime to a fixed squarefree odd
`m`; the limit over all primes is NOT proved, the tail of the product being uncontrolled. The agreement at
`n = 26` is not yet convergence: below `2^26` the 3-divisible and the 7-divisible states still have `1.035` and
`0.777` times their limiting densities `1/9` and `1/49` (same audit output); treated as independent conditions,
these two deviations change the proportion by about `-0.4 %` and `+0.5 %`, and cancel. A proof needs the
equidistribution of THM-5128 (B) uniformly in `p` up to a power of `2^n`, i.e. (B) above with its constants, plus
a large-sieve step.
