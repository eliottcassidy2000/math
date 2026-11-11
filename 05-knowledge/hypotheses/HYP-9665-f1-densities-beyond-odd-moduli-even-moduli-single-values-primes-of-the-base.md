---
id: HYP-9665
title: "F1 local densities beyond THM-5134: (A) even moduli 2^j m, j >= 2, carry the product of the flat 2-adic law and L_m (every ingredient is available; the proof is not written out); (B) question: does the multiplicity of a single odd N follow 2 kappa_sigma (8/pi^2) sigma_1(n)/n ln N on average over short ranges; (C) in every even base b > 2 the second base-b digit of N_b is equidistributed by value (b = 2 is THM-5128 (E), odd b is THM-5134 (G)); (D) question: are the first-order constants A_(m,d) computable from the phase clocks of the one-apex count"
status: >
  OPEN (all four parts); audited 2026-11-11 with THM-5134 / THM-5135 (PASS WITH CORRECTIONS, MISTAKE-763;
  independent audit by eight auditors, report
  05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/AUDIT_REPORT.md).
  (A): OPEN (not written out); every ingredient is available (audit). The sketch of the audit is not a written or
  audited proof; the case j = 1 is THM-5134 (D).
  (B): OPEN (a question). A proof of it needs at least a version of THM-5134 uniform in the modulus (the uniform gap of
  HYP-9660 (B) with a norm bound polynomial in the modulus, OPEN) and, beyond it, a second-moment bound for r (the
  collision count sum r(N)^2); neither is available, and (B) does not follow from the first alone.
  (C): OPEN (not written out) in even bases b > 2, where every ingredient is available as for (A) (audit); in base
  2 the statement is THM-5128 (E) with j = 2, and in odd bases it is PROVED in THM-5134 (G) (found by the audit).
  (D): OPEN (a question, not tested); the facts about fans that it quotes are KNOWN (THM-5128 (F), THM-5126 (A)),
  and the lists of the primes below 400 whose fan is the first state of type p are FINITE-EXACT.
  NUMERICAL evidence only where stated; nothing here is pre-registered (the flatness test of (A) was proposed by
  the audit and evaluated on its own census).
session: opus-2026-11-11-S92
date: 2026-11-11
source: 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/README.md
scripts:
  - 04-computation/experiments/s92_prime_powers_20261111/every_base.py, 04-computation/experiments/s92_prime_powers_20261111/every_base.out (part (C); sheet s >= 1, complete census regenerated after the audit)
  - 04-computation/experiments/s92_prime_powers_20261111/series_vs_counts.py, 04-computation/experiments/s92_prime_powers_20261111/series_vs_counts_25_24.out (part (B))
  - 04-computation/experiments/s92_prime_powers_20261111/bridges_exact.py, 04-computation/experiments/s92_prime_powers_20261111/bridges_exact.out (part (D); section S8)
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/proof_sketches_and_tables.txt (audit; sketches for (A) and (C), tables for (A), (B), (C))
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/even_moduli_census.py, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/even_moduli_census_22.out, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/even_moduli_census_24.out, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/even_moduli_census_26.out, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/even_moduli_census_28.out (audit, part (A); exact census by value)
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/rflat_detail.py, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/rflat_detail.out, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/recheck_12.out, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/class_factors_from_hist.out (audit, part (A); conditional flatness, the remainder at 12, class factors)
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/mixed_gap.py, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/mixed_gap.out, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/mixed_cells.py, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/mixed_cells.out (audit, part (A); characters of even order)
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/series_counts.py, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/series_counts_25.out, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/series_counts_trend.out, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/poisson_pair_floor.py, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/poisson_pair_floor.out (audit, part (B))
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/proof-5134-AE/series_vs_counts_audit_25_24.out, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/numerics-prereg/my_series_25_24.out (audit, part (B))
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/second_digit.py, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/second_digit_b6.out, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/second_digit_other_bases.out, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/second_digit_more.out, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/second_digit_b3.out (audit, part (C))
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/fan_primes.py, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/fan_primes.out, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/numerics-prereg/fan_primes.out (audit, part (D))
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/verify-HYP-9665/verify_hyp9665.py, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/verify-HYP-9665/verify_hyp9665.out, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/verify-HYP-9665/verify_census.py, 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/verify-HYP-9665/verify_census.out (verification of the corrections; the bounds on R_prod in (A), the algebra of (A), (C), (D), independent recomputation of the censuses quoted in (A), (B), (C))
related:
  - 01-canon/theorems/THM-5134-f1-local-densities-at-every-odd-modulus-prime-power-factor-and-the-four-square-singular-series.md (the proved law at odd moduli, (A)-(F); (G) every base and the second digit in odd bases)
  - 01-canon/theorems/THM-5128-f1-by-value-local-densities-first-order-class-law.md ((B) digit chain, cells and regions; (D2), (D4) first-order constants; (E) powers of 2; (F) fans; Remark 3)
  - 01-canon/theorems/THM-5126-unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law.md ((A) the balanced states are the fans)
  - 01-canon/theorems/THM-4792-two-sheets-mod-odd-m-merged-centre-reflections-unit-shift-and-alternation-classification.md ((c) alternation condition, Moree)
  - 05-knowledge/hypotheses/HYP-9660-local-densities-beyond-one-root-change-two-changes-uniform-gap-sign-of-ap-two-apex-bias-at-three.md ((B) uniform gap, OPEN)
  - 05-knowledge/results/session193_triangulated_proof_20261105/modular_sync/PROOF.md (Theorems 1, 3, 4; the count at every modulus)
  - 05-knowledge/results/s91_f1_local_densities_20261110/README.md (S91 note, section 3; fans and primes of even order)
  - 05-knowledge/results/one_apex_join_representation_mass_pi2_over_12_20261106.md (S82 note, section 7; rates for all a/q with q <= 120)
  - 05-knowledge/results/odd_tables_signed_carries_20261111/prime_layers/PHASE_CLOCKS.md (the phase data of (D))
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/AUDIT_REPORT.md (independent audit, findings F001-F111)
---

# HYP-9665 - what THM-5134 leaves open about the local densities of F1

Notation of THM-5134 (sheets `sigma`: sheet 0 is `s = 0` with `c_0 = -1`, sheet 1 is `s >= 1` with `c_1 = 0`;
`T_sigma(X)`, `#T = #T_sigma(X)`; types and their counts `D^(d) = D_sigma^(d)(X)`; `c_d`, `L_m`, `kappa_sigma`,
`A_(m,d)^sigma`) and of THM-5128 (`q(h) = h - psi(h)`, `g(h) = gcd(x(h), y(h))`, cells, rates `rho_m`).
`sigma_1(n)` is the sum of the divisors of `n`; `sigma` alone is the sheet index. Numbers written `a / b` are
sheet 0 / sheet 1. Session outputs are in `04-computation/experiments/s92_prime_powers_20261111/`; `audit/...`
stands for `05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/...`.

**(A) Even moduli (OPEN: not written out; every ingredient is available).** For `j >= 2` and odd `m >= 3`, the
limiting class factor of an odd class `c` modulo `2^j m` is `2 L_m(c - c_sigma)` on each sheet: the flat law on the
odd classes modulo `2^j` (THM-5128 (E)) times the odd law (THM-5134 (D)). The class factor is `2^j m` times the
share of the class; the even classes are empty, so the mean over all classes is 1. For `j = 1` the statement is
THM-5134 (D) itself, `N` being odd.

*Ingredients (all available; supplied by the audit).* A character `q -> e(k q/(2^(j-1) m))` of `Z/(2^(j-1) m)` of
even order (`2^(j-1)` does not divide `k`) is `q -> e(k' q/M'')` with `k'` odd and `M''` even. On the states `w`
modulo `2 M''` its twisted operator `T[w -> 2w] = 1`, `T[w -> 2w + 1] = e(k' c(w)/M'')`, `c(w) = ceil(w/2)`,
satisfies `T T^* = 2I`: the rows `w` and `w + M''` have the same two columns, with entries `(1, om)` and `(1, -om)`,
`om = e(k' c(w)/M'')`, because `c(w + M'') = c(w) + M''/2`. This is the row computation of THM-5128 (E) with an
even modulus `M''` in place of `2^j`; the rate is exactly `2^(-1/2)`. A character of odd order (`2^(j-1)` divides
`k`) is a character of `Z/m`, covered by THM-5128 (B1), (B2). So `(q(h) mod 2^(j-1) m, h mod m)` is equidistributed
on cells, up to `O_(j,m)((2 rho)^r)` on a cell with `r` free digits, `rho = max(rho_m, 2^(-1/2))`: the rate of
`2^(j-1) m` is the larger of `2^(-1/2)` and the rate of `m`, and it is the rate of `m` at every modulus computed.
For the count alone, equidistribution on the admissible residues (all of them for an odd modulus, the odd ones for
an even modulus) at every modulus `2^k m`, for arbitrary exterior words and with an exponential rate, is KNOWN
(Session193 modular_sync, Theorems 1, 3, 4; the S82 scan of all `a/q` with `q <= 120` includes even `q`). With
this ingredient the by-value proofs of THM-5128 (E) and THM-5134 (A)-(D) combine: for a fixed smaller state, on a
cell of the larger word with at least `j - 1` bottom digits prescribed, and up to the error of the cell,
`N mod 2^j` is uniform on the odd classes (`x1` is odd and `x2 = y2 - 2 q2`, with `y2 mod 2^j` constant on the
cell), `N mod m` is distributed as in THM-5134 (A), and the two are independent by the Chinese remainder theorem.
The cut-cell count of THM-5128 (B4) and the summation over the smaller state would then give, for every odd `c`,

    #{t in T_sigma(X) : N = c (mod 2^j m)} = 2^(1-j) sum over d | gcd(m, c - c_sigma) of (d/m) D_sigma^(d)(X) + O_(j,m)(X^(1 - theta'))

with an exponent `theta' > 0`. This is a sketch (audit/hypotheses/proof_sketches_and_tables.txt, section 1), not a
written or audited proof: the audit proved and tested the ingredient on characters of even order and did not
re-derive the summation line by line.

NUMERICAL (audit; mixed_gap.out, mixed_cells.out in audit/hypotheses/): `rho(T_k)/2 = 0.7071` for every character
of even order of `Z/M`, `M = 6, 12, 24, 10, 20, 40, 18, 36, 14, 28, 30, 60, 90`; the rate of `M` is the rate of its
odd part (`0.8660` at 6, 12, 24, 18, 36; `0.7418` at 10, 20, 40; `0.9320` at 14, 28; `0.8922` at 30, 60, 90); the
cell counts of `(q mod 2^(j-1) m, h mod m)` approach equidistribution at `0.74` to `0.84` per free digit.

*Test (conditional flatness; proposed and evaluated by the audit, not pre-registered).* For `2^j m` in
`{12, 20, 36, 72}` and each sheet let

    R_flat(X) = max over classes e modulo m and odd classes a modulo 2^j of | 2^(j-1) #{N = a (mod 2^j), N = e (mod m)} / #{N = e (mod m)} - 1 |.

(A) is refuted if `R_flat` fails to fall by a factor 2 from `X = 2^24` to `X = 2^28` at one of these moduli. The
remainder of the first-order product law,
`R_prod = (2^j m) max_{c odd} |#{N = c} - 2^(1-j) sum_{d | gcd(m, c - c_sigma)} (d/m) D^(d)| / #T`, is not a test of
(A): it lies between `2 R_m` and `2 R_m + 2 F R_flat`, where `R_m` is the remainder of THM-5134 (B) at `m` and `F`
the largest class factor modulo `m` (`m` times the share of a class), so it oscillates with the odd part (12,
sheet 0: `0.060, 0.028, 0.036, 0.034` at `X = 2^22, 2^24, 2^26, 2^28`, about twice the remainder
`0.028, 0.014, 0.018, 0.017` at 3). (Proof of the bounds: for the odd class with residues `a` modulo `2^j` and `e`
modulo `m`, the signed deviation counted in `R_prod` is `2 u_e + 2 f_e dl`, where `u_e` is the signed deviation of
the class `e` counted in `R_m`, `f_e` its class factor and `dl` the flatness defect of `(a, e)` counted in
`R_flat`; the defects of a class `e` have mean 0 over the odd `a`, so one of them does not have the sign opposite
to `u_e`. Checked in exact arithmetic on the audit's censuses, 13 moduli at `X = 2^20, ..., 2^28`:
audit/verify-HYP-9665/verify_hyp9665.out.) The falsifier written here before the audit ("a remainder that does not
decrease from `2^24` to `2^28` at 12 would refute it") was already met on sheet 0 (`0.028 -> 0.034`), although the
product structure holds to sampling accuracy; applied to the PROVED odd law, the same criterion would "refute"
THM-5134 at 9, 27, 81 (its table, sheet 0). It is withdrawn.

NUMERICAL (audit, exact census by value; even_moduli_census_22.out, ..., even_moduli_census_28.out,
rflat_detail.out, recheck_12.out, class_factors_from_hist.out in audit/hypotheses/). `R_flat` at 12:
`0.0021 / 0.0014`, `0.0009 / 0.0005`, `0.0005 / 0.0001`, `0.0002 / 0.0002` at `X = 2^22, 2^24, 2^26, 2^28`; at 20:
`0.0028 / 0.0018`, `0.0009 / 0.0007`, `0.0006 / 0.0003`, `0.0001 / 0.0002`; at 36: `0.0040 / 0.0024`,
`0.0020 / 0.0011`, `0.0005 / 0.0004`, `0.0005 / 0.0005`; at 72: `0.0053 / 0.0048`, `0.0036 / 0.0019`,
`0.0013 / 0.0009`, `0.0006 / 0.0006`. This is between 0.6 and 4.0 sampling scales `sqrt(2^(j-1) m/#T)` at these
sizes (4.3 at 72 at `2^20`). From `2^24` to `2^28` `R_flat` falls by factors between `2.1` (36, sheet 1:
`0.00108 -> 0.00051`) and about 12 (20, sheet 0: `0.00088 -> 0.00007`), where a pure sampling term would fall by
about `4.3`: the factor 2 is reached at the four moduli on both sheets, and (A) is not refuted (at 36 `R_flat` is
stationary over the last two octaves: `0.00052 / 0.00045` at `2^26`, `0.00050 / 0.00051` at `2^28`, that is 1.7 / 1.7
and 3.5 / 4.0 sampling scales; rflat_detail.out). The class factors at 12 at `X = 2^28` are on average
`1.6996 / 1.6936` on the four classes with 3 not dividing `N - c_sigma` (limit `2 L_3 = 1.7778`) and
`2.6008 / 2.6128` on the two others (limit `2.4444`): twice the factors at 3, as far from the limits as at odd
moduli.

**(B) Single values (OPEN: a question).** Let `r_sigma(N)` be the number of triples of the sheet `sigma` with
value `N`. Question: for `X^eps <= K <= X`, is

    sum_{N in (X, X + K], N odd} | r_sigma(N) - 2 kappa_sigma S_sigma(N) ln N | = o(K ln X),      S_sigma(N) = (8/pi^2) sigma_1(n)/n ?

Here `n` is the odd part of `N + [sigma = 0]`, as in THM-5134 (E). The statement is the F1 analogue of Ingham's
asymptotic for a sum of two products (THM-5134 (F)), asked on average only: F1 has misses at every scale tested
(THM-5128 Remark 3).

*What it needs.* Summed over a class modulo a fixed `m`, (B) with `K = X` gives back the limit law of
THM-5134 (D). It does not follow from class laws: these give the class means of `r` only (a count taking the
values 0 and twice its mean at random has the same class statistics). A proof of (B) needs at least a version of THM-5134
uniform in the modulus (the uniform gap of HYP-9660 (B) with a norm bound polynomial in the modulus, OPEN;
THM-5128 Remark 3) and, beyond it, a second-moment (minor-arc) input: a bound for the collision count
`sum r_sigma(N)^2` showing that `r/ln N` concentrates around its class mean. Neither is available; the repository
has only the dispersion index of THM-5128 Remark 3. In the other direction, (B) with an unquantified `o()` gives
back the class law at each fixed `m` only, so a remainder of THM-5134 (B) that is not uniform for `m` up to a
power of `ln X` would not refute it (the falsifier written here before the audit is withdrawn).

NUMERICAL (exploratory; series_vs_counts_25_24.out; odd `N` in `[2^24, 2^25)`, that is `K = X`): mean
multiplicities `5.918` (sheet 0) and `7.482` (sheet 1), against the means `5.50` and `7.23` of
`2 kappa_sigma S_sigma(N) ln N` (audit/numerics-prereg/my_series_25_24.out); correlation of `r` with
`sigma_1(n)/n`: `0.4016` and `0.4275` (with the squarefree-level product of THM-5128 (D3): `0.3929`, `0.4181`);
mean of `r` by quintile of `sigma_1(n)/n`: from `4.30` to `8.55` (sheet 0) and from `5.44` to `10.77` (sheet 1),
ratios `1.99` and `1.98`, while the quintile means of `sigma_1(n)/n` run from `1.001` to `1.609` (ratio `1.61`): the
multiplicity varies more steeply than the series at this size. Likewise the ratios by `v_3(n)` and `v_5(n)` are
above their limits (sheet 0: `1.455, 1.722` at 3 against `1.333, 1.444`; `1.318, 1.323` at 5 against
`1.200, 1.240`), as the positive `1/ln X` terms of THM-5134 (D) at 3 and 5 predict.

NUMERICAL (audit; series_counts_trend.out, series_counts_25.out, poisson_pair_floor.out in audit/hypotheses/):
the L1 remainder of (B), `mean |r - 2 kappa S ln N| / mean r`, is `0.605 / 0.570`, `0.540 / 0.508`,
`0.497 / 0.456` on the odd `N` of `[2^16, 2^17)`, `[2^20, 2^21)`, `[2^24, 2^25)`, against `0.532 / 0.505` for a
constant predictor and `0.454 / 0.414` for the Poisson floor of an ordered-pair count (`r = 2 Poisson(mean/2)`) on
the last window: there the remainder is 9 % / 10 % above the counting-noise floor, compatible with (B) and not a
test of it. On that window the series explains `16 % / 18 %` of the variance of `r` (squared correlation),
`Var(r)/mean(r)` is `2.79 / 3.06`, and `8.07 % / 4.86 %` of the odd `N` have `r_sigma(N) = 0`
(audit/proof-5134-AE/series_vs_counts_audit_25_24.out). There are no data on short ranges `K < X`.

**(C) The primes of the base: the second digit in even bases `b > 2` (OPEN: not written out).**
`N_b = 1 + alpha_1 b + alpha_2 b^2` (THM-5134 (G)), so `N_b = 1 (mod b)` and the second base-`b` digit of `N_b` is
`alpha_1 = psi1 + psi2 + (2^s y1 y2 - [s = 0])/2 (mod b)`. Statement: it is equidistributed modulo `b` by value on
each sheet. By base:

- `b = 2`: PROVED. The second binary digit of `N` is `N mod 4`: THM-5128 (E) with `j = 2`.
- odd `b`: PROVED, THM-5134 (G) (found by the audit). Modulo `b`, `alpha_1` is an affine form in the residue vector
  `(x2, y2)` of the larger state with a unit coefficient at `x2`, so it takes each class on exactly `b` residue
  vectors, and THM-5128 (B4)(iv) applies.
- even `b > 2`: OPEN (not written out); this is what remains of (C). Every ingredient is available (audit sketch,
  audit/hypotheses/proof_sketches_and_tables.txt, section 2). For `b = 2^a` the by-value proof of THM-5128 (E)
  applies: `psi2 = h2 - q2`, and on a cell with `a` bottom digits prescribed `h2` is constant modulo `2^a` while
  `q2` is uniform modulo `2^a` (unitarity of THM-5128 (E)). For `b = 2^a m` with `m > 1` odd it is the mixed
  equidistribution of (A).

NUMERICAL (every_base.out; sheet `s >= 1`, `N_b <= 2^24`, complete census regenerated after the audit): the
deviation from flatness `max_c |b share(c) - 1|` is `0.001` at `b = 4`, `0.052` at `b = 6` and `0.003` at `b = 10`
(`0.0034` to four decimals, audit); for comparison `0.000` at `b = 2`, `0.040` at `b = 3`, `0.002` at `b = 5`. So it
is at most `0.0034` at `b = 2, 4, 5, 10`. The census behind the pre-audit figures missed 3.3 % of the triples at
`b = 2` and at most 0.14 % in the other bases; the figures did not change at the printed precision.

NUMERICAL (audit, complete census; second_digit_b6.out, second_digit_other_bases.out, second_digit_more.out in
audit/hypotheses/). `b = 6`: `0.184 / 0.163`, `0.090 / 0.098`, `0.048 / 0.052`, `0.028 / 0.029` at
`X = 2^16, 2^20, 2^24, 2^28`, a decay of `0.855 / 0.865` per binary digit; the digit modulo 2 is flat at the noise
level (`0.00007 / 0.00004` at `2^28`) and the digit modulo 3 carries the deviation, which decays as in base 3
(`0.867` per binary digit there, THM-5134 (G); the rate of the digit chain at 3 is `sqrt(3)/2 = 0.8660`). `b = 4`:
`0.0007 / 0.0009` at `2^24`, `0.0001 / 0.0003` at `2^28`. `b = 10`: `0.0051 / 0.0034` at `2^24`, `0.0015 / 0.0007` at
`2^28`. `b = 12` at `2^24`: `0.064 / 0.062`, with the digit modulo 4 flat to `0.0015`.

Also open: for odd `b` the 2-adic law of `N_b` is not known.

**(D) First-order constants (OPEN: a question, not tested).** The constants of THM-5134 (C),
`A_(m,d)^sigma = lim_{H -> infinity} sum_{h < H} chi_d(h) w^sigma(h)`, are conditionally convergent sums over all
states of the centred indicator `chi_d = [type d] - c_d`: the series is not absolutely convergent (THM-5128 (D2)),
and the sum over the states of type `d` alone diverges like `c_d kappa_sigma ln H`. For a prime modulus `p`,
`A_(p,p)` is the constant `A_p` of THM-5128 (D2), known numerically only (THM-5128 (D4); `A_7` is not determined,
sign included).

KNOWN (THM-5128 (F), THM-5126 (A); S91 note, section 3): a prime `p` has a balanced state (`x = y`) of type `p`
(for the modulus `p`, that is with `p | g(h)`: a `p`-divisible state of THM-5128) if and only if `p | 2^t + 1` for
some `t`, that is if and only if `ord_p(2)` is even; the first one is the fan `h = 2^(t-1)`, `t = ord_p(2)/2`; a
prime with `ord_p(2)` odd has none. Added here: for a prime modulus this is the alternation condition of
THM-4792 (c) (Moree); and, FINITE-EXACT, the fan is the first state of type `p` exactly for
`p = 3, 5, 11, 13, 17, 43, 113` among the primes below 200, then for `241, 257, 331` among the primes below 400
(the S91 audit had recorded it for `3, 5, 11, 13, 17, 43` and for `257`, `65537`; added are 113, 241, 331 and the
completeness of the lists). It is not the first one for `p = 19, 29, 37, 41, ...` (first states of type `p`:
`218`, `1435`, `1350`, `266`, as in THM-5128 (F)). Receipts: bridges_exact.out S8 for the primes below 200 (it
searches the states below `2^16`), confirmed by the audit's searches to `2^24`
(audit/numerics-prereg/fan_primes.out) and to `2^26` for the primes below 400 (audit/hypotheses/fan_primes.out).

The states of type `p` are the `h` with `2h + 1 = 0` and `q(h) = 0 (mod p)`. This is a condition on
`(q mod p, h mod 2p)`, the state of the digit chain of THM-5128 (B) taken together with its additive functional `q`
modulo `p`; so they form a 2-automatic set (an automaton with at most `2 p^2` states reads the binary digits of `h`
from the top) and their positions are determined without further data: the question as first written here ("are
the positions of the states of type `p`, hence `A_(p,p)`, determined by the phase data?") has a trivial positive
answer. What is open is a closed form or an effective evaluation of `A_(p,p)`. Question (from the sweep; not
tested, the audit read the definitions only): do the phase data `(i, a_i)` of the odd-tables packet
(PHASE_CLOCKS.md), which describe the one-apex count along the doubling orbits modulo `p`, give a closed form for
`A_(p,p)`?
