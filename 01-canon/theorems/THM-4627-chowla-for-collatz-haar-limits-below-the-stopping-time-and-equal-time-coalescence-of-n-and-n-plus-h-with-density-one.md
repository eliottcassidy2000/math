---
id: THM-4627
title: "Chowla for Collatz: (1) for orbit functionals measurable with respect to the parity word up to the stopping time, ordinary two-point averages along any two non-proportional affine forms converge to the Haar limit of the full-period averages modulo 2^k, with the finite exceptional-height error retained; for the stopping-time parity s the limit correlation is -(1 + 2E[s]) at every odd lag, +(1 + 2E[s]) at every lag 2 mod 4, and 3/4 + (1/4)E[s(n)s(n+h) | n = 3 mod 4] at lags 0 mod 4; (2) for every h >= 1 the integer orbits of n and n+h coincide from a common time on, with equal odd-step counts, for a set of n of natural density one, so tau(n) = tau(n+h) and omega(n) = omega(n+h) whenever either side is finite, and the two-point correlations of (-1)^tau and (-1)^omega along (n, n+h) tend to +1 whenever the orbits reach 1; the integer correlation at lag h equals the fraction of pairs merging before reaching 1 to three decimals"
status: >
  CORRECTED 2026-10-13: part 1 uses coefficient-stopping truncations and a finite exceptional-height term; see incoming_cylinder_prime_repair_20261013.md.
  PROVED: statement 1 (Terras periodicity of parity words modulo 2^k and Terras's density theorem for the stopping
  time) and the exact lag laws; statement 2 given THM-4581 (mac-mini, PROVED: equal-time merging of y and y+h for
  almost every 2-adic y, u = v + e with e in Z[1/3], k = 0) and the periodicity argument; statement 2(d) is
  conditional on the orbits of n <= N reaching 1 (true for N <= 2^68, CITED Barina; in general the Collatz conjecture).
  FINITE-EXACT: the periodic censuses (depth <= 26, 2^26 residues) and the integer censuses to 2^28 (scripts below).
  NUMERICAL: the rates (HYP-9258). Self-audited; the narrow corrective audit is recorded in the results note (section 10).
session: opus-2026-10-13-S25
source: 05-knowledge/results/five_manuscripts_chowla_collatz_residue_bsd_20261013.md (section 3)
scripts:
  - 04-computation/experiments/collatz_chowla_20261013.py (+ .out, _periodic.npy, _int.npy)
  - 04-computation/experiments/collatz_merge_lags_20261013.py (+ .out)
  - 04-computation/experiments/collatz_tau_runs_20261013.py (+ .out)
related:
  - THM-4581 (Haar coalescence at equal Terras time; the input of statement 2)
  - THM-4565 (two readers of one tape; the pair chain used for the merge masses)
  - HYP-9240 (rivers of children: the depth-river form of the same coalescence)
  - HYP-9258 (rates, convergent denominators of log_2 3 at unequal-time coincidences, equidistribution of tau mod 2)
  - Terras (1976), Everett (1977): finite stopping time for almost all n; periodicity of parity words
  - the two-point Chowla manuscript (OpenAI, September 2026; PREPRINT CLAIM) whose Lemma 3.1 this theorem renders unnecessary for Collatz
---

# THM-4627 — Chowla for Collatz

## Setting

`T(n) = n/2` (`n` even), `(3n+1)/2` (`n` odd). `sigma(n)` is the stopping time (least `k` with `T^k(n) < n`), `kappa(n)` Terras's coefficient stopping time (least `k` with `3^(a_k) < 2^k`, `a_k` the number of odd steps among the first `k`), which depends only on `n mod 2^k` and equals `sigma(n)` outside a finite set in each residue class; `tau(n)` the number of `T`-steps to `1`, `omega(n)` the number of odd steps to `1`; `s = (-1)^sigma`, `t = (-1)^tau`, `lam = (-1)^omega`. A functional `f` of the orbit is *stopping-time measurable* if `f(n) = F(w_1(n), ..., w_sigma(n)(n))` for a bounded `F` on finite parity words (any bounded value when `sigma = infinity`).

## Statement

1. Let `f, g` be bounded stopping-time-measurable functionals and `a_i>=1,b_i>=0` with non-proportional affine forms. Use the coefficient-stopping models `f_k^*,g_k^*` defined in the repair below, and their full-period product mean `R_k`. With `H_k` the maximum affine descent threshold through depth k and `P_i(k)=P(kappa(a_i y+b_i)>k)`,
   `|N^-1 sum f(a1*n+b1)g(a2*n+b2)-R_k| <= ||f||||g||[P_1(k)+P_2(k)+(4*2^k+4*H_k)/N]`.
   Thus the ordinary averages converge to the Haar-model limit, first taking N to infinity and then k. Actual-sigma truncations are not exactly periodic at small sources. For `f=g=s`, the limiting lag laws survive: with `m=lim E[s]`, `m_odd=1+2m`, the correlation is `-m_odd` at odd h, `+m_odd` at h=2 mod4, and `3/4+(1/4)lim E[s(n)s(n+h)|n=3 mod4]` at h=0 mod4.
2. Fix `h >= 1`. (a) `M_h(K) = {n : exists k <= K with T^k(n) = T^k(n+h) and equal numbers of odd steps up to k}` has natural density equal to the Haar probability `mu_h(K)` of the same event for `(y, y+h)`, `y` Haar on `Z_2`. (b) `mu_h(K) -> 1` as `K -> infinity`. (c) Hence for `n` in a set of density one the orbits of `n` and `n+h` coincide from a common time on with equal odd-step counts; for such `n > 2^K` the merge precedes either orbit reaching `1`, so `tau(n) = tau(n+h)` and `omega(n) = omega(n+h)` whenever either side is finite, and `tau(n) < infinity iff tau(n+h) < infinity`. (d) If the orbits of all `n <= N` reach `1`, then `liminf_N (1/N) sum_(n <= N) t(n)t(n+h) >= 2 mu_h(K) - 1` for every `K`; so the correlations of `t` and of `lam` along `(n, n+h)` tend to `+1`.
3. (FINITE-EXACT.) To `N = 2^28`: `E[s] = -0.234140`, `m_odd = 0.53172`; the integer correlations of `s` equal `-0.53172` at every odd lag `<= 39`, `+0.53172` at every lag `2 mod 4`, and `0.724, 0.696, 0.744, 0.805, 0.746, 0.685, 0.726, 0.868, 0.718, 0.677` at lags `4, 8, ..., 40`; the depth-`26` periodic truncation gives `-0.5390` at odd lags with `P(kappa > 26) = 0.0155`. To `N = 2^26`, stopping when either orbit reaches `1`: `P(merge before 1) = 0.510, 0.419, 0.330, 0.360, 0.341, 0.302, 0.317, 0.284` and `C_t = 0.508, 0.417, 0.327, 0.357, 0.338, 0.298, 0.314, 0.281` at lags `1, 2, 3, 4, 5, 6, 8, 16`; merged pairs have equal `tau` and correlation exactly `1`, unmerged pairs correlation `-0.003` to `-0.005`; `P(merge)` at lag `1` rises `0.478, 0.495, 0.510` over `2^20, 2^23, 2^26`. Class-wide merge by depth `26` (exact over `2^26` residues): `0.421, 0.306, 0.198, 0.229, 0.211, 0.162, 0.174, 0.131`; at lag `1` the first-merge masses are `1/8` (depth 3), `1/16, 1/16, 1/32, 3/128, 1/64, ...`. The break density `P(tau(n) != tau(n+1))` falls `0.555 -> 0.481` from `2^16` to `2^28`, parity flips are half the breaks, the longest run of equal `tau` below `2^28` has length `120`, and the commonest nonzero differences are `+-19, +-38, +-27, +-46, +-8, +-57`.

## Proof

**1 (repaired).** The parity word through k, and its coefficient stopping time, depend only on the source modulo `2^k`. The actual first descent also depends on the affine carry and the source height. For `T^j(n)=(3^a n+c)/2^j`, coefficient descent becomes actual descent once `n>c/(2^j-3^a)`. Retain the maximum `H_k` of those finitely many thresholds. The periodic coefficient-stopping models then agree with the actual functionals except on their Haar tail or at sources at most `H_k`. Period remainders and the exceptional source count give exactly the stated bound; details and the hostile pair1/5 are in [the correction](../../05-knowledge/results/incoming_cylinder_prime_repair_20261013.md). Terras's tail theorem yields the limit. The lag identities follow by splitting n into even,1 mod4 and3 mod4 classes; the exceptional n=1 changes no limiting average. ∎

**2.** (a) Let the parity words of `n` and `n+h` up to step `k` have `a` and `a'` odd steps and constants `c, c'`, so `T^k(n) = (3^a n + c)/2^k`, `T^k(n+h) = (3^(a')(n+h) + c')/2^k`. Equality with `a = a'` is the condition `c = 3^a h + c'` on the two words, which holds for the whole class of `n mod 2^k` or not at all. Equality with `a != a'` is a nontrivial linear equation in `n`, satisfied by at most one `n` per class and step; these form a set of density zero for fixed `K`. Hence the density of `M_h(K)` exists and equals the Haar probability of the class-wide event, which is `mu_h(K)` for the pair `(y, y+h)`. (b) THM-4581 with `u = v + h`: `e = h in Z[1/3]`, `k = 0`; for almost every `y` the orbits of `y` and `y + h` meet at equal Terras time with `y + h` making exactly `0` fewer odd steps than `y`, i.e. with equal odd-step counts; the events `M_h(K)` increase to this almost-sure event, so `mu_h(K) -> 1`. (c) For `n in M_h(K)` with `n > 2^K`, `T^j(n) >= n/2^j > 1` and likewise for `n+h` for all `j <= K`, so neither orbit has reached `1` when they coincide at time `k <= K`; from then on the orbits are identical, so `tau(n) = k + tau(T^k n) = tau(n+h)` and the odd-step counts agree before `k` (by the merge condition) and after `k` (identical orbits). The density of `M_h(K) cap (2^K, infinity)` is `mu_h(K)`, and `K` is arbitrary. (d) On `M_h(K) cap (2^K, N]` the product `t(n)t(n+h)` equals `1`; elsewhere it is at least `-1`; so the average is at least `mu_h(K) - (1 - mu_h(K)) - o(1)`. The same for `lam`. ∎

**3.** Computation; see the scripts. The integer merge census stops when either orbit reaches `1` (a first version counted two orbits resting at `1` as a merge and was discarded).

## Scope and boundary

- Statement 1 concerns ordinary averages at every cutoff; no logarithmic weighting is involved. The limits are typically nonzero: the stopping-time parity is `2`-adically locally constant, which is the opposite of the situation for the Liouville function.
- Statement 2(d) is conditional on reaching `1`; 2(a)-(c) are unconditional. The merge inside the final cycle (equal phase) would give equal parity but not equal `tau`; density-one merges happen at bounded depth, before that.
- The rate at which `mu_h(K) -> 1` and the rate at which the integer break density decays are not proved here (HYP-9258); THM-4581 gives `P(no merge by T) >= c T^(-1/2)` with the matching upper bound at sketch level.
- Nothing here bears on whether orbits reach `1`.

## Reading

Below the stopping time, Collatz statistics along affine forms are exact Haar computations because the relevant events are periodic; the Chowla manuscript's central lemma (bounded independence fooling residue tests beyond the joint period) is the price paid for the absence of periodicity in multiplicative number theory, and Collatz does not pay it. Beyond the stopping time the natural functionals are locally constant along lags because the orbits of `n` and `n+h` merge, a theorem of the repository's own thread; the only Chowla-flavoured question left is the one-point equidistribution of `tau mod 2`, which is as open as the conjecture.
