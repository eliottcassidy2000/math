---
id: HYP-9400
title: "The untouchable density and the Pollack-Pomerance model at computable scales. (1) PP's model constant is Delta = 0.1786569 (+- 3e-7), not 'about 0.17': the quoted 0.171822 is the partial sum at y = 2e10, and S(y)/log y = Delta - 0.162108/log y. (2) PP's Conjecture 1.4 with this constant: the natural density of untouchables exists and equals 0.1786569. (3) PP's key assumption (s maps T_a = {gcd(n, A_y) = a} into T_a) fails at each small prime l for a share ~ C_l (log X)^(-1/(l-1)) of inputs (52% at l = 5 and 80% at l = 13 near 2e9); its signature is the target residue structure given by the kill/serve rule of the prime-extension families, which the class model cannot see. (4) Falsifiable predictions to 2^40 and for the decade (10^12, 10^13]."
status: >
  OPEN (opus-2026-10-26-S57). (1) is VERIFIED numerically (three methods agree: window averages to 1e-8,
  independent-valuation Monte Carlo, partial sums reproducing PP's printed 0.171822 to 7 digits). (3)'s rate is CITED
  (Landau-Wirsing / Narkiewicz type) and measured (fitted exponents -0.525, -0.258, -0.167, -0.101, -0.083 against
  -1/2, -1/4, -1/6, -1/10, -1/12 for l = 3, 5, 7, 11, 13 on X = 2^10..2^31). The kill/serve rule is PROVED (elementary).
  The census to 1e10 is FINITE-EXACT (= A057978). Pre-registered tests passed: U-P2..U-P5 (to 1e9) and U-C,
  U-P6..U-P11 (1e10, committed before the analysis).
session: opus-2026-10-26-S57
source: 05-knowledge/results/s57_untouchables_root27_20261026/README.md
related:
  - 05-knowledge/results/sheets_root27_rank_untouchables_20261025.md
  - 05-knowledge/results/session66_root27_aliquot_graph_interfaces_20261026/session66-aliquot/REPORT.md
  - 05-knowledge/results/collatz_aliquot_lehmer_five_20260927.md
scripts:
  - 04-computation/untouchables_s57/ (untouchable_sieve, pp_delta, pp_truncated, pp_classwise, pp_class_drift, pp_xbin_drift, pp_residue_drift, pp_assumption_failure, core_decomposition, odd_goldbach_dependence; outputs beside them)
---

# HYP-9400: the untouchable density against the Pollack–Pomerance model

Let `U(X)` count the untouchable numbers in `[1, X]` (`m >= 2` with `m != s(n)` for all `n`). Known:
`liminf U(X)/X > 0.0604529` (Chen–Zhao, certified instantiation in session66), `limsup < 1/2` (Luca–Pomerance), and
5 is the only odd untouchable below `4·10^18` (Goldbach). Existence of the density is open.

1. **The model constant (numerical).** `Delta = lim (1/log y) sum_(a <= y even) e^(-a/s(a))/a = 0.1786569`. Window
   averages of `(1/2) e^(-a/s(a))` over even `a` in `[10^k, 10^(k+1))` equal `0.17865689` for `k = 5..8`. The
   independent-valuation Monte Carlo gives `0.178644 +- 0.000011`, plus about `1.4e-5` for the primes above 2000.
   The partial sums satisfy `S(y)/log y = Delta - 0.162108/log y` at `y = 10^6..10^9` to `3e-8`, and at `y = 2·10^10`
   this gives `0.1718224`. Pollack–Pomerance print `0.171822` there.
2. **Conjecture (PP 2016, Conjecture 1.4, with the evaluated constant).** `U(X)/X -> 0.1786569`. Evidence:
   - the local density per decade is `0.1662, 0.1684, 0.1701, 0.1713` for `10^8..10^12`, increasing and decelerating;
   - `gap = 5.35/log^2 N` to the model value fits within `0.0006`;
   - a pure `1/log^2` extrapolation of the last decades gives `0.1774` and a `1/log` one gives about `0.184`.
   The data do not decide.
3. **The model's assumption fails at computable scales, with an exact signature.**
   - PP assume `s(T_a) ⊂ T_a`, true once `n > e^(e^y)`. It fails at the prime `l` whenever `l` does not divide
     `sigma` of the cofactor. The share of such inputs is `0.135, 0.521, 0.569, 0.816, 0.797` for
     `l = 3, 5, 7, 11, 13` near `2^31`, decaying like `(log X)^(-1/(l-1))`.
   - The class model sees a target only through `gcd(k, A_y)`, so it gives equal shares to the nonzero residues
     mod `l`. The kill/serve rule says instead: a core `a` with `l` not dividing `s(a)` cannot serve
     `k = sigma(a) (mod l)`, and a core with `l | s(a)` serves only that residue.
   - This rule creates the measured spreads of the share over nonzero residues at `W_33 = (2^32, 2^33]`: `1.453`,
     `1.290`, `1.281`, `1.201` (mod 5, 7, 11, 13; from `1.92`, `1.73`, `1.90`, `1.43` at `W_14`). The mean
     multiplicity by residue (`[10^7, 10^8]`, exact core decomposition) orders the shares inversely.
   - The same failure shows in the per-number law. For `x = k/s(k) >= 1.4` the share is flat at `0.220–0.226` from
     `2^24` to `2^34`, where PP predicts `e^(-x)` falling to `0.141`. The abundant end rises toward PP from below.
4. **Predictions** (each falsifiable by a census with residue and abundancy tallies):
   - (a) The nonzero-residue spreads decrease monotonically in the dyadic windows. At `W_40` the mod-5 spread lies
     in `[1.33, 1.42]` and the mod-13 spread in `[1.14, 1.20]`. They tend to 1 (PP's equal shares) only like
     `(log k)^(-1/(l-1))`, so the mod-13 spread stays above `1.10` for `k < 10^30`.
   - (b) The share at `x` in `[1.8, 2.0)` stays in `[0.215, 0.230]` through `2^40`.
   - (c) The local density in the decade `(10^12, 10^13]` lies in `[0.1718, 0.1728]`. Every decade density stays
     below `Delta = 0.1786569` for `N <= 10^15`.
   - (d) The argmax/argmin residues are stable: mod 5 at 1/4, mod 7 at 3/2, mod 11 at 1/5, mod 13 at 3/7.

Non-consequences: (1) evaluates the model, not the density. A census confirming (4) supports the mechanism of (3),
not the value in (2). Nothing here bears on Goldbach or the odd untouchables beyond sec. 3 of the source note.
