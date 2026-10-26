---
id: HYP-9400
title: "The untouchable density and the Pollack-Pomerance model at computable scales. (1) PP's model constant evaluates to Delta = 0.17865689; PP's printed 0.171822 (y = 2e10) and CGJM's partial sums to 2^40 are intermediate values, S(y)/log y = Delta - 0.1621075/log y. (2) PP's Conjecture 1.4 with this constant: the natural density of untouchables exists and equals 0.17865689. (3) At computable sizes the real map departs from PP's random map in two measured ways: the class assumption gcd(s(n), A_y) = gcd(n, A_y) fails jointly for 62% (y = 13) to 82% (y = 53) of n near 2^27; residue neutrality fails at each prime l for a share ~ C_l (log X)^(-1/(l-1)) (Ramanujan-Watson), with the residue signature given by the classical kill/serve rule of the prime-extension cores, which the class model cannot see. (4) Falsifiable predictions to 2^40 and for the decade (10^12, 10^13]."
status: >
  OPEN (opus-2026-10-26-S57; independently audited, 05-knowledge/results/s57_untouchables_root27_20261026/AUDIT.md,
  MISTAKE-672, corrections applied). (1) is VERIFIED numerically by three evaluations: window averages to 1e-8, the
  audit's deterministic convolution of the independent-valuation model (0.17865689139), and Monte Carlo to 1e-5. The
  printed partial sums of PP and of CGJM's Table 2 are reproduced to 6 digits. The limit is not stated in PP, CGJM,
  OEIS A005114/A057978, erdosproblems #418 or Booker (searched, also by the audit). (3a) is FINITE-EXACT. (3b)'s rate
  is CITED (Ramanujan's manuscript sec. 19; Watson 1935; Rankin 1961; Scourfield 1964; Ciolan-Languasco-Moree 2023)
  and measured: the fitted exponents -0.525, -0.258, -0.167, -0.101, -0.083 match -1/(l-1) for l = 3, 5, 7, 11, 13
  on X = 2^10..2^31. The kill/serve rule is classical: the residue exclusion of PP 2016 sec. 2 and Pomerance-Yang
  2014 Prop. 2.3. The census to 1e10 is FINITE-EXACT (= A057978). Pre-registered tests passed: U-P2..U-P5 to 1e9, and
  U-C plus U-P6..U-P11 at 1e10, committed before the analysis.
session: opus-2026-10-26-S57
source: 05-knowledge/results/s57_untouchables_root27_20261026/README.md
related:
  - 05-knowledge/results/sheets_root27_rank_untouchables_20261025.md
  - 05-knowledge/results/session66_root27_aliquot_graph_interfaces_20261026/session66-aliquot/REPORT.md
  - 05-knowledge/results/collatz_aliquot_lehmer_five_20260927.md
scripts:
  - 04-computation/untouchables_s57/ (untouchable_sieve, pp_delta, pp_truncated, pp_classwise, pp_class_drift, pp_xbin_drift, pp_residue_drift, pp_assumption_failure, pp_class_preservation, core_decomposition, odd_goldbach_dependence; outputs beside them)
---

# HYP-9400: the untouchable density against the Pollack–Pomerance model

Let `U(X)` count the untouchable numbers in `[1, X]`, i.e. the `m >= 2` with `m != s(n)` for all `n`. Known:
- `liminf U(X)/X > 0.0604529` (Chen–Zhao, certified instantiation in session66);
- `limsup U(X)/X < 1/2` (Luca–Pomerance);
- 5 is the only odd untouchable below `4·10^18` (Goldbach).

Existence of the density is open.

1. **The model constant (numerical).** `Delta = lim (1/log y) sum_(a <= y even) e^(-a/s(a))/a = 0.17865689`.
   - PP (sec. 3.2) note it is effectively computable and that their printed numbers are intermediate values; they did
     not evaluate it.
   - Window averages of `(1/2) e^(-a/s(a))` over even `a` in `[10^k, 10^(k+1))` equal `0.17865689` for `k = 5..8`.
   - The audit's deterministic convolution gives `0.17865689139`.
   - The partial sums satisfy `S(y)/log y = Delta - 0.1621075/log y`. This gives PP's `0.171822` at `y = 2·10^10` and
     CGJM's `0.171617, 0.172257, 0.172790, 0.172810` at `10^10, 10^11, 10^12, 2^40`.
2. **Conjecture (PP 2016, Conjecture 1.4, with the evaluated constant).** `U(X)/X -> 0.17865689`. Evidence:
   - the local density per decade is `0.1662, 0.1684, 0.1701, 0.1713` for `10^8..10^12`, increasing and decelerating;
   - `gap = 5.35/log^2 N` to the model value fits the last three decades within `0.0005`;
   - a pure `1/log^2` extrapolation gives `0.1773–0.1778`, and a `1/log` one gives `0.1845–0.1866`.

   The data do not decide.
3. **How the real map departs from the model at computable sizes.**
   - (a) PP assume `s(T_a) ⊂ T_a`, true once `n > e^(e^y)`. Near `2^27` the gcd class changes, at `y = 13`, for
     `62%` of all `n`, and at `y = 53` for `82%`. Per prime (uncapped) the failure is `0.18, 0.24, 0.18, 0.15, 0.12` for
     `l = 3, 5, 7, 11, 13`.
   - (b) The random map treats the nonzero residues of `s(n) (mod l)` alike. The real map does so only on inputs with
     `l | sigma(n)`, where `s(n) = -n`. The complementary share is `0.135, 0.521, 0.569, 0.816, 0.797` near `2^31`,
     decaying like `(log X)^(-1/(l-1))`.
   - Signature of (b) in the targets: the class model gives equal shares to the nonzero residues mod `l`. The classical
     kill/serve rule says a core `a` with `l` not dividing `s(a)` cannot serve `k = sigma(a) (mod l)`, and a core with
     `l | s(a)` serves only that residue. The measured spreads of the share over nonzero residues at
     `W_33 = (2^32, 2^33]` are `1.453`, `1.290`, `1.281`, `1.201` (mod 5, 7, 11, 13), against `1.92`, `1.73`, `1.90`,
     `1.43` at `W_14`. Over `W_24..W_33` the spreads decay like `(log k)^(-0.7..-1.1)`, which is small-core dominated.
     Mod 3 is non-monotone (minimum `1.0095` at `W_18`, where the argmax flips).
   - The mean multiplicity by residue, from the exact core decomposition on `[10^7, 10^8]`, orders the shares inversely
     mod 5 and broadly mod 7 (rank correlation 0.77).
   - Separately, and not explained by (b): for `x = k/s(k) >= 1.4` the share is flat at `0.220–0.226` from `2^24` to
     `10^10`, where PP predicts `e^(-x)` falling to `0.141`.
4. **Predictions** (each falsifiable by a census with residue and abundancy tallies):
   - (a) For `l = 5, 7, 11, 13` the nonzero-residue spreads keep decreasing in the dyadic windows. At `W_40` the mod-5
     spread lies in `[1.33, 1.42]` and the mod-13 spread in `[1.14, 1.20]`. The audit's extrapolations of the fitted
     decay are `1.392` and `1.174`.
   - (b) The share at `x` in `[1.8, 2.0)` stays in `[0.215, 0.230]` through `2^40`.
   - (c) The local density in the decade `(10^12, 10^13]` lies in `[0.1718, 0.1728]`. Every decade density stays
     below `Delta = 0.17865689` for `N <= 10^15`.
   - (d) The argmax/argmin residues are stable: mod 5 at 1/4, mod 7 at 3/2, mod 11 at 1/5, mod 13 at 3/7.

Non-consequences: (1) evaluates the model, not the density. A census confirming (4) supports the description in (3),
not the value in (2). Nothing here bears on Goldbach or the odd untouchables beyond sec. 3 of the source note.
