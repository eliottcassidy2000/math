---
id: THM-4487
title: "The dip spectrum of 3n+-1: #{n <= X : T^i(n) >= n^gamma for i <= log_2 n} = X^(h(gamma/log_2 3)+o(1)) for gamma in (log_4 3, 1] (Terras at gamma = 1, Korec at gamma = log_4 3), and the counting lemma of THM-4476 is sharp, so the thin-divergence exponent h(log_3 2) cannot be lowered by that lemma"
status: >
  PROVED (elementary; self-audited; two-sided exponent on both sheets) +
  FINITE-EXACT controls to 2^24. Let b = +-1, T_b(x) = x/2 (x even),
  (3x+b)/2 (x odd), alpha = log_2 3, h the binary entropy, and
  D_b(X, gamma) = #{n <= X : T_b^i(n) >= n^gamma for all 0 <= i <= floor(log_2 n)}.
  (1) For every gamma in (log_4 3, 1], c X^(h(gamma/alpha))/log^3 X <= D_b(X, gamma)
  <= C X^(h(gamma/alpha)) log^2 X, so log D_b(X, gamma)/log X -> h(gamma/alpha);
  for gamma in (log_2(3/2), log_4 3] the limit is 1. So the exponent is
  h(log_3 2) = 0.949956 at gamma = 1 (Terras's undecided count) and rises to
  1 exactly at Korec's exponent log_4 3 = 0.792481; the constants
  0.949956, 0.050044, 0.792481, 0.207519, 0.48796 (the Chernoff tilt =
  -E'(1)) and 1.052681 (= 1/E(1)) of the procgen, crossroads and opus lanes
  are values, slopes and reciprocals of the one function
  E(gamma) = h(max(1/2, gamma/alpha)).
  (2) The no-dip set F_b(X, theta) of THM-4476 satisfies
  c X^(h(rho))/log^2 X <= #F_b(X, theta) <= C X^(h(rho)) log^2 X, rho = (1-theta)/alpha,
  theta in (0, 1 - alpha/2). Hence the recursion N(X) <= k N(X^(1-theta)) +
  #F_b(X, theta) cannot yield an exponent below h(log_3 2): THM-4476's
  exponent is optimal for its method. The exponent of an actual divergent
  orbit remains OPEN (conjecturally no polynomially dense divergent orbit
  exists). Upper bounds: Terras bijection + entropy count of words with
  prescribed odd density; lower bounds: the cycle-lemma rotation of THM-4478
  section 4 applied to no-dip prefixes. The q-analogue has rho = gamma/log_2 q
  < 1/2 for all gamma <= 1 once q >= 5, so no dip theorem of this kind exists
  for 5n+1. Sheet-blind, defect-blind; not a divergence exclusion.
source: collatz-exponent-atlas-20260926 session (opus), 2026-09-26; the owner asked to sharpen the THM-4476 exponent or prove it optimal, and to find recurrent numbers across threads. Mechanism: Terras count (both directions) with the crossroads rotation for the lower bound.
depends_on:
  - 01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md (Terras bijection, carry bound, the no-dip set F_b)
  - 01-canon/theorems/THM-4478-critical-tube-affine-capacity-sharp-provability-price.md (section 4: cycle-lemma rotation, exponentially many sources in a growth band)
related:
  - 05-knowledge/results/collatz_dipspectrum_20260926_entropy_curve.md (full note)
  - 05-knowledge/results/constants_atlas_20260926_recurrent_numbers.md (the atlas that motivated it)
  - 05-knowledge/results/collatz_procgen_20260922_barrier_atlas.md (Korec 1994 as cited there)
  - 05-knowledge/results/collatz_procgen_20260922_synthesis.md (section 2k: one exponent in three settings)
script: 04-computation/experiments/collatz_dipspectrum_20260926.py
output: 05-knowledge/results/collatz_dipspectrum_20260926.out
script_sha256: 7db80f5077f9fe641b790a8a4a6455aa9f187db091cca907b61e692542da2244
output_sha256: 9b0440da365859dc28a6bd49970be6e9f133c55d967fa95d3fddaf38b626caa2
hash_basis: raw LF bytes
audit: >
  Self-audited: both directions re-derived; the lower-bound construction
  checked against the direct counts (which carry the expected 1/k prefactor:
  the observed four-doubling slopes sit near h(rho) - 1/(k ln 2)); a
  counting off-by-one in the first control script was found and fixed
  before filing (the earlier rows counted n < 2^(t+1) as n <= 2^t).
  Independent audit not yet performed.
---

# THM-4487 -- the dip spectrum and the sharpness of the thin-divergence exponent

**PROVED.** Full note:
[collatz_dipspectrum_20260926_entropy_curve](../../05-knowledge/results/collatz_dipspectrum_20260926_entropy_curve.md).

## 1. Statements

For `b = +-1`, `alpha = log_2 3`, `h` the binary entropy:

* **Dip spectrum.** For `gamma in (log_4 3, 1]`,
  `#{n <= X : T_b^i(n) >= n^gamma, 0 <= i <= floor(log_2 n)} = X^(h(gamma/alpha) + o(1))`;
  for `gamma in (log_2(3/2), log_4 3]` the exponent is `1`.
* **Sharpness.** `#F_b(X, theta) = X^(h((1-theta)/alpha) + o(1))` for
  `theta in (0, 1 - alpha/2)`; the thin-divergence exponent `h(log_3 2)` is
  the best the dichotomy of THM-4476 can give.

## 2. Proof in four lines

1. **Upper bound.** No dip below `n^gamma` at step `t = floor(log_2 n)` forces
   `3^o/2^t >= n^(gamma-1)/2` (the carry is `O(n^0.585)` and
   `gamma > log_2(3/2)`), i.e. `o >= (gamma/alpha) t - 2` odd letters; such
   words number `2^(t h(gamma/alpha)) poly(t)`, and each is one residue class
   modulo `2^t` (Terras), with one representative per dyadic block.
2. **Lower bound.** Among the `C(t, o)` words with `o = ceil((gamma/alpha) t) + 1`
   odd letters, at least `C(t, o)/t` have all partial sums
   `S_i = i - o_i alpha` at most `max(0, S_t)` (rotate after the maximum
   partial sum: THM-4478 section 4). Their representatives in `[2^t, 2^(t+1))`
   satisfy `T^i(n) >= n 2^(-max(0,S_t)) - n^0.585 >= 3 n^gamma - n^0.585 >= n^gamma`;
   the count is `C(t, o)/t = X^(h(gamma/alpha))/poly(log X)`.
3. **Both sheets.** The carry is nonnegative on the plus sheet and
   nonpositive on the minus sheet; both directions only use `|carry| <= (3/2)^t`.
4. **Sharpness.** The same construction with `theta` in place of `1 - gamma`
   on the window `[1, floor(log_2 X)]` and representatives in `(X/2, X]`
   gives `#F_b(X, theta) >= X^(h(rho) - o(1))`; the upper bound is the lemma
   of THM-4476.

## 3. Consequences

* **One curve.** `E(gamma) = h(max(1/2, gamma/alpha))` carries `0.949956`
  (`E(1)`), `0.050044` (`1 - E(1)`: the sharp price exponent of THM-4475/4478/4479
  and the Chernoff rate), `0.792481` (where `E = 1`: Korec), `0.207519`
  (`1 - alpha/2`: the drift), `0.48796` (`-E'(1)`: the Chernoff tilt) and
  `1.052681` (`1/E(1)`: THM-4476's growth threshold).
* **THM-4476 is optimal for its method.** Any lower exponent for a
  divergent orbit needs a constraint on no-dip points beyond one
  `log_2 X`-window, i.e. small integers in residue classes modulo more than
  `X` (the transversality barrier).
* **`5n+1`.** `gamma/log_2 5 < 1/2` for every `gamma <= 1`: the dip count has
  exponent `1` and there is no Korec-type theorem, consistent with the
  expected divergent orbits.

## 4. Controls

`collatz_dipspectrum_20260926.py 24`: exact counts `D_b(2^t, gamma)` for
`t <= 24`, eight values of `gamma`, both sheets (agreeing to within a
handful); four-doubling slopes and the prediction `h(rho) - 1/(k ln 2)` for
the `1/k` prefactor of ballot-type counts. See the note's table.

## 5. Non-consequences

Nothing here excludes a divergent orbit or bounds the density of divergent
integers; the theorem is a statement about residue classes modulo
`2^(log_2 n)`, blind to the sheet and to defects.
