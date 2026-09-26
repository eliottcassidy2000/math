---
id: THM-4487
title: "The dip spectrum of 3n+-1: #{n <= X : T^i(n) >= n^gamma for i <= log_2 n} = X^(h(gamma/log_2 3)+o(1)) for gamma in (log_4 3, 1] (Terras at gamma = 1, Korec at gamma = log_4 3), and the counting lemma of THM-4476 is sharp, so the thin-divergence exponent h(log_3 2) cannot be lowered by that lemma"
status: >
  PROVED + INDEPENDENTLY AUDITED (elementary; two-sided exponent on both sheets) +
  FINITE-EXACT controls to 2^24. Let b = +-1, T_b(x) = x/2 (x even),
  (3x+b)/2 (x odd), alpha = log_2 3, h the binary entropy, and
  D_b(X, gamma) = #{n <= X : T_b^i(n) >= n^gamma for all 0 <= i <= floor(log_2 n)}.
  (1) For every gamma in (log_4 3, 1], c X^(h(gamma/alpha))/log^(3/2) X <=
  D_b(X, gamma) <= C X^(h(gamma/alpha)) log X on both sheets (at gamma = 1 the
  rotation alone gives only S_i <= 0 and the negative minus-sheet carry
  needs a margin: the proof of record prepends one odd step to a rotated
  tail, the 223 crossroads repair, independently audited there; a tilted
  cycle lemma and a block construction are recorded alternatives); so
  log D_b/log X -> h(gamma/alpha);
  for gamma in (log_2(3/2), log_4 3] the limit is 1. So the exponent is
  h(log_3 2) = 0.949956 at gamma = 1 (Terras's undecided count) and rises to
  1 exactly at Korec's exponent log_4 3 = 0.792481; the constants
  0.949956, 0.050044, 0.792481, 0.207519, 0.48808 (the Chernoff tilt
  lambda* = log_3(1/log_2(3/2)) = -E'(1) exactly) and 1.052681 (= 1/E(1)) of the procgen, crossroads and opus lanes
  are values, slopes and reciprocals of the one function
  E(gamma) = h(max(1/2, gamma/alpha)).
  (2) The no-dip set F_b(X, theta) of THM-4476 satisfies
  c X^(h(rho))/log^(3/2) X <= #F_b(X, theta) <= C X^(h(rho)) log^2 X, rho = (1-theta)/alpha,
  theta in (0, 1 - alpha/2). Hence the recursion N(X) <= k N(X^(1-theta)) +
  #F_b(X, theta) cannot yield an exponent below h(log_3 2): THM-4476's
  exponent is optimal for its method. The exponent of an actual divergent
  orbit remains OPEN (conjecturally no polynomially dense divergent orbit
  exists). Upper bounds: Terras bijection + entropy count of words with
  prescribed odd density; lower bounds: the cycle-lemma rotation of THM-4478
  section 4 applied to no-dip prefixes. The q-analogue has rho = gamma/log_2 q
  < 1/2 for all gamma <= 1 once q >= 5, so no dip theorem of this kind exists
  for 5n+1. (3) General form (Theorem 4 of the note, PROVED, added after the
  audit and not covered by it): for every Conway map g(x) = (p_i x + q_i)/m
  (p_i coprime to m, a_i = p_i/m, max a_i > 1) and every 0 < gamma <= 1, the
  exponent of the dip count is the constrained maximum entropy
  E_g(gamma) = max{H_m(pi) : sum pi_i log_m a_i >= gamma - 1}, attained by
  the tilted law pi_i ~ a_i^lambda; E_g = 1 iff the uniform law meets the
  constraint, which for 3n+-1 is exactly Korec's log_4 3, and E_g(1) =
  1 - I(g) is the thin-divergence exponent of THM-4476's general form. Its
  multiplicative carry estimate (y_j = M_j n (1 + O(n^(-gamma) log n)) along
  a no-dip orbit) removes the hypothesis gamma > log_2(3/2) from Theorem 1:
  the statement holds for all gamma in (0, 1]. Sheet-blind, defect-blind; not a divergence exclusion.
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
  Targeted boundary audit (223 crossroads session, 2026-09-26): the gamma = 1
  lower bound's first-draft margin (n >= 3n) was invalid; repaired by
  prepending an odd step to a rotated tail (every positive prefix
  multiplier >= 3/2 on both sheets), checked on 3,142 rotated inputs and
  2,244 sheet realizations through length 16
  (crossroads223_20260926_boundary_audit.py; MISTAKES 2026-09-26). Full
  independent adversarial audit (separate agent, 2026-09-26; script
  04-computation/experiments/collatz_dipspectrum_20260926_audit.py, output
  05-knowledge/results/collatz_dipspectrum_20260926_audit.out): verdict SOUND.
  Independent adversarial audit (separate agent, 2026-09-26; script and output written blind
  to the control script): verdict SOUND after the author's same-day repair. Terras bijection,
  carry bounds, cycle lemma (exhaustive to length 12; good(t,o) >= C(t,o)/t to t = 16), the
  upper bound (0 violations of o >= rho t - 2 among all counted n <= 2^24, both sheets), the
  lower bound for gamma < 1 and for gamma = 1 on the plus sheet, and Prop. 2 (both bounds;
  #F reproduced exactly against the THM-4476 audit; construction 0 failures) all re-derived.
  A first-draft GAP at gamma = 1 on the minus sheet was confirmed (rotation gives only
  S_i <= 0; negative carry); the block repair was verified (c_L > 0 by irrationality of
  log_2 3; 0 failures on all tested concatenations) and the auditor's tilted cycle lemma
  (steps + alpha/t, o = ceil(t/alpha) + 2, checked t = 8..20) was adopted to restore the
  polynomial-log lower bound uniformly. Cor. 3 identities exact (lambda* = log_3(1/log_2(3/2))
  = -E'(1) = 0.488077; 1 - h* is the Chernoff rate). All 64 published counts per sheet
  reproduced (the script omits n = 1; the "0.792" column is gamma = 0.7925). Cosmetic items
  (a false crude binomial bound in the tools list, log powers, slack factors, rounding) were
  corrected in the text after the audit.
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
   satisfy `T^i(n) >= n 2^(-max(0,S_t)) - n^0.585 >= 3 n^gamma - n^0.585 >= n^gamma`
   for `gamma < 1`; the count is `C(t, o)/t = X^(h(gamma/alpha))/poly(log X)`.
   At `gamma = 1` the rotation gives only `S_i <= 0`; prepend one odd letter
   to a rotated tail of length `t - 1` (223 crossroads): every prefix
   multiplier is then `>= 3/2` and `T_b^i(n) >= (3/2) n - (3/2)^t >= n` on both
   sheets, with the same count `C(t-1, o)/(t-1)`. (Alternatives: tilted steps
   `x_j + alpha/t`, or concatenated rotated blocks at the cost of `X^(-o(1))`.)
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
  (`1 - alpha/2`: the drift), `0.48808` (`-E'(1) = log_3(1/log_2(3/2))`: the Chernoff tilt, exactly) and
  `1.052681` (`1/E(1)`: THM-4476's growth threshold).
* **THM-4476 is optimal for its method.** Any lower exponent for a
  divergent orbit needs a constraint on no-dip points beyond one
  `log_2 X`-window, i.e. small integers in residue classes modulo more than
  `X` (the transversality barrier).
* **`5n+1`.** `gamma/log_2 5 < 1/2` for every `gamma <= 1`: the dip count has
  exponent `1` and there is no Korec-type theorem, consistent with the
  expected divergent orbits.
* **Exact order at `gamma = 1` (THM-4495, post-audit).** `D_b(X, 1) = Theta(X^h (log X)^(-3/2))`
  on both sheets, and the residue count `|Bad_k|` obeys the exact identity
  `k W_k = sum_n B_n W_(k-n)` with binomial tails `B_n`; the bracket
  `log^(-3/2) .. log^(+1)` of (1) closes at its lower end for `gamma = 1`.
* **Sharper upper bound for every `gamma` (note 1.2b, post-audit).** A
  geometric binomial tail replaces `(t+1) max_o C(t,o)`:
  `D_b(X, gamma) <= C X^(h(rho)) (log X)^(-1/2)` for `gamma in (log_4 3, 1]`, so the
  bracket of (1) is `[log^(-3/2), log^(-1/2)]`, and (2)'s `log^2 X` becomes
  `(log X)^(-1/2)`; the same tail turns THM-4476's `eps` into `(log X)^(0.014+eps)`.
* **General form.** For every Conway map with `max p_i > m` and every
  `gamma in (0, 1]` the dip exponent is the constrained maximum entropy of
  the multiplier law (Theorem 4 of the note, post-audit); Korec's threshold
  is where the uniform law meets the constraint, the thin-divergence
  exponent is `E_g(1)`, and Theorem 1's carry condition is unnecessary.

## 4. Controls

`collatz_dipspectrum_20260926.py 24`: exact counts `D_b(2^t, gamma)` for
`t <= 24`, eight values of `gamma`, both sheets (agreeing to within a
handful); four-doubling slopes and the prediction `h(rho) - 1/(k ln 2)` for
the `1/k` prefactor of ballot-type counts. See the note's table.

## 5. Non-consequences

Nothing here excludes a divergent orbit or bounds the density of divergent
integers; the theorem is a statement about residue classes modulo
`2^(log_2 n)`, blind to the sheet and to defects.
