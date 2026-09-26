---
id: THM-4476
title: "Thin divergence: every orbit of x -> x/2, (3x+b)/2 (b odd) that is not eventually periodic has at most C X^(h*+eps) elements below X, h* = h(log_3 2) = 0.94996; hence its reciprocal sum converges (HYP-9160 PROVED), every divergent Collatz orbit is full-rate, and the real value of the Bernstein series of a non-periodic word is strictly below its 2-adic integer value"
status: >
  PROVED + INDEPENDENTLY AUDITED (elementary, full proof in the
  companion note; audit details in the audit field) + FINITE-EXACT controls of both
  counting lemmas. For odd b let T_b(x) = x/2 (x even), (3x+b)/2 (x odd) on
  Z, and let (x_i) be a T_b-orbit in Z\{0} that is not eventually periodic.
  Then for every eps > 0, #{i : |x_i| <= X} <= C(eps,|b|) X^(h*+eps) for all
  X >= 1, with h* = h(log_3 2) the binary entropy of log_3 2 and C
  independent of the orbit. Corollaries: (1) sum_i 1/|x_i| < infinity, so
  every divergent positive orbit of 3n+1 or 3n-1 has summable reciprocals
  over its odd iterates (HYP-9160); (2) for a non-eventually-periodic
  halving word d whose 2-adic Bernstein value R_2(d) is a rational with odd
  denominator, the real series R(d) = sum 2^(d_l)/3^(l+1) converges, and if
  R_2(d) = n is a positive integer then R(d) < n strictly, with the 3n-1
  orbit of n diverging at full rate c = n prod(1 - 1/(3m_l)) > 0; (3) no
  such orbit has m_j <= K j^a for all j with a < 1/h* = 1.05268, hence no
  plus-sheet log-band |Delta_j| <= C log_2 j with C < 0.02634, no
  minus-sheet log-band with C < 1.05268, and no bounded strip (the in-house
  no-bounded-strip theorem is recovered); (4) the union of all cycles of T_b, and more
  generally every T_b-invariant set on which T_b is injective, is thin in
  the same sense (a cycle with maximum M has at most C M^(h*+eps)
  elements); not proved for non-constant sign strategies by this method
  (the parity map is not a bijection there; witness -chi_{-4}); (5) for every divergent positive 3n+1
  orbit sum_l 2^(d_l)/3^l < infinity, so its discrepancy d_l - l log_2 3
  tends to -infinity; (6) all constants are uniform over orbits: sum 1/|x_i|
  <= K(b), hence c <= e^(K/3) n on the plus sheet and c >= kappa n on the
  minus sheet, i.e. |R(d)| <= K' |R_2(d)| with an absolute K'; (7) every
  infinite branch of the inverse tree is thin in the same sense; (8) an
  orbit is eventually periodic iff its reciprocal sum diverges, so the
  no-divergence half of Collatz (not the cycle half) is equivalent to
  sum_i 1/T^i(n) = infinity for every n >= 1, and the Dirichlet series of a
  divergent orbit has abscissa <= h*. The same proof gives, for every Conway/Matthews-Watts map
  g(x) = (p_i x + q_i)/m with p_i coprime to m, prod p_i < m^m and
  max p_i < m^2, the exponent max(1 - I(g), log_m max p_i/m) with I(g) the
  Chernoff rate of the multipliers (section 1.8 of the note). It does NOT
  exclude divergent orbits
  (expected to have N(X) of order log X) and does not bound the density of
  the union of divergent orbits.
  UPDATE 2026-09-26 (Lemma 1.4b and addendum 1.6b of the note, opus,
  post-audit and not covered by the audit): the eps is a small power of the
  logarithm, N(X) <= K(a,|b|) X^(h*) (log_2 X)^a for every
  a > lambda*/h* - 1/2 = 0.0138 (lambda* = 0.488077 the tilt of THM-4487),
  by the recursion (R) with theta = (1+eta) log_2 log_2 X/(h* log_2 X), the
  concavity bound h(rho) <= h* + lambda* theta, and the counting lemma with
  a geometric binomial tail (factor (log_2 X)^(-1/2) instead of log_2 X + 1);
  same for injective invariant sets (Cor. 4).
  UPDATE 2026-09-26 (THM-4499, opus): superseded by N(X) <= K X^(h*) (log_2 X)^a
  for every a > lambda*/h* - 3/2 = -0.9862, i.e. N(X) = o(X^(h*)): the
  no-dip count at the moving barrier carries THM-4495's ballot factor
  (log X)^(-3/2) (weighted Spitzer identity), which the recursion turns
  into the exponent a*. Divergent orbits remain not excluded.
  AUDIT of Lemma 1.4b and addendum 1.6b (separate agent, 2026-09-26;
  collatz_dipspectrum_20260926_orders_audit.py/.out sections (6)-(7)): every
  constant of 1.4b re-derived (ratio bound 0.8003 <= 0.81, 1/0.19 <= 5.3,
  p(1-p) >= 0.232858, factor 2, k >= 0.995 log_2 X, 25.66 <= 27; exact tail
  sums <= 5.82 against the claimed 12.8 for theta in (0, theta_1] and every
  k in [200, 3000]); two misrounded decimals fixed (rho >= 0.5654,
  (rho/(1-rho))^2 <= 2.93); 1.6b's concavity step, (R'), the identities and
  the induction CONFIRMED; presentation gaps recorded in the note
  (induction on floor(X); Cor. 4's recursion has #F_b + #F_(-b), so 27 -> 54;
  K is not effective). Verdict: SOUND.
  UPDATE 2026-09-26 (one-bit band lemma, found by the audit of the landing
  reassessment): every dipper of a landing point j has y_i in
  (2^(theta L) y_j, 2^(theta L + 1) y_j], and at most two of three consecutive
  orbit values lie in a 1-bit band, so the landing multiplicity in (D) is at
  most 2 ceil(k/3) rather than k; constants only.
source: collatz-squares-doubles-20260925 session (opus), 2026-09-25; the owner asked to prove HYP-9160 for discrepancy O(log l). Mechanism: the in-house no-bounded-strip theorem's Terras stopping-time count (collatz_guards_20260921_discrepancy.md, section 2a) plus a pigeonhole on landing points that replaces the strip hypothesis. No priority claimed.
depends_on:
  - 05-knowledge/results/collatz_sqdbl_20260925_squares_doubles_foundry.md (Proposition 6: the orbitwise two-place identity)
related:
  - 05-knowledge/hypotheses/HYP-9160-no-slow-divergence-real-2adic-coincidence.md (SETTLED by this theorem)
  - 05-knowledge/results/collatz_guards_20260921_discrepancy.md (no bounded strip; recovered as Corollary 3)
  - 05-knowledge/results/collatz_procgen_20260922_hard_class.md (Proposition T, two places for bounded-discrepancy words)
  - 05-knowledge/results/collatz_procgen_20260922_choice_ladder.md (dimension h(log_3 2) of the exceptional set)
  - 01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md (its strategy -chi_{-4} is the witness that the METHOD does not extend to non-constant shifts, and that unions of thin orbits can have positive density)
script: 04-computation/experiments/collatz_thin_20260925_counts.py
output: 05-knowledge/results/collatz_thin_20260925_counts.out
script_sha256: 721b2213ac70e3c2dd8d1e6fb109fd7c2c03a8d7918dcb1464d460c164d133b1
output_sha256: 9e3dd585e2fc46be83eb6ec0e709682b50009b0fb5f22838d404e4310fafd9a5
script_controls2: 04-computation/experiments/collatz_thin_20260925_controls2.py
output_controls2: 05-knowledge/results/collatz_thin_20260925_controls2.out
script_controls2_sha256: b7a2dbba3d4e9e95560734394f49d29ab90079497e214e0e438897dda6720cfc
output_controls2_sha256: 31a284957f46959d72c8c7ebd93716f8b568af658e1faca3b0073597d9f638a4
output_counts_to_2_24: 05-knowledge/results/collatz_thin_20260925_counts_to_2_24.out (same script with argument 24)
output_counts_to_2_24_sha256: e344d992f821fe2f801eba960b25992859ba587821ccef613ac8d2477036e822
hash_basis: raw LF bytes
audit: >
  Self-audit complete (every step re-derived; both counting lemmas
  checked numerically; a first draft's extension to non-constant sign
  strategies withdrawn, see MISTAKES 2026-09-25).
  Independent adversarial audit (separate agent, same day; script
  04-computation/experiments/collatz_thin_20260925_audit.py, output
  05-knowledge/results/collatz_thin_20260925_audit.out): verdict SOUND.
  Independent adversarial audit (separate agent, 2026-09-25; collatz_thin_20260925_audit.py/.out,
  written blind to counts.py). Theorem (1.1-1.6) and Cor. 1-3: SOUND. Every step re-derived; Terras
  bijection, carry bound, entropy step (exact binomial sums k<=400), E/D/ND partition,
  #(D)<=k N(X^(1-theta)), #(ND)<=#F checked on real segments (b=1,-1,5,-7); #F_{+-1}(2^20,theta)
  reproduced exactly; Prop. 6 identities exact for odd n<=300, both sheets. Cosmetic: (E) has at most
  k indices; distinctness is used in (ND) and 1.1, not (D). First-draft Cor. 4 (sign strategies) was
  WRONG and is withdrawn: the parity map is a bijection only for constant b (level-2 witnesses (+,-):
  word 11 unattained; (-,+): #F >= X/2). New Cor. 4 (injective invariant sets), 5, 6, 7, 9: sound
  (Cor. 4 needs #F_b + #F_{-b} and discards the <=(k+1)|S| points whose window meets S). Cor. 8
  OVERCLAIMED: sum_i 1/T^i(n)=infinity for all n is equivalent to "no divergent orbit", not to
  Collatz (a nontrivial cycle also has a divergent reciprocal sum). Controls: "growth exponents
  0.775/0.863" were log F/log X, not slopes; local slopes are 0.92-0.96 and 0.99-1.02, AT the bound.
  All three items were corrected in the text after the audit (same day); section 1.8 (general
  Matthews-Watts form) was spot-checked only, not audited line by line.
---

# THM-4476 -- thin divergence

**PROVED.** Full proof, corollaries and controls:
[collatz_thin_20260925_thin_divergent_orbits](../../05-knowledge/results/collatz_thin_20260925_thin_divergent_orbits.md).

## 1. Statement

Let `b` be odd, `T_b(x) = x/2` (`x` even), `(3x+b)/2` (`x` odd), and let
`(x_i)` be a `T_b`-orbit in `Z \ {0}` with pairwise distinct terms (i.e. not
eventually periodic). With `h* = h(log_3 2) = 0.949956...`:

```text
#{ i : |x_i| <= X } <= C(eps, |b|) X^(h* + eps)     for all X >= 1, every eps > 0.
```

## 2. Proof in five lines

1. **Terras.** The parity word of length `k` is a bijection `Z/2^k -> {0,1}^k`,
   and `T_b^i(y) = 3^o y/2^i + beta_i` with `|beta_i| <= |b|((3/2)^i - 1)`.
2. **Counting lemma.** For `theta in (0, 1 - (log_2 3)/2)` and
   `rho = (1-theta)/log_2 3 > 1/2`, the integers `y <= X` whose next
   `k = floor(log_2 X)` iterates never fall below `y X^(-theta)` have a word
   with at least `rho k - 2` odd letters (the carry is `O(X^0.585)`), so they
   number at most `2|b| X^(0.585+theta) + A(theta)(log_2 X + 1) X^(h(rho))`.
3. **Dichotomy.** A point `y_i <= X` of the orbit either dips below
   `y_i X^(-theta)` within `k` steps or lies in the set of step 2.
4. **Pigeonhole.** A dipper lands on a point of the *same* orbit below
   `X^(1-theta)`, and each landing point serves at most `k` dippers, so the
   dippers number at most `k N(X^(1-theta))`. Injectivity within the set
   is what this pigeonhole needs, which is why the argument does not extend
   to unions of orbits that merge.
5. **Bootstrap.** `N(X) <= k + k N(X^(1-theta)) + O(X^(h(rho)) log X)` gives
   `N(X) <= C X^gamma` for every `gamma > h(rho(theta))` by strong induction;
   `theta -> 0` gives `h*`. Sign changes happen only at odd `|x| <= |b|/3`,
   hence at most `|b|/3 + 1` times, so the orbit splits into boundedly many
   one-signed segments, each handled by conjugating `b -> -b`.

## 3. Corollaries

* **HYP-9160 (PROVED).** `sum 1/|x_i| < infinity` (Abel summation against
  `N(X) <= C X^gamma`, `gamma < 1`). Every divergent `3n+-1` orbit is
  full-rate: `m_L ~ c 3^L/2^(d_L)` with `0 < c < infinity`.
* **Two places.** If a non-eventually-periodic word `d` has 2-adic value
  `R_2(d) = n in Z_(>0)`, its real value satisfies `R(d) < n`. If
  `R_2(d) = -n`, then `R(d) < infinity`. For rational values with odd
  denominators `R(d)` converges.
* **Growth.** `m_j > j^a` infinitely often for every `a < 1.05268`; no
  bounded strip; no plus-sheet log-band with `C < 0.02634`; no minus-sheet
  log-band with `C < 1.05268`.
* **Periodic points.** The union of all cycles of `T_b` is thin, and so
  is every invariant set on which `T_b` is injective. Not proved for
  non-constant sign strategies: for `-chi_(-4)` all odd numbers share one
  word, so the counting lemma fails there.
* **Irrationality.** `R_2(d)` is irrational for every non-periodic word with
  `Delta_j >= -a log_2 j - O(1)`, `a < 1.05268` (bounded discrepancy and
  Sturmian words included).
* **Word of a divergent orbit.** `sum_l 2^(d_l)/3^l = 3R(d) < infinity` on
  the plus sheet: the leading factors are summable and the discrepancy tends
  to `-infinity`.
* **Uniformity.** `sum_i 1/|x_i| <= K(b)` for every orbit with distinct
  terms; so `n < c <= e^(K/3) n` (plus) and `kappa n <= c < n` (minus) with
  absolute constants, i.e. `|R(d)| <= K'|R_2(d)|`.
* **Inverse tree.** Every infinite branch of the inverse tree is thin; the
  proof only needs an injective `T_b`-chain.
* **Contracting Collatz-like maps.** Same statement and proof for every
  `g(x) = (p_i x + q_i)/m` with `p_i` coprime to `m`, `prod p_i < m^m`,
  `max p_i < m^2`; exponent `max(1 - I(g), log_m(max p_i/m))`.
* **Harmonic form.** Eventually periodic iff `sum 1/|x_i| = infinity`: the
  no-divergence half of Collatz says `sum_i 1/T^i(n) = infinity` for every
  `n >= 1` (a nontrivial cycle would satisfy it too). The union
  of thin orbits can have positive density (strategy `-chi_(-4)`), so the
  theorem does not bound the density of divergent integers; it does bound
  the density of periodic points.

## 4. Controls

`collatz_thin_20260925_counts.py`: exact binomial sums against `2^(k h(rho))`
for `k <= 320`; direct enumeration of the no-dip sets `F_(+-1)(X, theta)`
for `X <= 2^20` and `2^24`: `log F/log X = 0.775` (`theta = 0.03`) and `0.863`
(`theta = 0.10`) at `2^20`, with local doubling slopes `0.92`-`0.96` and
`0.99`-`1.02`, so the lemma's exponents `0.9635` and `0.9867` are essentially
attained by `F`; the two sheets
differ by at most four elements at every size. `collatz_thin_20260925_controls2.py`: the union of all
cycles of twelve maps `3n+b` (at most 86 periodic points below `2*10^5`); a
contracting `m = 3` Conway map (bijection mod `3^k`, rate `0.2513`, counts
below the bound); the parity map is a bijection mod `2^k` for `T_(+-1)` and
fails for `-chi_(-4)` (`k + 1` words among `2^k` residues).

## 5. Non-consequences

Divergent orbits are not excluded; they are expected to be exponentially
thin. The union of all divergent orbits is not bounded. `5n+1` is untouched
(`1/log_2 5 < 1/2`). Nothing here bears on cycles.
