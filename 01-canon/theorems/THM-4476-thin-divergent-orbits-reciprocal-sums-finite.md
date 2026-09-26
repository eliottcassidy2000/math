---
id: THM-4476
title: "Thin divergence: every orbit of x -> x/2, (3x+b)/2 (b odd) that is not eventually periodic has at most C X^(h*+eps) elements below X, h* = h(log_3 2) = 0.94996; hence its reciprocal sum converges (HYP-9160 PROVED), every divergent Collatz orbit is full-rate, and the real value of the Bernstein series of a non-periodic word is strictly below its 2-adic integer value"
status: >
  PROVED (elementary, full proof in the companion note; self-audited;
  independent audit status in the audit field) + FINITE-EXACT controls of both
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
  no-bounded-strip theorem is recovered); (4) the same holds for every
  level-k sign strategy of THM-4474. It does NOT exclude divergent orbits
  (expected to have N(X) of order log X) and does not bound the density of
  the union of divergent orbits.
source: collatz-squares-doubles-20260925 session (opus), 2026-09-25; the owner asked to prove HYP-9160 for discrepancy O(log l). Mechanism: the in-house no-bounded-strip theorem's Terras stopping-time count (collatz_guards_20260921_discrepancy.md, section 2a) plus a pigeonhole on landing points that replaces the strip hypothesis. No priority claimed.
depends_on:
  - 05-knowledge/results/collatz_sqdbl_20260925_squares_doubles_foundry.md (Proposition 6: the orbitwise two-place identity)
related:
  - 05-knowledge/hypotheses/HYP-9160-no-slow-divergence-real-2adic-coincidence.md (SETTLED by this theorem)
  - 05-knowledge/results/collatz_guards_20260921_discrepancy.md (no bounded strip; recovered as Corollary 3)
  - 05-knowledge/results/collatz_procgen_20260922_hard_class.md (Proposition T, two places for bounded-discrepancy words)
  - 05-knowledge/results/collatz_procgen_20260922_choice_ladder.md (dimension h(log_3 2) of the exceptional set)
  - 01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md (Corollary 4 applies to its strategies)
script: 04-computation/experiments/collatz_thin_20260925_counts.py
output: 05-knowledge/results/collatz_thin_20260925_counts.out
script_sha256: 721b2213ac70e3c2dd8d1e6fb109fd7c2c03a8d7918dcb1464d460c164d133b1
output_sha256: 9e3dd585e2fc46be83eb6ec0e709682b50009b0fb5f22838d404e4310fafd9a5
hash_basis: raw LF bytes
audit: >
  Self-audit complete (every step re-derived; both counting lemmas
  checked numerically). An independent adversarial audit by a separate
  agent was launched at the time of this commit; its result is recorded
  here in the following commit.
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
   dippers number at most `k N(X^(1-theta))`. This is where distinctness of
   the orbit is used, and why the argument does not extend to unions of
   orbits.
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
* **Strategies.** Same for every residue-determined odd shift `b(x)`, in
  particular all strategies of THM-4474.

## 4. Controls

`collatz_thin_20260925_counts.py`: exact binomial sums against `2^(k h(rho))`
for `k <= 320`; direct enumeration of the no-dip sets `F_(+-1)(X, theta)`
for `X <= 2^20`, growth exponents `0.775` (`theta = 0.03`) and `0.863`
(`theta = 0.10`) against the bounds `0.9635` and `0.9867`; the two sheets
differ by at most four elements at every size.

## 5. Non-consequences

Divergent orbits are not excluded; they are expected to be exponentially
thin. The union of all divergent orbits is not bounded. `5n+1` is untouched
(`1/log_2 5 < 1/2`). Nothing here bears on cycles.
