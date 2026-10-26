---
id: THM-4800
title: "Pepin's test as the valuation-two Collatz branch. Let A(x) = (3x+1)/4, the Syracuse step on x = 1 (mod 8); A fixes 1 (the trivial cycle) and is conjugate to y -> (3/4)y, y = x - 1. (1) On integers the Fermat line of THM-4793 is the A-orbit of F_m: S^k(F_m) = A^k(F_m) = 3^k 2^(2^m - 2k) + 1 for k < 2^(m-1). (2) For m >= 1, F_m is prime iff A has exactly two cycles on Z/F_m, the trivial fixed point 1 and one cycle of length 2^(2^m) through 0 = F_m, iff the half-period iterate A^((F_m - 1)/2) is the reflection x -> 2 - x about the trivial cycle; the integer Fermat line is the initial arc of that cycle. (3) The reflection x -> 2 - x commutes with A on Z and sends F_m to 1 - 2^(2^m), whose valuation-two run is the mirror image of the Fermat line. (4) For prime F_m, A is a cyclic self-complementing permutation of the Paley graph P(F_m) fixing 1. (5) Mersenne counterpart: B(x) = (3x+1)/2 fixes -1 (the negative trivial cycle) and carries the Mersenne line of THM-4556; modulo a Mersenne prime, B^((M_p - 1)/2) is the reflection about -1. This is necessary for primality and not known to be sufficient: no composite M_p with p < 1400 passes it. Among the 11 Mersenne primes with 3 <= p <= 127, 3/2 is a primitive root for 7. (6) F_(2^n) = 2 (mod F_n), so the Fermat lines of F_n and of F_(2^n) are always mirror images under x -> 2 - x. They lie on one A-cycle iff F_n is prime, and then the line of F_(2^n) is exactly the second half of the Pepin cycle, which is the dynamical form of kind-pasteur's concurrent landing-point test F_n | L_(2^n). For F_5, F_6 and F_7 they lie on different cycles."
status: >
  PROVED (elementary; proofs below). Builds on klein's THM-4793, which gives the odd-branch form x -> 3x+1 of Pepin's
  test and the Fermat line. The new parts are:
  - the valuation-two carrier A, with its fixed point at the trivial cycle;
  - the identification of the integer Fermat line with the initial arc of the Pepin cycle, which places THM-4793's
    prime-blind line inside the prime-sensitive cycle;
  - the half-turn form and the mirror line;
  - the Paley reading and the Mersenne counterpart;
  - (6), the mirror and antipodal-arc reading of kind-pasteur-S5's concurrent Theorem L (CITED; landing-point
    divisibility).
  INDEPENDENTLY AUDITED (05-knowledge/results/s58_fermat_pepin_collatz_20261026/AUDIT.md, MISTAKE-673): (1)-(5)
  confirmed and every number recomputed. (6) is corrected: "antipodal arcs of one cycle" holds only for prime F_n,
  while the mirror holds always. (5)'s "not sufficient" is corrected to "not known to be sufficient".
  FINITE-EXACT: collatz_pepin_valuation_two_20261026.py.
  - The line and mirror identities for m = 2..11.
  - The two-cycle structure and the half-turn for F_1..F_4.
  - Exact cycle counts for the composite F_5, F_6 and F_7 from their complete factorisations: 410, 1822 and 518
    cycles; ord(3/4) is short of 2^(2^m) by 8.59, 10.81 and 9.00 bits.
  - The Mersenne half-turn and ord(3/2) for the 11 Mersenne primes with 3 <= p <= 127.
  - (6) for n = 1..7: the mirror; the same-cycle test through the 2-adic valuations of ord_q(3/4); for n <= 5, a walk
    of the whole cycle through 0.
  CITED: Pepin 1877; Euler's criterion; quadratic reciprocity; Sachs 1962 and Ringel 1963 (self-complementing
  permutations).
session: opus-2026-10-26-S58
source: 05-knowledge/results/s58_fermat_pepin_collatz_20261026/README.md
depends_on:
  - 01-canon/theorems/THM-4793-the-fermat-line-under-collatz-valuation-two-descent-exact-merge-of-consecutive-lines-and-pepin-as-an-affine-cycle-type.md
related:
  - 01-canon/theorems/THM-4556-the-mersenne-line-is-a-chain-of-debt-states-odd-shift-distance-2-adic-periodicity.md
  - 01-canon/theorems/THM-4792-two-sheets-mod-odd-m-merged-centre-reflections-unit-shift-and-alternation-classification.md
  - 05-knowledge/results/fermat_five_pepin_collatz_20261025/PEPIN_COLLATZ.md
scripts:
  - 04-computation/fermat_pepin_s58/collatz_pepin_valuation_two_20261026.py (+ collatz_pepin_valuation_two.json)
---

# THM-4800: Pépin's test as the valuation-two Collatz branch

**Setting.** `S(x) = (3x+1)/2^(v_2(3x+1))` on odd `x`. On `x = 1 (mod 8)` the valuation is exactly 2 and
`S = A`, `A(x) = (3x+1)/4`. Since `A(1) = 1`, the trivial cycle `1 -> 4 -> 2 -> 1` is a valuation-two fixed point. The
map `A` is affine with multiplier `3/4` and is conjugate, by `y = x - 1`, to `y -> (3/4) y`. Write
`F_m = 2^(2^m) + 1` and `K = 2^m`.

**(1) The Fermat line is an A-orbit.** THM-4793's run law gives `S^k(F_m) = A^k(F_m) = 3^k 2^(K-2k) + 1` for
`0 <= k < K/2`. Indeed `A^k(x) = 1 + (3/4)^k (x - 1)`, and `x - 1 = 2^K`.

**(2) Pépin as a cycle type and as a half-turn.** Take `m >= 1`, so `F_m` is odd and prime to 3, and `A` is a
permutation of `Z/F_m`.
- If `F_m` is prime, then `(3/F_m) = -1` by reciprocity (`F_m = 2 (mod 3)`), and 4 is a square. So `3/4` is a
  quadratic non-residue, hence a generator of the cyclic group `F_(F_m)^*` of order `2^K`. Then `y -> (3/4)y` fixes 0
  and is one cycle on the `2^K` units, so `A` has exactly two cycles: `{1}` and a cycle of length `2^K` through every
  other residue.
- If `F_m` is composite, the nonzero nonunits are never in the orbit of a unit, so there are at least three cycles.
- `(3/4)^((F_m-1)/2) = (3/F_m)` by Euler's criterion. So for prime `F_m` the half-period iterate is `y -> -y`, i.e.
  `A^((F_m-1)/2)(x) = 2 - x`. Conversely, if that iterate is the reflection, then `3/4` has order exactly `2^K` modulo
  `F_m`; since `2^K` divides `ord_q(3/4)`, which divides `q - 1`, for every prime `q | F_m`, it follows that
  `q > 2^K`, so `q = F_m`.
- `2^(F_m - 1) = 1 (mod F_m)` because `ord_(F_m)(2) = 2^(m+1)` divides `2^K`. So `(3/4)^((F_m-1)/2) = 3^((F_m-1)/2)`, and
  the half-turn test is literally Pépin's congruence.
- `A^k(F_m) = A^k(0) (mod F_m)`. So the integer Fermat line, the first `K/2` points of the orbit of `F_m`, is the
  initial arc of the length-`2^K` cycle through `0`.

THM-4793 found the line *prime-blind*. That remains true of the line, but the line is the visible arc of a cycle whose
length decides primality. The integer orbit leaves the branch after `K/2 - 1` steps, at the valuation-three step,
because `|3/4|_2 = 4`: `A` repels 2-adically from its fixed point 1, and `F_m` is `2^(-K)`-close to 1 in `Z_2`.
Meanwhile `|3/4|_inf = 3/4` contracts the real values toward 1 and `|3/4|_3 = 1/3` contracts 3-adically (product
formula). Modulo `F_m`, `A` is a permutation whose period is Pépin's question.

**(3) The mirror.** `A(2 - x) = 2 - A(x)` identically. Over `Z`, `2 - F_m = 1 - 2^K` is again `1 (mod 8)`, so its
valuation-two run is `S^k(1 - 2^K) = 2 - S^k(F_m) = 1 - 3^k 2^(K-2k)`, the Fermat line reflected through the trivial
cycle. Modulo `F_m` the half-turn of (2) sends `F_m` (that is, 0) to this mirror point (that is, 2). This is the
valuation-two analogue of the reflection centres in THM-4792: the one-step maps reflect about `-+1/2`, and `A`
reflects about its fixed point.

**(4) Paley.** For prime `F_m`, `x -> ax + b` with `a` a non-residue maps the Paley graph `P(F_m)` onto its
complement. So `A` is a self-complementing permutation fixing the vertex 1. Such permutations have one fixed point and
cycles of length divisible by 4 (Sachs, Ringel). Pépin says this one is **cyclic**: a single cycle of length `2^K`.

**(5) The Mersenne counterpart.**
- On `x = -1 (mod 4)` the Syracuse step is `B(x) = (3x+1)/2`, which fixes `-1`, the negative trivial cycle.
  THM-4556's Mersenne line is `S^k(M_p) = B^k(M_p) = 3^k 2^(p-k) - 1`, and modulo `M_p = 2^p - 1` it is the initial
  arc of the `B`-orbit of 0.
- For a Mersenne prime with `p >= 3`: `M_p = 7 (mod 24)`, so `(3/M_p) = -1` and `(2/M_p) = 1`, and `3/2` is a
  non-residue. Euler's criterion then makes `B^((M_p-1)/2)` the reflection `x -> -2 - x` about `-1`.
- Since `2^((M_p - 1)/2) = 1`, the half-turn is the congruence `3^((M_p - 1)/2) = -1 (mod M_p)`, an Euler-type
  condition. It is necessary. Whether it is sufficient for `M_p` with `p` prime is not known: `M_p - 1` is not a power
  of two, so Pépin's argument does not apply, but no composite `M_p` with `p < 1400` (207 cases) passes.
- Data for the 11 Mersenne primes with `3 <= p <= 127`: `3/2` is a primitive root for `p = 3, 5, 7, 17, 19, 89, 107`,
  and `B` has 10, 4, 10 and 4 cycles for `p = 13, 31, 61, 127`. These data concern prime `M_p` and say nothing about
  sufficiency.
- The genuine Mersenne test, Lucas–Lehmer, lives instead in the norm-one torus of `F_(M_p^2)`, the cyclic group of
  order `M_p + 1 = 2^p`. That torus is the non-split analogue of `F_(F_m)^*`, of order `F_m - 1 = 2^K`.
- So the Fermat line (valuation two, contraction `3/4`) sees the split 2-power torus, while the Mersenne line
  (valuation one, expansion `3/2`) sits in the split group `F_(M_p)^*`, which is not a 2-group. That asymmetry is why
  only the Fermat half-turn is a *proved* primality test.

**(6) The two Fermat lines are mirror images, and antipodal arcs of one cycle exactly when `F_n` is prime.**
Concurrently, kind-pasteur-2026-10-26-S5 proved (Theorem L in
`05-knowledge/results/fermat_landing_points_pepin_20261026.md`) that `F_n` is prime iff `F_n` divides the Syracuse
landing point `L_(2^n) = (3^M + 1)/2` of the Fermat line of `F_(2^n)`, where `M = (F_n - 1)/2`. The present carrier
explains this dynamically.
- **The mirror, always.** `F_(2^n) = 4^M + 1`, and `4^M = 2^(F_n - 1) = 1 (mod F_n)` because `ord_(F_n)(2) = 2^(n+1)`
  divides `F_n - 1`. So `F_(2^n) = 2 (mod F_n)`. The Fermat line of `F_(2^n)`, which is
  `A^j(4^M + 1) = 3^j 4^(M-j) + 1` for `j < M`, reduces to `2 - A^j(F_n) (mod F_n)`. So the two lines are mirror
  images under `x -> 2 - x`, for every `n`, prime or not.
- **One cycle iff prime.** In `y = x - 1` coordinates the two starting points are `y = -1` and `y = +1`. They lie in
  one cycle iff `-1` is a power of `3/4` modulo `F_n`.
  - For prime `F_n` this holds by (2). The line of `F_(2^n)` (`M` points) is then exactly the **second half** of the
    Pépin cycle, from `A^M(0) = 2` to `A^(2M-1)(0)`.
  - The line of `F_n` covers only the first `2^(n-1)` points of the first half; at `n = 2` that is 2 of 8.
  - For composite `F_n` it can fail. At `F_5`, `F_6` and `F_7` the 2-adic valuations of `ord_q(3/4)` over the two prime
    factors are `(7, 4)`, `(6, 8)` and `(9, 7)`, so `-1` is not a power of `3/4` and the two lines lie on **different**
    cycles.
- **Theorem L.** The landing value of the `F_(2^n)` line before the final halving is `A^M(4^M + 1) = 3^M + 1`. Then
  `F_n | 3^M + 1` iff `A^M(2) = 0 (mod F_n)` iff `(3/4)^M = -1`. That is the half-turn of (2), i.e. the closing of the
  single Pépin cycle.
- **Checks.** `F_(2^n) = 2 (mod F_n)` for `n = 1..7`. The landing divisibility holds for `n <= 4` and fails at
  `n = 5, 6, 7`.
  - The pre-audit script tested only the mirror identity. That identity holds for every modulus, so it did not check
    the same-cycle claim (MISTAKE-673).
  - The corrected script tests that claim through the 2-adic valuations above, the audit's values.
  - For `n <= 5` it also walks the whole cycle through 0. At `n = 1..4`, 2 sits at position `M`. At `F_5` the cycle
    has length 11,167,360 and does not contain 2.

**Scope.** Everything here is elementary. The content is the identification of the carrier, not a new primality
test: (2) is equivalent to Pépin's congruence. Nothing is claimed about the existence of further Fermat primes.
kind-pasteur-S5 also independently audited THM-4793 (confirmed) and states, concurrently with the S58 note, the
algebraic and Aurifeuillean no-go for Fermat numbers.
