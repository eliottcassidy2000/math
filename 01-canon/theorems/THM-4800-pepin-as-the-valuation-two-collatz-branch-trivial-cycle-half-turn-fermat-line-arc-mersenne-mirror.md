---
id: THM-4800
title: "Pepin's test as the valuation-two Collatz branch. Let A(x) = (3x+1)/4, the Syracuse step on x = 1 (mod 8); A fixes 1 (the trivial cycle) and is conjugate to y -> (3/4)y, y = x - 1. (1) On integers the Fermat line of THM-4793 is the A-orbit of F_m: S^k(F_m) = A^k(F_m) = 3^k 2^(2^m - 2k) + 1 for k < 2^(m-1). (2) For m >= 1, F_m is prime iff A has exactly two cycles on Z/F_m, the trivial fixed point 1 and one cycle of length 2^(2^m) through 0 = F_m, iff the half-period iterate A^((F_m - 1)/2) is the reflection x -> 2 - x about the trivial cycle; the integer Fermat line is the initial arc of that cycle. (3) The reflection x -> 2 - x commutes with A on Z and sends F_m to 1 - 2^(2^m), whose valuation-two run is the mirror image of the Fermat line. (4) For prime F_m, A is a cyclic self-complementing permutation of the Paley graph P(F_m) fixing 1. (5) Mersenne counterpart: B(x) = (3x+1)/2 fixes -1 (the negative trivial cycle) and carries the Mersenne line of THM-4556; modulo a Mersenne prime, B^((M_p - 1)/2) is the reflection about -1, which is necessary but not sufficient for primality (3/2 is a primitive root for 7 of the 11 Mersenne primes with p <= 127)."
status: >
  PROVED (elementary; proofs below). Builds on klein's THM-4793, which gives the odd-branch form x -> 3x+1 of Pepin's
  test and the Fermat line. The new parts are:
  - the valuation-two carrier A, with its fixed point at the trivial cycle;
  - the identification of the integer Fermat line with the initial arc of the Pepin cycle, which places THM-4793's
    prime-blind line inside the prime-sensitive cycle;
  - the half-turn form and the mirror line;
  - the Paley reading and the Mersenne counterpart.
  FINITE-EXACT: collatz_pepin_valuation_two_20261026.py.
  - The line and mirror identities for m = 2..11.
  - The two-cycle structure and the half-turn for F_1..F_4.
  - Exact cycle counts for the composite F_5, F_6 and F_7 from their complete factorisations: 410, 1822 and 518
    cycles; ord(3/4) is short of 2^(2^m) by 8.59, 10.81 and 9.00 bits.
  - The Mersenne half-turn and ord(3/2) for the 11 Mersenne primes with p <= 127.
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
- `M_p - 1 = 2(2^(p-1) - 1)` is not a power of two, so this half-turn is only an Euler-pseudoprime condition. Checked
  data for `p <= 127`: `3/2` is a primitive root for `p = 3, 5, 7, 17, 19, 89, 107`; `B` has 10, 4, 10 and 4 cycles for
  `p = 13, 31, 61, 127`.
- The genuine Mersenne test, Lucas–Lehmer, lives instead in the norm-one torus of `F_(M_p^2)`, the cyclic group of
  order `M_p + 1 = 2^p`. That torus is the non-split analogue of `F_(F_m)^*`, of order `F_m - 1 = 2^K`.
- So the Fermat line (valuation two, contraction `3/4`) sees the split 2-power torus, while the Mersenne line
  (valuation one, expansion `3/2`) sits in the split group `F_(M_p)^*`, which is not a 2-group. That asymmetry is why
  only the Fermat half-turn is a primality test.

**Scope.** Everything here is elementary. The content is the identification of the carrier, not a new primality
test: (2) is equivalent to Pépin's congruence. Nothing is claimed about the existence of further Fermat primes.
