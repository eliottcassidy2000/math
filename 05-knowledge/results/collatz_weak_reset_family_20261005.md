# The excess-rate family F_q: weak resets paid by a 1/q deadline

2026-10-05, opus session `opus-2026-10-05-S8` (weak-reset family). Owner's seed:
continue the proved-family program whose "next precise target is a cost bound
admitting weaker resets such as 55 -> 47" (Codex,
[run-block family, section 6](collatz_run_block_family_20261005.md)).

**Status.** PROVED (elementary): the positive-cycle-point characterization of
first descents (Proposition R), the segment inequality and the two source-only
deadlines (Theorem D), the compiled floor on `F_q` through the backward
compiler, nesting and the inheritance of `q` along running minima, and
containment of the run-block tree in `F_3`. FINITE-EXACT: all declared
checks (1,602,929), including the exhaustive segment audit below `2^16`.
VERIFIED: the census of `q(n)` below `2^20` (float excess, exactly re-verified
at the reported sources). **OPEN:** universal entry; no single `F_q` covers
every source, and the union over `q` is exactly the rooted sources. Not a canon
promotion.

## 0. Answer

A weak reset is paid by its excess rate. For a first-descent segment of length
`l` (odd steps) and total valuation `A`, define the excess `A - l log_2 3 > 0`
and call the segment `q`-admissible when `2^(qA-l) >= 3^(ql)`, i.e. the excess
is at least `l/q`. Let `F_q` be the odd sources all of whose chained first-
descent segments (the running minima of the orbit, down to ROOT) are
`q`-admissible, and `q(n)` the least such `q`. Then:

- every `n` in `F_q` has the source-only deadline
  `tau(n) <= floor(1.051 q log_2 n) + max{tau(y): y odd, y < 10q}`, hence,
  through Codex's backward compiler, the positive floor
  `W(n) >= eta(n, T'_q(n)) = 2/((C+2) binom(C+1, floor((C+1)/2)))`;
- the Codex run-block tree (strong guard `a >= r+2`) lies inside `F_3`, and
  `F_3` is eight times larger below `2^20` (31,920 members against 3,859);
- `55 -> (1,1,3) -> 47` is admissible at `q = 13`, but 55 belongs to `F_31`,
  not earlier, because its next running minimum 47 opens the 34-step
  excursion to 23 whose excess rate is `1/31`; `7` enters at `q = 7`
  (segment `(1,1,2,3)`), `27` at `q = 104` (its 37-step excursion), `739` at
  `q = 1`;
- the families are nested, closed under passing to running minima, and
  exhaust the rooted sources as `q -> infinity`; their complement at any
  finite `q` is the set of sources whose running-minimum chain contains a
  segment of excess rate below `1/q`, and these are 2-adic shadows of large
  positive rational cycle points (Proposition R): weak resets are
  near-cycles;
- the price of coverage is the floor: at `q = 31` the compiled floor of 47 is
  `10^(-94.5)` against the actual weight `10^(-9.8)`, at `q = 104` the floor
  of 27 is `10^(-242)`; the deadline is tight only up to the factor
  `q`-at-the-worst-segment, applied to all of `log_2 n`.

So the requested cost bound exists, is explicit, and admits every exact
first-reset descent at some finite price; what it cannot do is keep the price
bounded, because the excess rate of a first-descent segment can be
arbitrarily small (below `2^20` the smallest is `1/6951`, at `n = 432923`,
a 94-step segment). This is the same obstruction as before in a sharper
coordinate: `q(n)` is a source-intrinsic measure of how close its worst
excursion is to a positive rational cycle.

## 1. Inheritance and the board

- Closest proved mechanism: the run-block family's block `(1^r, a)` with its
  strong guard `a >= r+2`, expense `E = sum(a_i - 1)` and ceiling
  `4(4/3)^E <= n-1` ([run-block family](collatz_run_block_family_20261005.md),
  (1)-(9)); the deadline-to-floor compiler with the hidden valuation budget
  `2^A <= n (10/3)^tau`
  ([backward measurement compiler](collatz_backward_measurement_compiler_20261005.md),
  sections 2-3); the first-descent induction of
  [inductive floor receipts](collatz_inductive_floor_receipts_20261005.md).
- Canonical hostiles: 7 (two rises, a valuation-two reset to 13 still above 7)
  and 27 (the 37-step excursion), both excluded from the run-block tree.
- Corrected near miss: this author's claim, corrected by Codex's audit of the
  previous session, that a deadline alone gives a floor; the sibling budget
  `K` is a hidden coordinate, which is why the floor below is compiled from
  the deadline through the valuation budget and not from `2/((T+1)(T+2))`.
- Least-used sidecar: the positive rational cycle point of a non-rising word,
  dual to the negative points of the shadow theorem.
- Board: **first-descent segment / excess rate / running minimum / positive
  cycle point / source-only deadline / compiled floor / coverage proportion.**

The audit of the previous note (Codex, mistakes ledger 2026-10-05,
measurement-independence audit) is accepted in full: the original
Proposition F quantified over models it did not construct, the `-8u` bound is
a certified lower estimate rather than the readout, the deadline converse
reversed an upper inequality, and coefficient-cone exit is not a ROOT
deadline. Nothing in this note depends on the withdrawn statements.

## 2. Proposition R: first descents are shadows of positive rational cycle points

Let `w = (a_1, ..., a_l)` be a valuation word with `A = sum a_i`, carry
`c_w = sum_i 3^(l-i) 2^(a_1+...+a_(i-1))`, and suppose `w` is non-rising,
`2^A > 3^l`. Put `x_w = c_w/(2^A - 3^l) > 0`.

**Proposition R (PROVED).** For an odd `x` that follows `w`, the endpoint
`y = (3^l x + c_w)/2^A` satisfies `y < x` iff `x > x_w`. The sources following
`w` are the cylinder `x = x_w mod 2^(A+1)`; its members below `x_w` do not
descend along `w` (these are the exceptional members of THM-4512), and the
rational point `x_w` is a positive cycle of the `3x+1` map on `Q` with period
word `w`.

*Proof.* `y < x` is `c_w < (2^A - 3^l) x`. The cylinder statement is the
shadow theorem's forward half, valid for any word; `x_w` is the fixed point
of the affine map `x -> (3^l x + c_w)/2^A`. QED.

The exhaustive audit below `2^16` (321,236 chained segments) confirms the
non-rising words, the descent criterion, the forward formula and the cylinder
residue, and tests the non-descent of cylinder members below `x_w`. Example:
`55 -> 83 -> 125 -> 47` has `w = (1,1,3)`, `c_w = 19`, `2^5 - 3^3 = 5`,
`x_w = 19/5 = 3.8`; the block of Codex's section 6 is the shadow of the
rational cycle `19/5`. The weaker the reset (the smaller `2^A - 3^l` relative
to `c_w`), the larger the cycle point and the more cylinder members fail to
descend. Together with the shadow theorem this gives one picture: rising
words are shadows of negative rational cycles (ancestors from below),
non-rising words are shadows of positive rational cycles (descents from
above), and the integers sit between the two.

## 3. Theorem D: the family and its deadline

**Definition.** For an odd `x > 1` its first-descent segment is the orbit
prefix up to the first value `y < x`. Its word `w`, length `l` and valuation
`A` are as above; it is `q`-admissible iff `2^(qA - l) >= 3^(ql)`. The chain
of `n` is the sequence of segments started at the successive running minima
`n = x_0 > x_1 > ... > 1`. `F_q` is the set of odd `n` whose chain consists
of `q`-admissible segments, and `q(n)` is the least admissible `q` over the
chain (`q(1) = 0`).

**Lemma D1 (segment inequality; PROVED).** If a segment from `x` is
`q`-admissible and `x >= q`, then `y^(2q) 2^l <= x^(2q)`, i.e.
`log_2 y <= log_2 x - l/(2q)`. If `x >= 10q`, then
`log_2 y <= log_2 x - 0.952 l/q`.

*Proof.* Along the segment every state satisfies `x_i >= x`, so
`2^(a_i) x_i = 3x_(i-1) + 1 <= 3x_(i-1)(1 + 1/(3x))` and
`y <= x (3^l/2^A)(1 + 1/(3x))^l <= x 2^(-l/q) e^(l/(3x))`. For `x >= q`,
`e^(l/(3x)) <= 2^(l/(2q))`; for `x >= 10q`, `<= 2^(0.048 l/q)`. QED. (Checked
exactly as `y^(2q) 2^l <= x^(2q)` on all 317,954 segments below `2^16` with
`x >= q`.)

**Theorem D (source-only deadline; PROVED).** Let `n` be in `F_q` and let
`t_q = max{tau(y): y odd, y < q}`, `t_(10q)` likewise below `10q` (finite
tables; every odd source below `2^20` is rooted, which covers `q <= 10^5`).
Then

    tau(n) <= floor(2 q log_2 n) + t_q,
    tau(n) <= floor(1.051 q log_2 n) + t_(10q).

*Proof.* Sum Lemma D1 over the segments whose start is at least the threshold:
the total length is at most `2q log_2 n`, respectively `1.0504 q log_2 n`,
because the logarithms telescope from `log_2 n` down to a positive number.
Once a running minimum falls below the threshold, the remaining odd steps are
those of that source's own orbit, bounded by the table. QED. (Both deadlines
hold for all 524,287 odd sources below `2^20` with their actual `q(n)`; the
sharp one is attained to a ratio `0.692` at worst.)

**Corollary D2 (compiled floor; PROVED given Theorem D).** With
`T = T'_q(n)`, Codex's budget gives `2^A <= n (10/3)^T`, `N = L + K <= C(n,T)`
and `W(n) >= eta(n,T) = 2/((C+2) binom(C+1, floor((C+1)/2))) > 0`; for the
designated leaf `z = rho(n)` the same at `T+1`, and the localized selector
degree `d` with `8(16/25)^d <= eta_z/2`. No ROOT word of `n` enters the
formula; membership in `F_q` is decided by at most `T'_q(n)` forward steps.

**Proposition D3 (structure; PROVED).** `F_q` is contained in `F_(q+1)`; if
`y` is a running minimum of `n` then `q(n) >= q(y)`, so `F_q` is closed under
passing to running minima; `n` is in `F_q` iff its segments above `y` are
admissible and `y` is in `F_q`; the union of all `F_q` is the set of rooted
sources, and `q(n)` is finite exactly when `n` is rooted. The run-block tree
is contained in `F_3`: its blocks are first-descent segments with excess
`>= 0.415 (r+1) >= (r+1)/3`, and its base `5 -> 1` is `1`-admissible.

## 4. Numbers

Examples (chain of segments with per-segment least `q`):

| `n` | `q(n)` | `tau` | chain |
|---:|---:|---:|---|
| 5 | 1 | 1 | `5 -(4)-> 1` |
| 7 | 7 | 5 | `7 -(1,1,2,3)-> 5`, then `5 -> 1` |
| 23 | 2 | 4 | `23 -(1,1,5)-> 5` |
| 739 | 1 | 4 | `739 -(1,8)-> 13 -(3)-> 5` |
| 55 | 31 | 41 | `55 -(1,1,3)-> 47` [q 13], `47 -(34 steps)-> 23` [q 31], `23 -> 5` |
| 47 | 31 | 38 | the same from 47 |
| 27 | 104 | 41 | `27 -(37 steps)-> 23` [q 104], then as above |
| 1161 | 169 | 66 | a 66-step chain with one segment of rate `1/169` |

Deadlines and compiled floors (`T'_q` is the sharp deadline):

| `n` | `q` | `tau` | `T'_q` | `C` | `log10 eta` | `log10 W(n)` | leaf | degree |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 1 | 1 | 8 | 10 | -3.4 | -0.5 | 3 | 26 |
| 23 | 2 | 4 | 15 | 21 | -6.9 | -2.1 | 15 | 44 |
| 739 | 1 | 4 | 16 | 25 | -8.1 | -2.8 | 15765 | 53 |
| 7 | 7 | 5 | 61 | 83 | -25.9 | -1.9 | 9 | 142 |
| 47 | 31 | 38 | 226 | 310 | -94.5 | -9.8 | 501 | 499 |
| 55 | 31 | 41 | 234 | 321 | -97.8 | -10.8 | 1173 | 517 |
| 27 | 104 | 41 | 584 | 800 | -242.2 | -10.1 | 27 | 1256 |
| 1161 | 169 | 66 | 1874 | 2568 | -774.7 | -16.8 | 1161 | 4003 |

For comparison, the run-block expense ceiling gives 23 the floor `w(5,3) =
1/420` and degree 26; the excess-rate floor at `q = 2` is `10^(-6.9)` with
degree 44. The strong guard is sharper where it applies; `F_q` applies where
the guard does not.

Coverage below `2^20` (proportion of odd sources with `q(n) <= q`):

| `q` | 1 | 2 | 3 | 4 | 8 | 16 | 32 | 64 | 100 | 200 | 500 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| share | 0.0016 | 0.0032 | 0.0609 | 0.0665 | 0.1487 | 0.3656 | 0.5105 | 0.5770 | 0.9508 | 0.9734 | 0.9961 |

The staircase is the inheritance of Proposition D3: a long excursion with a
small excess rate (the 34-step `47 -> 23`, rate `1/31`; the 37-step
`27 -> 23`, rate `1/104`) is a running minimum for every source whose chain
passes through it, so `q(n)` jumps at those values. Below `2^20`, 42.7% of the
sources have `q(n) >= 50`; the largest value is `q = 6951` at `n = 432923`,
whose critical segment has 94 steps and excess `94/6951`.

## 5. What this does and does not complete

- It completes the stated next target: a source-bounded expense that pays
  every exact first-reset descent, with the price `1/q` equal to the inverse
  excess rate, an explicit deadline, an explicit compiled floor, and a
  decidable recognizer. The run-block tree is the `q <= 3` sub-forest with
  single-reset blocks and unit parents.
- It does not approach universal coverage: for every `q` the complement is
  nonempty, and by Proposition R it consists of the 2-adic shadows of
  positive rational cycle points with small `2^A - 3^l`, together with
  everything that inherits them through running minima. Exhaustion of the
  rooted sources is only in the limit `q -> infinity`, where the deadline and
  the floor degenerate.
- The honest quantitative content is the coverage curve against `q` and the
  price curve `log10 eta` against `q`; both are explicit and reproducible.
- Next admissible targets: (i) a segment-adaptive deadline that charges each
  segment at its own rate without reading the word, which would require a
  source-only bound on the excess of the *next* running minimum, i.e. the
  deadline of the orbit itself; (ii) the distribution of `q(n)`, in
  particular whether the proportion with `q(n) > q` decays like a power of
  `q` (the census suggests a slow decay with plateaus); (iii) whether the
  popular hard segments (`47 -> 23`, `27 -> 23`) can be certified by a
  separate finite table so that `F_q` is replaced by `F_q` relative to a bank
  of excursions, which is the recycling the program allows for finitely many
  sources.

## 6. Reproduction

[Script](../../04-computation/experiments/collatz_weak_reset_family_20261005.py),
[output](collatz_weak_reset_family_20261005.out),
[JSON](collatz_weak_reset_family_20261005.json):

```text
python3 04-computation/experiments/collatz_weak_reset_family_20261005.py --census-bits 20 --json 05-knowledge/results/collatz_weak_reset_family_20261005.json
python3 -O 04-computation/experiments/collatz_weak_reset_family_20261005.py --census-bits 16
```

1,602,929 checks in 100 s. Exact integer and Fraction arithmetic for
Proposition R, Lemma D1 and the floors; numba float excess for the census,
with exact re-verification of `q(n)` at 7, 27, 55 and the maximizer. The
run-block recognizer is re-implemented from the published grammar (blocks
`(1^r, a)`, `a >= r+2`, unit parents, base 5), not imported. Hostiles: the
cylinder members below `x_w` must not descend; the deadlines must hold at every
source with its actual `q(n)`; the floors must not exceed the actual weights.
