---
id: THM-387
name: lrc-arc-menu-structure
status: PROVED
date: 2026-06-01
session: claudebox-2026-06-01-S521
depends_on:
  - THM-384
  - THM-381
related:
  - HYP-1987
  - HYP-1993
---

# THM-387: Arc-confined half-turn tournaments are transitive backbones with a flipped up-set

## Setting

By THM-384, the marked observer of an `n`-runner Lonely Runner system is lonely
at time `t` iff all `m = n-1` moving runners lie in the closed safe arc
`[1/n, 1-1/n]` of length `L = (n-2)/n`.  At such a `t` the runner-runner
half-turn sub-tournament (the object S511/S512 track inside `A000568(m)`) is a
*half-turn tournament of `m` points confined to a circular arc of length `L`*.

This theorem characterizes those tournaments, independently of any speed box.

Let `m` distinct points lie inside a closed arc of length `L < 1` on `R/Z`.
Lift them to reals `x_1 < x_2 < ... < x_m` in an interval of length `L`
(no wraparound, since `L<1`).  The project's half-turn rule orients
`a -> b` iff `frac(x_a - x_b) in (0, 1/2)`.  Assume genericity: no two points are
exactly antipodal (`x_j - x_i != 1/2`), so every pair is oriented (a genuine
tournament; the non-generic antipodal case is the degenerate "tie wall" that
pollutes the S512 boundary-witness counts).

## Statement

1. **(Short/long rule.)** For `i < j`:
   `j` beats `i` if `x_j - x_i < 1/2` ("short" pair), and
   `i` beats `j` if `x_j - x_i > 1/2` ("long" pair).

2. **(Up-set structure.)** Hence the tournament is the transitive
   "later-beats-earlier" backbone with exactly the long pairs reversed.  The set
   `S = { (i,j) : i<j, x_j - x_i > 1/2 }` of reversed pairs is an **up-set of the
   interval-containment order** on pairs (`(i,j) <= (i',j')` iff `i' <= i` and
   `j' >= j`).  Up-sets of this (type-A root) poset number the Catalan number
   `C_m`.

3. **(Sub-`1/2` arcs are transitive.)** If `L <= 1/2` then `S = empty` and the
   tournament is transitive (`H = 1`).  In particular the `n=4` LRC arc
   (`L = 1/2`) admits only the transitive class; the `score=(0,1,1)` "class"
   reported by the S512 boundary scan is a degenerate antipodal tie, not a
   tournament.

## Proof

**(1)** Since `0 < x_j - x_i <= L < 1`, the real `x_j - x_i` is already in
`(0,1)`, so `frac(x_j - x_i) = x_j - x_i`.  Thus `j` beats `i`
(`frac(x_j - x_i) in (0,1/2)`) iff `x_j - x_i < 1/2`.  For the reverse,
`frac(x_i - x_j) = 1 - (x_j - x_i) in (0,1)`, so `i` beats `j`
(`frac(x_i - x_j) in (0,1/2)`) iff `1 - (x_j - x_i) < 1/2` iff `x_j - x_i > 1/2`.
By genericity exactly one of `<1/2`, `>1/2` holds, so the pair is oriented.

**(2)** If `(i,j) in S` and `i' <= i`, `j' >= j` with `i' < j'`, then
`x_{j'} - x_{i'} >= x_j - x_i > 1/2`, so `(i',j') in S`.  Thus `S` is an up-set.
The interval-containment order on `{(i,j): 1<=i<j<=m}` is the type-`A_{m-1}` root
poset; its order ideals (equivalently up-sets) are counted by `C_m`.

**(3)** If `L <= 1/2` then `x_j - x_i <= x_m - x_1 <= L <= 1/2`, and under
genericity `x_j - x_i < 1/2` for every pair, so `S = empty`: the backbone is the
transitive tournament, `H = 1`. ∎

## Significance

This turns the S511/S512 "reachable class inside `A000568(m)`" object into a
purely order-theoretic one: every LRC source class is a *flip-up-set* tournament.
The realizable up-sets at the LRC arc `L=(n-2)/n` (`> 1/2` for `n >= 5`), and the
resulting iso-class count (the **box-independent geometric menu**), are studied in
HYP-1993, which gives the exact closed form `menu(m) = A000016(m)` for `m >= 4`
and the labelled count `2^{m-1}`.

## Verification Record

`04-computation/lrc_arc_menu_geometry_s521.py` (+ `..._confirm_s521.py`) build
these tournaments by exact difference-constraint feasibility and reproduce the
S512 *non-degenerate* reachable classes exactly for `n = 5,6,7`, correcting the
`n=4` degenerate. Outputs: `05-knowledge/results/lrc_arc_menu_geometry_s521.out`,
`05-knowledge/results/lrc_arc_menu_confirm_s521.out`.
