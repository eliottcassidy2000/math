---
id: THM-4511
title: "Wall theorem for Gilbreath's automaton: a lone defect of size 4 at distance F from the leading 1, in a 0/2 sea, never crosses the first sea 2 on its left, whatever lies to its right; it destroys the leading 1 iff the F-1 cells between are all 0. Hence its extinction probability in a uniform random sea is exactly 2^(1-F), re-emissions from the stationary copy included."
status: >
  PROVED (left-looking rule; column sequences are {0,4}* then 2 then {0,2}^omega,
  propagated from the defect column leftward; the wall column stays 2 while its
  right neighbour is in {0,4}) + FINITE-EXACT (exhaustive over all 2^(F-1+10)
  sea contexts for F = 3..12: minimum column reached by a 4 equals F minus the
  trailing zeros, and p_4(F) = 2^(1-F) exactly; Monte Carlo at F = 14, 16
  consistent). SCOPE: one defect of size 4 with 0/2 everywhere else; several
  defects or sizes >= 6 are NOT covered (for a lone 6 the wall is breached to a
  4 and the exact extinction exceeds the front-only tail by 6-10 per cent).
source: opus-2026-09-26 session gilbreath-fermat-platonic-20260926
depends_on: [S9 note collatz_oscillation_gilbreath_20260926.md (fronts move left at speed 1, the 0/2 sea is XOR), Odlyzko 1993 (persistence heuristic, cited), Lucas (kernel support)]
verification: 04-computation/experiments/gilbreath_extinction_20260926.py -> gilbreath_extinction_20260926.out (wall test, exact table, Monte Carlo); 04-computation/experiments/gilbreath_fermat_tower_20260926.py part 5 (front path is a unit-triangular image of the sea bits)
---

# THM-4511 -- the wall theorem for a lone size-4 defect

## Statement

Run the absolute-difference automaton `a_(r+1)(i) = |a_r(i) - a_r(i+1)|` from
a row with `a_0(0) = 1`, `a_0(F) = 4`, and `a_0(i) in {0, 2}` for every other
column (finitely or infinitely many columns to the right of `F`). Let `c` be
the largest column with `1 <= c < F` and `a_0(c) = 2`, if any.

1. If `c` exists, then for every row `r` all columns `i <= c` are in `{0, 2}`
   and `a_r(0) = 1`: no entry `>= 4` ever reaches column `c` or anything left
   of it, whatever the cells to the right of `F` are.
2. If `c` does not exist, `a_(F-1)(1) = 4` and `a_F(0) = 3`.

Consequently a lone size-`4` defect destroys the leading `1` iff the `F - 1`
cells between them are all `0`; if those cells are uniform random in `{0, 2}`
the probability is exactly `2^(1-F)`, and the stationary copy's re-emitted
fronts contribute nothing.

## Proof

Column `i` depends only on columns `i` and `i + 1`. Columns `> F` start in
`{0, 2}` and stay in `{0, 2}` (`|x - y|` for `x, y in {0, 2}`). Column `F`
starts at `4`: while its right neighbour is `0` it stays `4`; the first time
the neighbour is `2` it becomes `2`; afterwards it stays in `{0, 2}`.

Call a column sequence *good* if it takes values in `{0, 4}` up to some row,
then possibly the value `2` once, then values in `{0, 2}` for ever (or values
in `{0, 4}` for ever). Column `F` is good. Suppose column `i + 1` is good and
column `i` starts at `0` (`c < i < F`). While column `i + 1` is in `{0, 4}`,
column `i` stays in `{0, 4}`. If column `i + 1` shows `2` at row `s`, column
`i` is in `{0, 4}` at row `s` and becomes `|x - 2| = 2` at row `s + 1`, after
which both columns are in `{0, 2}` and column `i` stays there. So column `i`
is good, and by induction columns `F - 1, ..., c + 1` are good.

Column `c` starts at `2`. While column `c + 1` is in `{0, 4}`, `|2 - x| = 2`
keeps it at `2`. When column `c + 1` shows `2`, column `c` becomes `0`, and
afterwards both are in `{0, 2}`, so column `c` stays in `{0, 2}` for ever.
Columns `< c` start in `{0, 2}` and their right neighbours stay in `{0, 2}`,
so by induction leftward they stay in `{0, 2}`; the leading entry is
`|1 - a(1)| = 1`. This proves 1.

For 2, on the alphabet `{0, 4}` the rule is XOR in units of `4`
(`|4 - 0| = 4`, `|4 - 4| = 0`), so over the zeros the defect's cone is
Pascal's triangle mod `2` and its left edge `(s, F - s)` equals `4` for every
`s <= F - 1`; the cone of `(F - 1, 1)` is `[1, F]`, so the right side is
irrelevant, and `a_F(0) = |1 - 4| = 3`. ∎

## Remarks

* The first front alone: the cells it meets are `L_s(b) = XOR_(j subset s)
  b(F - 1 - s + j)`, `s = 0..F-2`, a unit-triangular hence invertible
  `F_2`-linear image of the sea bits, so exactly `C(F-1, w)` sea words give a
  path of weight `w` (checked `F <= 12`); for size `4` the theorem shows this
  first-front count is the whole story.
* Several defects are not covered: a second `4` arriving at column `F` after
  the first has died can find column `c` at `0` and cross. Sizes `>= 6` are
  not covered: `|6 - 2| = 4` breaches the wall to a `4`, and the exact
  extinction probabilities (`gilbreath_extinction_20260926.out`) exceed the
  front-only tail `sum_(t<=j-2) C(F-1,t) 2^-(F-1)` by `6`-`10` per cent for
  size `6` and `4`-`8` per cent for size `8`.
* Consequence for the primes below `200000`: of the `29` fresh fronts in rows
  `1`-`64`, `23` have size `4`; the whole random-model risk of the triangle is
  concentrated in row `1` (a `4` at distance `3`, risk `1/4`) and row `2`
  (distance `8`, `1/128`), total `0.258`.
