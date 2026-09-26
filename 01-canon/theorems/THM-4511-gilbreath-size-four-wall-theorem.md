---
id: THM-4511
title: "Wall theorem for Gilbreath's automaton: in a row of 0s, 2s and 4s after the leading 1, no 4 ever crosses the first sea 2 on the left of the first 4, whatever 0/2/4 pattern lies to its right, and the leading 1 is destroyed iff the first 4 is preceded by zeros only. For a lone defect of size 4 at distance F the extinction probability in a uniform random 0/2 sea is exactly 2^(1-F), re-emissions from the stationary copy included; a finite Gilbreath triangle is decided at its first all-{0,2,4} row."
status: >
  PROVED (columns look only right; every column sequence has a {0,4} phase
  followed by a {0,2} phase, propagated from the right; the wall column stays 2
  while its right neighbour is in {0,4}) + INDEPENDENTLY AUDITED (the audit
  removed the first draft's restriction to a single 4, which the proof never
  used, and verified the multi-4 form exhaustively on all 3^L rows of length
  L <= 13 and on 200000 random rows of length 40) + FINITE-EXACT (lone defect:
  all 2^(F-1+R) sea contexts, R = 10 (9 at F = 12), F = 3..12: the minimum column
  reached by a 4 is F minus the trailing zeros and p_4(F) = 2^(1-F) exactly;
  Monte Carlo at F = 14, 16 consistent; certificates for the primes below 2e5,
  1e6, 1e7). SCOPE: values 0, 2, 4 only; entries >= 6 are NOT covered (for a
  lone 6 the wall is breached to a 4 and the finite-context extinction exceeds
  the front-only tail by 6-10 per cent). Novelty not checked against Odlyzko
  1993 beyond memory.
source: opus-2026-09-26 session gilbreath-fermat-platonic-20260926
depends_on: [S9 note collatz_oscillation_gilbreath_20260926.md (fronts move left at speed 1, the 0/2 sea is XOR), Odlyzko 1993 (persistence heuristic, cited), Lucas (kernel support)]
verification: 04-computation/experiments/gilbreath_extinction_20260926.py -> .out (wall test, exact table, Monte Carlo); gilbreath_certificate_20260926.py -> .out (first all-{0,2,4} row certificates); gilbreath_fermat_platonic_20260926_audit.py -> 05-knowledge/results/gilbreath_fermat_platonic_20260926_audit.out (multi-4 exhaustive check); gilbreath_fermat_tower_20260926.py part 5 (front path is a unit-triangular image of the sea bits)
---

# THM-4511 -- the wall theorem for size-4 defects

## Statement

Run the absolute-difference automaton `a_(r+1)(i) = |a_r(i) - a_r(i+1)|` from
a row with `a_0(0) = 1` and `a_0(i) in {0, 2, 4}` for every `i >= 1` (finitely
or infinitely many columns), containing at least one `4`. Let `F` be the first
column with `a_0(F) = 4` and let `c` be the largest column with `1 <= c < F`
and `a_0(c) = 2`, if any.

1. If `c` exists, then for every row `r` all columns `1 <= i <= c` are in
   `{0, 2}` and `a_r(0) = 1`: no entry `>= 4` ever reaches column `c` or
   anything left of it, whatever `0/2/4` pattern lies to the right of `c`.
2. If `c` does not exist, `a_(F-1)(1) = 4` and `a_F(0) = 3`.

Consequently a row of `0`s, `2`s and `4`s destroys the leading `1` iff its
first `4` is preceded by zeros only. For a lone size-`4` defect at distance
`F` whose `F - 1` left cells are uniform random in `{0, 2}` the probability is
exactly `2^(1-F)`, and the stationary copy's re-emitted fronts contribute
nothing.

**Corollary (light cone).** If row `r` has `a_r(1..t)` in `{0, 2, 4}`, the
leading `1` survives at least until row `r + t` unless the first `4` among
those cells is preceded by zeros only, in which case it is destroyed at row
`r + F`. In particular a finite Gilbreath triangle is decided at its first
all-`{0,2,4}` row `r_4`: the leading `1` survives every later row iff the
first `4` of row `r_4` is preceded by a `2`, and all danger to the leading `1`
comes from entries `>= 6`.

## Proof

Column `i` depends only on columns `i` and `i + 1`, and cell `(r, i)` only on
`a_0(i), ..., a_0(i + r)`; so for a statement about rows `<= R` and columns
`<= C` the row may be truncated at column `C + R`, and we may assume finitely
many `4`s. Call a column sequence *good* if it takes values in `{0, 4}` up to
some row and values in `{0, 2}` from the next row on (either phase may be
empty). Every column right of the last `4` is in `{0, 2}` for ever, hence
good.

If column `i + 1` is good, column `i` is good whatever its initial value.
While column `i + 1` is in its `{0, 4}` phase, column `i` stays in `{0, 4}` if
it started there (`|x - y|` for `x, y in {0, 4}`) and stays `2` if it started
at `2` (`|2 - 0| = |2 - 4| = 2`). From the first row where column `i + 1` is in
`{0, 2}`: if column `i` is `2` or `0` there, both columns are in `{0, 2}` from
then on; if column `i` is `4`, it stays `4` while its neighbour is `0`,
becomes `|4 - 2| = 2` at the neighbour's first `2`, and then both are in
`{0, 2}`. So the sequence of column `i` is a `{0, 4}` phase followed by a
`{0, 2}` phase, and by induction from the right every column is good.

Columns `c + 1, ..., F - 1` start at `0`, so their `{0, 4}` phases begin at
row `0`; column `c` starts at `2`, stays `2` while column `c + 1` is in
`{0, 4}`, becomes `0` at the first row where column `c + 1` is `2`, and stays
in `{0, 2}` for ever. Columns `1, ..., c - 1` start in `{0, 2}` and their right
neighbours stay in `{0, 2}`, so by induction leftward they stay in `{0, 2}`;
the leading entry is `|1 - a(1)| = 1`. This proves 1.

For 2, on the alphabet `{0, 4}` the rule is XOR in units of `4`
(`|4 - 0| = 4`, `|4 - 4| = 0`), so over the zeros the cone of the first `4` is
Pascal's triangle mod `2` and its left edge `(s, F - s)` equals `4` for every
`s <= F - 1`; the cone of `(F - 1, 1)` is `[1, F]`, so nothing to the right of
`F` matters, and `a_F(0) = |1 - 4| = 3`. The corollary follows by truncating
row `r` at column `t`. ∎

## Remarks

* The first front alone: the cells it meets are `L_s(b) = XOR_(j subset s)
  b(F - 1 - s + j)`, `s = 0..F-2`, a unit-triangular hence invertible
  `F_2`-linear image of the sea bits, so exactly `C(F-1, w)` sea words give a
  path of weight `w` (checked `F <= 12`); for size `4` the theorem shows this
  first-front count is the whole story.
* The audit's structural reading: on `{0, 2, 4}`, `|a - b| = a + b mod 4`, so
  the `2`-pattern is an autonomous XOR automaton and `4`s propagate by XOR only
  through cells whose two source cells are not `2`.
* Sizes `>= 6` are not covered: `|6 - 2| = 4` breaches the wall to a `4`, and
  the finite-context extinction probabilities (`gilbreath_extinction_20260926.out`,
  `R = 10`) exceed the front-only tail `sum_(t<=j-2) C(F-1,t) 2^-(F-1)` by
  `6`-`10` per cent for size `6` (values `R`-stable to seven digits) and
  `4`-`8` per cent for size `8` (values still moving in the fifth digit at
  `R = 14`).
* Certificates: for the primes below `200000` the first all-`{0,2,4}` row is
  `r_4 = 59` (first all-`0/2` row `65`), below `10^6` it is `95` (`= r_2`),
  below `10^7` it is `132` (`r_2 = 135`); in each case the first `4` of row
  `r_4` is preceded by a `2`, and no computed row has a first `4` preceded by
  zeros only. The random-model risk numbers of the S9 heuristic (front-only
  tails, sum `0.258` for the fresh fronts of rows `1`-`64`) are not
  consequences of this theorem: those rows contain entries `>= 6`.
