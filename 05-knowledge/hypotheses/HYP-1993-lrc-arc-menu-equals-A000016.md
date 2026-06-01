# HYP-1993: the box-independent LRC arc menu is A000016(m) = (1/2m) Σ_{d|m odd} φ(d) 2^{m/d}

**Session:** claudebox-2026-06-01-S521
**Status:** OPEN (closed form VERIFIED exactly for m=4..10; structure THM-387 PROVED;
supersedes the menu/growth part of HYP-1987).

## Setup

By THM-384 + THM-387, at a lonely time the `m = n-1` movers form a half-turn
tournament that is the transitive backbone with a flipped **up-set** of "long"
pairs (separation `> 1/2`), realized by points in a circular arc of length
`L = (n-2)/n`.  Define the **geometric menu** `menu(m)` = number of
iso-classes of such tournaments realizable in an arc of length `L`.  This is the
*box-independent* version of S512's arithmetic "reachable menu" (HYP-1987): it
does not depend on any finite speed box, only on the arc.

## Claims

**(A) Closed form.** For every `m >= 4` (i.e. `n >= 5`),
```
menu(m) = A000016(m) = (1/(2m)) * Σ_{d | m, d odd} φ(d) * 2^{m/d}.
```
`A000016` = "number of distinct output sequences from a binary m-stage shift
register feeding back the complement of the last stage" (a classical
complemented-necklace count).
VERIFIED exactly, `menu(m) = A000016(m)`:

| m  | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|----|---|---|---|---|---|---|----|
|menu| 2 | 4 | 6 |10 |16 |30 | 52 |
|A000016|2|4|6|10|16|30|52|

Predicts `menu(11)=94`, `menu(12)=172`.  The boundary case `m=3` (`n=4`,
`L=1/2`) is exceptional: `menu(3)=1` (transitive only, THM-387.3), vs
`A000016(3)=2`.

**(B) Labelled count.** The number of geometrically **realizable flip-up-sets**
(labelled, before isomorphism) is exactly `2^{m-1}` for `m >= 4` (and `1` at
`m=3`, `L=1/2`).  Verified m=4..11: 8,16,32,64,128,256,512,1024.

**(C) L-independence — REFUTES HYP-1987's "menu grows with L".** Both `menu(m)`
and the labelled count `2^{m-1}` are **constant for all `L in (1/2,1)`**: they
jump from `1` (transitive) to their full value the instant `L` exceeds `1/2`,
and do not grow further as `L -> 1`.  Verified at `L = 0.505` and `L = 0.995`
for m up to 9, plus a fine sweep.  So the menu depends only on the predicate
`L > 1/2`, not on the value of `L`.  HYP-1987 predicted the reachable menu grows
with `L`; that is false for the geometric (true) menu.

## Relation to HYP-1987 / S512

The geometric menu **equals the non-degenerate arithmetic-reachable menu** of
S512 exactly for `n = 5,6,7` (`menu = 2,4,6` matches the genuine-tournament
classes), and corrects `n=4`.  S512's larger counts (2,2,6,6) included
degenerate antipodal "tie wall" near-tournaments (score-sums `< C(m,2)`), which
THM-387 excludes.  So HYP-1987's core equivalence (the true target is the
arc-confined menu) is CONFIRMED and sharpened; only its "grows with L" growth
claim and its n=4 count are corrected.

## Why A000016? (open structural question — the real lead)

`A000016(m)` counts cyclic equivalence classes under shift+complement.  The
half-turn rule carries an antipodal (complement) symmetry, and the menu count is
a *cyclic/necklace* count even though the points lie in an arc, not the full
circle.  **Open:** find an explicit bijection between flip-up-set tournament
iso-classes and complemented binary necklaces of length `m`.  Such a bijection
would (i) prove (A) with a closed form, (ii) explain why `2^{m-1}` labelled
objects collapse to `A000016(m)` classes, and (iii) hand the formalizer a clean
finite target over THM-369 (cf. HYP-1986 source-gap route).  Ekhad-Zeilberger
(arXiv:1112.6207) give the g.f./recurrence for related shift-register sequences.

## Evidence / scripts

- `04-computation/lrc_arc_menu_geometry_s521.py` (+ `.out`): menu, labelled
  counts, full class lists (H/score), L-sweep.
- `04-computation/lrc_arc_menu_confirm_s521.py` (+ `.out`): exact
  refinement-canonical recount (validated vs brute permutation canon for
  m<=8, and brute-cross-checked at m=9), Fibonacci/A164142 guesses refuted,
  A000016 closed form confirmed, L-independence at L=0.505 vs 0.995.

## See

THM-387, THM-384, THM-381, HYP-1987, HYP-1981, HYP-1951, THM-369;
OEIS A000016.
