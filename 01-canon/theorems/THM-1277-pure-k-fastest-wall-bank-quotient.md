---
id: THM-1277
title: PURE-K FASTEST-WALL BANK QUOTIENT -- complete fastest teeth are whole floods or two distinct-owner chronological seams
status: PROVED (exact pure-K wall-bank extraction; sharp alternating-wall count; first-wall private fastest tooth; complete-tooth selected/unselected dichotomy; two distinct boundary-owner seams for every selected tooth; arbitrary prefix/suffix whole-tooth floods layered pointwise with the full chronological seam invoice; exact lcm, four-owner, functional, harmonic, and localized regular-run forms; c=140 two-wall/no-complete-tooth guardrail; optimization-safe exact referee; sorry-free Lean arithmetic consumer).  The complementary near-top tax can decay under primitive scaling and does not close LRC(14)
source: codex-2026-07-19 pure-K wall-bank selection audit
depends_on: [THM-1198, THM-1253, THM-1266, THM-1273, THM-1275]
related: [THM-1156, THM-1241, THM-1260, THM-1272]
script: 04-computation/lrc14_pure_k_fastest_wall_bank_quotient_thm1277.py
output: 05-knowledge/results/lrc14_pure_k_fastest_wall_bank_quotient_thm1277.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCPureKFastestWallBankQuotient.lean
script_sha256: HASH_PENDING
output_sha256: HASH_PENDING
formalization_sha256: HASH_PENDING
---

# THM-1277 -- pure-`K` fastest-wall bank quotient

## 1. Statement

Use THM-1273's notation.  In particular, a complete `c`-safe gap `G` is
covered by the strict danger combs of

```text
c<d1<d2<d3<d4<d5<h,                    h=d6,          (1)
```

`S` is the centered complete closed `d1`-safe component, and after reflection

```text
K=S intersect G,             E=S minus G             (2)
```

occur in that order.  Normalize `S` to `[0,1]` and write `b` for the `K/E`
interface.  In THM-1273's small-protrusion branch one may choose a regular
point

```text
x in int(K),       b-x>1/12,                          (3)
```

which is closed-safe for `d2,d3,d4,d5` and strictly dangerous for `h`.

Let `H_0` be the open `h`-tooth containing `x`.  Moving from `x` toward `b`,
let `delta` be one if the right wall of `H_0` occurs before `b`, and zero
otherwise.  Let `P` be the number of later complete `h`-teeth whose **two**
walls both lie in `(x,b)`.  If `W_K` is the number of `h` walls in `(x,b)`,
then the alternating wall order gives exactly

```text
delta=1_(W_K>0),
P=floor((W_K-1)/2)                    when W_K>0.     (4)
```

More importantly, (3) gives the phase-free lower bank

```text
h>=2d1  => delta=1,
P>=floor((h-2d1)/(14d1)).                             (5)
```

Choose any deletion-minimal tooth subcover of `G`.  Among the `P` complete
fastest teeth, let

```text
F=# unselected teeth,       A=# selected teeth,      F+A=P.  (6)
```

Then:

* every unselected complete fastest tooth is wholly covered by selected low
  teeth and supplies one whole-tooth **flood**;
* if `delta=1`, `H_0` is selected and its right wall supplies one selected
  `h`--low chronological seam;
* every selected complete fastest tooth supplies two selected chronological
  seams, one at each wall, and their two lower owners are distinct members of
  `{d2,d3,d4,d5}`.

All these forced seam occurrences are distinct.  The `F` whole teeth are
pairwise disjoint and can be added pointwise to **all** selected seams, even
when they occur before the first or after the last selected fastest address.
If `j_0` is the first-wall owner when `delta=1`, and `(j_r^-,j_r^+)` are the
two boundary owners of selected complete tooth `r`, the exact weighted bank
invoice is

```text
F6 >= cF/(8h)
      +(c/16)[delta/lcm(h,j_0)
              +sum_(selected r){1/lcm(h,j_r^-)+1/lcm(h,j_r^+)}].  (7)
```

Here `F6=sum_i Pbar(6d_i/(7c))-1`.  In particular, distinctness of the two
owners gives the phase-free consequences

```text
F6 >= cF/(8h)+c delta/(16h d5)
                 +cA/(16h)(1/d4+1/d5)               (8)
   >= c delta/(16h d5)
                 +cP/(16h)(1/d4+1/d5).              (9)
```

The harmonic companion is

```text
sum_i 1/d_i-1/c
 >=7F/(6h)+(7/12)[delta/lcm(h,j_0)
              +sum_(selected r){1/lcm(h,j_r^-)+1/lcm(h,j_r^+)}]  (10)
 >=7 delta/(12h d5)
              +7P/(12h)(1/d4+1/d5).                 (11)
```

Thus THM-1273's wall bank is not merely a bank of possibly duplicated point
events.  After quotienting the two walls of every complete fastest tooth, it
is a bank of disjoint proof obligations: a whole flood or a two-sided seam
fork with distinct lower-owner colours.

## 2. Exact alternating-wall count

In normalized `S` coordinates an `h` tooth has length

```text
a=d1/(6h),                                           (12)
```

and the following `h`-safe gap has length `6a=d1/h`.  Since `x` is in the
interior of `H_0`, the first wall in the positive direction is a right wall.
The successive wall-free cell lengths are therefore

```text
at most a, then 6a,a,6a,a,... .                      (13)
```

Suppose there are exactly `P` complete fastest teeth after `H_0` and before
`b`.  The longest possible interval starts arbitrarily near the left wall of
`H_0`, crosses `P` complete tooth/gap periods, crosses one further safe gap,
and stops arbitrarily near the far wall of the next tooth.  Consequently

```text
b-x <=(7P+8)a.                                       (14)
```

Equations (3), (12), and (14) give

```text
h<(14P+16)d1.                                        (15)
```

If `h>=2d1`, then `a<=1/12<b-x`, so the first wall occurs and `delta=1`.
If an integer `Q` satisfies `(14Q+2)d1<=h`, (15) forbids `Q>=P+1`; hence
`Q<=P`.  Taking `Q=floor((h-2d1)/(14d1))` proves (5).  The constants are
sharp for phase-free information: the two unfinished endpoint tooth pieces
in (14) really can approach total length `2a`.

This improves the coarse `W_K>=floor(h/(12d1))` obtained by bounding every
wall-free cell by the longest cell.  It retains the alternating `1:6` wall
geometry instead of erasing it.

## 3. Why every wall in this bank is crossed

Every wall in `(x,b)` lies in `int(K)`, hence in `int(G)`.  The owner `d1`
is closed-safe throughout `S`; `h` is exactly safe at its own wall.  Since
the strict fast combs cover `G`, at least one of `d2,d3,d4,d5` is strictly
dangerous there.  The same remains true after passing to a deletion-minimal
individual-tooth subcover: some **selected** lower tooth contains the wall.

The crosser label has an exact finite residue form.  Write a fastest wall as

```text
z=a/(14h),                 a=14n+-1.                 (16a)
```

For an integer `u`, let `|u|_m` denote the least absolute residue modulo
`m`.  Then a lower owner `j` crosses this wall exactly when

```text
||jz||<1/14
 iff ||ja/(14h)||<1/14
 iff |ja|_(14h)<h.                                  (16b)
```

Thus every selected/flood cell can retain the wall numerator `a` and the
four-bit crosser mask using only multiplication and centered reduction modulo
`14h`.  This is the exact endpoint digit which a runner-order tournament
forgets; the referee checks (16b) on both wall signs and negative as well as
positive addresses.

This is the useful strengthening of THM-1273.  Its longer `67/540` needle can
enter `E` and encounter bare walls.  The shorter left-hand needle `(x,b)` is
purely inside `K`, so its wall bank has no bare branch at all.

The point `x` itself is covered only by `h`: it is safe for `d1,...,d5` and
strictly inside `H_0`.  Hence every minimal subcover selects `H_0`.  If its
right wall lies before `b`, the selected low tooth crossing that wall overlaps
`H_0` positively.  THM-1253's minimal-chain separation makes this pair a
literal chronological seam.  This is the `delta` seam in (7).

## 4. Pair the two walls of a complete fastest tooth

Fix a complete fastest tooth `H` counted by `P`.

If `H` is not selected, the selected subcover still covers every point of
`H subset int(G)`.  Since `d1` is safe on `S`, selected teeth of the four
middle owners cover all of `H`.  The full `h` comb is simultaneously active,
so the whole interval `H` lies at multiplicity at least two.  It has physical
length `1/(7h)`.

Now suppose `H` is selected.  At each of its walls choose a selected lower
tooth which crosses that wall.  The two choices cannot be the same tooth: an
open interval containing both walls contains the closure of `H`, making `H`
deletable.  They cannot even be distinct teeth of the same owner `j`.  Two
`j` teeth are separated by a safe gap of length `6/(7j)`, whereas

```text
|H|=1/(7h)<1/(7j)<6/(7j).                            (16c)
```

Thus the two lower owners are distinct.  Each lower tooth overlaps `H`
positively.  In a deletion-minimal interval cover, overlapping selected
intervals are consecutive: a third selected interval between them in endpoint
order would be contained in their union.  The two intersections are therefore
raw THM-1253 seams.  They are distinct and, by minimality, disjoint; otherwise
the two lower teeth would cover all of `H` and make it redundant.

This proves the promised quotient:

```text
complete h tooth
       |-- unselected --> whole h flood
       `-- selected   --> ordered pair (j_left,j_right), j_left!=j_right,
                          and two chronological seams.              (17)
```

## 5. The prefix/suffix flood extension of the layered invoice

THM-1275 added skipped fastest teeth between successive selected fastest
addresses to the complete chronological seam family.  Here an unselected
complete tooth in `(x,b)` may lie before the first or after the last selected
fastest address, so it need not occur in THM-1275's variable `X`.  The same
pointwise proof nevertheless applies.

Let `C(t)` be the number of active full fast combs.  At a point of one of the
unselected teeth, the selected subcover supplies a low tooth, so `C(t)>=2`.
If the point also belongs to a raw seam, that seam uses two selected low teeth:
no selected `h` tooth can overlap the unselected tooth.  The full unselected
`h` tooth is then a third active owner and `C(t)>=3`.  Since both the flood
teeth and the seam intervals are internally disjoint,

```text
sum_(unselected H) 1_H(t)+sum_(all raw seams W)1_W(t) <=C(t)-1.   (18)
```

Normalization by `7c/6` and the density floor `3/4` charges a whole fastest
tooth by `c/(8h)` and a seam by `(7c/8)|W|`.  The endpoint numerator law gives

```text
|W|>=1/[14 lcm(h,j)],                                (19)
```

which proves (7).  THM-1253's singleton-discrepancy conversion gives (10).

For a selected complete tooth, its two distinct owners are at most `d5`,
and the smaller one is at most `d4`.  Therefore

```text
1/lcm(h,j^-)+1/lcm(h,j^+)
 >=1/(hj^-)+1/(hj^+)
 >=(1/h)(1/d4+1/d5).                                 (20)
```

Also `1/lcm(h,j_0)>=1/(h d5)`.  Finally a whole-tooth flood is stronger than
the right side of (20), because `1/d4+1/d5<=2`.  This proves (8)--(11).

## 6. Localized star completion

The quotient retains more than scalar mass.  Put

```text
e=#{i<=5:h>6d_i}.                                    (21)
```

The block `H_0,H_1,...,H_P` consists of consecutive fastest addresses;
`H_0` is selected and `F` of the remaining teeth are not.  Among adjacent
selected pairs, let `T` count transitions having at least two internal low
teeth.  The other adjacent selected transitions are regular rungs.

The `F` missing teeth split the selected teeth into at most `F+1` address
blocks.  The `T` turns split those into at most `F+T+1` regular vertex runs.
THM-1266 bounds each such run by `e+1` selected fastest teeth.  Since there
are `A+1=P-F+1` selected teeth including `H_0`,

```text
A+1 <=(e+1)(F+T+1),                                  (22)
P <=(e+2)F+(e+1)T+e,                                 (23)
F+T >=ceil_+((P-e)/(e+2)).                            (24)
```

In particular `P>e` forces a located flood or turn.  Under the top-dominated
condition `h>=(7/2)d5`, every such turn pays at least `c/(8h)` by THM-1275,
so independently of (9),

```text
F6 >= c/(8h) ceil_+((P-e)/(e+2)).                    (25)
```

This is a genuine local star-completion statement.  It does not improve
THM-1272's `d6/c<798` cut: THM-1275's global private count supplies at least
as many fastest occurrences earlier and with the sharper denominator `e+1`.

## 7. Why the complementary branch remains open

The intended target was the complementary near-top branch

```text
h<(7/2)d5.                                           (26)
```

Here `e<=4`, and (5), (24) force many located floods or turns once `h/d1`
is large.  But a turn packet has only the exact positive reciprocal invoice

```text
sum_q 1/j_q-6/h>0;                                   (27)
```

condition (26) does not lower it by `1/h`.  The universally forced credit is
therefore (7), or its gcd-forgetting form (9).  The latter can decay: for
primitive near-coprime `h,j`, `1/lcm(h,j)` is of order `1/(hj)`.  With only
`P=O(h/d1)` cells, (9) is of order `c/(d1 d5)`, which tends to zero along
large primitive carriers with fixed ratios.

This is the honest obstruction.  THM-1277 improves the structural state and
creates an exact located lcm bank, but it does not produce a ratio-only
projective terminal in (26).  A closing consumer must retain gcd/endpoint
digits, correlate the distinct boundary-owner pairs across several teeth,
or obtain a non-arithmetic lower bound for the multi-low turn bracket.

## 8. The exact `c=140` guardrail

In THM-1266/1273's sharp cell, take the displayed control point `x`, the
interface `b`, and the displayed exterior point `y`:

```text
c=140, d1=254, h=1805,
x=7476011/12938240, b=1133/1960, y=7425603/12837160. (28)
```

Exactly one fastest wall lies in the pure-`K` interval `(x,b)`:

```text
14603/25270 = right wall of 1805@1043.                (29)
```

Thus `delta=1` and `P=0`.  On the longer mixed `K/E` needle `(x,y)`, the
wall word is

```text
14603/25270 = right wall of 1805@1043,
2923/5054   = left  wall of 1805@1044,
14617/25270 = right wall of 1805@1044.                (30)
```

The first two are crossed by the single tooth `256@148`.  They bound the
complete **fastest-safe gap** between different fastest teeth and straddle
the `K/E` interface; they are not the two walls of one fastest tooth.  This
shows why the mixed wall count cannot be inserted unchanged into the pure-`K`
invoice.  The correct quotient retains one pure-`K` boundary seam and no
complete tooth.  The ratio corollary agrees:

```text
floor((1805-2*254)/(14*254))=0.                       (31)
```

The row is not a six-cover, and its concrete `x` was not chosen to realize the
deeper `b-x>1/12` option in the abstract proof, so it is only a geometry
control.  It shows why both the alternating wall types and the `K/E` location,
rather than an untyped matching of wall events, must survive the quotient.

## 9. Tournament, Fano, Kakeya, and carrier audit

The normalized centered component is the Kakeya needle.  Individual walls
are too fine a vertex set: the first two mixed walls in (30) would be falsely
counted as one pure-`K` complete-tooth obligation.  Pair walls by both the
alternating fastest-cell type and their `K/E` location, and use the `P`
complete pure-`K` fastest teeth as vertices.  The binary switch is

```text
unselected/flood  versus  selected/two-seam fork.    (32)
```

Chronological order is the tie Hamiltonian path.  On `P` vertices its
tournament is transitive, with score histogram `(0,1,...,P-1)`, no directed
cycles, `P` singleton SCCs, no edge flips, and one Hamiltonian path.

The selected vertices carry an additional ordered boundary-owner edge

```text
j_left -> j_right,             j_left!=j_right,       (33)
```

on the four vertices `{d2,d3,d4,d5}`.  This is a directed multigraph, not a
tournament: pairs may be absent, repeated, or occur in both orientations.
There are six unordered and twelve ordered colours, so seven selected teeth
force a repeated unordered colour and thirteen force a repeated oriented
colour.  Those repetitions are a concrete placed input for the still-open
Fano/`chi_7` probe, but no `chi_7` restriction follows locally; THM-1260's
surjectivity guardrail remains intact.

We challenged runners, gaps, individual walls, wall pairs, complete fastest
teeth, lower tooth occurrences, selected seams, residues, Fano lines, and
proof obligations as vertices.  The faithful carrier is

```text
(oriented S and K; x and b; alternating h-wall word; complete h-tooth cells;
 selected/flood bit; ordered distinct lower-owner pair; exact seam intervals
 and lcm endpoint digits).                            (34)
```

It preserves the cover predicate, paired-wall topology, minimal-subcover
status, spatial disjointness, and the layered invoice.  The runner quotient
destroys all of them; the untyped wall quotient destroys which two walls bound
a tooth rather than a safe gap; the owner-pair quotient destroys multiplicity
and chronology.

## 10. Verification and scope

The dependency-free exact referee has no Python `assert` nodes.  It exhausts
the alternating wall bank over rational phases, checks the sharp bound (15)
and all integer thresholds in (5), enumerates the selected/flood/turn block
capacity (22)--(24), verifies the pointwise layered truth table and every
normalization coefficient, rechecks the lcm endpoint quantum and the
distinct-owner improvement, checks the exact wall-residue crosser identity,
and reconstructs the one-pure-wall/three-mixed-wall `c=140` guardrail.  Normal
and optimized outputs are byte-identical.

The sorry-free Lean module kernel-checks the wall-length consumer, threshold
rounding interface, layered multiplicity rule, functional and harmonic bank
coefficients, distinct-owner improvement, localized run capacity, and all
displayed `c=140` wall identities.  Selection of regular points, minimal
subcover extraction, interval-chain adjacency, and the geometric wall pairing
are the explicit paper providers.

THM-1277 does not prove that a complementary-branch turn has scale-free mass,
close `d6<(7/2)d5`, prove uniform six-comb noncoverage, empty the sporadic
branch, or prove LRC(14).  It turns THM-1273's wall bank into the strongest
currently justified disjoint selected/flood invoice and identifies the exact
remaining loss.  Frozen hashes are recorded above after verification.  ∎
