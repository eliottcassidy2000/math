# LRC body `H`: upper-median two-star matching scout

Date: 2026-08-02
Status: **FINITE-EXACT SCOUT + SAMPLED INDEPENDENT ARCHITECTURE AUDIT.**
This is not a theorem, not a proved dependency, and not a closure of a
THM-2941 reflected cone or of `LRC(14)`.

## Inheritance and target

The hostile body left by the corrected THM-2941/THM-3135 reflected analysis is

```text
H=(1,2,3,4,6,12),                 L=lcm(H,14)=168.            (1)
```

THM-3135 proves that every standard single-pair uniform envelope for `H`
orients to a DAG.  The least-used surviving operation is therefore genuinely
multi-pair: keep one safe upper-median cell and allow a small graph of overlap
channels to answer each hostile level assignment.

The scout fixes safe cell `j=90`, two pivot labels `0,1`, and the nine edges

```text
01,02,03,04,05,12,13,14,15.                                 (2)
```

Thus the certificate graph is the union of two stars.  For six **distinct**
integer levels injected into `[m,m+D]`, it tests whether at least one edge in
`(2)` has exact intersection gain larger than the total inherited singleton
debt

```text
sum_(i=0)^5 H_i/[7(q_i L-H_i)].                              (3)
```

## Finite-exact universes

The exact certificate is positive in both declared universes:

```text
6<=D<=50,       m=D;                                         (4)
6<=D<=15,       D<=m<=3D.                                   (5)
```

These are respectively `45` boundary cells and `220` two-parameter cells.
The quantifier inside each cell is over every injection of the six labels into
`[m,m+D]`; the implementation does not sample assignments.

The weakest certified scaled floor in `(4)` occurs at `D=m=8`:

```text
pivots=(8,12),
threshold floor=3248767010812/10^15,
relaxed debt ceiling=2623186501005/10^15,
certified gap=625580509807/10^15>0.                           (6)
```

Before directed rounding, the weakest leaf gain and relaxed debt in this
certificate are

```text
2928/901265,
147749412194406/56324402454058375,                           (7)
```

whose exact difference is

```text
6307145263970526/10082068039276449125>0.                     (8)
```

The floor/ceiling calculation at scale `10^15` is one-sided: gains are rounded
down and debts up, so positivity of the integer certificate implies the exact
rational inequality.

## Why the exhaustive quantifier is tractable

After fixing pivot levels `q_0,q_1`, assign to each leaf `i=2,3,4,5` the
larger of its two pivot-edge gains.  For a trial gain threshold, the hostile
problem is a maximum-weight injective matching from four leaves to the
remaining levels, with allowed edges defined by being below that threshold.

The singleton debt is strictly decreasing with the level.  Consequently an
optimal leaf uses one of its first four allowed levels: if it used a later
one, the four earlier and better levels could not all be occupied by the other
three leaves.  This turns the apparent permutation census into a four-left
matching problem.  The independent audit replaces that truncation by a
level-by-level subset dynamic program and agrees on the six hostile controls

```text
(D,m)=(6,6),(8,8),(10,10),(12,12),(15,45),(20,20).           (9)
```

The two scripts load the same canonical THM-2941 interval geometry by separate
routes; independence is in the subset-DP matching architecture, not in the
underlying arc engine.

## Mechanism suggested by the scout

The positive signal is not a universal winning edge.  The winning edge moves
among the two pivots and four leaves as the interval phases change.  The useful
structural object is instead a low-overlap **channel graph**: pivot blindness
can occur only at a short list of ratios, while simultaneous blindness of all
leaf channels appears incompatible with the injective matching constraint.

This points to an analytic Hall-type statement:

> For every threshold below total debt, the four leaf neighborhoods in the
> two-star channel graph cannot admit a matching that keeps every pivot/leaf
> gain below threshold.

Proving that statement would explain the computation and could survive beyond
the scanned phase boxes.  A raw monotonicity extrapolation will not: exact arc
overlap oscillates with phase and gcd data.

## Scope and stopping boundary

The scout does **not** cover repeated levels, `D>50` on the boundary, or most
cells with `D>15` and `m>D`.  It also does not replay every inherited filter
needed to promote an audited THM-2941 cone, and it does not repair the global
DAG obstruction of THM-3135.  The next decisive task is to classify the
two-star low-overlap neighborhoods symbolically and prove the proposed Hall
deficit, then replay the resulting statement with the complete inherited
filters and canonical repeated-level discharge.

## Reproduction

Run

```text
python3 04-computation/lrc_H_upper_median_two_star_matching_scout.py --max-D 50
python3 -O 04-computation/lrc_H_upper_median_two_star_matching_scout.py --max-D 50
python3 04-computation/lrc_H_upper_median_two_star_matching_scout.py --max-D 15 --all-m --m-multiples 3
python3 -O 04-computation/lrc_H_upper_median_two_star_matching_scout.py --max-D 15 --all-m --m-multiples 3
python3 04-computation/lrc_H_upper_median_two_star_matching_audit.py
python3 -O 04-computation/lrc_H_upper_median_two_star_matching_audit.py
```

Normal, optimized, and stored transcripts match byte for byte.  The script
and output hashes are

```text
primary script   928083a7e81da8491e54468d121cb9fb404349255d5dd32c87b2d52700a26218
audit script     6cff95e576184558018458545cd7bb46d7472bb6f5275dfd850389b76de6d3be
boundary output  1bd7962badf57fd5d1ec7eb344841fbd22625ecc54ab20a637d17d0e8a721d0c
box output       fb815e95914cc64a5ecded206a2380bfc91877298e7df628b3542e2ed1487d28
audit output     348d063a97fe24de2e3bc1570f1370947345a951d9045b88ffc08f080b8dede6
```
