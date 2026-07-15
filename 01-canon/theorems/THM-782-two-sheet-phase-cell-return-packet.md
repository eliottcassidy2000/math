---
id: THM-782
title: Uniform phase-cell return packets inside every ten-speed two-sheet core
status: PROVED (elementary corollary of the quantitative phase-pigeonhole theorem and settled LRC(<=13))
source: codex-2026-07-14-S9
depends_on:
  - THM-780
  - LRC(<=13)
related:
  - THM-772
  - THM-774
  - THM-776
  - HYP-6820
---

# THM-782 — A uniform phase-cell return packet in every two-sheet core

## Statement

Let `U` be any set of ten distinct positive integers, put `H=max(U)`, and
write

```text
G_U={tau in R/Z : ||u tau||>1/13 for every u in U}.       (1)
```

There are `t_0,s_0 in R/Z` and a measurable set `A_0 subset R/Z` such that,
for

```text
B=A_0-s_0,
```

all of the following hold:

```text
0 in B,                         |B|=|A_0|>=72^(-10),      (2)
||u b||<1/72                    (u in U, b in B),         (3)
||u(t_0+b)||>1/13+1/10296       (u in U, b in B).         (4)
```

In particular,

```text
|G_U|>=72^(-10).                                           (5)
```

Moreover, `G_U` has at most `2 sum_(u in U) u` connected components.  Hence
it has a component `J` satisfying

```text
|J| >= 72^(-10)/(2 sum_(u in U)u)
    >= 72^(-10)/(20H).                                    (6)
```

Now suppose that a two-sheet packet

```text
S=2U union {x,y},             x<y odd,                    (7)
```

is tight at `1/13`.  With the folded diamond

```text
H_(x,y)={tau:
  ||(x+y)tau/2||+||(y-x)tau/2|| >= 11/13},                (8)
```

THM-774's exact criterion forces

```text
t_0+B subset G_U subset H_(x,y),                           (9)
```

and the component `J` in (6) is also contained in `H_(x,y)`.  Thus a
hypothetical tight packet must contain an anchored simultaneous-return packet,
not merely a scalar amount of safe measure.

## Proof

Settled `LRC(<=13)`, applied to the ten speeds in `U`, gives a time `t_0`
with

```text
||u t_0||>=1/11             for every u in U.             (10)
```

Apply the finite phase-cell proof of THM-780 with `d=10`, `beta=1/11`, and
`alpha=1/13`.  Since

```text
1/(beta-alpha)=1/(1/11-1/13)=143/2,
```

the phase partition has

```text
M=ceil(143/2)=72                                               (11)
```

intervals in each coordinate and `72^10` joint cells.  For the orbit map

```text
phi:R/Z -> (R/Z)^10,       phi(s)=(us)_(u in U),            (12)
```

one joint-cell preimage `A_0` has measure at least `72^(-10)`.  Choose
`s_0 in A_0`.  Any two points of `A_0` have all their phase coordinates in
the same half-open interval of length `1/72`; consequently `B=A_0-s_0`
contains zero, has the same measure as `A_0`, and satisfies (3).

The circle triangle inequality and (10) now give, for every `u in U` and
`b in B`,

```text
||u(t_0+b)|| >= ||u t_0||-||u b||
              > 1/11-1/72
              = 1/13+1/10296.                              (13)
```

This proves (2)--(5).

Every boundary point of `G_U` lies on a wall

```text
||u tau||=1/13.
```

For a fixed `u` this wall has exactly `2u` points on the circle.  Therefore
the open set `G_U` has at most `2 sum u` components.  Combining that count
with (5), and then using `sum u<=10H`, proves (6).

Finally, THM-774/776 prove the lossless two-sheet equivalence

```text
M(2U union {x,y})<=1/13  iff  G_U subset H_(x,y).          (14)
```

For a tight packet (14), (4), and (6) give (9) and the asserted component
containment.  This completes the proof.

## What this changes — and what it does not

THM-780 supplies the uniform measure floor that the earlier two-sheet
discussion lacked.  Its value is deliberately crude:

```text
72^(-10) << 8/117,                                          (15)
```

where `8/117` is THM-774's sharp upper bound for the measure of a folded
diamond.  Thus (5) cannot by itself contradict `G_U subset H_(x,y)`.  This is
not just a constants problem: the exact core-19 closure audit contains
primitive divisor-complete cores with `|G_U|<=8/117`.  A proof based only on
the two scalar measures is therefore structurally insufficient on those
cores.

The new theorem-facing object is the containment in (9): a translate of one
heavy joint-phase atom, anchored at zero after subtraction, whose every point
is a simultaneous `1/72` return for all ten core frequencies.  The most direct
uniform residual can now be asked as a **structured noncontainment theorem**:

> Can an admissible folded diamond contain `t_0+B` for a primitive,
> divisor-complete ten-core satisfying THM-772's speed bounds?

This is strictly richer than comparing `|G_U|` and `|H_(x,y)|`, but it is not
yet a proof that such containment is impossible.

## Assumption challenge and Tournament Analysis

The useful vertices in the proof are the `72^10` joint phase cells, not the
ten runners, their gaps, or their pairwise arcs.  Subtracting an anchor from a
heavy cell preserves the simultaneous-return predicate (3) and, after adding
the deep time, the LRC-safe predicate (4).  It destroys the cell's absolute
phase address and most endpoint incidence data; those must be restored when
testing containment in the folded diamond.

There is no clean antisymmetric switch on these cells.  Same-cell membership
is an equivalence relation, and anchored subtraction is a group operation.
A forced tournament would discard precisely the common-cell information used
in (3).  The honest tournament fingerprint is therefore “no predicate-
preserving orientation”; runner tournaments remain telemetry, while the
phase-cell/return-packet incidence is the proof carrier.
