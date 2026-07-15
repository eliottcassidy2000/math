---
id: THM-782
title: Uniform phase-cell return packets inside every ten-speed two-sheet core
status: PROVED (elementary corollary of the quantitative phase-pigeonhole theorem and settled LRC(<=13)); quantitatively strengthened and locally delimited by THM-789
source: codex-2026-07-14-S9
depends_on:
  - THM-780
  - LRC(<=13)
related:
  - THM-772
  - THM-774
  - THM-776
  - THM-789
  - HYP-6820
---

# THM-782 — A uniform phase-cell return packet in every two-sheet core

> **Strengthening (THM-789).** Passing from one anchored difference `A-s_0` to
> the symmetric difference set `A-A` doubles the measure floor to
> `2*72^(-10)` and improves the component-width floor to
> `72^(-10)/(5 max(U))`. Tightness also forces the pointwise thickness tax (3)
> and the erosion `E_U subset H minus R_U`. An exact admissible triple shows
> the full natural return set can remain inside the diamond at one deep time,
> so the residual is global selection among deep components, not further local
> refinement of an arbitrary anchor.

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

Finally, THM-774 proves the lossless two-sheet equivalence (restated and used
in the bounded theorem THM-776)

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
not just a constants problem.  The exact core-19 closure audit gives

```text
U={1,2,3,5,7,8,9,10,11,12},       (x,y)=(1,9),
|G_U|=41/858 < 8/117=|H_(1,9)|.                            (16)
```

Thus even comparison of the two exact scalar measures does not refute
containment for this admissible measure profile; the incidence geometry is
essential.

The new theorem-facing object is the containment in (9): a translate of one
heavy joint-phase atom, anchored at zero after subtraction, whose every point
is a simultaneous `1/72` return for all ten core frequencies.  The most direct
uniform residual can now be asked as a **structured noncontainment theorem**:

> For every admissible `(U,x,y)`, can one choose a `1/11`-deep time `t_0`, a
> phase cell `Q` with `|phi^(-1)(Q)|>=72^(-10)`, and an anchor
> `s_0 in phi^(-1)(Q)` so that
> `t_0+(phi^(-1)(Q)-s_0)` is not contained in `H_(x,y)`?

This is strictly richer than comparing `|G_U|` and `|H_(x,y)|`, but it is not
yet a proof that such containment is impossible.

## Assumption challenge and Tournament Analysis

The useful vertices in the proof are the `72^10` joint phase cells, not the
ten runners, their gaps, or their pairwise arcs.  Subtracting an anchor from a
heavy cell preserves the simultaneous-return predicate (3) and, after adding
the deep time, the LRC-safe predicate (4).  It destroys the cell's absolute
phase address and most endpoint incidence data; those must be restored when
testing containment in the folded diamond.

This proof exhibits no natural useful antisymmetric switch on these cells.
Same-cell membership is an equivalence relation, and anchored subtraction is
a group operation.  An unlabelled runner or cell tournament, without the
cell/anchor incidence labels, does not encode (3).  Runner tournaments remain
telemetry here, while the labelled phase-cell/return-packet incidence is the
proof carrier.
