# LRC(14): the canonical uniform-unit atlas is not a `D`-invariant state space

**Status: VERIFIED-EXACT SCOUT / THM-2680 INPUT.  Not a theorem dependency.**

The closest positive mechanism is THM-2635: after choosing its canonical rail
in each base cell, the same four fixed-half unit labels work everywhere,

```text
U^0={9},       U^1={3,8,10}.
```

The closest hostile is THM-2670: its formal clock matrices do not compose
physical witnesses, while the natural digit map `D(x)={13x}` reverses the
clock bookkeeping and imposes only the numerical equality `j(Dx)=h(x)`.
The least-used sidecar is the literal successor half/carry geometry retained
before integration in THM-2630/2670.

## The exact connection

The source is a current canonical unit label `(epsilon,h)`.  The target is a
possible child label `(epsilon',h')`.  The numerical map forgets the source
rail, base cell, component, and physical point and retains only

```text
j'=h,
r'=2j'+epsilon'+kappa' mod 13,
h'=-r'-1 mod 13.
```

For the four source labels and four child bit pairs, this gives sixteen
rows.  Only three rows land back in the canonical uniform-unit atlas:

```text
(source epsilon,h; child epsilon',kappa',h')
(1,8;  0,0,9),
(1,8;  1,0,8),
(1,10; 1,1,3).
```

But a literal successor half can meet the `h'` digit only when

```text
kappa'=0 => h'<=6,
kappa'=1 => h'>=6.
```

Each of the three rows violates this condition.  Hence every possible
state-independent uniform-unit child is structurally empty.  The exact
carrier reconstruction checks all `15,552` rail-labelled candidates,
retains the safe/danger/guard-free cospan, and finds eight positive child rows,
all outside the canonical uniform atlas.  Two counts distinguish why the
quotient matters: `367` and `345` positive rail occurrences project to only
`347` and `325` distinct `(d,s,ell)` cells.

## What this destroys, and what survives

Destroyed: the hope that the four globally uniform THM-2635 labels form a
closed finite state space for iterating `D` without extra memory.

Survives: positive nonunit child atoms, cell-dependent unit sets `U_c^epsilon`,
and every physical two-edge fibre product.  No map on rail, source, base-cell,
or connected-component labels has been constructed, so the computation does
not say that a positive current atom lacks a positive child atom at `Dx`.

The next decisive object should therefore retain the missing coordinate rather
than enlarge the same quotient.  The smallest honest test is a physical
two-edge set

```text
E_(b,a;s0,j,h) intersect D^(-1) E_(c,b;s1,h,k),
```

with both source labels and the current/child rail components recorded.  If
that object is positive, its induced cell-dependent unit label is precisely
the sidecar needed to repair closure.  If it is empty, the first failed factor
will identify whether the obstruction is the base carrier, delayed word,
source transport, or literal half geometry.

Reproduce with

```bash
python3 04-computation/lrc14_dilation_canonical_unit_atlas_no_go_20260728.py
python3 -O 04-computation/lrc14_dilation_canonical_unit_atlas_no_go_20260728.py
```

Both executions must byte-match
`05-knowledge/results/lrc14_dilation_canonical_unit_atlas_no_go_20260728.out`.
