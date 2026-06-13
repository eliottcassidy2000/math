# Fourteen-Runner Scalar-Gauge Quotient

Session: `codex-2026-05-31-S367`

The scalar-ramp obstruction from S364 was the right dead end to stare at.  In
the `n=14`, `k=13` micro-staircase cell system, every affine ramp

```text
v_i = m i mod 14
```

blocks every `(s, alpha-cell)` candidate.  This is not a random obstruction and
not evidence against the route.  It is a gauge symmetry.

Indeed, replacing `v_i` by `v_i + m i` changes

```text
s v_i + floor(14 {i alpha})
```

to

```text
s v_i + floor(14 {i alpha}) + s m i       mod 14.
```

Modulo `14`, this is the same floor-vector one gets after shifting the cell
parameter by `s m / 14`.  Since the full micro-staircase arrangement contains
all alpha cells, scalar addition only reindexes candidates.  Thus the affine
line should be quotiented out, not fought coordinate-by-coordinate.

The normalized representative is

```text
v -> v - v_1 (1,2,...,13),  so v_1 = 0.
```

All scalar ramps collapse to the single zero vector.  The proof target becomes:

```text
Every nonzero normalized vector v in (Z/14Z)^13 has some shift s and some
micro-staircase cell alpha such that

    s v_i + floor(14 {i alpha}) notin {0,13} mod 14

for every i=1,...,13.
```

This is the quotient version of the S363 micro-staircase lemma.

## Computation

`04-computation/lonely_runner_k13_scalar_gauge_s367.py` rebuilds the exact
cell system:

```text
n=14, k=13
patterns=812
candidates=11368
```

The script verifies, using raw unnormalized scores, that random gauge shifts
preserve coverage:

```text
Gauge invariance random check: trials=200, failures=0
```

It also confirms that every scalar ramp is a full blocker and that all of them
normalize to zero.

After quotienting, local search again finds no full non-scalar blocker.  The
best local optimum is the very simple half-turn perturbation

```text
v = (0,0,0,0,0,7,0,0,0,0,0,0,0),
```

which covers `11312/11368` candidates and misses `56`.

## Exact Small Families

The new useful evidence is not just local search.

Single-coordinate half-turn missed counts are:

```text
coord:  2  3  4  5  6  7  8  9  10 11 12 13
miss: 196 224 126 154 56 294 280 336 168 350 84 168
```

Exact normalized support scans found:

```text
support 1: scanned 156,    best missed 56  at support (6,)
support 2: scanned 11154,  best missed 112 at support (8,12)
support 3: scanned 483340, best missed 126 at support (4,6,8)
```

The complete normalized `2`-torsion cube is small enough to exhaust:

```text
scanned=4095
best_missed=56
best_count=1
```

So inside the binary half-turn quotient, the coordinate-6 bump is the unique
best near-blocker and still leaves many open cells.

## The 56 Misses

For the extremal vector, all misses occur at odd shifts:

```text
shift_hist=[(1,8),(3,8),(5,8),(7,8),(9,8),(11,8),(13,8)]
gcd_shift_hist=[(1,48),(7,8)]
min_margin_hist=[(1,56)]
```

There are eight alpha patterns, each repeated for all seven odd shifts.  Four
widths occur, each with count `14`:

```text
1/728, 1/882, 1/1176, 1/1386.
```

The first missed interval is

```text
alpha in [29/182, 9/56)
bins=(2,4,6,8,11,13,1,3,6,8,10,12,1).
```

At any odd shift, the coordinate-6 half-turn changes the sixth residue from
`13` to `6`, turning a boundary-blocking cell into a genuine safe cell.  This
is the cleanest small witness family found so far.

## Proof Shape

I now think the `k=13` proof should split into three finite lemmas.

1. **Gauge lemma.**  Prove for all `n` that scalar addition `v_i -> v_i+m i`
   is candidate-reindexing, hence it preserves full-cell coverage.

2. **Zero-class lemma.**  For `n=14`, prove that the zero normalized class is
   the only full blocker.  Equivalently, every nonzero normalized vector has a
   safe micro-cell.

3. **Extremal half-turn lemma.**  Treat `2`-torsion vectors first.  The exact
   scan suggests the worst case is uniquely the coordinate-6 half-turn, whose
   eight explicit mirror intervals are already a certificate.  A proof could
   handle this cube by a parity/interval argument before addressing nonbinary
   residues.

The tempting next move is a shift-splitting argument.  Even shifts can hide
half-turn perturbations; odd unit shifts expose them.  Nonunit shifts, especially
`s=7`, should belong to quotient descent.  The remaining unit shifts should
force a positive-margin cell unless the vector is already scalar.

## Dead Ends And New Routes

The full normalized vector space is still far too large for direct exhaustion:
`14^12` representatives.  The small-support scans do not prove the whole
lemma, but they sharpen the terrain.  Adding support does not improve the
near-blocker through support `3`, and the exact binary cube has a unique
extremal.  That points away from a chaotic SAT obstruction and toward a
structured inequality over shifts and cell intervals.

One creative route is to define a "distance from scalar line" by the first
coordinate where `v_i` is not compatible with an affine ramp, then use odd
shifts to move that defect into a half-turn or near-half-turn window.  The
coordinate-6 certificate shows that a half-turn defect can be enough by itself.

Another route is to treat the eight missed intervals as a target stencil.  For
a general nonzero normalized vector, try to use the scalar gauge plus a unit
shift to make the residues agree with one of these stencils on enough
coordinates that the remaining coordinates have positive margin automatically.
