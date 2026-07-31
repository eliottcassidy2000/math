---
id: THM-2984
title: "Projected k3 signed-ray attainment and primitive-unit phase gate"
status: >
  PROVED.  On an exact periodic label ray with contribution A/z, existence
  of a scalar-admissible height is decided by a three-way sign/attainment
  law.  Retaining the primitive unit also makes every fixed-cell phase
  independent of height, giving a finite unit-by-unit strict-open exclusion
  gate.  Its cardinality shadow is sharp: more than
  beta(d)=2 floor((d-1)/14)+1 fixed-safe residues clear every unit, while
  beta(d) residues need not.  This is a reusable refinement of the projected
  k=3 denominator quotient; it does not assert that any new atlas row is
  empty and does not improve the current proved cap by itself.
source: codex-lrc14-k3-signed-ray-phase-gate-2026-07-30
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
related:
  - THM-2979-projected-k3-z275-ten-body-status-and-located-torsion-closure
  - THM-2981-projected-k3-z270-to-z247-cardinality-torsion-descent
---

# THM-2984 -- signed-ray attainment and primitive-unit phase gate

**PROVED.**

This theorem records the coordinate discarded by the current projected
`k=3` denominator envelope.  A denominator remembers the finite cyclic
group in which a high label moves, but not its primitive direction.  Taking
the maximum over those directions also conflates an attained value zero with
a negative ray whose supremum is zero and is never attained.  Both losses can
be avoided at finite cost.

Nothing below proves that a particular atlas row is empty.  To obtain such a
closure one must still supply the exact scalar residuals and a fixed-safe cell
certificate for every surviving primitive direction.

## 1. Exact signed-ray attainment

Let `L,z_0` be positive integers and put

```text
z_m=z_0+mL,                    m>=0.                     (1)
```

Let `A` and `sigma` be rational numbers.  Suppose the contribution on this
ray is exactly

```text
Delta(z_m)=A/z_m.                                         (2)
```

Here `sigma` is the contribution of all fixed labels minus the required
scalar floor.  Then

```text
there exists m>=0 with sigma+A/z_m>=0                    (3)
```

if and only if precisely the corresponding sign condition holds:

```text
A>0:    sigma+A/z_0>=0;
A=0:    sigma>=0;
A<0:    sigma>0.                                         (4)
```

For `A>0`, the sequence `A/z_m` is nonincreasing, so its maximum is attained
at `m=0`.  For `A=0`, it is identically zero.  For `A<0`, it is strictly
negative and increases to zero.  Thus `sigma<=0` makes `(3)` impossible,
including the boundary `sigma=0`.  If `sigma>0`, choose `m` so large that
`z_m>-A/sigma`; then `sigma+A/z_m>0`.  This proves `(4)`.

The strict inequality in the last line is essential.  Replacing it by
`sigma>=0` silently promotes an unattained supremum into a label.

For the THM-2941 projected rays, the inherited identity

```text
(z+L) delta(z+L)=z delta(z)                              (5)
```

gives `(2)` with `A=z_0 delta(z_0)` on every fixed nonzero residue ray.
Consequently `(4)` is an exact infinite-tail decision, not a search horizon.

## 2. Primitive-unit phase is height-free

Let `d|L`, let `1<=u<d` with `gcd(u,d)=1`, and suppose

```text
z_0 == (L/d)u  (mod L).                                  (6)
```

For an integer cell address `c`, evaluate the label at the grid point
`t=c/L`.  Equations `(1)` and `(6)` give

```text
z_m t == uc/d  (mod 1)                                  (7)
```

for every `m>=0`: the term `mc` is integral.  Thus the unbounded height has
disappeared, while the primitive unit `u` remains.

Write `r` for the least nonnegative residue of `uc` modulo `d`.  The weak
integer inequality

```text
14 min(r,d-r) >= d                                      (8)
```

is equivalent to

```text
||z_m t|| >= 1/14.                                      (9)
```

The lonely-runner danger comb is strict-open, so equality in `(8)` is safe.
If `t` is already outside the danger combs of all fixed labels, `(8)` proves
that `t` is outside the high-label comb for every height on this primitive
ray.

## 3. Finite ray-resolved closure gate

Fix a projected one-high case with high denominator `d`.  For every primitive
unit `u mod d`:

1. compute its exact amplitude `A_u` and first high label `z_0(u)`;
2. use `(4)` to discard the unit if no height can meet the scalar floor;
3. for each remaining unit, locate a fixed-safe cell address `c` satisfying
   `(8)`.

If step 3 succeeds for every scalar-reachable unit, the one-high case is
empty uniformly over the infinite label tail.  The witness cell may depend on
`u`; requiring one pair of cells to work for all units is sufficient but not
necessary.

This gate is strictly more informative than first maximizing over primitive
units and then asking for one unit-blind torsion pair.  It retains two
sidecars that the coarser quotient destroys:

```text
primitive direction u,          zero attained versus zero only approached.
```

Strictness already appears at `d=8`.  Take a fixed-safe residue set
`S={1}`.  It contains no pair, hence no pair-difference torsion certificate
of any order.  But every primitive unit is one of `1,3,5,7`, and

```text
min(u,8-u)>=1,             so 14 min(u,8-u)>=14>8.
```

Thus the absolute-cell test `(8)` closes all primitive directions.  The
ray-resolved gate can therefore certify a case that every pair-only gate
necessarily leaves undecided.

The three boundary controls are:

```text
A>0, sigma=-A/z_0:  admissible, attained at the first high point;
A=0, sigma=0:       admissible at every height;
A<0, sigma=0:       inadmissible at every height although the supremum is 0.
```

They separate all weak/strict directions in `(4)`.  The next computational
use must still pin the exact universe, upstream hashes, ordinary and optimized
transcripts, and positive and hostile controls before claiming an atlas
closure.

## 4. Sharp cardinality shadow of the unit gate

The unit-resolved test also has a closed cardinality threshold which is
different from the pair-collision threshold.  For `d>=2`, let

```text
B_d={r in {0,...,d-1}:14 min(r,d-r)<d}.
```

Put

```text
beta(d)=2 floor((d-1)/14)+1.                            (10)
```

Then

```text
|B_d|=beta(d).                                         (11)
```

Indeed, if `b=floor((d-1)/14)`, the strict inequality defining `B_d`
is equivalent to

```text
r in {0,1,...,b} or r in {d-b,...,d-1}.
```

The displayed intervals are disjoint because `b<d/14`; they have `b+1` and
`b` elements, respectively, so their union has `2b+1` elements.  Strict
openness is visible in the `d-1` in `(10)`: residues at exact distance
`d/14` are safe.

Let `S` be any set of fixed-safe cell residues modulo `d`.  Multiplication by
a primitive unit `u` permutes the residue classes, so the set of cells bad for
that unit is `u^{-1}B_d` and also has `beta(d)` elements.  Consequently

```text
|S|>beta(d)  implies that for every primitive u mod d
             some c in S satisfies 14 min(uc mod d,d-(uc mod d))>=d.    (12)
```

This is sharp among arguments using only `|S|`: at equality take
`S=u^{-1}B_d`, for which every cell is bad for that unit.  Thus `(12)` is the
exact independence-number analogue for the bipartite unit--cell incidence
graph.  It complements, rather than subsumes, the Cayley pair threshold
`|S|>d/R` in THM-2979:

```text
pair view:  S must contain a short-order difference;
unit view:  S must escape every multiplicative translate of B_d.
```

There is also a useful cell-count form.  If `d|L` and `C` complete grid cells
map to `S`, each residue modulo `d` has at most `L/d` preimages.  Hence

```text
C>beta(d)L/d  implies |S|>beta(d),                     (13)
```

and the unit gate closes all primitive directions.  When `(13)` fails, the
actual shape of `S` can still make every intersection `S ∩ u^{-1}B_d`
proper; this is exactly the information retained by the finite ray-resolved
test and lost by a count-only quotient.
