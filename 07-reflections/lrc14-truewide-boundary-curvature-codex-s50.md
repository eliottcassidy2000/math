# LRC14 true-wide boundary curvature

**Source:** codex-2026-06-20-S50.  User asked for another long LRC(14) proof
session, to pull/push frequently, integrate other agents, and think about
boundary functions, curvilinear convergence, bounded harmonic boundary values,
and maps of the Riemann sphere.

## The analogy that survived contact with the arithmetic

The useful part of the boundary-function analogy is the separation between an
endpoint boundary value and the path by which it is approached.  In Fatou-type
language, controlled nontangential approach and wild curved approach are not
the same object.  In the repo's Cayley-monad language, the Riemann sphere is
already the compactified home for rational transform data and boundary winding.

For true-wide LRC rows, the translation is:

```text
core B                      = boundary datum / state word
far speed u                 = one approach arc
two far speeds u,v          = two approach arcs interacting
p0(B union {u,v})           = endpoint value
I_B(u,v)                    = curvature of the approach
```

where

```text
I_B(u,v)=p0(B union {u,v})-p0(B union {u})-p0(B union {v})+p0(B).
```

Positive curvature means the two far speeds are complementary.  Negative
curvature means the second far speed overlaps the first.

## Exact result

The script `04-computation/lrc14_truewide_boundary_curvature_codex_s50.py`
stores `05-knowledge/results/lrc14_truewide_boundary_curvature_codex_s50.out`.
It scanned exact k=9 rows

```text
core=(0)+6-subsets of [1,14], far pair in [15,24].
```

Counts:

```text
raw rows       135135
primitive rows 135065
positive I     131003
negative I       3961
zero I            101
```

The direct-risk leader stayed exactly where HYP-2675 put it:

```text
E=(0,4,6,8,10,12,14,15,16)
p0=321/980
cap-p0=11681/70070
core=(0,4,6,8,10,12,14)=2*(0,2,3,4,5,6,7)
I_B(15,16)=-13/1470
```

That negative curvature is the main information.  The top true-wide row is not
a positive two-far synergy.  It is a dilated-core overlap row: the `16` peel is
large, the `15` peel is smaller, and together they overlap slightly.

The strongest positive-curvature row was:

```text
E=(0,1,4,8,10,12,14,16,20)
p0=89/336
cap-p0=11021/48048
I_B(16,20)=307/1960
```

So positive curvature is real, but the rows where it is largest are not the
cap-threatening rows in this box.

## What this does to HYP-2678

HYP-2678 says true-wide should split into a `d=1` scale-invariant finite route
and a `d>=2` signed dimension-penalty route.  The curvature scan agrees, but
sharpens the first half:

```text
d=1 / dilated core:
  finite curvature ledger, not a new analytic two-far theorem.

d>=2 / high-growth:
  positive curvature may happen, but observed direct p0 is lower;
  prove an excess/curvature cap or signed dimension penalty.
```

The KPS third-pocket row also clarifies the taxonomy.  Over core `(0,3,5)`, the
sorted far path `(16,28,30,33,35)` keeps `p0=0` until the last two additions;
all two-far curvatures over that tiny core are zero.  It is a multi-far
threshold/packet-cancellation row, not a two-far curvature exception.

After the atlas checkpoint, mac-mini S3 claimed THM-548 as the analytic
counterpart.  Its proposed decorrelated limit is

```text
Phi_2(B)=(2*p2(B)-p1(B))/49,
```

and the deviation is supported by resonance frequencies `m*u+n*v`.  That is
the exact analytic version of the boundary-function analogy: off resonance the
two approach arcs see independent boundary data; on resonance they share a
small integer relation, which is precisely the Freiman/scale-invariance finite
case.

## Next proof obligation

Prove a finite theorem for the d=1 branch:

```text
If the true-wide row has a dilated AP/AP-subset core plus O(1) dilation-breaking
extras, scale it to bounded height and verify the finite curvature ledger.
```

For the complement, use the exact curvature atlas as a gate:

```text
high excess or d>=2 + large positive curvature
  => direct p0 already has a comfortable cap gap

high excess or d>=2 + small/negative curvature
  => signed dimension-penalty / one-far Abel bounds should close it
```

THM-548 turns the second line into a concrete analytic problem: bound
`I_B(u,v)-Phi_2(B)` off the resonance set, then pass the resonant set to the
finite curvature ledger.

No LRC(14) proof is claimed here.  The gain is structural: the live k=9 leader
is now a finite d=1 overlap target, not an unbounded two-far synergy target.
