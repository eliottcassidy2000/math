# Pivot capacity is scale-dependent, and full-X still needs left covariance

**Status: FINITE-EXACT + VERIFIED in normal and optimized modes.**  The
half-cycle admits a unique maximin support selector, and every pivot has full
moving `C_13` character support on the full THM-2707 packet bank.  Neither
fact constructs a THM-2334 endpoint pair.  No scalar row is excluded and
LRC(14) remains open.

## Inheritance pass

The closest proved mechanism is the paired-cylinder addendum of
`THM-2701-literal-singleton-word-one-step-dilation-nilpotence.md`: on both
THM-2698 event cylinders, each provisional pivot

```text
k in {27,40,53,66}
```

has all twelve moving `C_13` characters.  The corrected near miss is
`THM-2563-paired-dipole-deep-target-corner-and-partial-bare-boundary.md`:
the moving term `-eta.v` is not the endpoint difference `eta.(u-v)` because
the left Fourier multiindex is absent.  The canonical hostile is
`THM-2568-full-x-transition-annihilation-and-refined-pair-drift-boundary.md`:
a lawful common danger-to-safe endpoint completion has zero full-`X` current
in every coarse target character.  The least-used sidecar is the residue-part
structure of `THM-2707-full-physical-lift-fibre-common-simplex-and-packet-scc.md`,
now sharpened by the finite-exact `304`-endpoint / `3042`-midpoint atom macro.

The live concepts were therefore:

```text
right moving-shift support;
owner-pivot choice;
full physical packet SCC;
atom-based two-step macro;
left/right spectral covariance;
full-X pair-to-target pushforward.
```

The computation below separates them instead of calling every cyclic label a
target residue.

## 1. The prime-cyclotomic rank test

For a packet centre `x`, private root `r`, and provisional pivot `k`, define

```text
S_(x,k)={s in F_13:
  d_1(2*13^5 x-r/13)=1,
  u_1(13^3 x-s/13)=1,
  u_1(kx+s/13)=1}.                                  (1)
```

On every cylinder used below the three factors are constant, so the physical
integral is the cylinder length times `1_(s in S_(x,k))`.

The exact Fourier test is especially rigid at the prime `13`.  If
`c=(c_0,...,c_12)` is rational, then

```text
sum_s c_s zeta_13^(hs)=0 for one h!=0
iff c_0=c_1=...=c_12.                               (2)
```

Indeed the degree-at-most-twelve polynomial `sum c_s X^s` then has a
primitive thirteenth root, hence is a rational multiple of
`Phi_13=1+X+...+X^12`.  Thus rank after deleting the constant direction is
exactly moving-Fourier rank, and a nonconstant aggregate has **all** twelve
moving characters, not merely one.

## 2. The half-cycle has a unique maximin pivot

For the two THM-2698 centres, the eight supports are exactly

```text
event 0:
  27: {0,3,4,5,8,9,10,11,12}
  40: {0,1,2,3,4,5,8,9,10,11}
  53: {0,1,2,3,4,5,8,11,12}
  66: {0,1,2,3,4,5,8,9,10,11,12}

event 1:
  27: {1,2,3,4,7,8,9,10,11}
  40: {1,2,3,4,5,6,7,8,9,10}
  53: {1,2,3,6,7,8,9,10,11}
  66: {1,2,3,4,5,6,7,8,9}.                         (3)
```

At either event the four support rows have rational rank four and moving
Fourier rank four.  No nonzero rational weighting of the four pivots makes
the aggregate constant; the same kernel remains zero when both events are
imposed simultaneously.  Equal-pivot aggregation is nonconstant at each
event, so it preserves all twelve moving colours.

Use the event-symmetric capacity triple

```text
Cap(k)=(min_i |S_(i,k)|,
        sum_i |S_(i,k)|,
        product_i |S_(i,k)|).
```

Then

```text
Cap(27)=(9,18,81),   Cap(40)=(10,20,100),
Cap(53)=(9,18,81),   Cap(66)=(9,20,99).             (4)
```

Hence `k=40` is the unique lexicographic maximin selector on this particular
two-event carrier.  This is not a hidden symmetry choice: the simultaneous
`AGL(1,13)` automorphism group of the labelled pair of support systems is the
identity, and no affine shift reindexing exchanges the two events.

Equation `(4)` is a useful positive selection rule, but it is **extrinsic**
to THM-2309/2350.  Those theorems permit any distinct graft units after an
owner pivot is chosen; they do not assert that physical paired-support
capacity chooses the relation allocation.

## 3. The full packet orbit has every colour but no stable capacity selector

The second universe is all `13^6=4,826,809` physical addresses, with exactly
the THM-2707 filters:

```text
private root/unit residue;
literal rail support;
present support.
```

This recovers `3346` packets.  At each packet `n`, equation `(1)` is evaluated
at

```text
q_n={z+7n/13^6},                 root(n)=6+n mod13. (5)
```

Every inserted factor is constant on the whole inherited open interval `I`;
the smallest exact inserted-wall radius is `56447` times the pulled-back
radius of `I`.  Thus this is a whole-open-cylinder result, not a midpoint
sample.

All `3346*4=13384` support sets are nonempty and proper.  By `(2)`, every
one has all twelve moving characters.  The `13384 x 13` support matrix has
full rank `13`; its quotient by constants has rank `12`.  No fixed nonzero
weighting of the four pivots makes the support aggregate constant at every
packet.  At every individual address the equal-pivot aggregate is also
nonconstant.

The exact size census is frozen in the companion output.  Its sharp selector
features are:

```text
unique local winners:
  k=27 at 467 packets,  k=40 at 122,
  k=53 at 244,          k=66 at 394;

all-four tie: 1375 packets;

full-orbit support masses:
  27:33256, 40:32867, 53:32856, 66:32581.           (6)
```

Thus every candidate uniquely wins somewhere, while more than forty percent
of the packets do not select among them at all.  Total-mass aggregation
uniquely chooses `k=27`, not the half-cycle maximin choice `k=40`.  The
displayed THM-2707 triangle and eleven-cycle instead tie `27` and `40` at
worst support `11`.  This is an exact scale/carrier conflict:

```text
THM-2698 two-event maximin:       k=40;
full THM-2707 mean capacity:      k=27;
displayed THM-2707 cycle floor:   tie {27,40}.       (7)
```

Capacity can define a pivot only after the carrier, aggregation law, and
objective have been retained.  It is not a carrier-independent canonical
allocation.

## 4. The atom macro is a base carrier, not a pivot action

The following-atom macro splits the same packet bank into

```text
A={304 residue-zero atom endpoints},
T={3042 physical transit packets}.                       (8)
```

Every `A -> T` edge is physical, so there are

```text
304*3042=924768                                          (9)
```

atom-to-transit edges and every such edge can be returned to an atom
endpoint.  The pivot profiles refine this promising based carrier, but do
not transform covariantly on the whole relation.

For each edge the scout asks whether the four support profiles at its two
vertices agree exactly, or become equal after one **common** reindexing of
the shift coordinate.  The exact joint-profile counts are

```text
exact equality:                         26857 / 924768;
translation / dihedral / AGL equality:  29299 / 924768. (10)
```

Even allowing all of `AGL(1,13)`, more than `96%` of the physical edges have
no joint four-pivot support transport.  No individual pivot works uniformly
either; the AGL-compatible edge counts are

```text
k=27:153847,  k=40:163227,
k=53:171929,  k=66:156987,                               (11)
```

all strictly below `924768`.

Two exact controls locate the boundary:

```text
0 -> 1 -> 0,       signed numerators (7,-7):
  preserves the exact four-pivot profile;

0 -> 106 -> 0,     signed numerators (742,-742):
  fails joint profile matching even up to common AGL(1,13). (12)
```

Therefore the `304` endpoints and `3042` midpoints are a real physical base
carrier on which a spectral cospan might be built.  But the packet relation
itself does not supply a pivot-covariant target action.  Agreement on the
small loop in `(12)` is a positive seed, not covariance on the macro bank.
The concurrent all-depth phase-torus probe gives the complementary reason:
physical diagonal translation and the target anti-diagonal dipole intersect
only at the identity for these blocker/graft pairs.

## 5. Moving aggregation is not full-X endpoint aggregation

The data in Sections 2--4 is a right moving-shift profile `K(t)`.  It has no
left twist variable.  The missing datum is visible without a numerical
counterexample.  Two joint tables have the same observed fixed-left row:

```text
A_frozen(s,t)=K(t),
A_covariant(s,t)=K(t-s),

A_frozen(0,t)=A_covariant(0,t)=K(t).                 (13)
```

Their common-twist diagonals are opposite:

```text
A_frozen(s,s)=K(s),       nonconstant, hence all moving colours;
A_covariant(s,s)=K(0),    constant, hence no moving colour.       (14)
```

Thus the right profile, even with every primitive character nonzero and even
after pivot/SCC aggregation, cannot determine the coarse endpoint verdict.
The distinction is exactly the missing covariance of the left factor.

For an actual THM-2334 current the indices `u,v` are global Fourier
multiindices, and the target is `eta.(u-v)`.  A packet address, lift residue,
root, owner, or shift label is not either index.  If one lawfully completes
the old dangerous head and the repaired safe endpoint by a common target
twist, THM-2568 applies and gives

```text
sum_X Jhat_X(q)=0             for every q in F_13.   (15)
```

Consequently full-pivot or full-packet aggregation preserves the **moving
right spectrum**, but full-`X` danger-to-safe recombination does not restore
an endpoint amplitude.  The high-leverage next object is a factorwise
THM-2350-covariant left/right spectral cospan on one atom-macro ancestry at a
fixed marked triangle `(X,m,Y)`, or a proved normal/jet/orientation which
weights `X` before the annihilating sum.

## 6. Exact reproduction

Run

```bash
python3 04-computation/lrc14_paired_pivot_full_packet_aggregation_probe_20260728.py
python3 -O 04-computation/lrc14_paired_pivot_full_packet_aggregation_probe_20260728.py
```

Both modes byte-match

```text
05-knowledge/results/lrc14_paired_pivot_full_packet_aggregation_probe_20260728.out
```

with LF-normalized SHA-256 values

```text
script  24a5f2c087b2a376c0eb7311f0dcd7519b8e6deb9b5585d6dd3c199b3a7d1251
output  cfebce49986323322139d51f39693a9a0f633c0b7dcf7ad7773c5c52328f41a6
```

The output records the complete support-size and argmax censuses, every
aggregate profile, exact rational ranks and kernel dimensions, orbit counts,
atom/transit splits, individual and joint covariance counts, both macro
controls, and the full-`X` scope boundary.
