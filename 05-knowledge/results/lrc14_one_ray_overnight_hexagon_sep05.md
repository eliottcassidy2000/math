# Complete one-ray closure and a sharp gap for every non-norm-four projection

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The `nc2_seed` referee checked the all-height argument, complete support and
strict cutoff, equality scope, short-relation corollary, and every finite-head
row against the independent physical-sheet representation. This report does
not use a scarce theorem identifier or treat reserved THM-4423 as a proved
dependency.

Let `w=(a,b,c)` be primitive, distinct, positive, odd, and ternary-unit, with
`a<b<c`. Let `Lambda(w)` be the **complete** owner-live raw carrier set of
[THM-4422 — projection deficit and Beatty-row reduction](../../01-canon/theorems/THM-4422-lrc14-projection-deficit-and-beatty-row-reduction.md),
and let `E_i` denote its three degree-zero network projections in coordinate
order. Suppose the real span of `Lambda(w)` has dimension at most one.

Then `min_i E_i<=6/77` at every height. If the support is nonempty and its
primitive direction `v` does not have `l1` norm four, the stronger sharp bound
is

```text
E_i <=12/161 <6/77 for every i.                         (1)
```

Equality in (1) occurs exactly for the primitive triple `(1,19,23)` and
coordinates `i=2,3`. Its physical mass is only `12/437`: the theorem controls
the actual network upper certificates, not just the physical comb mass.
Signed norm-four directions use the already proved THM-4422 theorem.
Empty support gives three zero capacities.

Consequently every proposed counterexample to the universal network target
must have at least two distinct primitive live directions. In particular,
all triples admitting a full-support ternary-unit relation of `l1` norm at
most fourteen now satisfy the network target; the older shell results
controlled their physical mass. No bound for arbitrary multi-ray support,
chart entry, synchronization, or LRC(14) is claimed.

## Inheritance and the paper-inspired probe

The closest proved mechanism is THM-4422's residue-deleted ordinal count for
the signed norm-four families, combined with
[THM-4414 — six-separated contact capacity collapse](../../01-canon/theorems/THM-4414-lrc14-six-separated-contact-capacity-collapse.md),
which identifies these exact raw sums with the component networks. The
equality hostile is `(1,5,11)`. The corrected near miss is THM-4422's
refutation of every fixed convex average: a selected projection is a maximum
deficit consumer, not an average. The least-used coordinate here is the
largest **primitive** carrier coefficient; it bounds the entire ray's
eligible length without dropping the raw multiplier.

The [empty-hexagon paper, Section 4.1](https://arxiv.org/html/2403.00737v1#S4.SS1)
uses an extremal interior point to replace a vertex while preserving a
smaller witness from which an empty hexagon follows. The transferable research
move is to identify exactly what a local witness allows one to discard.
It is not a theorem about LRC carriers: central symmetry, collinear
multiples, a missing zero vector, and residue-deleted lattice points violate
the paper's native point-set setting. No SAT relaxation is used here.

The concept board is:

| Concept | Question / operation | Retained coordinate or obstruction |
|---|---|---|
| Primitive live direction | Choose a direction as local witness | Must first prove it contains the complete support |
| Maximum coefficient `M` | Select the narrowest coordinate strip | Raw multiplier `k` remains present |
| Mod-three ordinal list | Delete every third multiplier | Both signs and strict endpoint remain |
| Network projections | Bound each summand by its common cap | Does not identify a capacity with physical mass |
| Multi-ray replacement | Replace all carriers by one minimal witness | Fails at `(17,23,25)`; another live direction is lost |

The connection contract is therefore scoped: complete rank-one raw support
maps to a residue-deleted interval of integer multipliers; this preserves
every summand and owner predicate. Its count controls every network capacity.
Replacing an arbitrary support by one direction is not part of the map.

## 1. Complete one-ray address and the narrow coordinate

Write

```text
r=3/14, q=3/(7c),
Lambda(w)={C in Z^3 : w.C=0, C_i!=0 mod3,
                         |C_i|<r(sum(w)-w_i) for every i}.
```

The set is centrally symmetric. If nonempty and contained in one real line,
choose its primitive integer direction `v`, oriented by making its first
nonzero coordinate positive. Every coordinate of `v` is a ternary unit:
if a coordinate vanished modulo three, no multiple could be owner-live.
Every integer point on the line is `k v`. Conversely, every such point
satisfying the displayed strict inequalities and residue gate belongs to the
complete support. Thus, with

```text
B=min_i 3(sum(w)-w_i)/(14|v_i|),
K=strict_floor(B),
Lambda(w)={+/-k v : 1<=k<=K, 3 does not divide k},
N=|Lambda(w)|=2(K-floor(K/3)).                           (2)
```

This retains the full raw address `k`, not just its projective direction.
The strict floor is the greatest integer strictly less than `B`; ordinary
floor would wrongly add zero-roof endpoints when `B` is integral.

Since every speed is odd, `w.v=0` implies that `sum_i |v_i|` is even.
If `max_i |v_i|<4`, the unit condition restricts every absolute coefficient
to `1` or `2`. Primitivity excludes `(2,2,2)`; parity excludes the patterns
with odd total. The only remaining pattern is `(1,1,2)`.
Therefore a direction outside signed norm four has

```text
M=max_i |v_i| >=4.                                     (3)
```

Choose a coordinate attaining `M`. The other two speeds sum to strictly less
than `2c`, so

```text
B<3c/(7M)<=3c/28.                                      (4)
```

This is the extremal selection that is legal: one narrow coordinate bounds
the whole **already proved one-dimensional** support. It does not discard
any component.

## 2. Counting pays every projection

Writing `K=3t+j`, with `j=0,1,2`, gives

```text
K-floor(K/3)=2t+j <=(2K+2)/3.
```

Since `K<B`, (2) gives the strict envelope

```text
N <=(4K+4)/3 <4B/3+4/3.
```

Each summand in every THM-4414 projection is at most `q`. Consequently,

```text
E_i <=qN <12/(49M)+4/(7c)
             <=3/49+4/(7c),                  every i.  (5)
```

For `c>=43`, the rightmost expression is at most

```text
3/49+4/(7*43)=157/2107 <12/161.                         (6)
```

Thus only the complete finite head `c<43` remains for the sharp bound.
This is an analytic height reduction within the complete one-ray class;
it is not an extrapolation of a height census.

## 3. Exact finite head and sharpness

There are precisely `363` primitive sorted distinct positive odd ternary-unit
triples with `c<43`. The verifier exhausts all of them, before applying the
support filter. Their complete dictionaries split as follows:

| support class | count |
|---|---:|
| empty | 30 |
| signed norm-four ray | 41 |
| one ray outside norm four | 236 |
| two or more directions | 56 |

Every one of the `236` relevant rows satisfies (1). The only row attaining
the bound in any coordinate is

```text
w=(1,19,23), v=(4,1,-1), B=9/4,
Lambda(w)={+/-v,+/-2v},
(E_1,E_2,E_3)=(12/437,12/161,12/161),
physical mass=12/437.                                  (7)
```

Equations (5)--(6) make every larger-height row strict, establishing the
claimed all-height equality classification. For the norm-four branch,
THM-4422 proves `min_i E_i<=6/77` and gives equality only at `(1,5,11)`.
Therefore the whole one-ray theorem has that same sole target equality.

The stronger claim for **every** projection genuinely needs the norm-four
exception. For example

```text
w=(1,5,7), v=(2,1,-1),
(E_1,E_2,E_3)=(8/245,6/49,4/35),                       (8)
```

so two projections exceed the target even though the selected one closes it.

## 4. Short-relation corollary and an unbounded family

[THM-4386 — canonical component relation and zero-defect incidence,
Lemma 2.1](../../01-canon/theorems/THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence.md)
proves that every full-support ternary-unit relation `u.w=0` with
`||u||_1<=14` has zero defect at every physical component.
Explicitly, if `C=w cross n` is its raw address, distinct owners imply
`u.n=0 mod3`, while strict eligibility implies
`|u.n|<(3/14)||u||_1<=3`. Hence `u.n=0`.

Both `u` and `C` annihilate `w,n`. When the owner gate is live, `w,n` are
independent; their common orthogonal is a line. Every live `C` is therefore
parallel to `u`. Applying the theorem above upgrades this entire short-unit
relation sector to a **network** closure, including the non-norm-four
every-projection bound (1). The previously known physical shell ceilings do
not by themselves prove this stronger consumer.

The support condition is much broader than any bounded count. For example,
the family

```text
w_t=(1,6t+1,6t+5), t>=1,
v=(4,1,-1),
B=(18t+9)/28,
N=2(strict_floor(B)-floor(strict_floor(B)/3))             (9)
```

has complete one-ray support by its norm-six relation. Its carrier count is
unbounded, while the theorem bounds every projection. The exact control
`(1,1201,1205)` has `172` live carriers and

```text
(E_1,E_2,E_3)=(310116/10130435,516/8435,516/8435).
```

Thus this analytic result does not duplicate the separate reserved
THM-4423 bounded-atlas/ninety-carrier program. At the startup revision
`62e1b34d82a4`, that theorem is an honest empty reserved stub and asserts
no mathematical result.

## 5. Exact obstruction to an unlicensed extremal replacement

The smallest multi-ray triple, ordered first by maximum speed and then
lexicographically, is

```text
w=(17,23,25),
Lambda(w)={+/- (8,-7,1), +/- (5,5,-8)}.                 (10)
```

The shortest live primitive direction in `l1` norm is `v=(8,-7,1)`.
Its own cutoff is `B=9/7`, so its exact ordinal list contains only `+/-v`.
Treating that local witness as if it carried the whole support would predict
two components. Even the loose count envelope would give

```text
N<4(B+1)/3=64/21,
```

whereas the actual support has `N=4`. The first failed implication is
existence of a minimal live direction to completeness of that direction.
The omitted coordinate is the second owner-live ray `(5,5,-8)`.
The strongest survivor is the exact per-ray address/count formula; the
repair is to keep every ray or separately prove collinearity.

For reference the true projections in this hostile are

```text
(106/4025,12/425,2546/68425),
physical mass=744/68425.
```

This does not refute the universal network target: all three projections
here are below `6/77`. It refutes the proposed replacement step. The later
`(19,23,29)` configuration from THM-4422 remains the first height-79 dense
multi-ray control. Both hostiles are retained in the verifier.

## 6. Reproduction, independence, and open scope

The script
[lrc14_one_ray_overnight_hexagon_sep05.py](../../04-computation/lrc14_one_ray_overnight_hexagon_sep05.py)
imports no repository mathematics. One implementation enumerates the complete
integer relation box and computes exact raw roof sums. A second clears
denominators and constructs all six literal shifted sheet configurations;
it computes every projection from contact-edge minimum lengths and physical
mass from literal overlaps. The representations agree on every one of the
363 finite-head rows and all wide controls. Canonical equality, the
norm-four every-projection hostile, empty support, strict endpoints,
owner-residue deletions, and multi-ray hostiles are explicit controls.

```bash
python3 -B 04-computation/lrc14_one_ray_overnight_hexagon_sep05.py
python3 -B -O 04-computation/lrc14_one_ray_overnight_hexagon_sep05.py
```

The finite head is part of the proof; the wide controls are only verification.
All checks are explicit and remain active under optimization. The open
consumer is now entirely multi-ray after combining this theorem with
THM-4422. Neither that reduction nor the short-relation corollary supplies
entry, synchronization, or LRC(14).

The normal and optimized output streams are byte-identical, recording `7,478`
explicit checks. Frozen raw-byte SHA256:

```text
source 6b41a879700632aa934650f27dafe9d076c051ddcee3262fabc07556a7aaf875
output c098ee8a43644f349d21f16257596c84991b732596059bb603db1769cbb73a2f
```
