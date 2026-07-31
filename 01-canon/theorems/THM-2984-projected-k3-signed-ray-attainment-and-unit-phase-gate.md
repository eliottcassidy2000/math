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
  beta(d) residues need not.  Gcd-stratum capacities and the exact
  multiplicative transporter refine this threshold whenever residue shape is
  retained.  The resulting transport complex is flag through d=42 and first
  becomes non-flag at d=43, with an explicit irreducible three-cell
  certificate.  This is a reusable refinement of the projected k=3
  denominator quotient; it does not assert that any new atlas row is empty
  and does not improve the current proved cap by itself.
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
graph.  Its representation differs from the Cayley pair threshold
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

## 5. Gcd-stratum refinement and the exact obstruction

The total count in `(12)` is not the end of the finite reduction.  Let

```text
O_g={r mod d:gcd(r,d)=g},
```

where `g|d`.  Every primitive unit preserves every `O_g`, because

```text
gcd(ur mod d,d)=gcd(ur,d)=gcd(r,d).                    (14)
```

Therefore the same injection argument can be applied inside one orbit
stratum:

```text
|S intersect O_g|>|B_d intersect O_g|
   for some g|d
implies that every primitive u has a c in S with uc notin B_d.          (15)
```

This can be stronger than `(12)`: a large part of `B_d` may lie in gcd
classes which `S` never uses.  The capacities in `(15)` are themselves
explicit.  With `b=floor((d-1)/14)`, the zero stratum has

```text
|B_d intersect O_d|=1.
```

For a proper divisor `g<d`, reflection `r -> d-r` pairs the two bad flanks
and gives

```text
|B_d intersect O_g|
 =2 #{a:1<=a<=floor(b/g), gcd(a,d/g)=1}.                (16)
```

Indeed every nonzero bad residue on the lower flank is uniquely `r=ga` with
the displayed conditions, and the two flanks are disjoint since `b<d/14`.
Thus the shape-sensitive refinement costs only partial-totient counts.  For
prime `d`, for example, the unit stratum has capacity `beta(d)-1`: the extra
bad residue counted by `beta(d)` is zero and cannot receive a unit cell.

The intrinsic remaining object is the multiplicative transporter

```text
T_d(S,B_d)={u in (Z/dZ)^*:uS subset B_d}.               (17)
```

If `U_active` is the set of primitive directions surviving the signed
attainment law `(4)`, the absolute-cell gate is unresolved exactly when

```text
T_d(S,B_d) intersect U_active is nonempty.              (18)
```

This separates scalar reachability from residue geometry.  It also gives two
cheap compressions below the stratum threshold.  If `S` contains a unit
`s`, then `us` determines `u`, so only unit residues of `B_d` can index a
transporter candidate.  If `|S|=|B_d|` and one `u_0` transports `S` into
`B_d`, cardinality forces `u_0S=B_d`.  Comparing any other transported image
with this equality, and conversely applying any stabilizer to it, gives the
exact identity

```text
T_d(S,B_d)=Stab(B_d) u_0.                               (19)
```

The stabilizer in `(19)` is explicit.  For `2<=d<=14`, `B_d={0}`, so every
unit stabilizes it.  For every `d>=15`,

```text
Stab(B_d)={+1,-1}.                                      (20)
```

To prove this, note first that `1 in B_d`; hence a stabilizing unit has a
least absolute representative of size `1<=a<=b`.  Since reflection already
stabilizes `B_d`, replace the unit by its negative if necessary and take this
representative to be `+a`.  If `a>=2`, put
`k=floor(b/a)+1`.  Then `k<=b`, while

```text
b<ak<=b+a<=2b<d-b.
```

Thus `k` is bad but `ak` is not, a contradiction.  Reflection shows that
both signs do stabilize the symmetric band.  In particular, at the sharp
equality boundary and `d>=15`, a nonempty transporter has exactly two
elements, `u_0` and `-u_0`; there is no large hidden multiplicative symmetry.

The same stabilizer acts on every transporter, without the equality-size
hypothesis:

```text
u in T_d(S,B_d), s in Stab(B_d)  implies  su in T_d(S,B_d).             (20a)
```

Hence every transporter is a union of stabilizer cosets.  For `d>=15` it is
in particular closed under `u -> -u`, and every forced stabilizer coset has
two elements.
The mask can therefore be stored on the projective unit set

```text
PU_d=(Z/dZ)^*/{±1}.                                     (20b)
```

This remains exact even when scalar reachability distinguishes the two
signs.  Mark a projective class active when either representative belongs to
`U_active`; then `(18)` is nonempty exactly when the projectivized transporter
contains an active class.  Thus reversal-paired rays halve the finite unit
search without assuming their amplitudes agree.

Equations `(15)--(20b)` identify the correct recursive search object: a
bipartite unit--cell incidence relation with multiplicative orbit sidecars,
not a tournament.  Only after this transporter is nonempty is a
unit-indexed carrier-overlap or spanning-tree calculation warranted.

## 6. Every short-order pair is already a universal unit certificate

The pair and unit views are not merely parallel.  A located pair certificate
feeds the unit gate directly.  Suppose `c_1,c_2 in S` are distinct and their
difference modulo `d` has additive order `2<=r<=7`.  For every primitive
unit `u`, the difference `u(c_2-c_1)` still has order `r`.  Writing an
order-`r` residue as `(d/r)a` with `gcd(a,r)=1` shows that its circular
distance from zero is at least

```text
d/r >= d/7.                                             (21)
```

If both `uc_1` and `uc_2` lay in `B_d`, their circular separation would be
strictly less than `d/14+d/14=d/7`, contradicting `(21)`.  Thus for every
primitive `u`, at least one member of this same pair satisfies the absolute
cell test `(8)`.  This extends THM-2979's exact order-seven argument verbatim
to every order at most seven.

Consequently the logical hierarchy is one-way:

```text
located order-at-most-seven pair
  => universal primitive-unit cell escape
  => closure of every scalar-active primitive ray.       (22)
```

The converse fails: the `d=8`, `S={1}` control in Section 3 has no pair but
escapes the bad band for every primitive unit.  The unit-resolved gate is
therefore a genuine extension of pair-only closure, while the pair witness
remains a cheaper unit-blind certificate when present.  This also explains
why a terminal already closed by THM-2979's located torsion test must pass the
THM-2984 absolute-cell test: agreement of those two audits is a theorem, not
independent numerical luck.

The count-only comparison has the same direction.  If `R` is THM-2979's
largest divisor of `d` at most seven, then

```text
beta(d) <= ceil(d/7) <= d/R.                            (23)
```

The first inequality follows immediately from `(10)`, and the second from
`R<=7` and `R|d`.  Thus the pair-density hypothesis `|S|>d/R` always implies
the absolute-unit count hypothesis `|S|>beta(d)`.  The gain can be strict:
for `d=14m`, `beta(d)=d/7-1`.  Retaining absolute cell addresses can therefore
save a whole residue fibre before any transporter-shape refinement is used.

This admits a useful graph-theoretic rephrasing without inventing an
orientation.  The bad band `B_d` is itself an independent set in
THM-2979's Cayley graph `G(d,7)`: two bad residues have circular separation
strictly below `d/7`, whereas every nonzero residue of additive order at most
seven has circular size at least `d/7`.  Hence

```text
pair gate:       S is not independent in G(d,7);
all-unit obstruction: S is contained in one multiplicative-unit translate
                      of the particular independent set B_d.          (24)
```

For a restricted scalar-active set, `(24)` is intersected with `U_active`
exactly as in `(18)`.  The second condition is far more rigid than mere
pairlessness.  Thus after a pair gate fails, the next object is not an
arbitrary independent set: it is the recognition problem for the unit orbit
of one short symmetric band, refined first by gcd strata and then by its
two-element stabilizer.  This is the exact additive--multiplicative interface
which the denominator-only quotient had erased.

## 7. The transport complex is genuinely higher-order

There is a still more faithful finite object.  Define

```text
K_d={S subset Z/dZ:T_d(S,B_d) is nonempty}.              (25)
```

This is a simplicial complex: it is downward closed, and its maximal faces
are the unit translates `u^{-1}B_d`.  Its one-skeleton remembers whether two
cells can be bad for one common unit.  It does **not** determine the whole
complex.  Since `0` belongs to every translate,

```text
T_d(S union {0},B_d)=T_d(S,B_d).                        (25a)
```

Thus `K_d` is always a cone with apex `0`, more precisely the join of that
vertex with its deconed nonzero complex.  The full complex is contractible;
the decone carries all minimal nonfaces and flag failures and is the
meaningful object for any nontrivial topology.

A clean hostile control occurs at the prime modulus `d=43`.  Here

```text
B_43={0,±1,±2,±3}
S={1,2,15}.                                              (26)
```

Every two-element subset of `S` is a face of `K_43`, with the following
explicit primitive multipliers:

```text
u=1:   {1,2}   -> {1,2};
u=3:   {1,15}  -> {3,2};
u=23:  {2,15}  -> {3,1}                 modulo 43.       (27)
```

But `S` itself is not a face.  Indeed, if `uS subset B_43`, the image of `1`
forces `u in {±1,±2,±3}`.  Requiring `2u in B_43` leaves only `u=±1`,
and then `15u=±15` is not in `B_43`.  Hence

```text
T_43({1,2,15},B_43)=empty.                              (28)
```

Thus `S` is a three-clique in the one-skeleton of `K_43` but is not a
two-simplex: `K_43` is non-flag.  Moreover all nonzero differences in `S`
have additive order `43`, since `43` is prime, so `S` is independent in
`G(43,7)` and no short-order pair certificate is present.

This is the first non-flag modulus.  For `2<=d<=42`, the band radius
`b=floor((d-1)/14)` is at most two, and the three possible ranges have a
uniform description.

* For `b=0`, the decone is empty.
* For `b=1`, its facets are the disjoint antipodal pairs `{±a}`, so it is
  flag.
* For `b=2`, its facets are `{±a,±2a}`.  Projectivize each antipodal pair.  If
  `d` is even, the classes `[a]` and `[2a]` lie in the distinct gcd strata
  `O_1` and `O_2`; the projective compatibility graph is bipartite and has no
  triangle.  If `d` is odd, multiplication by `[2]` permutes the projective
  unit classes and the compatibility graph is a disjoint union of its orbit
  cycles.  A triangle would force `[2]^3=[1]`, equivalently
  `8 == ±1 (mod d)`, so `d` would divide `7` or `9`, impossible for
  `29<=d<=42`.

In either `b=2` case, a clique uses at most two projective classes and is
contained in the edge witnessing those classes.  Every projective edge lifts
to its entire sign-saturated four-vertex facet `{±a,±2a}`; restoring the
antipodal pairs therefore replaces each projective vertex by a simplex and
preserves flagness.
Together with `(26)--(28)`, this proves the exact threshold

```text
K_d is flag for 2<=d<=42, while K_43 is not flag.        (28a)
```

The exact relation with the short-order graph is

```text
E(G(d,7)) subset {two-element nonfaces of K_d},
and this inclusion can be strict.                                  (29)
```

Indeed a short-order pair cannot lie in one translate of `B_d` by Section 6,
while the converse need not hold: for `d=15`, the pair `{1,2}` has difference
of order `15` but cannot be carried into `B_15={0,±1}`.  Even the resulting
unit-compatibility graph does not determine the full complex, by `(26)--(28)`.

The missing face `(26)` is an irreducible three-cell unit certificate: no
singleton or pair closes the universal unit gate, while the triple does.  A
pairwise orientation would discard still more information.  The natural
recursive sidecar is therefore the full intersection mask

```text
M_d(S)=intersection_{c in S}{u:uc in B_d}=T_d(S,B_d).    (30)
```

Its nonemptiness pattern is the face poset of `K_d`; the mask itself retains
which units survive, information still needed when only `U_active` is
scalar-reachable.  Small minimal nonfaces, and active-unit versions obtained
by intersecting `(30)` with `U_active`, are the proof certificates to search
for.
