---
id: THM-2984
title: "Projected k3 signed-ray attainment and primitive-unit phase gate"
status: >
  PROVED + FINITE-EXACT + HOSTILE-AUDITED + LEAN-CORE-FORMALIZED.  On an
  exact periodic label ray with contribution A/z, existence
  of a scalar-admissible height is decided by a three-way sign/attainment
  law.  Retaining the primitive unit also makes every fixed-cell phase
  independent of height, giving a finite unit-by-unit strict-open exclusion
  gate.  Its centered cardinality shadow is sharp: more than
  beta(d)=2 floor((d-1)/14)+1 fixed-safe residues clear every unit, while
  beta(d) residues need not.  An arbitrarily translated open danger band has
  the different sharp capacity ceil(d/7), and its continuous center compresses
  exactly to a finite affine interval-orbit complex.  The complement of that
  complex's one-skeleton is exactly the short-order graph.  The affine
  complex is flag exactly for d in {2,3,4,5,6,7,8,10,12,14,15}; in particular
  its first failure is d=9 and every d>=16 fails.  The centered and translated
  bounds must not be conflated.  Gcd-stratum capacities and the exact
  multiplicative transporter refine the centered threshold whenever residue
  shape is retained.  The resulting centered transport complex is flag
  through d=42 and first becomes non-flag at d=43, with an explicit
  irreducible three-cell certificate.  This is a reusable refinement of the
  projected k=3 denominator quotient; it does not assert that any new atlas
  row is empty and does not improve the current proved cap by itself.
source: codex-lrc14-k3-signed-ray-phase-gate-2026-07-30
audit: >
  The exact referee uses explicit require gates, including regression guards
  for both the 499 formula rows and the six finite flag-gap rows.  Ordinary
  and optimized Python are byte-identical to the frozen transcript.  The
  dedicated Lean module elaborates directly and its displayed axiom audit has
  only propext, Classical.choice, and Quot.sound.  This does not claim that the
  aggregate project root was rebuilt in the same audit.
depends_on:
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-2979-projected-k3-z275-ten-body-status-and-located-torsion-closure
related:
  - THM-2981-projected-k3-z270-to-z247-cardinality-torsion-descent
  - THM-2992-signed-quartic-edge-block-discriminant-parity-and-keller-owner-line-boundary
  - HYP-2935-lrc14-bigraded-relation-signature
  - HYP-3003-summand-multiplicand-farey-basis-merge
  - MISTAKE-340
script: 04-computation/lrc14_projected_k3_signed_ray_gate_thm2984.py
output: 05-knowledge/results/lrc14_projected_k3_signed_ray_gate_thm2984.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCSignedRayGate.lean
script_sha256: 23be12885b44efb1e23acfb6a58b9dd285ef57cba666c60f9235bde05fc27af5
output_sha256: 9e33756dd352042ba6a4571a0af1af626e70e0991e94624c0dfe18711a229dbc
semantic_sha256: 756ae7ed291fc438648ed76f57b07eae1cde4a1879b184f9564b6ca8adc94f48
formalization_sha256: 6f3ab284fc1936ebc0b565bf812fcdc38e6940d17e93d607b802f80d6261bc91
hash_basis: raw file bytes
---

# THM-2984 -- signed-ray attainment and primitive-unit phase gate

**PROVED + FINITE-EXACT + HOSTILE-AUDITED + LEAN-CORE-FORMALIZED.**

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

### 4.1 The sharp translated-band capacity

The exact centered count `beta(d)` must not be used after an arbitrary local
phase translation.  For a real phase `theta` and a primitive unit `u mod d`,
put

```text
B_d(theta,u)={r mod d:||theta+ur/d||<1/14}.             (13a)
```

Then the sharp uniform capacity is

```text
|B_d(theta,u)|<=kappa(d):=ceil(d/7).                    (13b)
```

Indeed multiplication by `u` permutes the residues.  The strict danger arc
has length `1/7` on the unit circle, hence length `d/7` after scaling the
residue lattice to unit spacing.  If there are no occupied residues the bound
is trivial.  Otherwise cut the circle outside this arc and lift its `k`
occupied lattice points to integers

```text
n_1<...<n_k
```

in one open real interval of length `d/7`.  Unit spacing and strict openness
give

```text
k-1<=n_k-n_1<d/7,
```

so the integer `k` is at most `ceil(d/7)`.  Conversely a suitably shifted
open interval of length `d/7` contains exactly `ceil(d/7)` lattice points;
thus `(13b)` is sharp among translation-uniform count bounds.

Consequently

```text
|S|>ceil(d/7)
implies that for every primitive u and every theta
some c in S lies outside B_d(theta,u).                  (13c)
```

The distinction from `(12)` is real.  At `d=28`,

```text
beta(28)=3,                 kappa(28)=4,
u=1, theta=-3/56:           B_28(theta,u)={0,1,2,3}.

Equivalently, after scaling by 28, all four lattice points lie in the open
interval (-1/2,7/2) of length 4.                        (13d)
```

Thus four residues can all lie in one translated open danger band although
no centered band contains more than three.  The `beta(d)` gate is exact for
the absolute phase `(8)`; `(13b)--(13c)` is the required gate when aligned
coordinates leave an arbitrary translation `theta`.  This scope guard is
especially important in projected terminal arguments: a local shift cannot
be silently recentered unless the fixed-cell sidecar proves that operation.

### 4.2 The affine transport complex removes the continuous phase

The translated gate has an exact finite shape-sensitive form.  Put
`kappa=ceil(d/7)` and let

```text
C_a={a,a+1,...,a+kappa-1} mod d.                        (13e)
```

Since

```text
kappa-1<d/7,                                            (13f)
```

every `C_a` lies in some open circular interval of length `d/7`.  Conversely,
if lattice residues lie in such an interval, cut and lift as in `(13b)`.
Their integer span is at most `kappa-1`, so they extend to a consecutive
`kappa`-block.  Therefore the translated obstruction is the simplicial
complex

```text
K_d^tr={S:there exist u in (Z/dZ)^* and a in Z/dZ
             with uS subset C_a}.                       (13g)
```

Its maximal faces are precisely the affine images `u^{-1}C_a`, all of size
`kappa`.  Equivalently the continuous parameter `theta` in `(13a)` is
compressed without loss to the finite affine transporter

```text
T_d^tr(S)={(u,a):uS subset C_a}.                         (13h)
```

Thus `(13c)` is exactly the facet-cardinality gate for the orbit of one short
cyclic interval under the affine group `(Z/dZ) semidirect (Z/dZ)^*`.  At or
below `kappa`, the correct next search is membership in this affine orbit,
or small minimal nonfaces of `K_d^tr`, rather than another count.

The centered complex studied below is a subcomplex of `K_d^tr`: its symmetric
band `B_d` is one cyclic consecutive block of length `beta(d)<=kappa`, possibly
properly contained in a facet.  The inclusion can be strict, as `(13d)`
shows.  Translation also destroys the gcd strata anchored at zero.  Hence the
gcd refinements in Section 5 apply to the centered transporter `(17)`, not
automatically to `(13h)`; an affine application must retain the translation
coordinate or replace gcd by a genuinely affine invariant.
The first such survivor is the multiset of pair-difference strata

```text
{gcd(c-c',d):c,c' in S},                               (13i)
```

equivalently the additive orders `d/gcd(c-c',d)`.  Translation cancels in
each difference and multiplication by a unit preserves its gcd with `d`.
This is why THM-2979's short-order pair gate survives arbitrary local
translation even though the absolute point-stratum capacities `(16)` do not.

In fact the pair relation is exact.  For distinct `x,y mod d`, put
`g=gcd(x-y,d)`.  Then

```text
{x,y} in K_d^tr
  iff g<=kappa-1
  iff the additive order of x-y is greater than 7.       (13j)
```

For the forward direction, two points in one `kappa`-block have a nonzero
integer difference of absolute value at most `kappa-1`; unit multiplication
preserves the gcd with `d`.  Conversely, reduction of units
`(Z/dZ)^* -> (Z/(d/g)Z)^*` is onto (prime-power lifting and CRT), so some unit
`u` has `u(y-x)=g mod d`.  If `g<=kappa-1`, the two images lie together in
the block beginning at `ux`.  Finally

```text
g<=ceil(d/7)-1  iff  g<d/7  iff  d/g>7,
```

which proves `(13j)`.  Therefore the two-element nonfaces of `K_d^tr` are
**exactly** the edges of THM-2979's short-order Cayley graph:

```text
{two-element nonfaces of K_d^tr}=E(G(d,7)).              (13k)
```

This also gives the first affine flag failure without enumeration.  For
`2<=d<=7`, `kappa=1`, so `K_d^tr` is discrete and flag.  For `d=8`, its
two-point facets have unit difference; every such difference is odd, so the
one-skeleton is bipartite by parity and has no triangle.  At `d=9`, however,
`kappa=2` and every pair in `{0,1,2}` has unit difference.  The triple is a
three-clique but cannot be a face of a complex whose facets have size two.
Thus

```text
K_d^tr is flag for 2<=d<=8, while K_9^tr is not flag.    (13l)
```

There is also an exact gap formula which removes the translation coordinate
before any search over `a`.  For nonempty `S` and a fixed unit `u`, sort the
distinct residues of `uS` cyclically and let `Delta_1,...,Delta_s` be their
positive forward gaps, summing to `d`; for a singleton take `Delta_1=d`.
Then

```text
#{a:uS subset C_a}
  =sum_i max(0,Delta_i-(d-kappa)).                       (13m)
```

Indeed the complement of `C_a` is a consecutive block of `d-kappa`
residues.  It avoids `uS` precisely when it lies among the `Delta_i-1`
unselected residues inside one cyclic gap, and that gap admits
`max(0,(Delta_i-1)-(d-kappa)+1)` placements.  In particular,

```text
there exists a with uS subset C_a
  iff max_i Delta_i>=d-kappa+1.                         (13n)
```

Thus the continuous-center obstruction is equivalently a unit-indexed
maximum-gap profile.  Below the cardinality gate, a terminal computation need
not enumerate centers: it records the cyclic gaps of each `uS`; the affine
transporter is empty exactly when every unit has maximum gap at most
`d-kappa`.

### 4.3 Additive compositions under diagonal unit multiplication

The whole complex has a second exact normal form which exposes the operation
split hidden by `(13g)`.  Let `S` have `s>=1` elements and let

```text
Delta_(s-1)(kappa)
 ={(a_1,...,a_(s-1)): every a_i>=1 and sum_i a_i<=kappa-1},   (13o)
```

with the unique empty composition when `s=1`.  Then `S` is a face of
`K_d^tr` if and only if there are an ordering `(x_1,...,x_s)` of `S`, a unit
`u mod d`, and `a in Delta_(s-1)(kappa)` such that

```text
u(x_(i+1)-x_i)=a_i mod d,       1<=i<s.                (13p)
```

If `uS subset C_A`, order the distinct images by their offsets
`0<=t_1<...<t_s<=kappa-1` and take `a_i=t_(i+1)-t_i`.  Conversely, telescope
`(13p)`: the images are `ux_1,ux_1+a_1,...,ux_1+sum_i a_i`, all in the block
beginning at `ux_1`.  This proves the equivalence, including all wraparound
cases because `kappa<d`.

Thus addition supplies a short positive-composition simplex, while
multiplication by a unit tests the diagonal orbit of the ordered difference
chain.  Translation has disappeared.  In this precise finite setting, this
realizes the operation split proposed in HYP-2935 and HYP-3003; it does not
promote their broader synthesis claims.

At a prime modulus `p` and for `s>=2`, every adjacent difference is
invertible.  For a fixed ordering, projectivize its difference chain:

```text
[x_2-x_1:...:x_s-x_(s-1)] in P^(s-2)(F_p).             (13q)
```

The ordering contributes a face exactly when this point lies in the
projective image of `Delta_(s-1)(kappa)`: equality of projective points
supplies one common nonzero scalar, hence the single unit required in
`(13p)`, rather than separate coordinatewise units.

Distinctness has also become a standard arrangement complement.  Every
consecutive subsum of the difference chain is nonzero, since

```text
sum_(k=i)^j (x_(k+1)-x_k)=x_(j+1)-x_i != 0.            (13q1)
```

Thus ordered `s`-point configurations modulo affine translation and unit
scaling lie in the projectivized complement of the type-`A_(s-1)` braid
hyperplanes; reordering acts by `S_s`.  The face predicate asks whether this
reordering orbit meets the projective short-composition image.

For a triple this is the finite short-summand ratio set

```text
R_(p,kappa)={b/a mod p:a,b>=1 and a+b<=kappa-1};        (13r)
```

one tests the six orderings of the triple.  If `r=(z-y)/(y-x)`, their ratio
orbit is

```text
{r,1/r,-1-r,-1/(1+r),-r/(1+r),-(1+r)/r},              (13s)
```

with repetitions at special orbits.  The set `(13r)` is the modular image
of a Farey-like positive triangle; injectivity of distinct short fractions
is not asserted.  This is the first genuinely higher-order sidecar after the
pair graph `(13k)`, and gives a small projective lookup table for triple
certificates.

The six transformations in `(13s)` are exactly an anharmonic `S_3` action,
not an accidental list.  It is generated by

```text
s(r)=1/r,                       c(r)=-1/(1+r),
s^2=c^3=1,                      scs=c^(-1).             (13s1)
```

Indeed `c^2(r)=-(1+r)/r`, and applying `s` to `r,c(r),c^2(r)` supplies the
other three entries of `(13s)`.  This is the same abstract quotient
`S_4/V_4~=S_3` that acts on the three quartic matchings in THM-2992, but the
shared group must not be mistaken for a cross-problem proof map.  Here it
records information destroyed by reordering a translated-and-scaled
three-point residue configuration; there it records information destroyed by
passing from four labelled sheets to their three unordered matchings.  No
loneliness predicate, quartic inertia, or owner sheet is transported between
the two objects without an additional sidecar.

At composite `d`, a difference may be a zero divisor and this
normalization is unlawful; the diagonal unit orbit together with the
difference-gcd strata `(13i)` is the correct nonuniform replacement.

### 4.4 Complete flag classification of the affine complex

The pair graph, the exact facet size, and the gap normal form together give a
complete classification.  With `flag` meaning that every clique in the
one-skeleton is a face,

```text
K_d^tr is flag
  iff d is in {2,3,4,5,6,7,8,10,12,14,15}.             (13t)
```

Equivalently, the affine complex is non-flag at `d=9,11,13` and at every
`d>=16`.  This is a uniform statement, not a finite-range census.

Let

```text
R=max{r:r divides d and 1<=r<=7}.
```

By THM-2979, the independence number of the short-order graph is
`alpha(G(d,7))=d/R`.  By `(13k)`, every independent set is a clique in the
one-skeleton of `K_d^tr`, while every face has at most
`kappa=ceil(d/7)` vertices.  Therefore

```text
d/R>kappa  implies K_d^tr is non-flag.                 (13u)
```

The reverse weak inequality `kappa<=d/R` is `(23a)`, so only equality needs
further analysis.  Write `d=Rq`.  Then

```text
q=ceil(Rq/7)  iff  (7-R)q<7.
```

Keeping the defining condition that `R` is the largest divisor of `d` at
most seven gives the exact equality list.  Explicitly, `R=1` yields no
`d>=2`; `R=2,3` force `q=1`; and `R=4,5,6` allow respectively
`q<=2,3,6`, after which the largest-divisor condition removes no further
listed value.  Thus

```text
d/R=kappa
iff d is a multiple of 7, or
    d in {2,3,4,5,6,8,10,12,15,18,24,30,36}.           (13v)
```

For `d<=7`, `kappa=1` and the complex is discrete.  For
`d=8,10,12,14`, `kappa=2`; compatibility in `(13j)` means that the
difference is a unit, hence odd.  The one-skeleton is bipartite by parity and
has no triangle, so its edge complex is flag.

At `d=15`, compatibility again means unit difference.  A clique has at most
three vertices, since its residues modulo three must be distinct.  Normalize
any triangle by translation and unit multiplication to `{0,1,t}`.  Requiring
both `t` and `t-1` to be units modulo `15` leaves exactly

```text
t in {2,8,14}.
```

Each resulting triple is an affine image of `{0,1,2}`: `t=2` is literal,
multiplication by `2` sends the `t=8` triple to `{0,1,2}`, and the `t=14`
triple is the cyclic block `{14,0,1}`.  Hence every clique is
a face and `K_15^tr` is flag.  In contrast, `{0,1,2}` is already a three-clique
nonface at each of `d=9,11,13`, since `kappa=2` there.  These arguments settle
all equality cases below `16` and also display the first three failures.

Four isolated equality moduli and two exceptional multiples of seven have
explicit three-clique nonfaces.  In the following table, the listed units are
the positive representatives of `U(d)/{+-1}`.  An entry records the positive
cyclic gap word of `uS` and its maximum; replacing `u` by `-u` reverses the
word and preserves the maximum.  Every pair difference in `S` has additive
order greater than seven, so `(13j)` makes `S` a clique.  Every displayed
maximum is strictly below the face threshold `d-kappa+1` in `(13n)`, so it is
not a face.

| `d` | `kappa` | clique `S` | threshold | exact `+/-u` gap rows |
|---:|---:|---|---:|---|
| 18 | 3 | `{0,1,5}` | 16 | `1:(1,4,13;13)`, `5:(5,2,11;11)`, `7:(7,10,1;10)` |
| 21 | 3 | `{0,1,5}` | 19 | `1:(1,4,16;16)`, `2:(2,8,11;11)`, `4:(4,16,1;16)`, `5:(4,1,16;16)`, `8:(8,11,2;11)`, `10:(8,2,11;11)` |
| 24 | 4 | `{0,1,10}` | 21 | `1:(1,9,14;14)`, `5:(2,3,19;19)`, `7:(7,15,2;15)`, `11:(11,3,10;11)` |
| 30 | 5 | `{0,1,8}` | 26 | `1:(1,7,22;22)`, `7:(7,19,4;19)`, `11:(11,17,2;17)`, `13:(13,1,16;16)` |
| 36 | 6 | `{0,1,11}` | 31 | `1:(1,10,25;25)`, `5:(5,14,17;17)`, `7:(5,2,29;29)`, `11:(11,2,23;23)`, `13:(13,22,1;22)`, `17:(7,10,19;19)` |
| 42 | 6 | `{0,1,10}` | 37 | `1:(1,9,32;32)`, `5:(5,3,34;34)`, `11:(11,15,16;16)`, `13:(4,9,29;29)`, `17:(2,15,25;25)`, `19:(19,3,20;20)` |

It remains to treat the equality family `d=7m`.  The cases `m=1,2` are the
flag moduli `7,14`; the cases `m=3,6` are the rows `21,42` above.  For every
`m>=4`, `m!=6`, put

```text
S_m=({0,1,...,m-1} minus {1}) union {m+1}.              (13w)
```

This set has `m=kappa` vertices and is a one-skeleton clique.  Differences
among its old interval points have absolute representatives at most `m-1`
and therefore additive order greater than seven.  The differences from its
new point are `{2,3,...,m-1,m+1}`.  Only the last requires a separate check,
and

```text
ord_(7m)(m+1)=7m/gcd(m+1,7m)
              =7m if 7 does not divide m+1, and m otherwise.           (13x)
```

This is greater than seven in the stated range; `m=6` is exactly the failed
small construction.

Finally `S_m` cannot be a facet.  Let `N_q(S)` count unordered pairs at
circular difference `+/-q`.  Direct inspection of `(13w)` gives

```text
N_1(S_m)=m-3,                 N_2(S_m)=m-2.             (13y)
```

If `S_m` were a facet, it would be a unit-step progression
`{A+jv:0<=j<m}` for a unit `v`.  An index separation `h` contributes `m-h`
pairs whenever `hv=+/-q`.  For fixed `q` there is at most one such
`1<=h<m`: equal signs would force `7m` to divide `h-h'`, while opposite signs
would force it to divide `h+h'`, impossible because `h+h'<=2m-2<7m`.
Thus `(13y)` forces

```text
3v=epsilon_1,                 2v=2 epsilon_2  (mod 7m),
epsilon_1,epsilon_2 in {+1,-1}.                         (13z)
```

After multiplying these equations by `2` and `3`, respectively, `7m` would
divide one of `+/-4,+/-8`, impossible for `m>=4`.  Hence `(13w)` is a clique
nonface.  Together `(13u)--(13z)` prove `(13t)`.

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

The conclusion is translation-robust.  If both
`theta+uc_1/d` and `theta+uc_2/d` lay in the translated band `(13a)`, their
difference would again have circular distance strictly below `1/7`; the
common `theta` cancels, contradicting `(21)` after division by `d`.  Hence
every facet `u^{-1}C_a` of `K_d^tr` is an independent set in THM-2979's
short-order Cayley graph `G(d,7)`.  The unit-lifting converse in `(13j)`
shows that this is exact, not merely a necessary condition:

```text
E(G(d,7))={two-element nonfaces of K_d^tr}.              (21a)
```

This is the precise affine survival mechanism behind the difference strata
`(13i)`.

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

For the translated gate `(13c)`, the correct comparison is instead

```text
kappa(d)=ceil(d/7)<=d/R.                               (23a)
```

Thus the same pair-density count `|S|>d/R` also implies universal escape from
every translated band, without selecting a pair.  What cannot be done is to
replace `kappa(d)` by the smaller centered value `beta(d)` in that translated
argument; `(13d)` is a canonical boundary mechanism to test whenever such a
recentring is proposed.

### 6.1 Complete-cell projection corollary

The translated-band gate becomes the exact missing sidecar when an argument
keeps whole complete cells rather than only their distinguished endpoints.
Let `C` be a set of complete-cell addresses throughout which the body and all
fixed drift labels are safe, put `S=C mod d`, and choose once a representative
`c_r in C` for each `r in S`.  Write a high label as

```text
z=(L/d)u+mL,                 m in Z,
```

and parametrize the complete cell with address `c` by `x=(c+y)/L`,
`0<=y<1`.  Then

```text
zx == uc/d+(u/d+m)y  (mod 1).                           (23b)
```

Apply `(13b)` with `theta=(u/d+m)y`.  At every fixed local coordinate `y`, at
most `ceil(d/7)` selected representatives are high-danger.  The pointwise
high-safe count

```text
M_z(y)=#{r in S:x_r=(c_r+y)/L is high-safe}
```

satisfies

```text
M_z(y)>=|S|-ceil(d/7).                                  (23c)
```

In particular the inherited hypothesis `|S|>d/R>=ceil(d/7)` gives
`M_z(y)>=1` for every `y`.  For each `y`, choose one such `r`.  The point
`x_r` is safe for the body, every fixed drift, and the high drift, while
`phi_L(x_r)=y`.  In the notation of THM-2941 `(25g)`, this proves

```text
P_(E,Z)=T.                                               (23d)
```

An aligned completion would instead force `P_(E,Z) subset U_A`.  This is
impossible whenever `mu(U_A)<1`; in the three-comb applications the proved
THM-1166 cap is `mu(U_A)<=36/91`.  This is the pair-free complete-cell
implication used by the THM-2981 descent, without choosing a
collision pair or a pair-selection policy.  Notice that `(23c)` is a
pointwise section argument, not a claim that the raw Haar mass of one
complete cell is one.

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

This is the first non-flag modulus for the centered complex `K_d`.  For
`2<=d<=42`, the band radius
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

This threshold is centered only.  It does not transfer to the affine
translated complex: for the same `d=43` witness,

```text
3*{1,2,15}={2,3,6} subset C_0={0,1,...,6},
```

so `{1,2,15}` is a face of `K_43^tr`.  No first-nonflag threshold for
`K_d^tr` follows from the centered witness: its separate complete
classification is `(13t)`, with first failure already at `d=9`.

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

## 8. Reproducibility and formalization audit

The exact standard-library referee is
`04-computation/lrc14_projected_k3_signed_ray_gate_thm2984.py`; its frozen
transcript is
`05-knowledge/results/lrc14_projected_k3_signed_ray_gate_thm2984.out`.
It is a hostile finite audit of the uniform proofs above, not a proof
dependency.  Among its controls are all subsets through `d=15`, every
cardinality-equality split through `d=5000`, all `58` signed units in the six
sporadic affine-gap tables, and the symbolic multiple-seven obstruction for
every `4<=m<=1000`, `m!=6`.  Every load-bearing check uses explicit
`require`; Python optimization cannot erase the audit.  From the repository
root, ordinary and optimized runs are byte-identical to the stored output:

```bash
PYTHONHASHSEED=0 python3 04-computation/lrc14_projected_k3_signed_ray_gate_thm2984.py
PYTHONHASHSEED=0 python3 -O 04-computation/lrc14_projected_k3_signed_ray_gate_thm2984.py
```

The root-imported Lean module
`04-computation/lean/TournamentH7/TournamentH7/LRCSignedRayGate.lean`
formalizes the exact sign/attainment trichotomy, primitive-unit height-free
phase, strict-open centered band cardinality and sharp count gate, generic
translated lattice-capacity arithmetic, and the equality-split/pair-count
arithmetic used in the affine flag classification.  It does not claim a Lean
formalization of every transporter, gap-word, or composition-orbit argument
in this document.  Direct elaboration and the dedicated Lake target pass, and
the module is imported by the project root.  This audit does not claim a fresh
aggregate `lake build TournamentH7`.  Its previously exposed
`LRCCoherentBlockerChronology` obstruction has been repaired separately by
restoring the intended integer types of its tooth-address variables; see
MISTAKE-340.  Every displayed `#print axioms` report for `LRCSignedRayGate` is
contained in
`{propext, Classical.choice, Quot.sound}`; there is no `sorryAx`, theorem-
specific axiom, or `native_decide` dependency.

Frozen SHA-256 hashes are

```text
source    23be12885b44efb1e23acfb6a58b9dd285ef57cba666c60f9235bde05fc27af5
output    9e33756dd352042ba6a4571a0af1af626e70e0991e94624c0dfe18711a229dbc
semantic  756ae7ed291fc438648ed76f57b07eae1cde4a1879b184f9564b6ca8adc94f48
Lean      6f3ab284fc1936ebc0b565bf812fcdc38e6940d17e93d607b802f80d6261bc91
```
