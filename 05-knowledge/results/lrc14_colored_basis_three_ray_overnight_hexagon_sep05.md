# A native live lattice basis, two three-direction circuits, and the LRC14 network closure

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The `observer_collision` referee independently checked the general
colored-basis lemma, including primitive reduction, centered-boundary cases,
dead-color repair, strict convex feasibility, the complete three-direction
classification, the cap-six refinement, and the all-height network envelope.
Root independently checked the basis, classification, and cap refinement.
The `nc2_seed` referee separately passed the analytic network bound and
replayed the complete finite head using carrier and literal sheet engines.
No scarce theorem ID or shared navigation is edited by this note.

There are three conclusions, with different scopes:

1. A centrally symmetric convex body in a rank-two lattice, after deleting
   an additive subgroup, contains a basis made of surviving points whenever
   its surviving points span the plane.
2. For an index-three deleted subgroup, exactly three primitive surviving
   directions have only two circuit types: `(1,1,1)` and `(1,1,2)`.
3. Every primitive sorted distinct positive odd ternary-unit LRC14 triple
   whose **complete** raw support has exactly three primitive directions
   satisfies `E_i<6/77` for **each** of the three degree-zero networks, at
   every height. This uses an analytic tail and a complete finite head.

Together with the proved one- and two-direction notes, a hypothetical failure
of the universal network target must have at least four primitive live
directions. This does not prove entry, synchronization, or LRC(14).

## Inheritance, incoming work, and the actual geometric transfer

The closest proved mechanism is the exact network roof formula in
[THM-4414 — six-separated contact capacity collapse](../../01-canon/theorems/THM-4414-lrc14-six-separated-contact-capacity-collapse.md).
The determinant and complete ordinal-list bounds come from the
[two-direction note](lrc14_two_ray_overnight_hexagon_sep05.md), now routed by
[THM-4428 — two-direction closure and sharp one-direction gap](../../01-canon/theorems/THM-4428-lrc14-two-direction-network-closure-and-sharp-one-direction-gap.md).
The canonical dense hostile is `(19,23,29)`, first isolated by
[THM-4422 — projection deficit and Beatty-row reduction](../../01-canon/theorems/THM-4422-lrc14-projection-deficit-and-beatty-row-reduction.md).
The corrected near miss is discarding all but a chosen ray; the least-used
sidecar is the index-three owner lattice, not just the ambient plane.

Incoming commit `48aaa19246` independently tested all `7,646` exactly-three-
direction rows through height 499, conjectured precisely the two circuit
types, and found the smallest non-additive hostile `(41,47,49)`. We credit
that independent finite signal and hostile in the
[empty-core research board](empty_core_overnight_sep06_board.md); the present
proof settles its stated circuit-exhaustion question. The complete height-499
atlas is not rerun or repackaged here. Reserved incoming claims are not used
as dependencies.

The live concept board is: complete raw support; primitive rays versus raw
multiples; owner subgroup; extremal determinant descent; colored circuits;
the three exact network roofs. The source inspiration is Heule--Scheucher,
[arXiv:2403.00737, Section 4.1](https://arxiv.org/html/2403.00737v1#S4.SS1):
an extremal replacement preserves the native witness and improves an integer
or geometric objective. The map here replaces a pair of owner-live lattice
points by a smaller-index pair while keeping them in the **actual** strict
carrier body. It preserves nonzero owner residue and full support; it loses
the other raw points and roof weights, so those must be restored separately
for the network estimate. No SAT relaxation, abstract order-type realization,
or empty-point-set theorem is asserted to imply an LRC statement.

## 1. A live basis in a centrally symmetric convex body

Let `L` be a rank-two lattice, `H` an additive subgroup of `L`, and `B` a
convex subset of `L tensor R` satisfying `B=-B`. Neither boundedness nor
closedness is required. Set `S=(L\H) intersect B`. If `S` contains two
linearly independent points, then it contains a lattice basis of `L`.

Choose independent `u,v in S` minimizing the positive integer
`d=|det_L(u,v)|`. This minimum exists by well-ordering. Both vectors are
primitive in `L`: if `u=k u0`, `k>=2`, then convexity and `0 in B` put `u0`
in `B`; and `u0 in H` would imply `u in H`. Replacing `u` by `u0` would
therefore lower `d` while keeping a live independent pair.

Suppose `d>1`. A nontrivial coset of `L/(Z u+Z v)` has a representative

```text
x=t u+s v,             |t|,|s|<=1/2.                    (1)
```

Both coefficients are nonzero: if, for example, `s=0`, primitiveness of `u`
forces `t` to be an integer, and the centered bounds force `x=0`, contrary
to the chosen nontrivial coset. Since `|t|+|s|<=1`,

```text
x in conv{u,-u,v,-v} subset B.                          (2)
```

If `x notin H`, use `(u,x)`: its determinant magnitude is `|s|d`, positive
and at most `d/2`, a contradiction. If `x in H` and `|s|<=|t|`, let

```text
y=u-sign(t)x.
```

The coefficient sum relative to `u,v` is
`1-|t|+|s|<=1`, so `y` lies in the same convex hull. Also `y notin H`,
because otherwise `u=y+sign(t)x` would belong to the subgroup. The live pair
`(u,y)` again has determinant magnitude `|s|d` in `(0,d)`. If `|t|<|s|`,
interchange `u,v` and use `v-sign(s)x`. Every case contradicts minimality;
therefore `d=1`.

The cases `|t|=|s|=1/2` cause no problem. If `B` is defined by strict linear
inequalities, every point of the finite convex hull in (2) still satisfies
each inequality strictly. No closure operation or endpoint carrier is added.

### Both structural hypotheses matter

Deleting an arbitrary set instead of a subgroup can leave only
`+/-(1,1), +/-(1,-1)` in the symmetric convex square `(-11/10,11/10)^2`.
Their minimum nonzero determinant is two, not one. Subgroup closure is used
exactly in the dead-color repair.

Central symmetry about zero cannot be dropped: in the convex rectangle
`9/10<x<11/10`, `-1/10<y<21/10`, delete the parity subgroup
`H={(x,y):x-y is even}`. The surviving lattice points are `(1,0),(1,2)`;
their determinant is two. The strongest survivor is the theorem with both
convexity about zero and subgroup deletion retained.

## 2. The LRC owner body is exactly a deleted-subgroup lattice body

Let `w=(a,b,c)` have the hypotheses in conclusion 3, and write

```text
Lambda={C in Z^3 : w.C=0},
L0={C in Lambda : a C1=b C2=c C3 (mod 3)},
H=3 Lambda,
B={C in w-perp : 14|C_i|<3(sum(w)-w_i), all i}.          (3)
```

Primitivity of `w` makes `Lambda` saturated of rank two. Reduction modulo
three identifies `Lambda/3 Lambda` with the two-dimensional kernel of `w`
over `F3`: surjectivity follows by correcting an integer lift with three
times a Bezout vector for `w`. In that kernel, equal weighted coordinates
form a line of three elements. Hence

```text
[Lambda:L0]=3,          [L0:H]=3.                        (4)
```

The two nonzero elements of that line are exactly the vectors whose three
coordinates are all nonzero modulo three. Thus the complete owner-live raw
carrier set of THM-4414 is **exactly**

```text
Omega=(L0\H) intersect B.                              (5)
```

In particular, this is not an arbitrary colored point set. Dividing a live
vector by its positive integer content preserves liveness and the strict
body, so its primitive direction representative belongs to `Omega`.
The lemma therefore gives, for every multi-direction support, primitive live
`u,v` with

```text
u cross v=+/-3w.                                       (6)
```

Indeed a basis of the saturated kernel has cross product `+/-w`, and `L0`
has index three. This is stronger than the older statement that every live
pair has cross-product multiplier divisible by three: at least one pair
attains the smallest possible multiplier.

## 3. Exactly three directions admit only two primitive circuits

The following argument applies in any rank-two `L` with an index-three
subgroup `H`, a centrally symmetric convex body, and exactly three primitive
directions in `(L\H) intersect B`.

Take the live basis `u,v` from Section 1. Choose their signs independently
so that a quotient character `chi:L->F3`, with kernel `H`, satisfies
`chi(u)=chi(v)=1`. Let the third primitive direction be

```text
z=m u+n v,  gcd(|m|,|n|)=1,  mn!=0,  m+n!=0 (mod 3).   (7)
```

We first show that one of `|m|,|n|` equals one. Suppose both are at least two.
If they have the same sign, negate `z` so that `m,n>=2`. Then `u+v` lies in
`conv{u,v,z}` with barycentric coefficients

```text
(n-1)/(m+n-1), (m-1)/(m+n-1), 1/(m+n-1).
```

It is live and on a fourth direction: primitivity rules out `m=n>=2`.
If the signs differ, negate `z` as necessary and write `z=m u-ell v`,
`m,ell>=2`. Again `m!=ell` by primitivity. When `m>=ell+1`, the live point
`2u-v` lies in `conv{u,-v,z}`, with coefficients

```text
2(ell-1)/(m+ell-1), (m-ell-1)/(m+ell-1), 2/(m+ell-1).
```

These are nonnegative and sum to one. It lies on a fourth direction because
`ell>=2` and `gcd(m,ell)=1`. When `ell>=m+1`, use the symmetric live point
`u-2v`. Every case
contradicts the complete three-direction hypothesis.

After interchanging `u,v` and possibly negating `z`, write `n=1`. The live
condition excludes `m=2,-1`. If `m>=3`, the segment from `v` to `z` contains
the fourth live direction `u+v`; if `m<=-3`, it contains the fourth live
direction `v-2u`. Finally `m=0` is not a third direction. Consequently

```text
{directions} = +/-{u,v,u+v},
          or = +/-{u,v,v-2u}.                            (8)
```

Their primitive circuit coefficient multisets are `(1,1,1)` and `(1,1,2)`.
For LRC, the three absolute pairwise cross-product multipliers relative to
`w` are respectively `(3,3,3)` and `(3,3,6)`.

The incoming hostile `(41,47,49)` has complete support

```text
+/-{(11,5,-14),(14,-7,-5),(17,-19,4)}.                   (9)
```

These displayed representatives form an arithmetic progression. It realizes
the second type and refutes the universal additive/A2-only claim. The type
classification forgets raw multiplicities; it does **not** say that there
are only six or eight raw carriers. For example `(5,191,199)` has exactly
three directions and sixteen live carriers. The first failed implication
in the false `N<=8` extrapolation is primitive-ray count implying a bounded
ordinal multiplier count.

### A full cap has at most six carriers, and a three-direction cap is additive

Incoming commit `e10dff3181` proves in the
[full-cap carrier note](overnight_20260906_lrc_cap_carriers.md) that a full
owner dictionary with no three collinear points has at most eight carriers:
there are four parity buckets in each of two owner colors, and two points
in one bucket force a live midpoint. The census finds no eight-point cap.
Central symmetry and primitive reduction explain that absence analytically.

In the abstract index-three setting, a live point `C=2D` would force `D`
live and in the body: subgroup closure proves liveness, and convexity about
zero proves feasibility. The three distinct points `-C,D,C` would be
collinear. Thus a full cap contains **no zero parity class** in `L/2L`.
Each owner color can use at most one of each of the three remaining parity
classes, by the same midpoint argument. Therefore

```text
N<=2(2^2-1)=6.                                         (9a)
```

The argument is exact and all-height, not a conclusion from the missing
eight-point rows. A noncollinear full cap therefore has two or three
primitive directions, each contributing precisely its positive and negative
primitive representatives. In the three-direction case, the AP circuit in
(8) is forbidden: `v,u,-(v-2u)` are three distinct collinear live points.
Consequently every three-direction full cap has the additive/A2 circuit.
The example `(23,29,37)` from the incoming cap note realizes six points,
so (9a) is sharp. This does not replace that note's sharper network constant
`204/5957` or its equality analysis; it strengthens its structural bound and
embeds its entire noncollinear class in the two/three-direction closure.

More generally, let `L` be rank `d`, let `H` have index three, let `B=-B`
be convex, and let the **complete** set `(L\H) intersect B` have no three
collinear points. The identical zero-parity exclusion and same-color
midpoint argument give `N<=2(2^d-1)`. No assertion of sharpness in higher
dimension is made. Arbitrary selected subsets, a nonsymmetric body, or an
arbitrary list of admitted cosets are outside this corollary's hypotheses.

## 4. An all-height three-direction envelope for every projection

For a primitive live direction `d`, let `M_d=max_i|d_i|`. Every direction
in a multi-direction support has `M_d>=7`: the short-relation zero-defect
lemma in
[THM-4386 — canonical component relation, Lemma 2.1](../../01-canon/theorems/THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence.md)
would otherwise force every carrier parallel to a relation of norm at most
fourteen. The parity and ternary-unit details are proved in the two-ray note.
For independent directions `d,e`, the owner residue gives

```text
M_d M_e >= P,             P=3c/2.                       (10)
```

Let the three maxima be `M1,M2,M3`, and `m=min Mj`. If `m>=sqrt(P)`, then
`sum 1/Mj<=3/sqrt(P)`. Otherwise the other two maxima are at least `P/m`,
so

```text
sum 1/Mj <=1/m+2m/P,           7<=m<sqrt(P).
```

The right side is convex in `m`, so its maximum on `[7,sqrt(P)]` is at an
endpoint. (If `sqrt(P)<7`, only the first case is possible.) Therefore

```text
sum_(j=1)^3 1/Mj <= max(3/sqrt(P), 1/7+14/P).           (11)
```

Crucially, keep the complete positive and negative multiplier lists:

```text
B_d=min_i 3(sum(w)-w_i)/(14|d_i|),
K_d=max{k in Z : k<B_d},
N_d=2(K_d-floor(K_d/3)).                                (12)
```

The surviving multipliers are exactly `+/-k` with `1<=k<B_d` and `3` not
dividing `k`. Since `B_d<3c/(7M_d)` and `N_d<4B_d/3+4/3`, the cap
`q=3/(7c)` in each THM-4414 summand gives, separately for each `i`,

```text
E_i < (12/49) sum_(j=1)^3 1/Mj +12/(7c)
    <= max( (36/49)sqrt(2/(3c))+12/(7c),
            12/343+4/c ).                              (13)
```

This bound uses all three directions and all actual multiples; neither a
physical-mass inequality nor the minimum of the three networks replaces
the claimed conclusion. The circuit classification is additional structure,
not a logically necessary input to this particular reciprocal bound.

Both branches in (13) decrease with `c` and are below `6/77` at `c=99`.
For the radical branch, `6/77-12/(7*99)=2/33>0`, and the squared radical is
`96/26411<4/1089`. For the rational branch,

```text
12/343+4/99=2560/33957 <2646/33957=6/77.                 (14)
```

Thus (13) proves the claim for all `c>=99`. Height 99 itself is excluded by
the ternary-unit condition; the first eligible tail height is 101. The only
remaining proof obligation is the complete finite head `c<99`.

## 5. Exact finite head, abstract hostile tests, and reproduction

The finite universe is all `5,409` primitive sorted distinct positive odd
ternary-unit triples with `c<99`, without a shell, rank, or count prefilter.
Exactly `3,500` are multi-direction; each admits the basis (6). Exactly
`1,791` have three directions: `1,107` additive and `684` AP-circuit type.
Every projection in this complete three-direction head is strictly below
`6/77`; the head maximum across rows and coordinates is

```text
w=(5,37,43),
(E1,E2,E3)=(240/11137, 18/301, 2822/55685).              (15)
```

The value `18/301` is a **finite-head maximum**, not a claimed sharp global
constant. There are no equality cases at `6/77` in this three-direction
family, because both the finite head and analytic tail are strict.

The new
[exact script](../../04-computation/lrc14_colored_basis_three_ray_overnight_hexagon_sep05.py)
transparently reuses the previous congruence-row engine and a separately
implemented literal shifted-sheet network engine. It checks complete raw
dictionaries before direction filtering at wide controls, and compares all
three raw projection sums against physical sheet networks for every selected
head row. Proof tests use explicit exceptions, not optimized-away assertions.

The independent abstract lattice tests use centrally symmetric intersections
of three open integer strips in `Z^2`, all nonzero linear characters modulo
`2,3,4,5`, and the noncyclic subgroup `2Z^2`. There are `5,115` spanning
body/subgroup tests, including `2,526` live-centered-representative descents
and `810` dead-representative repairs; both three-direction forms are tested
independently of the LRC specialization. All `242` noncollinear index-three
cap controls have at most six points and no zero-parity point. Explicit
hostile tests retain both
the failed subgroup and failed central-symmetry relaxations from Section 1.

Reproduce from the repository root:

```bash
python3 -B 04-computation/lrc14_colored_basis_three_ray_overnight_hexagon_sep05.py
python3 -B -O 04-computation/lrc14_colored_basis_three_ray_overnight_hexagon_sep05.py
```

The saved [output](lrc14_colored_basis_three_ray_overnight_hexagon_sep05.out)
records exact check counts, the two circuit types, finite-head fractions,
and controls `(19,23,29)`, `(41,47,49)`, `(5,191,199)`, `(7,611,613)`.
The last has thirteen directions and is a live-basis control only; it is not
claimed as a three-direction instance or as proof of arbitrary-support
network closure. A semantic digest covers complete selected carrier sets,
circuit types, and all three rational projection sums.

The normal and optimized transcripts byte-match. There are `50,286` primary
checks, `37,194` reused row-engine checks, and `1,794` literal-sheet checks.
Frozen raw-LF SHA-256 values:

```text
source 9b33a16dc69ee61e4bc260c6f7e7ccc77a79b9d39231c8369979546e0d400700
output 1f1434c644397758ee4ee9a7f5c5f7d9999a56846126b29f55f6659f27069f4c
```

## 6. What remains missing for arbitrary support

The native basis lemma removes a genuine geometric obstruction: one can
legally perform extremal replacement while retaining the owner gate and
strict carrier body. It does not license deleting the rest of the support
from a weighted network sum. The next object is the residue-deleted convex
body with all three projection roofs, not a cosmetic tournament on rays.

For each projection, layer cake shrinks one of the three strips while the
other two stay fixed. Consequently its roof can remain positive at a
transverse boundary, unlike the physical three-facet minimum. A prospective
continuous-area estimate needs an explicit owner-coset boundary discrepancy;
it cannot silently transfer a physical tent estimate to a network roof.
Four or more directions, universal comparable-region closure, entry,
synchronization, and LRC(14) remain **OPEN** here.
