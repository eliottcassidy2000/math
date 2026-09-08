# Independent audit of the infinity-eleven two-cusp braid bundle

**Status: independent analytic/source audit PASS; normal, optimized, and
frozen-output replay PASS. Audit frozen.** This accepts the proposed
literal complement computation and its stated inherited whole-class
consumer. Primary status promotion remains with the producer/root.
No source, witness, older family note, or Git state was changed.

The primary is
[planar_jc48_sep08_infinity11_braid.md](planar_jc48_sep08_infinity11_braid.md).
I read its complete final proof, the full new source, the five raw
rational witnesses, and the proved
[two-cusp family transport](planar_jc48_sep08_two_cusp_family.md)
with its [independent audit](planar_jc48_sep08_two_cusp_family_audit.md).
The new certificate uses the already independently accepted
[moving-tube method](planar_jc48_sep08_moving_tubes_audit.md), but the
present acceptance includes the actual new polynomial, paths, and group
consumer; method acceptance alone is not used as their proof.

## 1. The actual polynomial and the good infinity-eleven point

The literal parametrization is

    U=t^4-(8/3)t^3-2t^2+8t,
    V=t^6-4t^5-(3/2)t^4+(548/27)t^3-(368/9)t.

The source independently reconstructs every coefficient of
`F(u,v)=Res_t(U-u,V-v)`, verifies `F(U,V)=0`, and verifies the complete
monic degree-four condition in `v`. There is no implicit approximate
quartic substituted for this resultant. The full nonzero discriminant
is the one displayed in the primary, with three distinct cubic-node
values and the separate values `13/3,8/3,-19/3`.

There is also an independent way to check the finite inventory from
this exact discriminant. My direct computation gives

    disc_t(U-u)=-256(3u-13)(3u-8)(3u+19)/27.

At a root of the cubic `H(u)` in the primary, all four roots of
`U(t)=u` are therefore simple analytic functions of `u`. In that
neighborhood the monic resultant is the product of
`v-V(t_i(u))`. A discriminant zero of order two forces exactly one
pair difference to vanish, to order exactly one. No triple can
occur, and the two corresponding graph slopes are distinct. Thus
each of the three simple roots of `H` gives one ordinary node.

The source checks that the only common zeros of `U'` and `V'` are
`t=1,-1`, that each cusp image has only its declared double parameter
preimage, and that both ordinary-cusp determinants are nonzero.
Their discriminant orders three account for the entire discriminant
at their projections. The remaining simple projection fold at
`t=2` has discriminant order one and is a smooth point of the image.
There is no additional finite singularity at a value with nonzero
vertical discriminant.

The monic equation for `t` makes the parametrization finite onto its
image. Its image is irreducible. The nonzero vertical discriminant
excludes a repeated generic image of two parameter roots, hence the
parametrization has degree one onto its image. This recovers the
actual affine normalization without assuming it from a passport.

At infinity put `w=1/t`, `X=U/V`, and `Z=1/V`. I independently
expanded the literal rational functions through order eleven and
obtained

    Z-X^3+(7/3)X^4-(2237/108)X^5
        =-(800/81)w^11+O(w^12).

All earlier terms in this expression vanish, while `ord(X)=2` and
the coefficient shown is nonzero. This is the actual infinity type
`(2,11)`, with one infinity preimage. It agrees with the already
proved family formulas at `s=2,a=-10/3,c=184/27`. Thus the literal
curve is indeed the specified good point of `G11`, with two ordinary
cusps and three ordinary nodes.

## 2. The complete new moving certificates

The dense Taylor arrays include every base degree through six and
root degree through four. In the affine substitution, an entry
`b_jk` contributes

    b_jk du^j binom(k,l) dz^(k-l)

to the coefficient of `h^(j+k-l) w^l`. This is the same complete
expansion proved in the moving-method audit. The new source also
compares it with a direct symbolic affine substitution for this
new literal polynomial.

The Rouché comparison uses the lower norm of the constant linear
coefficient times the radius, against the sum of upper norms of
**every** other coefficient, including all varying linear terms.
The strict positive margin pays the whole real parameter interval.
The exact affine L-infinity distance examines every endpoint and
all coordinate/equal-absolute-coordinate breakpoints. Its positivity,
with the explicitly tested radius-sum inequalities, makes the four
moving disks disjoint. Four one-root disks exhaust the verified
monic quartic, so the actual roots are all simple throughout every
accepted path.

At adjacent endpoints, the disks of a given label have the identical
centre and are nested. Their unique roots coincide. The initial
fixed tiny disks supply an additional root isolation: their radii
are the initial pair lower distances divided by `32*64`. They are
pairwise disjoint and have one root each. They need not be smaller
than every later moving radius; concentric nesting at the first
endpoint already identifies the same root. The same reasoning
identifies each final centre with its corresponding initial root.

All five witnesses have exactly the same ordered initial centre row,
and each final unordered row equals it. I independently checked the
raw shape of all **741** rows: two coordinates per base point and
four two-coordinate root centres. The default conversion permits
some shorter tuples syntactically, so this independent raw-universe
check is recorded explicitly; no such tuple occurs in these frozen
witnesses.

The contraction of each actual root to its affine centre remains
inside a disjoint convex disk. It glues at every endpoint, including
the common basepoint. Every loop therefore uses the same actual-root
to rational-centre base identification. No independently chosen
conjugation is inserted between loops.

The positive common-translation and the collision/curvature hostiles
remain present and pass in the adapted source. Default verification
uses no numerical root solver: `numpy` is confined to the explicitly
separate proposal/scout functions. Their numerical choices do not
authorize any segment of the default exact replay.

## 3. All thirty base edges and an independent crossing computation

The exact common base is `4+3i`; every declared radius is `1/32`.
The source checks each witness point against its declared straight
edge, with strictly increasing rational edge parameter, and requires
all six edges to finish. I additionally reconstructed those parameters
using a determinant collinearity test and a dot-product parameter,
rather than the source's coordinate-division test. The exact endpoint
row indices were:

| Path | Six edge endpoint row indices |
|---|---|
| smooth | `33,42,51,60,69,102` |
| cusp_plus | `48,61,75,88,102,148` |
| cusp_minus | `60,73,86,99,112,172` |
| node_real | `66,75,84,93,102,169` |
| node_lower | `54,63,72,81,90,145` |

Thus no corner, final return, or section of a stem is omitted. The
rational node centres are not falsely treated as exact roots of
the cubic. Separate isolation of the enclosed critical values is
unnecessary: the root tubes certify that these actual paths avoid
the discriminant, and every such based loop imposes a valid relation.

I also recomputed all crossings without the source's current-order
update algorithm. At each exact crossing time I calculated the
instantaneous projected position of every strand; the generator
index is one plus the number of the other strands strictly to its
left. This gives the same adjacent crossing and sign. The raw
crossing counts are `3,3,7,10,10` in the order above, and after only
adjacent inverse cancellations the independently reconstructed words
are respectively

    [1,2,-1],
    [1,1,1],
    [2,1,3,3,3,-1,-2],
    [2,1,3,3,1,1,-3,-3,-1,-2],
    [2,-3,2,2,3,-2].

Endpoint projection ties, coincident crossing times, and collisions
were excluded exactly. Multiplication by `1+i/4` preserves complex
orientation. The sign convention is the inherited below-stem one:
an initially left strand passing below gives the positive
counterclockwise half twist. The exact tube homotopy, not agreement
with a stored word alone, relates these centre polygons to the
actual roots.

## 4. Common-access topology and the arbitrary-group implication

Monicity gives a uniform root bound over any compact base disk.
Choose a constant reference fibre point below that bound and use it
as a section. In particular the five selected loops have section
lifts contracting over their actual fillings. Their geometric
actions therefore fix the same positive meridian tuple exactly,
not merely up to separate simultaneous conjugacies.

For completeness, fibre-generation of the full affine complement
does not require the selected five loops to generate the punctured
base. Take a larger disk containing all six critical values, with
the same type of outside section. The regular-base bundle group is
generated by the four fibre meridians and lifts of a full set of
base generators. Restoring the critical fibres kills all those
section lifts. Any loop in the full complement can be perturbed
off the finitely many critical vertical fibres. Hence the four
fibre meridians surject onto the actual complement group. Retaining
only five certified relations yields a quotient bound in the
correct direction; omitted relations cannot enlarge that group.

In the chronological representation convention
`H(a,b)=(aba^-1,a)`, the positive geometric loop action is the
inverse-Artin map and contravariance gives this Hurwitz action.
A prefix followed by a core and the inverse prefix tests the core
on the prefix-transformed tuple. I independently checked the exact
five implications in an arbitrary group:

1. The smooth prefix `[1]` exposes `(a,c)` to one half twist, so
   `a=c`.
2. The first cusp gives `aba=bab`.
3. The real-node prefix `[2,1,3,3]` exposes
   `(abcb^-1a^-1,a)` to a square. With `c=a` and the braid relation,
   the first entry is `aba b^-1a^-1=b`. Thus `[a,b]=1`, and the
   braid relation now gives `a=b=c`.
4. The lower-node prefix `[2,-3]` exposes `(bcb^-1,d)`, now `(a,d)`,
   to a square; hence `[a,d]=1`.
5. The other cusp prefix `[2,1]` exposes `(b,d)` to a cube. With
   `b=a` and commutation, `bdb=dbd` forces `d=a`.

All four positive generators coincide. The implication has no
finite-group, degree, or transitivity hypothesis. Equal tuples
fix all five words, as the source's reverse control checks.

The actual group is consequently cyclic. The reduced irreducible
equation `F` maps a positive meridian to winding one in `C*`.
This proves that the cyclic group is infinite. Equivalently the
map from the selected-relation quotient `Z` to the actual group,
followed by winding, is the identity of `Z`; no omitted relation
can create torsion. The conclusion is the actual affine group,
not a semilocal passport or only a permutation quotient.

## 5. The whole-support consumer and the complete G11 transport

If the literal curve were the whole irreducible nonproperness curve
of a nonautomorphic planar Keller map, the complement would carry
its connected finite étale cover of generic degree `d>1`. The
birational degree-one exclusion and the retained-sheet argument
are the fully proved consumer in
[the earlier two-cusp braid, Section 5](planar_jc48_sep08_two_cusp_braid.md).
In particular `F` composed with the Keller map is a nonconstant
polynomial. A component of its zero divisor maps dominantly onto
the irreducible curve by quasifiniteness, and gives an actual
retained sheet at a generic smooth target point. A positive
meridian fixes it.

But a transitive action of a cyclic group generated by that
meridian is a single d-cycle. It has no fixed label for `d>1`.
This gives the literal whole-support exclusion in every mapping
degree. It does not apply to the curve merely as one component of
a larger nonproperness support.

The new primary's Section 6 correctly invokes the already proved
family theorem, rather than deriving a family from numerical
continuity. I read that theorem's full parameter and transport
argument. It retains both `s=1` and `s=-1` projection-critical
cases, and never divides by `s-1` or `s+1`. The derivative/common
zero conditions give a birational finite normalization. The
actual source and target affine operations reduce the declared
(4,6) class to its stated parameter spaces.

Within the infinity-eleven graph, the good locus is a nonempty
Zariski open subset of an irreducible line. Proper resolved
off-diagonal incidence pays its openness; both ordinary cusp
sections and the marked infinity branch have fixed resolutions.
The first three infinity blowups separate the original infinity
line from the branch point. Later even translations fix the
exceptional axis and its marked intersections, so vanishing an
intermediate even coefficient is not discarded. Resolving the
entire finite étale node multisection avoids a global node-label
assumption. The resulting proper normal-crossing pair transports
the actual affine complement along paths in the connected good
locus, preserving positive curve meridians up to access conjugacy.

The literal point audited here is `s=2` in this good locus. Thus
its cyclic group transports to **all of G11**. Together with the
previously proved G7 and G9 suppliers and the complete infinity
exhaustion, this pays the proposed whole-class consumer: every
irreducible curve with a birational polynomial normalization of
degrees (4,6), exactly two ordinary finite cusps and otherwise only
ordinary nodes, has the stated cyclic affine complement and is
excluded as whole Keller nonproperness support. The normalization
degree pair is not the generic degree of a Keller map.

No isotopy across different infinity strata is claimed. The older
family file's retained OPEN G11 endpoint is a historical proof
boundary paid by this new certificate after promotion, not a
dependency asserting its own conclusion. Other degree pairs,
additional finite singularities, and larger reducible supports
remain outside the consumer. JC(2) remains OPEN.

## 6. Exact reproduction and pins

Independent commands were

    python3 -B 04-computation/planar_jc48_sep08_infinity11_braid.py
    python3 -B -O 04-computation/planar_jc48_sep08_infinity11_braid.py

Both pass **6,793 always-active gates** over all **736 moving
segments**. Their complete outputs are byte-identical to the
frozen **2,461-byte** output. There are 39 exact polynomial/group/
method controls before the individual path verifications.

| Path | Segments | Path gates | Gzip bytes | Raw bytes |
|---|---:|---:|---:|---:|
| smooth | 102 | 936 | 6923 | 26744 |
| cusp_plus | 148 | 1351 | 9885 | 38692 |
| cusp_minus | 172 | 1575 | 11857 | 45075 |
| node_real | 169 | 1554 | 12621 | 46034 |
| node_lower | 145 | 1338 | 8592 | 41259 |

The source is **14,600 bytes**, SHA256
`2df6b5892ea04d0a771936c9ca21838bc727d131804cdefeae03bcfc09e2b26f`.
The frozen output and both independent replays have SHA256
`ca8153a95aa81dcfd06991d30bdabe03c11898fd7dc36ad28f55a394161f6549`.
All five raw/compressed witness pins were independently recomputed:

| Path | Gzip SHA256 | Raw SHA256 |
|---|---|---|
| smooth | `821f52e0199693af3222e6fe8cc506c6f9341194bcc38b2c079c9a8ae999e4a0` | `3e773ba9e2a041f8a51a5a45acc2f92e7d7d92f3d12fae57a443b0fe305b79e9` |
| cusp_plus | `31c8b46dfd0775e80790ab0188602d020579fe990b857c817418a58feee94cea` | `8e9ca689d19d2ba241e71147cf51b1d3d6028425a61afd8bc63c66980a8c8c92` |
| cusp_minus | `05332d61d143522b7ea2184b7dc964201dc66c91e5c3da155ab76d4cfb56678e` | `7c1dd3a237dc7db81683e0d7aca2a9f723de45f242e730ef14d1bc650317529f` |
| node_real | `f1097a644150953f07a40cf57b06e437a820358bbbfbc5707ce503d83f986e22` | `a1cb6797e007708a5cae21c9d03df76c6125cbfe46d0257e471d02ba6ce7a1b8` |
| node_lower | `e8051f7e47493fbe5ccf47a0693f8d849f4213df4a14ab848823343159380860` | `b3828eb89e8d9bccc2749fe41c29bf0b0012a2c43bcbb963444dbe6e48f2178f` |

The new Section 6 consumer was read after its addition. No source,
witness, or mathematical correction remains necessary. This audit
records acceptance; it does not itself edit the primary status or
the earlier family proof.
