# Six certified local lassos for the three-cusp sextic

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The [full independent audit](planar_jc48_sep08_three_cusp_braid_audit.md)
accepts the actual geometry, all six rational paths, the marked local
clusters and pair repair, the group quotient, and the degree-floor
consumer. Independent normal and optimized replays match the frozen output. The earlier numerical
scouts remain heuristic provenance, not proof inputs.

## 1. Actual object and precise conclusion

Let C be the curve parametrized by

    U(t)=t^4-(8/3)t^3-2t^2+8t,
    V(t)=t^6-4t^5-(3/2)t^4+(52/3)t^3-32t.                 (1)

It is a rational sextic with birational normalization A1, three ordinary
finite cusps, exactly three ordinary finite nodes, and one infinity
branch of type (2,9). The four positive meridians of a regular vertical
fibre generate its affine complement group. This note proves that their
six certified common-access relations give a surjection

    Artin(D4) -> pi_1(C^2 minus C),                       (2)

with the same four marked positive generators. No assertion that (2)
is an isomorphism is made: the six lassos need not be a free generating
system for the punctured projection base, and additional actual
relations may remain.

The local marking is retained. It pays the actual cusp/node-pair input
of the [three-cusp passport](planar_jc48_sep08_three_cusp_passport.md).
After that independently audited conditional theorem is invoked, C
cannot be the whole irreducible nonproperness curve of a polynomial
Keller map of mapping degree 5,6 or7. The inherited mapping-degree
2,3,4 exclusions then give a necessary degree of at least eight. This
is not an all-degree Keller exclusion or a realization theorem.

The closest mechanism is the [five-loop infinity-eleven
certificate](planar_jc48_sep08_infinity11_braid.md), using the
[audited moving-centre method](planar_jc48_sep08_moving_tubes.md).
The current six-word group is different: the [exact S4
hostile](planar_jc48_sep08_three_cusp_group.md) blocks cyclicity and
blocks generation by any one declared cusp subgroup. The corrected
near miss is treating an algebraic conjugating prefix as an actual
local access path. The new sidecar certifies the incoming stem, a
local two-root disk, and its identification with the marked pair.

The five live concepts are actual quartic fibres, the complete six-value
discriminant, exact root motion, marked local clusters, and actual sheet
counts. The map from root motion to its braid retains the common gauge;
the map to (2) discards additional possible global relations. The
cluster sidecar retains the local meridian subgroup and intersection
counts that a bare whole-loop word would lose.

## 2. Complete geometry and literal quartic

Set F(u,v)=Res_t(U(t)-u,V(t)-v). Every coefficient of F is a literal
Gaussian-real rational in the standalone source. A separate symbolic
resultant and substitution verify them. The polynomial is monic of
degree four in v and has total degree six. Its discriminant is

    disc_v(F) = -167772160000/5559060566555523
      * (3u-13)^3(3u-8)^3(3u+19)^3 H(u)^2,
    H(u)=729u^3+3483u^2+5547u-78359.                      (3)

The cubic H has three distinct roots and avoids the three displayed
rational values. Independently,

    disc_t(U-u)=-(256/27)(3u-13)(3u-8)(3u+19).            (4)

The nonzero discriminant (3) proves that, for generic u, the four
distinct preimages under U have four distinct V values. Thus the
parametrization is generically one-to-one onto its irreducible image;
the resultant has no repeated image factor. The homogeneous degree-six
parametrization extends to P1 without basepoints and maps its sole
infinite parameter to [0:1:0]. Its affine normalization is A1 since U
is monic and makes the parametrization finite.

The monic derivative gcd is (t-2)(t-1)(t+1). At t=-1,1,2, respectively,
the second/third derivative determinants U''V'''-V''U''' are
-5760,128,1152. Each is nonzero. At each such e, the monic gcd of
U(t)-U(e) and V(t)-V(e) is exactly (t-e)^2, so there is no other
preimage at the cusp image. These are three ordinary cusps with
projection values -19/3,13/3,8/3.

At a root u0 of H, (4) makes all four source parameters analytic
functions of u near u0. Hence all four image heights V_i(u) are analytic.
The discriminant product is the square of the product of their six
pairwise differences. Its exact order two at u0 forces exactly one
pairwise difference to vanish, to order one. Thus precisely two
immersed branches meet, with different slopes as graphs over u: one
ordinary node, with no triple image. This gives three distinct nodes.
At a cusp projection, the ordinary cusp alone contributes order three
to (3); any additional collision there would increase that order.
All possible affine singular images project into (3), so the inventory
is complete. This argument does not infer a reduced collision scheme
from a repeated projection factor.

For an independent infinity/genus control put z=1/t, X=U/V and Z=1/V.
Exact numerator/denominator valuations give

    X=z^2+O(z^3),
    Z-X^3+(7/3)X^4=-(160/27)z^9+O(z^10).                (5)

The local coordinate shear in (5) gives the single (2,9) branch.
Its delta invariant is four. The three ordinary cusps and three nodes
supply the remaining six units of the rational sextic arithmetic
genus ten. This is consistent with the direct finite inventory, rather
than a replacement for the complete discriminant check.

## 3. Rational loops, moving tubes and common based words

All six loops use u_*=4+3i and r=1/32. For the centre c in the table,
the exact vertices are

    u_*, c+r, c+ir, c-r, c-ir, c+r, u_*.                 (6)

The node centres approximate their critical values; the local
Rouché checks in §4, not their decimals, isolate the actual values.

| Name | Exact centre c | Exact reduced word |
|---|---|---|
| cusp_plus | 13/3 | `[-1,-2,1,1,1,2,1]` |
| cusp_minus | -19/3 | `[2,1,3,3,3,-1,-2]` |
| cusp_two | 8/3 | `[2,1,1,1,-2]` |
| node0 | 1688683/524288 | `[2,-1,2,2,1,-2]` |
| node1 | -2096807/524288 - (4371107/1048576)i | `[2,1,2,3,3,-2,-1,-2]` |
| node2 | -2096807/524288 + (4371107/1048576)i | `[2,3,2,2,-3,-2]` |

For a Gaussian rational w write A(w)=|Re(w)|+|Im(w)| and
B(w)=max(|Re(w)|,|Im(w)|), so B(w)<=|w|<=A(w). On a rational base
segment u(h)=u0+h*du, the four proposed centres are affine,
z_i(h)=z_i0+h*dz_i. Expand the complete polynomial exactly:

    F(u0+h*du,z_i0+h*dz_i+w)=sum C_jk h^j w^k.

The verifier accepts only if

    B(C_01) r_i > sum_(j,k)!=(0,1) A(C_jk) r_i^k.        (7)

For every real h in [0,1], this is a strict Rouché comparison with the
same nonzero linear polynomial C_01*w. Pairwise affine centre distances
are bounded below exactly in the B norm by testing endpoints and all
zeros of real part, imaginary part, their sum and their difference.
The chosen radii are one sixteenth of the minimum relevant pairwise
distance. The four disks stay disjoint, so (7) gives one continuously
labelled root in each at every h.

Adjacent segments share identical labelled centre endpoints; their
concentric isolating disks glue to the same actual roots. All six
witnesses share identical initial centres in the same order, and their
final unordered centre sets equal that initial set. Tiny initial root
disks pay a common actual-root-to-centre contraction. Interpolating
inside disjoint disks proves that the actual root path and the
polygonal centre path represent the same based configuration loop.

After projection by the orientation-preserving complex multiplier
1+i/4, every polygonal crossing is calculated with exact rational
arithmetic. Tied endpoint real coordinates, simultaneous crossings,
nonadjacent crossings and zero imaginary crossing separation are
rejected. A left strand passing below gives a positive letter. With
below-stem positive meridians the chronological representation
transport is H_i(a,b)=(aba^-1,a). The geometric loop action is the
inverse-Artin pullback; the representation transport is contravariant.
This is the convention of the audited predecessor, used unchanged.

Monicity in v gives an epimorphism from the four-generator vertical
meridian group to the affine complement group. Over a compact region
filling all six base loops choose a constant v_* above the root bound.
The section (u,v_*) contracts those base loops inside the actual
complement, so each exact Hurwitz word fixes the same generator tuple,
not just its simultaneous-conjugacy class. No assertion that the six
loops generate the punctured base group is needed for this necessity.

## 4. The new marked local-cluster sidecar

A whole-loop word of the form P,core,P^-1 is not by itself an
identification of P with a geometric local access path. Here the
incoming edge in (6) and the complete two-root cluster are certified
separately.

For each c, the source gives three literal affine holomorphic centres

    w_j(u)=w_j0+L_j(u-c),        j=0,1,2,

with rational complex w_j0,L_j and fixed rational radii R_j. They are
valid on the entire complex disk |u-c|<=r, not just the diamond boundary.
The first disk contains two roots and the other two contain one root
each. To prove this, expand exactly, with a complex parameter h,

    F(c+h,w_j0+L_j h+w)=sum D_ell,k h^ell w^k.

For m_0=2 and m_1=m_2=1 the verified strict inequalities are

    B(D_0,m_j) R_j^m_j
      > sum_(ell,k)!=(0,m_j) A(D_ell,k) r^ell R_j^k.       (8)

Rouché compares with D_0,m_j*w^m_j uniformly for every complex |h|<=r.
This uses the same full affine coefficient compiler as (7), with each
complex parameter power now bounded by r^ell. Separation is certified
by the elementary lower bound

    B(w_i0-w_j0)-A(L_i-L_j)r > R_i+R_j.                 (9)

The radii (cluster, other, other) are respectively (1/64,1/64,1/64)
at all three cusps; (5/64,1/64,1/64) at node0; and
(3/16,1/64,1/64) at node1 and node2. The full centres and slopes are
literal data in LOCAL in the source, with no numerical choice in replay.

At each node centre, linear Rouché for H at radii r and r/2 proves
exactly one root in each. The smaller disk lies strictly inside the
r-diamond. At a cusp centre, constant Rouché proves H has no zero in
the r-disk; its single rational cusp value is the centre. Direct
rational separation excludes every other cusp value. Thus each disk
has exactly one discriminant point, of the type proved in §2. The two
other root disks are unramified throughout, so the unique singular
interaction lies in the two-root disk.

At the incoming stem endpoint c+r, four additional tiny Rouché disks
identify the actual roots with the stored centres. Each tiny disk is
strictly contained in exactly one of the local disks, giving the
partition 2+1+1. Rational separation also proves the actual projected
root order agrees with the centre order. The two local roots occupy
adjacent positions; every other root lies outside the projected strip
of the two-root disk. Hence the two adjacent below-stem meridians are
the local two-puncture meridians, transported along the certified stem.
The other punctures do not intervene in their local access strip.

More explicitly, the two-root disk over the base disk is a proper
local two-sheet family with one ordinary cusp or node and no other
critical value. It identifies this two-puncture system with a nearby
Milnor-fibre meridian pair by transport inside that disk. The remaining
punctures lie in the separate one-root disks. This pays an actual local
pair, not just two letters whose product has a matching cycle type.
The outside-root section fixes the common basepoint throughout.

The exact incoming stem words and one-based adjacent pair positions
are as follows. Put e=bcb^-1 and f=e^-1ae.

| Name | Actual incoming stem | Pair positions | Actual stem pair |
|---|---|---|---|
| cusp_plus | `[-1,-2]` | 1,2 | (b,c) |
| cusp_minus | `[2,1]` | 3,4 | (b,d) |
| cusp_two | `[2,-1]` | 1,2 | (e,e^-1ae) |
| node0 | `[2,-1,-2]` | 2,3 | (b,b^-1fb) |
| node1 | `[2,1,2]` | 3,4 | (a,d) |
| node2 | `[2,3]` | 2,3 | (e,bdb^-1) |

At cusp_two the formal prefix [2] gives (a,e); the actual prefix
adds H_1^-1 inside this very pair. At node0 the formal prefix [2,-1]
gives (f,b); the actual prefix adds H_2^-1 inside that pair. Thus both
repairs are local distinguished-basis changes, not asserted u-base
loops or arbitrary conjugators. For an inside-pair inverse Hurwitz move

    (sigma,tau) -> (tau,tau^-1 sigma tau),
    (A,B) -> (B,tau^-1 A),

actual local retained sets are transported with their meridians. Since
B is fixed pointwise by tau, intersection cardinality is preserved:
|B intersect tau^-1 A|=|A intersect B|. The deleted-overlap cardinality
is preserved as well, because the complement of B is tau-invariant.
The joint local subgroup is unchanged. We apply the ordinary-cusp
injection to the directly certified actual pair and transport these
counts. We do not infer a new abstract re-access equation from the
algebraic constraints alone. When all retained sets are full fixed
sets, their intersections are simply the joint fixed sets of these
unchanged local subgroups.

Finally the node2 pair is simultaneously conjugated by b^-1 to (c,d),
transporting both sets together. These statements pay all original
access-pair count data in the passport note. In particular its crucial
cusp_plus, cusp_minus and node2 pairs require no inside-pair repair.

## 5. Actual group consequence and exact scope

The certified whole-loop words impose exactly the six literal relations
in the passport note. Its reversible free-group calculation gives
Artin D4 with central vertex b and leaves a,c,d. Combining this with the
actual four-positive-meridian epimorphism proves (2). The group could
be a further quotient. The S4 control satisfies these six relations
and shows that cyclicity is not a consequence of them.

The actual local marking in §4 also pays the declared three-cusp and
three-node pairs. If C were the whole irreducible nonproperness curve
of a nonautomorphic Keller map, actual retained-page constancy and
Euler integration would therefore satisfy the precise hypotheses of
the three-cusp passport, giving the stated exclusion of degrees5,6,7.
Its cited geometric degree2,3,4 inputs yield the necessary degree floor8.
This uses mapping degree, not the degrees4 and6 of (1).

Nothing here identifies a polynomial source from a finite group action,
excludes every higher-degree passport, treats an individual component
inside a larger nonproperness set, or transports the result to an
unproved parameter family. Those require separate suppliers.

## 6. Reproducibility and freeze

The standalone [source](../../04-computation/planar_jc48_sep08_three_cusp_braid.py)
contains the literal F coefficients, all six loops, complete moving
tube verification, the complex-parameter local clusters, exact marked
stem computations, and the group controls. The numerical proposal
routine is used only for --produce and --scout. Default replay consumes
rational data and checks every accepted segment and every local disk
without numerical root finding.

The independent symbolic geometry and affine-substitution controls are
also rerun. Hostiles include swapped affine strands that collide
between separated endpoints, a nonlinear moving-root remainder that
cannot be discarded, and arbitrarily long common translation that the
complete coefficient expansion cancels correctly. The local-disk
addition separately pays actual adjacency and the two inside-pair
basis changes; those conditions are not inferred from the old method.

```sh
python3 -B 04-computation/planar_jc48_sep08_three_cusp_braid.py
python3 -B -O 04-computation/planar_jc48_sep08_three_cusp_braid.py
```

The six witnesses contain **807 rational moving segments**. Normal and
optimized replays pass **9,534 always-active gates**, with byte-identical
frozen output. Exact source and output pins are:


- `planar_jc48_sep08_three_cusp_braid.py`: 20,428 bytes; SHA256 `369c55582cb7f53cbd7a6bfa1c28d267190e9f58dbc9c9d22291be2ddbf05e27`.
- `planar_jc48_sep08_three_cusp_braid.out`: 4,892 bytes; SHA256 `2a2fc8fe5d91ca032e249c703a96b19738b4363431f836933536ac006d11d252`.

Every compressed and uncompressed witness is pinned below. The source
independently recomputes these pins on replay.

| Path | Segments | Compressed SHA256 | Raw JSON SHA256 |
|---|---:|---|---|
| [cusp_plus](planar_jc48_sep08_three_cusp_braid_cusp_plus_certificate.json.gz) | 161 | `96328859a5811a761d67d6d471a61d97d418aabc5771ba6219a321bea80c570d` | `970f0e8a9c11a9199398c5a70213f45253ca4817adda8979840bf36fbe05ab2a` |
| [cusp_minus](planar_jc48_sep08_three_cusp_braid_cusp_minus_certificate.json.gz) | 174 | `c630b15221dc82d808ce8dacfed9f5bb468bb335974af7a82c939171160b2f6b` | `41c7915654f748b41b83a55ab369dcdc07e62f327c8e0abbcae32030c183b190` |
| [cusp_two](planar_jc48_sep08_three_cusp_braid_cusp_two_certificate.json.gz) | 135 | `8c27247e644f433f6881ec7976411fcec5f9b24fc300bbab8699dd52a376750f` | `230c49f542429e1be158ebeb0dc3bccc713999acef2251b751562850d17790ef` |
| [node0](planar_jc48_sep08_three_cusp_braid_node0_certificate.json.gz) | 99 | `5f2d88d58c7509996c79f1fc4608bb6167a4fb07858bfc9860f635affe5481a1` | `8710b8f189dc7806ea6cd8930f6a2a866d37dd57f6992224d997b70b5c2f2356` |
| [node1](planar_jc48_sep08_three_cusp_braid_node1_certificate.json.gz) | 141 | `e9321d7e7f0b2fa454afd09a1bd6d2ff5c55961b3979c61ad550c9f712964cca` | `8c82826f8039a318dbd9383d76697b43eaa097161c9f8af5bccbe9c2b440c023` |
| [node2](planar_jc48_sep08_three_cusp_braid_node2_certificate.json.gz) | 97 | `902768fdd1968fc697f25e9a975318f50094401974e952a1d67bdc39533999ab` | `fa10c1e76d49683f65e1660b4db4755e7bc1bfdd12db3a7befbad5b6999ac1bd` |

Source, output and all six witnesses are frozen. The independent full
analytic, local-marking, source and exact replay audit passes, including
a separate original-monomial reconstruction of all 18 complex-disk
Rouché comparisons. Its 16,841-byte audit file has SHA256
`75b1ead5c1c1291a18876b15d2c0bfe4052af82c102ee19fd5cf287c2f2476fa`.
Earlier frozen bundles are unchanged.
