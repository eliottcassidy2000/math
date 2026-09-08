# Independent audit of the three-cusp braid and marked local pairs

**Status: full independent analytic/source audit PASS; normal,
optimized, and frozen 9,534-gate output agree.** The actual result
is a marked Artin-D4 quotient and a mapping-degree floor of eight
for the declared whole nonproperness support. It is not an
all-degree exclusion or an identification of the whole complement
group with Artin D4.

Audited:

- [Primary proof](planar_jc48_sep08_three_cusp_braid.md).
- [Frozen source](../../04-computation/planar_jc48_sep08_three_cusp_braid.py),
  all six compressed rational witnesses, and the
  [frozen output](planar_jc48_sep08_three_cusp_braid.out).
- The [proved marked group and conditional degree-eight
  passport](planar_jc48_sep08_three_cusp_passport.md).
- The already audited moving-centre method and the
  [infinity-eleven predecessor](planar_jc48_sep08_infinity11_braid.md).

The new load-bearing input is the identification of the actual
local meridian pairs, including their retained sets. A whole-loop
word alone would not establish that identification. The proof and
source now supply it independently.

## 1. Literal geometry and completeness

The source retains the exact parametrization

    U=t^4-(8/3)t^3-2t^2+8t,
    V=t^6-4t^5-(3/2)t^4+(52/3)t^3-32t.

The entire literal coefficient table agrees with the independent
symbolic resultant `Res_t(U-u,V-v)`, and substitution of `(U,V)`
vanishes exactly. Its being monic quartic in `v` is used throughout;
no leading coefficient was discarded in either certifier.

The nonzero `v` discriminant is exactly

    -167772160000/5559060566555523
    *(3u-13)^3*(3u-8)^3*(3u+19)^3
    *(729u^3+3483u^2+5547u-78359)^2.

I independently recomputed

    disc_t(U-u)=-(256/27)(3u-13)(3u-8)(3u+19),
    disc(H)=-94850468536320000!=0.

The ordinary-cusp data are

| Parameter | `U` | `V` | `U'' V'''-V'' U'''` |
|---:|---:|---:|---:|
| `-1` | `-19/3` | `109/6` | `-5760` |
| `1` | `13/3` | `-115/6` | `128` |
| `2` | `8/3` | `-40/3` | `1152` |

The derivative gcd and each actual image-preimage gcd agree with
the primary. Thus these are ordinary cusps, each with exactly one
normalization preimage. At any root of `H`, the four parameter
branches are analytic and unramified under `U`. The discriminant
has exact order two, so precisely one height difference has a
simple zero. This is one transverse double point, with neither
three branches nor tangency. The cubic has three simple roots
outside all cusp projections. Each cusp accounts for the full
order-three discriminant factor at its own projection. These
facts exhaust the affine singularities.

Nonvanishing of the discriminant proves that the generic four
parameters over `u` have four different image heights. Hence the
parametrization is birational onto its irreducible image, rather
than a multiple parametrization hidden in the resultant. Monicity
of `U` makes the affine parametrization finite, so `A1` is its
affine normalization. Homogenization to degree six has no base
point, including the unique infinite parameter, whose image is
`[0:1:0]`.

The literal infinity shear has first surviving term

    Z-X^3+(7/3)X^4=-(160/27)t^(-9)+O(t^(-10)),
    X=t^(-2)+O(t^(-3)).

It gives the single `(2,9)` branch. Its delta four, together with
three cusp and three node contributions, totals ten as required
by the rational sextic. This genus check is consistent with the
independent discriminant inventory; it was not used to replace it.

## 2. Complete moving-root verification and common gauge

I read the full source, including both Taylor expansions, the
exact affine-distance minimum, every raw-data check, and the
chronological free-group action. For the moving path compiler,
an original translated coefficient `b_jk` contributes

    b_jk * du^j * binom(k,l) * dz^(k-l)

to `h^(j+k-l) w^l`. Thus every term from the monic quartic is
retained. The linear Rouché inequality bounds all powers of the
real path parameter by one. The exact minimum of the affine
`L-infinity` distance is attained at an endpoint or a change
among the four affine absolute-value branches; the implemented
real, imaginary, sum and difference zeros include all of them.
Four disjoint one-root disks exhaust the quartic at every point
of each accepted segment.

Shared labelled centre endpoints yield concentric isolating disks
on adjacent segments. Their common root is the same even when
the radii differ. The initial tiny disks identify the common
actual roots, and the same root-to-centre contraction works at
the endpoint of the whole unordered loop. Interpolating inside
the disjoint disks is a configuration homotopy; it does not
assume a numerical root matching is correct.

All six witnesses have the same actual base `4+3i` and the same
ordered initial centres. Their unordered final centres agree with
the initial set. I independently checked that every raw base
point has exactly two rational coordinates and every root row
has exactly four centres with two coordinates each: **813 rows,
807 segments** in total. This also pays the precise frozen-data
shape despite the complex-number constructor's optional defaults.

For a second exact path check I used determinant collinearity
and dot-product edge parameters, rather than the producer's
selected real/imaginary coordinate parameter. Every segment
lies on the declared directed edge, with strict progress. The
six edge-end indices are

| Path | Edge ends |
|---|---|
| cusp_plus | `51,66,81,95,110,161` |
| cusp_minus | `61,74,87,100,113,174` |
| cusp_two | `42,55,68,81,94,135` |
| node0 | `31,40,49,58,67,99` |
| node1 | `51,60,69,78,87,141` |
| node2 | `30,39,48,57,66,97` |

I separately extracted crossings from the projected affine
strands by their instantaneous ranks, without updating the
producer's running strand order. The projected real tie time is
exact rational; the left strand before that time and the sign
of its imaginary separation determine the crossing. Simultaneous
third-strand ties and zero imaginary separations are excluded.
This independent path reproduces every full word and every
incoming stem. The raw crossing counts are respectively
`7,7,7,8,8,6`; cancelling adjacent inverses gives the six declared
words. In particular the extra cancellations in cusp_two and
node0 do not erase their distinct actual incoming stems.

The orientation convention is unchanged from the audited
predecessor: multiplication by `1+i/4` preserves orientation;
a left strand going below gives the positive letter; below-stem
positive meridians have inverse-Artin geometric pullback and
contravariant chronological tuple transport
`H_i(a,b)=(aba^-1,a)`. The proof does not mix a reversed word
with a different tuple convention.

Finally monicity pays the vertical free-group epimorphism. A
constant outside-root section over a compact disk containing
all critical values and the filling regions contracts the base
loops in the actual complement. Their words therefore fix the
same meridian tuple exactly, not only its simultaneous-conjugacy
class. The six chosen lassos need not be a basis of the punctured
base group for these necessary relations to hold.

## 3. The new complex-disk cluster proof is independently checked

The local certifier operates over the entire **complex** parameter
disk `|u-c|<=1/32`. It compares the translated quartic with a
constant nonzero monomial of degree two for the cluster and
degree one for each other disk. The factor `(1/32)^j` bounds
each complex parameter power. Uniform strict centre separation
ensures that these degree counts describe disjoint disks and
exhaust all four roots. The linear one-root disks cannot contain
a critical collision, since the Rouché count is with multiplicity.

As an independent calculation I parsed only the literal coefficient
and centre data, and re-expanded the original monomials directly.
For an original term `f_ij u^i v^j`, substituting
`u=c+h`, `v=z0+Lh+w` gives the coefficient contribution

    f_ij binom(i,a) binom(j,k) binom(j-k,b)
             c^(i-a) z0^(j-k-b) L^b

to `h^(a+b)w^k`. This route does not use the producer's staged
Taylor compiler. All **18** resulting whole-complex-disk Rouché
margins are strictly positive, and all three pair separations per
path, eighteen in total, pass. The cluster radii are exactly
those printed in the primary: `1/64` at the cusps, `5/64` at
the real node, and `3/16` at the conjugate nodes; every other
root disk has radius `1/64`.

For each node, linear Rouché applied to the exact cubic `H`
gives one root in the radius-`1/32` disk and one in the smaller
radius-`1/64` disk. The latter is strictly inside the diamond,
so the lasso encloses that critical value. At the rational cusp
centres, the constant comparison excludes any root of `H`.
Direct separation excludes the other rational cusp values.
Thus each base disk contains exactly one discriminant point,
with the type already established from the literal curve.

At the stem endpoint `c+1/32`, the additional tiny one-root
disks are strictly contained in exactly one local disk, with
partition `2+1+1`. Their projection-error bound is valid because
the modulus of `1+i/4` is at most `5/4`. The two cluster roots
are adjacent in the actual projected order, and every spectator
lies outside the projected strip of the cluster disk. Hence
the standard below-stem paths to the cluster punctures have a
common local access stem with no intervening puncture.

After translating its holomorphic centre, the two-root disk is
a proper local two-sheet family over the parameter disk, with
one ordinary cusp or node and no other critical value. Moving
toward a sufficiently small neighborhood of that singular value
within the disk identifies this pair with the ordinary local
two-puncture meridians. The two spectator roots remain in their
separate one-root disks. This establishes a geometric local
pair, rather than inferring one from the algebraic whole word.

## 4. Exact pair transport and the retained-set qualification

Let `e=bcb^-1`, `f=e^-1ae`. The independently recovered actual
stem pairs are

| Path | Actual stem | Actual pair |
|---|---|---|
| cusp_plus | `[-1,-2]` | `(b,c)` |
| cusp_minus | `[2,1]` | `(b,d)` |
| cusp_two | `[2,-1]` | `(e,e^-1ae)` |
| node0 | `[2,-1,-2]` | `(b,b^-1fb)` |
| node1 | `[2,1,2]` | `(a,d)` |
| node2 | `[2,3]` | `(e,bdb^-1)` |

The formal prefix in cusp_two differs from its actual stem by
one inverse Hurwitz move within the colliding pair. The same
is true for node0. The frozen free-word checks verify precisely
those identities. Their difference is not being interpreted as
a newly certified loop in the `u` base.

For an actual local distinguished-basis change

    (sigma,tau) -> (tau,tau^-1 sigma tau),
    (A,B) -> (B,tau^-1 A),

the retained sets travel with the access paths. Since `B` is
pointwise fixed by `tau`, the intersection cardinality is
unchanged. The complement of `B` is also `tau`-invariant, so the
deleted-overlap cardinality is unchanged as well. The two
generators give the same local subgroup. When retained sets are
the full fixed sets, their intersection equality also follows
immediately from this unchanged subgroup.

The primary correctly avoids a stronger unsupported algebraic
inference. An abstract initial equation `B=(sigma tau)A` alone
does not automatically supply that particular re-access equation
after changing the distinguished basis; one would additionally
need the relevant peripheral preservation statement. Here the
ordinary-cusp injection is applied to the directly certified
actual pair, and only its actual counts are transported. That
is sufficient for every use in the degree argument.

At node2 both meridians and both retained/deleted sets are
simultaneously conjugated by `b^-1`, giving `(c,d)` without
losing their intersection. The crucial pairs for the degree-six
and degree-seven support contradiction are cusp_plus, cusp_minus,
and node2; the first two are already literal actual stems, and
the third has only this common conjugation. Thus no rebased
abstract re-access assumption is used in the decisive step.

## 5. Actual group and mapping-degree consumer

The complete six words impose exactly the six relations in the
proved passport. I checked its reversible reduction once more:
from the first braid one gets `e b e^-1=c` and `c e c^-1=b`;
conjugating node0 gives `[a,c]=1`, and then conjugating the
third cusp gives `aba=bab`. The other leaf commutators and
central braids are retained. Each step reverses, giving Artin
D4 with the same marked letters and no involution relations.

The four-generator vertical epimorphism and these six actual
necessary relations therefore give a surjection

    Artin(D4) -> pi_1(C^2 minus C).

Omitted base relations can make the target a further quotient;
they cannot invalidate this surjection. The finite Weyl group
of order 192 is not inferred from this Artin presentation.

For the hypothetical whole irreducible Keller support, the
actual normalization and three-cusp/three-node inventory are
paid by Section 1. The local marked data are paid by Sections
3--4. The inherited actual Euler identity is

    1=-2k+n1+n2+n3+sum node overlaps.

The uniform cusp inequality and its equality consequences apply
to each certified actual cusp pair. The degree-five argument
uses only those counts and the transposition obstruction in
D4. In degree six, and in the surviving degree-seven case,
the meridians are single three-cycles and all retained sets are
their full fixed sets. The node2 supports of `c,d` are then
disjoint, while each must meet the support of `b` in at least
two letters by cusp_plus and cusp_minus. This is impossible
for the three-letter support of `b`. The degree-seven smaller
retained-count case is excluded by the actual saturation/node
parity argument, without requiring any formal third-pair
re-access equation.

Thus the already proved conditional passport is legitimately
consumed: mapping degrees five, six and seven are excluded.
The inherited classical geometric exclusions in degrees two
through four, together with the nonautomorphic degree-one
boundary, give the stated floor eight. These are mapping
degrees, not the parametrization's polynomial degrees four
and six.

No conclusion about all larger mapping degrees, polynomial
realization of a passport, an individual component of a larger
nonproperness set, or a whole parameter family is justified
or claimed. In particular this is not an all-degree JC(2)
exclusion for the curve.

## 6. Frozen replays and all witness pins

Independent commands:

```sh
python3 -B 04-computation/planar_jc48_sep08_three_cusp_braid.py
python3 -B -O 04-computation/planar_jc48_sep08_three_cusp_braid.py
```

Both pass **9,534 always-active gates** and reproduce all
**4,892 bytes** of the frozen output. There are 53 initial
geometry/algebra/method gates. The source is **20,428 bytes**,
SHA256 `369c55582cb7f53cbd7a6bfa1c28d267190e9f58dbc9c9d22291be2ddbf05e27`.
The frozen output and both independent replays have SHA256
`2a2fc8fe5d91ca032e249c703a96b19738b4363431f836933536ac006d11d252`.

| Path | Segments | Path gates | Gzip bytes | Raw bytes |
|---|---:|---:|---:|---:|
| cusp_plus | 161 | 1863 | 10474 | 41941 |
| cusp_minus | 174 | 2051 | 12021 | 45806 |
| cusp_two | 135 | 1567 | 8832 | 35230 |
| node0 | 99 | 1171 | 6490 | 26406 |
| node1 | 141 | 1689 | 11651 | 40570 |
| node2 | 97 | 1140 | 7986 | 27852 |

All raw and compressed SHA256 pins were independently recomputed:

| Path | Gzip SHA256 | Raw JSON SHA256 |
|---|---|---|
| cusp_plus | `96328859a5811a761d67d6d471a61d97d418aabc5771ba6219a321bea80c570d` | `970f0e8a9c11a9199398c5a70213f45253ca4817adda8979840bf36fbe05ab2a` |
| cusp_minus | `c630b15221dc82d808ce8dacfed9f5bb468bb335974af7a82c939171160b2f6b` | `41c7915654f748b41b83a55ab369dcdc07e62f327c8e0abbcae32030c183b190` |
| cusp_two | `8c27247e644f433f6881ec7976411fcec5f9b24fc300bbab8699dd52a376750f` | `230c49f542429e1be158ebeb0dc3bccc713999acef2251b751562850d17790ef` |
| node0 | `5f2d88d58c7509996c79f1fc4608bb6167a4fb07858bfc9860f635affe5481a1` | `8710b8f189dc7806ea6cd8930f6a2a866d37dd57f6992224d997b70b5c2f2356` |
| node1 | `e9321d7e7f0b2fa454afd09a1bd6d2ff5c55961b3979c61ad550c9f712964cca` | `8c82826f8039a318dbd9383d76697b43eaa097161c9f8af5bccbe9c2b440c023` |
| node2 | `902768fdd1968fc697f25e9a975318f50094401974e952a1d67bdc39533999ab` | `fa10c1e76d49683f65e1660b4db4755e7bc1bfdd12db3a7befbad5b6999ac1bd` |

The independent raw-path and direct-monomial local-cluster
calculation used only parsed literal data, not imported primary
mathematical functions. Its equations and complete finite
universe are recorded above. No source, witness, or mathematical
correction remains necessary. The owner may promote the primary
after accepting this frozen audit; earlier scout artifacts retain
their own historical heuristic status.
