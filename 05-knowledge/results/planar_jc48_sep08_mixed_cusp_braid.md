# An actual marked Artin H4 quotient for the mixed (5,3,3) cusp curve

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The complete exact verifier and independent analytic/source/witness audit
pass. The theorem uses the full rational paths and complex-disk marking.

## 1. Literal curve and precise conclusion

Let C be the image of the polynomial normalization

    U=t4-(8/3)t3-2t²+8t,
    V=t6-(3/2)t4-16t3+48t.                             (1)

It is a birational rational sextic with one finite (2,5) cusp, two
ordinary (2,3) cusps, three ordinary nodes, and one infinity (2,7)
branch. The exact geometry is supplied by the
[mixed-cusp boundary proof](planar_jc48_sep08_three_cusp_boundary.md),
§§2–4, and independently checked in the present source as summarized
in §2.

**Theorem.** Four positive curve meridians generate pi1(A2\C). In a
positive free-meridian basis x,y,z,w they satisfy

    yzy=zyz,       ywy=wyw,
    xwxwx=wxwxw,
    [x,y]=[x,z]=[z,w]=1.                               (2)

Thus the actual affine complement group is a quotient of the Artin
group with chain

    z --3-- y --3-- w --5-- x,

called Artin H4. No involution relations are imposed. The six fixed-word
relations certified here have exactly this presentation, but the
actual complement may have additional relations. We do not assert
that the six chosen lassos freely generate the punctured-base group,
or that the actual complement equals Artin H4.

The local marking is also supplied. The two ordinary cusp pairs are
simultaneously transported to (w,y) and (z,y), the (2,5) cusp pair to
(x,w), and the three node pairs to (x,z),(x,y),(z,w). The precise
actual access paths, within-pair corrections and simultaneous maps
are given in §5. They permit subsequent use of actual retained and
deleted subsets in a hypothetical Keller monodromy. No Keller map,
no all-degree exclusion, and no ordinary-cusp inequality at the
five-cusp are inferred in this note.

The inherited mechanism is the [audited moving-centre Rouché
method](planar_jc48_sep08_moving_tubes.md), including its
[independent audit](planar_jc48_sep08_moving_tubes_audit.md), together
with the marked local-cluster sidecar of the
[ordinary three-cusp certificate](planar_jc48_sep08_three_cusp_braid.md).
The literal polynomial, all paths, all clusters and the group
presentation are new and are separately checked here.

The named hostile is sigma=(123), tau=(345): it satisfies the
five-term braid but violates the ordinary half-support inequality.
The [boundary note](planar_jc48_sep08_three_cusp_boundary.md), §6,
records its exact scope. The corrected near miss is transporting D4
relations across a degeneration to a fifth cusp. A second cheap
hostile arose in root proposals: the nearest pair at the endpoint
of the node0 stem is not its collision pair. The exact complex-disk
test rejected that proposal. The final certificate identifies the
pair by the entire parameter disk and endpoint root assignment.

The five live concepts are the literal discriminant, common based
root transport, local collision clusters, positive free basis, and
actual retained-set marking. Root and orthogonal_returns independently
derived the reversible H4 algebra before the path freeze. Root is
exploring further H4 consequences independently; none is a dependency
of this certificate.

## 2. Exact singularity inventory

Put F(u,v)=Res_t(U(t)-u,V(t)-v). The source lists all rational
coefficients of F explicitly and reconstructs the resultant. It is
monic of degree four in v and has total degree six. Direct substitution
F(U(t),V(t)) vanishes. Let

    H(u)=(u-3)(2187u²-13962u+22139)
        =2187u³-20523u²+64025u-66417.

The complete discriminants are

    disc_t(U-u)=-(256/27)(3u-13)(3u-8)(3u+19),
    disc_v F=-(1048576/847288609443)
       *(3u-13)³*(3u-8)³*(3u+19)^5*H(u)².             (3)

The three roots of H are simple and avoid the three critical U-values.
The common derivative zeros are exactly t=-1,1,2. The fibre gcd at
each corresponding target is (t-e)², so that cusp has no additional
normalization preimage. The jets at e=1,2 are ordinary and the jet
at e=-1 vanishes. At e=-1, writing h=t+1 and u0=U(t)-U(-1), one has

    V(t)-V(-1)-(9/2)u0-(1/16)u0²=4h5+O(h6),

with all lower terms zero and u0 of order two. This is exactly a
(2,5) cusp. The other two are ordinary. Birationality is explicit:
the field degree divides gcd(4,6)=2. A nontrivial degree-two involution
would make both coordinates even about its fixed point. That point
would have to be the sole zero ordinary jet e=-1, but U'''(-1)=-40
is nonzero. Thus the field degree is one.

At a root of H all four U-branches are analytic and U is a local
coordinate. Order two of the vertical discriminant forces exactly
one pair of analytic V-values to have a simple difference. This
gives a transverse node and excludes a tangent pair or triple image.
The cusp exponents3,3,5 exhaust the discriminant at the other three
values. A singular image must come either from a critical parameter
or distinct parameters with equal image; either gives a multiple
vertical root and has been included in (3). Hence the inventory is
complete.

At infinity let z=1/t, X=U/V, Z=1/V. Exact cleared-numerator checks give

    X=z²+O(z³),       Z-X³=8z7+O(z8).

This is the single infinity branch of type (2,7). The independent
genus control is 10=2+1+1+3+3, agreeing with three finite nodes.

## 3. Six literal rational paths in one based gauge

The common basepoint is u*=4+3i. Each path has vertices

    u*, c+r, c+ir, c-r, c-ir, c+r, u*.                  (4)

The exact centres, radii and chronological crossing words are:

| Name | Centre c | r | Word |
|---|---|---|---|
| cusp_plus | 13/3 | 1/32 | [-1,-2,-3,2,2,2,3,2,1] |
| cusp_five | -19/3 | 1/32 | [3,2,1,3,3,3,3,3,-1,-2,-3] |
| cusp_two | 8/3 | 1/32 | [3,2,1,1,1,-2,-3] |
| node0 | 3077421/1048576 | 1/64 | [-1,3,-2,1,2,2,-1,2,-3,1] |
| node1 | 3 | 1/64 | [-1,3,-2,1,1,2,-3,1] |
| node2 | 904195/262144 | 1/32 | [-1,3,-2,3,3,2,-3,1] |

The last three centres are rational approximations, not claims that
the node projections are rational. Exact Rouché bounds in §4 locate
one H-root inside each node diamond. The larger first and last
node estimates used during proposal search are not proof inputs;
only the declared radii above occur in the final witnesses.

Each gzip witness contains a complete ordered rational subdivision
of the six edges and four rational complex root centres at every
vertex. The verifier checks edge membership and strict subdivision
order, four entries at every row, exact common initial centres/order,
and the same unordered final centres. There is a single actual-root
to-centre contraction at u*, shared by all six paths.

For z in C define A(z)=|Re z|+|Im z| and B(z)=max(|Re z|,|Im z|),
so B(z)<=|z|<=A(z). On a segment u(h)=u0+h du, 0<=h<=1, let
z_i(h)=z_i0+h dz_i. Expand the literal polynomial completely as

    F(u0+h du,z_i0+h dz_i+w)=sum_{j,k} C_jk h^j w^k.

For a positive rational radius R_i the accepted Rouché inequality is

    B(C_01)R_i > sum_{(j,k)!=(0,1)} A(C_jk)R_i^k.       (5)

All powers of h, including nonlinear residual terms, are included.
Since |h|<=1, (5) gives exactly one root in each moving disk. The
verifier also computes the exact minimum of the L-infinity distance
between every pair of affine centres on [0,1], considering endpoints,
coordinate zeros and |Re|=|Im| breakpoints. That minimum exceeds the
sum of the two radii, so the four disks are disjoint throughout.

At consecutive segments the disks are concentric at their common
vertex. Both contain their unique actual root, so the labelled root
branches glue. The actual configuration is isotopic to the polygonal
centre configuration through these disjoint disks, with coherent
endpoint contractions. This identifies the actual braid.

The projection used for crossing order is multiplication by 1+i/4
followed by real part. All endpoint orders, crossing times, adjacency
at each crossing and nonzero imaginary separation are checked with
rational arithmetic. Crossings are processed chronologically, with
the convention H_i(a,b)=(aba^-1,a) when the left strand passes below
in the positive half twist. This is representation transport; the
geometric map on loops with below-stem meridians is the inverse Artin
pullback. The two conventions are not interchanged silently.

## 4. Entire local disks and actual pair marking

To identify the local cusp/node pairs, a whole complex parameter
disk is paid in addition to its boundary braid. For each name use
|u-c|<=r, with r from (4). The source's literal LOCAL table gives
three rational affine holomorphic centres

    z_i(u)=z_i0+L_i(u-c),

one intended for the colliding pair and two for spectators. Expand

    F(c+h,z_i0+L_i h+w)=sum D_jk h^j w^k.

The exact full-disk tests are

    B(D_0m)R_i^m > sum_{(j,k)!=(0,m)} A(D_jk)r^j R_i^k,
    m=2 for the pair, m=1 for each spectator.            (6)

The uniform separation bound is

    B(z_i0-z_j0)-A(L_i-L_j)r > R_i+R_j.                 (7)

Thus the configuration is separated into a degree-two disk and two
degree-one disks over the entire parameter disk, not just sampled
points or the lasso boundary. The final vertical cluster radii are:

| Name | Pair radius | Spectator radii |
|---|---|---|
| cusp_plus | 1/64 | 1/64,1/64 |
| cusp_five | 1/64 | 1/64,1/64 |
| cusp_two | 1/32 | 1/64,1/64 |
| node0 | 5/32 | 1/64,1/64 |
| node1 | 7/64 | 1/64,1/64 |
| node2 | 15/64 | 1/64,1/64 |

The rational centres and slopes are explicit in the source; they are
rechecked against F by (6), not trusted as numerical root proposals.
At node0, choosing the nearest endpoint pair proposed the wrong local
cluster. The final centre-based choice, and the smaller node0/node1
parameter disks, pass the exact tests. No maximal-radius claim is made.

For a node, Taylor expansion of H and a linear Rouché comparison
prove exactly one root inside both radius r and radius r/2. The
smaller disk lies inside the diamond (4). At each cusp, a constant
Rouché comparison proves H has no root in its parameter disk. All
other cusp projections are excluded by direct rational distances.
Consequently each path encircles exactly its declared singular value.

At the end of the incoming stem, tiny individually isolated root
disks assign the transported roots to the three clusters in counts
2,1,1. Exact projected-strip separation places the two colliding
punctures in adjacent positions of the transported below-stem meridian
basis, with spectators outside the local disk's strip. This supplies
the actual local distinguished pair and its common access path.
The actual incoming stem words are:

| Name | Incoming stem | Pair positions, one based |
|---|---|---|
| cusp_plus | [-1,-2,-3] | 2,3 |
| cusp_five | [3,2,1] | 3,4 |
| cusp_two | [3,2,-1] | 1,2 |
| node0 | [-1,3,-2,1] | 2,3 |
| node1 | [-1,3,-2] | 1,2 |
| node2 | [-1,3,-2,-3] | 3,4 |

These local pairs, their singularity types and the exact global words
are separate certified obligations. A formal word decomposition alone
would not have identified the actual retained-page pairs.

## 5. Reversible H4 algebra and retention-safe maps

Use a,b,c,d for the original four fibre meridians. Each full word is
P,core,P^-1, with these exact free-prefix pairs:

| Name | Formal prefix P | Core | Core pair |
|---|---|---|---|
| cusp_plus | [-1,-2,-3] | H2³ | (c,d) |
| cusp_five | [3,2,1] | H3^5 | (b,c) |
| cusp_two | [3,2] | H1³ | (a,bcdc^-1b^-1) |
| node0 | [-1,3,-2,1] | H2² | (b,(cdc^-1)^-1b^-1ab(cdc^-1)) |
| node1 | [-1,3,-2] | H1² | (b,cdc^-1) |
| node2 | [-1,3,-2] | H3² | ((cdc^-1)^-1b^-1ab(cdc^-1),c) |

A Hurwitz square fixes its pair iff the two letters commute. For odd
m=3 or5, H^m fixes its pair iff the alternating words of length m
are equal. One direct check uses product invariance and the second
coordinate: with g=(ab)^((m-1)/2), that coordinate is gag^-1. Its
equality to b is precisely ga=bg, the stated alternating relation;
the first coordinate then follows from the fixed pair product.

Now make the free Hurwitz basis change P0=[-1,3,-2], defining

    x=b,
    y=cdc^-1,
    z=y^-1b^-1ab y,
    w=c.                                               (8)

This is an automorphism of the free fibre group. The exact inverse is

    a=xyz y^-1x^-1,       b=x,
    c=w,                 d=w^-1yw.                     (9)

Every new generator is a positive conjugate of an original curve
meridian. The change of basis is a fibre meridian-basis change, not
an asserted loop in the u-base. It preserves generation and positivity.

The three node pairs become directly (x,z),(x,y),(z,w). The first
cusp pair is (w,w^-1yw), simultaneously conjugated by w to (w,y).
The fifth-cusp pair is (x,w). The remaining cusp pair is

    ((xy)z(xy)^-1, (xy)y(xy)^-1),

simultaneously conjugated by (xy)^-1 to (z,y). These are exact free
identities; no node relation is needed to obtain them. They give
exactly (2), and every step reverses. Thus the group defined by the
six fixed-word relations is precisely the displayed Artin H4
presentation in the stated positive basis.

Compare the actual stems in §4 to this formal table. Only cusp_two
and node2 differ: each actual stem adds one inverse half twist
inside its colliding pair. For an actual pair and retained subsets,

    (sigma,tau) -> (tau,tau^-1 sigma tau),
    (A,B)       -> (B,tau^-1 A).                        (10)

Because B is pointwise tau-fixed, both |A intersect B| and the
deleted-overlap cardinality are unchanged. The joint local subgroup
is unchanged. These are the already paid local distinguished-basis
changes inside the certified cluster. Any cusp injection should be
applied to the directly accessed actual pair before transporting its
count. No re-access equation is inferred from algebra alone after
a basis change. The fifth cusp needs no inside-pair correction.

The simultaneous maps following (9) preserve actual subset cardinalities,
joint fixed counts, moved support intersections and actual deleted
overlaps. Thus the theorem retains the marked local data needed for
later Keller analysis, without identifying distinct singularities'
retained sets or silently replacing them by full inertia fixed sets.

## 6. From root braids to the actual affine complement

Monicity in v gives a vertical fibre consisting of four punctures
and a continuous outside-root section over the regular projection
base. Over a compact disk containing all six chosen lassos and their
fillings one may choose a constant value v* beyond a uniform root
bound. This section extends over all critical u-values in that disk
and stays outside C. Each base lasso is therefore null in the actual
complement when based through that section.

Van Kampen then imposes the *fixed tuple* relation for each actual
Hurwitz word, not merely equality up to simultaneous conjugation.
The same fibre-and-section argument gives a surjection from the
four-punctured vertical fibre group to pi1(A2\C): the regular-base
loops become null after critical fibres are restored. Thus the four
positive fibre meridians generate the actual complement. It suffices
that each of the six certified loops supplies a valid relation; no
free-basis claim about their stems is needed.

The common basepoint, identical rational starting centres and coherent
root contractions ensure all six necessary relations use the same
four initial meridians. Sections3–5 consequently prove the actual
marked quotient assertion of §1.

The map preserves a concrete finite root configuration and marked
local pairs. It loses possible additional global relations and any
polynomial source of a hypothetical covering. A transitive finite
representation of (2), if later found, would not be a Keller source.
Conversely a proved obstruction to all relevant representations may
consume this theorem only with its actual retained-sheet hypotheses.

## 7. Exact universe, reproduction and independent audit

Source: [planar_jc48_sep08_mixed_cusp_braid.py](../../04-computation/planar_jc48_sep08_mixed_cusp_braid.py).
Output: [planar_jc48_sep08_mixed_cusp_braid.out](planar_jc48_sep08_mixed_cusp_braid.out).
The six matching `_*_certificate.json.gz` files retain all rational
subdivision vertices and root centres; the output records their raw
and compressed byte counts and hashes.

    python3 04-computation/planar_jc48_sep08_mixed_cusp_braid.py
    python3 -O 04-computation/planar_jc48_sep08_mixed_cusp_braid.py

The default consumer uses no floating-point root solver. It reconstructs
the literal resultant, discriminants and jets with exact SymPy
arithmetic, checks the free-word algebra, then verifies every stored
rational segment, crossing and entire complex-parameter cluster.
Always-active checks include the coherent endpoint contraction,
all six common starting fibres, and the complete local marking.

The inherited moving-method controls remain present: nonlinear
coefficient residuals cannot be discarded, exact interior separation
minima are retained, and an arbitrarily long common root translation
cancels correctly. The `--produce NAME` and `--scout STEPS` modes use
floating-point proposals; only an accepted rational witness and its
default exact replay are mathematical evidence.

Normal and optimized runs both pass 17,100 always-active gates over
1,517 moving segments, with byte-identical output. Source, output and
all six witnesses are frozen. The [independent complete audit](planar_jc48_sep08_mixed_cusp_braid_audit.md)
accepts every raw edge and crossing, all whole-disk clusters and endpoint
assignments, the free H4 basis and the actual pair transports. It retains
the generation/quotient-only scope and uses no ordinary bound at exponent five.

| Artifact | Bytes | SHA256 |
|---|---:|---|
| Source `.py` | 21,702 | `b6f884b7bf5871bbf3ba1d12bc94a1f1344c0add07f9832770a7d53ad536a0ff` |
| Output `.out` | 5,138 | `94fc466fea4d6bcad8d16c11f0e38a5191addb407533668cb8c9d5684fa883a8` |

Every witness uses the results-directory prefix
`planar_jc48_sep08_mixed_cusp_braid_` and suffix
`_certificate.json.gz`. Its two hashes pin both the stored compressed
file and the full literal decompressed JSON, including every root row:

| Name | Segments | Gzip bytes / SHA256 | Raw bytes / SHA256 |
|---|---:|---|---|
| cusp_plus | 144 | 9,828 / `437f88f34f79d41cac1058db65cd33e8fd7299c00e3f2aba862fab7c50be503d` | 37,461 / `130735471f691aec805ff5fe98b39ed751cc2a911d2a480f16979cf48eae0fc3` |
| cusp_five | 891 | 62,110 / `8b9a6e10554dd20fe382df44eeaf16df4af939312db9fbdcefceb341eb8baef4` | 241,772 / `9fc963f7b0b290ef6315e74d94d03b3f3a0591b7e54d83a0c4b8ca54e528e107` |
| cusp_two | 151 | 9,625 / `bcff0e7ec660bab00365806ebcff35c7d364451c5b9589787fce42d5c1e40c84` | 38,932 / `606536beaeb6301ad214ad0b63b1958c90b948fbd9f604582af67ec94233ad95` |
| node0 | 115 | 8,091 / `ed262265199c118e9fdee078fd8c026377bc82260df43bfef8e62fe99e77c278` | 30,614 / `8cb17de07488147a1a60b6a88f8a2107d3eabcb2e8cd9a427b31da2dec83b018` |
| node1 | 116 | 7,905 / `a1d90904c556445aaeeedf9719552c51875574ed2559c9f6818c72a4dc01381b` | 30,019 / `3e1f0ae0f5a261c964453b522bf44c516afa9e7028b7480f937b7ccf9976d674` |
| node2 | 100 | 6,852 / `71764a003a8e159f6ca395baa1616e91fade4ea85fde54a9a8fa10aca23e1a3a` | 26,463 / `2dab02efcc3b318e85e393999410528eddbf0a061e3aef29585a9400fe52836c` |
