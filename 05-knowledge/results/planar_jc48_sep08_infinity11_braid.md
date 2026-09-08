# Five common-access loops for the infinity-eleven two-cusp sextic

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.** All five
rational witnesses pass independent normal and optimized exact replays.
The [full bundle audit](planar_jc48_sep08_infinity11_braid_audit.md)
accepts the actual geometry, moving-tube realization, common access gauge,
arbitrary-group elimination, and whole-class consumer.

## 1. Actual target, inherited mechanism, and conclusion

The target is the actual polynomially parametrized curve
\[
 U=t^4-\frac83t^3-2t^2+8t,\qquad
 V=t^6-4t^5-\frac32t^4+\frac{548}{27}t^3-\frac{368}{9}t.
\tag{1}
\]
The complete geometry is reconstructed in
[the two-cusp family note](planar_jc48_sep08_two_cusp_family.md), §6,
and its exact source: two ordinary finite cusps, exactly three ordinary
nodes, birational normalization by \(\mathbb A^1\), and infinity
branch \((2,11)\). Its finite critical projection values are all
retained here. The result is
\[
 \pi_1(\mathbb C^2\setminus C)=\mathbb Z,
\tag{2}
\]
and consequently exclusion of C as the whole irreducible nonproperness
curve of a nonautomorphic polynomial Keller map. No claim about a
component inside a larger support is intended.

The closest proved mechanism is the
[literal two-cusp four-loop certificate](planar_jc48_sep08_two_cusp_braid.md),
whose common-access relations force four positive meridians to coincide.
The hostile is the exact local/Euler/three-generator A6 passport in
[the two-cusp passport](planar_jc48_sep08_two_cusp_passport.md): it pays
local singularity data but has no representation of the actual global
complement. The corrected near miss is transferring the type-nine group
to type eleven merely by specialization. This proof requires new actual
loops. The least-used sidecar is the
[moving-centre Rouché method](planar_jc48_sep08_moving_tubes.md), which
retains common root motion before taking absolute-value bounds.

The five live concepts are the actual quartic fibre, its six projection
critical values, rational root tubes, marked Hurwitz transport, and
arbitrary-group elimination. The map from a certified root path to a
word preserves the marked four-strand braid, not just endpoint
permutation. The passage from the word to the group preserves common
access paths and positive meridians. It discards detailed root motion,
whose sidecar is the frozen rational witness. No finite permutation
census is used as a surrogate for an arbitrary-group proof.

## 2. Literal quartic and explicit base loops

Let
\[
 F(u,v)=\operatorname{Res}_t(U(t)-u,V(t)-v).
\]
The source contains every rational coefficient of F explicitly and
independently compares it with this resultant. It is monic of degree
four in v, and substitution of (1) gives zero. Its full discriminant is
\[
\begin{split}
 \operatorname{disc}_vF={}&
 -\frac{104857600000000}{328256967394537077627}
 (3u-13)^3(3u-8)(3u+19)^3H(u)^2,\\
 H(u)={}&59049u^3+1526283u^2+13980187u+44759481.
\end{split}\tag{3}
\]
The three roots of H are distinct and avoid the other three values.
The smooth fold is \(u=8/3\); the cusp projections are \(13/3\)
and \(-19/3\). No projected cusp and node have been amalgamated.

Use common basepoint \(u_*=4+3i\) and radius \(r=1/32\). For
any centre c below the declared loop has vertices
\[
 u_*,\ c+r,\ c+ir,\ c-r,\ c-ir,\ c+r,\ u_*.
\tag{4}
\]
All vertices are exact Gaussian rationals. The three cubic-node
centres are rational approximations, not alleged exact critical values.
Only the actual root-path certification is needed: the proof does not
require that a loop surround exactly one critical value.

| Name | Exact centre | Exact reduced word |
|---|---|---|
| smooth | \(8/3\) | \([1,2,-1]\) |
| cusp plus | \(13/3\) | \([1,1,1]\) |
| cusp minus | \(-19/3\) | \([2,1,3,3,3,-1,-2]\) |
| real node | \(-8835419/2^{20}\) | \([2,1,3,3,1,1,-3,-3,-1,-2]\) |
| lower node | \((-9133949-3934729i)/2^{20}\) | \([2,-3,2,2,3,-2]\) |

The upper-node centre \((-9133949+3934729i)/2^{20}\) and its
word \([2,1,3,2,2,-3,-1,-2]\) remain in the scout universe but
are not required for the five-word proof. Meshes 256 and 512
agree on all six words; this is a cheap consistency check only.

## 3. Exact moving disks and the polygonal braid

For a Gaussian rational z, put
\[
 A(z)=|\Re z|+|\Im z|,\qquad
 B(z)=\max(|\Re z|,|\Im z|).
\]
Then \(B(z)\le |z|\le A(z)\). On each rational base segment,
write \(u(h)=u_0+h\Delta u\), \(0\le h\le1\), and use four
rational affine centre strands \(z_i(h)=z_{i0}+h\Delta z_i\).
Expand exactly
\[
 F(u_0+h\Delta u,z_{i0}+h\Delta z_i+w)
     =\sum_{j,k}C^{(i)}_{jk}h^jw^k.
\tag{5}
\]
All terms with the same pair (j,k) are combined before bounding. This
retains the cancellation of shared parameter/root motion. The sufficient
strict Rouché condition for radius \(r_i>0\) is
\[
 B(C^{(i)}_{01})r_i>
 \sum_{(j,k)\ne(0,1)}A(C^{(i)}_{jk})r_i^k.
\tag{6}
\]
For every h, (6) compares F on \(|w|=r_i\) with the same nonzero
linear polynomial \(C^{(i)}_{01}w\). Thus there is exactly one root
in the moving disk for every h. This is a uniform rational certificate,
not a bound on a first derivative alone.

Disk separation is paid uniformly as well. For each i,j compute exactly
\[
 d_{ij}=\min_{0\le h\le1}
 B(z_{i0}-z_{j0}+h(\Delta z_i-\Delta z_j)).
\tag{7}
\]
The max of the four affine functions giving this B-norm is piecewise
linear. Its minimum is attained at an endpoint or where the real part,
imaginary part, their difference, or their sum vanishes. The source
tests every such rational candidate in [0,1]. It requires \(d_{ij}>0\)
and uses
\[
 r_i=\frac1{16}\min_{j\ne i}d_{ij}.
\tag{8}
\]
Consequently \(r_i+r_j<d_{ij}\le|z_i(h)-z_j(h)|\), so all four
moving disks are disjoint.

Adjacent segments have the identical labelled endpoint centres. Their
endpoint disks may have different radii, but they are concentric and
nested. Each disk contains exactly one root of the same polynomial, so
the root labels agree. At the common basepoint all five witnesses use
identical initial labelled centres. Initial tiny-disk isolation supplies
a common endpoint displacement, and the final unordered centre set is
exactly the initial set. These facts give one common marked root fibre.

For each h the four roots can be joined to their four centres inside
their disjoint disks. Straight interpolation cannot cause a collision.
This homotopy is continuous through the shared endpoints and identifies
the actual root braid with the polygonal centre braid, using the same
initial and final short displacements for every loop. The verifier then
computes every signed crossing exactly after the fixed projection
\(z\mapsto(1+i/4)z\). Endpoint real-coordinate ties, coincident
crossing times, nonadjacent crossings, and zero imaginary crossing
separation are all rejected. This retains the whole braid, not merely
its root permutation.

The source includes independent symbolic substitution for (5), an
arbitrarily long common-translation positive control, an
endpoint-separated swapped-label collision hostile, and a curvature
hostile with vanishing first-order residual. The latter fails (6), as
it should: keeping only the linear motion would be insufficient.

## 4. Common access paths and arbitrary-group elimination

Use positive meridians in the common reference fibre, ordered by the
projection above, with stems below. Our chronological representation
transport convention is
\[
 H_i(a,b)=(aba^{-1},a),\qquad H_i^{-1}(a,b)=(b,b^{-1}ab).
\tag{9}
\]
The geometric transport on loops is the inverse-Artin pullback;
representation transport is contravariant and gives (9). The signed
crossing convention is the same as in the audited predecessor. Each
word is applied chronologically, so a prefix P, core K, and inverse
prefix \(P^{-1}\) impose the fixedness of K on the P-transformed
tuple. No arbitrary simultaneous conjugation is inserted between loops.

For a monic polynomial in v, the vertical four-meridian group surjects
onto the actual affine complement group. Over the compact union of the
five loop fillings choose a constant \(v_*\) beyond every root bound.
The section \((u,v_*)\) contracts every one of these base loops inside
the actual complement. Therefore their marked Hurwitz actions fix the
same tuple \((a,b,c,d)\) in the actual group. The five loops need not
generate the critical-value base group for these necessary relations.

The following implication holds in every group, without a transitivity
assumption. A half-twist core fixes a pair exactly when its entries are
equal; a square core exactly when they commute; a cube core exactly
when they satisfy the braid relation.

1. The smooth word has prefix [1], core [2]. Its relevant prefix pair
   is (a,c), hence \(a=c\).
2. The cusp-plus word gives \(aba=bab\).
3. The real-node word has prefix [2,1,3,3], core [1,1]. Its relevant
   pair is \((abcb^{-1}a^{-1},a)\). With c=a and the braid
   relation, its first entry is b. Thus \([a,b]=1\), and the
   braid relation now gives \(a=b=c\).
4. The lower-node word has prefix [2,-3], core [2,2]. Its relevant
   pair is \((bcb^{-1},d)\), now (a,d), so \([a,d]=1\).
5. The cusp-minus word has prefix [2,1], core [3,3,3]. Its relevant
   pair is (b,d), so \(bdb=dbd\). With a=b and commutation,
   this gives \(d=a\).

Thus all four generators coincide. The source checks the literal word
decompositions and prefix free words. The reverse control is immediate:
every equal tuple is fixed by all five words. This group implication
was independently checked by `orthogonal_returns` before any rational
path production; that audit does not itself certify the geometric words.

The certified witnesses therefore make the actual group cyclic. A reduced
irreducible equation of C maps a positive meridian to winding one in
\(\mathbb C^*\), so the cyclic group is infinite, giving (2).
The whole-support Keller exclusion then follows from the audited cyclic
consumer: a generic positive meridian retains at least one sheet, while
transitive cyclic degree-d monodromy has no fixed sheet for d>1.

## 5. Frozen finite-exact witnesses and independent audit

The [standalone source](../../04-computation/planar_jc48_sep08_infinity11_braid.py)
contains the literal coefficients, all six declared loops, the moving
root-tube and exact-crossing verifiers, and the group controls. Its
numerical proposal routine is used only by `--produce` and `--scout`;
a default replay consumes stored rational data and proves every
accepted path segment with exact arithmetic.

The moving method has an independent full analytic/source audit. The five
new witnesses contain **736 segments** in total. Both full exact replays
pass **6,793 always-active gates**, with byte-identical frozen output.
The fixed polynomial, common reference fibre, all thirty straight base
edges, and every accepted moving tube are rechecked by this new source.
The independent audit of the complete actual-curve bundle passes,
including fresh normal and optimized replays. The earlier frozen family
bundle has no backward dependency on this result.

| Path | Segments | Compressed bytes | Raw bytes |
|---|---:|---:|---:|
| smooth | 102 | 6923 | 26744 |
| cusp_plus | 148 | 9885 | 38692 |
| cusp_minus | 172 | 11857 | 45075 |
| node_real | 169 | 12621 | 46034 |
| node_lower | 145 | 8592 | 41259 |

Reproduction (no numerical root solver is used by these replays):

```bash
python3 04-computation/planar_jc48_sep08_infinity11_braid.py
python3 -O 04-computation/planar_jc48_sep08_infinity11_braid.py
```

The [frozen output](planar_jc48_sep08_infinity11_braid.out) contains the individual path gate counts,
raw and compressed hashes, and exact crossing words. Frozen file pins:

* Source SHA256: `2df6b5892ea04d0a771936c9ca21838bc727d131804cdefeae03bcfc09e2b26f`.
* Output SHA256: `ca8153a95aa81dcfd06991d30bdabe03c11898fd7dc36ad28f55a394161f6549` (2461 bytes).
* smooth witness: [planar_jc48_sep08_infinity11_braid_smooth_certificate.json.gz](planar_jc48_sep08_infinity11_braid_smooth_certificate.json.gz), compressed SHA256 `821f52e0199693af3222e6fe8cc506c6f9341194bcc38b2c079c9a8ae999e4a0`; raw SHA256 `3e773ba9e2a041f8a51a5a45acc2f92e7d7d92f3d12fae57a443b0fe305b79e9`.
* cusp_plus witness: [planar_jc48_sep08_infinity11_braid_cusp_plus_certificate.json.gz](planar_jc48_sep08_infinity11_braid_cusp_plus_certificate.json.gz), compressed SHA256 `31c8b46dfd0775e80790ab0188602d020579fe990b857c817418a58feee94cea`; raw SHA256 `8e9ca689d19d2ba241e71147cf51b1d3d6028425a61afd8bc63c66980a8c8c92`.
* cusp_minus witness: [planar_jc48_sep08_infinity11_braid_cusp_minus_certificate.json.gz](planar_jc48_sep08_infinity11_braid_cusp_minus_certificate.json.gz), compressed SHA256 `05332d61d143522b7ea2184b7dc964201dc66c91e5c3da155ab76d4cfb56678e`; raw SHA256 `7c1dd3a237dc7db81683e0d7aca2a9f723de45f242e730ef14d1bc650317529f`.
* node_real witness: [planar_jc48_sep08_infinity11_braid_node_real_certificate.json.gz](planar_jc48_sep08_infinity11_braid_node_real_certificate.json.gz), compressed SHA256 `f1097a644150953f07a40cf57b06e437a820358bbbfbc5707ce503d83f986e22`; raw SHA256 `a1cb6797e007708a5cae21c9d03df76c6125cbfe46d0257e471d02ba6ce7a1b8`.
* node_lower witness: [planar_jc48_sep08_infinity11_braid_node_lower_certificate.json.gz](planar_jc48_sep08_infinity11_braid_node_lower_certificate.json.gz), compressed SHA256 `e8051f7e47493fbe5ccf47a0693f8d849f4213df4a14ab848823343159380860`; raw SHA256 `b3828eb89e8d9bccc2749fe41c29bf0b0012a2c43bcbb963444dbe6e48f2178f`.


## 6. The inherited whole-class consumer

**PROVED, with the same independent audit.**
The proved [two-cusp family theorem](planar_jc48_sep08_two_cusp_family.md),
§§2–4, exhausts every birational polynomial normalization of degrees
(4,6), with exactly two ordinary finite cusps and otherwise only
ordinary nodes, into the three nonempty good loci \(G_7,G_9,G_{11}\).
It proves that each good locus is path connected and that proper marked
resolution transports the actual affine complement, preserving positive
meridians. Its §6 identifies (1) as the good point s=2 in \(G_{11}\).

Consequently the cyclic group established here for (1) transports to
**every member of \(G_{11}\)**. The earlier theorem already supplies
\(\pi_1=\mathbb Z\) for \(G_7\) and \(G_9\). The entire declared two-ordinary-cusp
(4,6) class therefore has infinite cyclic affine complement and is
excluded as whole irreducible Keller nonproperness support in every
mapping degree.

The transport is within \(G_{11}\), whose marked infinity type is
constant. It does not assert isotopy across its boundary into \(G_9\),
and it does not discard the repeated projection-critical cases s=±1.
Curves with extra finite singularities, other normalization degree
pairs, or occurrence merely as one component of a larger support are
outside this consumer. JC(2) remains OPEN.

Final independent audit: 16,194 bytes, SHA256
`01443bf46c6236cccbc33a9902c983f7935c7eb78beec228f60abb906b91ee39`.
All source, output, and witness pins above are unchanged by status promotion.
