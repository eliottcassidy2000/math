# Independent audit of the actual mixed-cusp H4 supplier

**Status: INDEPENDENT ANALYTIC / SOURCE / SIX-WITNESS AUDIT PASS.**
The [mixed-cusp braid proof](planar_jc48_sep08_mixed_cusp_braid.md), its complete standalone implementation, all six new rational witnesses, and both independent execution modes pass. No mathematical, source, or marking correction is required. The coordinating agent owns status promotion.

The accepted conclusion is an actual **marked quotient of Artin H4** for the stated literal affine complement. The group defined by the six certified fixed-word relations is exactly that Artin presentation in a changed free fibre basis. The actual complement may have further relations. Neither a free basis of the punctured projection base nor equality of its complement group with Artin H4 is asserted. No Coxeter involutions, Keller realization, or Keller exclusion are inferred here.

## 1. Exact geometry and inherited proof boundary

The literal parametrization is

    U=t⁴−(8/3)t³−2t²+8t,
    V=t⁶−(3/2)t⁴−16t³+48t.

The [boundary geometry](planar_jc48_sep08_three_cusp_boundary.md) and its [independent audit](planar_jc48_sep08_three_cusp_boundary_audit.md) already establish this exact good control. The new source nevertheless reconstructs its resultant from U and V and compares every rational coefficient with its literal table. The resulting F(u,v) is monic of v-degree four and total degree six, and F(U,V)=0 exactly.

The common critical parameters are precisely −1,1,2. Each critical target has fibre gcd (t−e)², so no additional normalization preimage shares it. The ordinary jets at 1 and 2 are nonzero. At −1 the exact tangent and even-term removal gives first odd term 4h⁵, with U of order two; this is a (2,5) branch. Since U'''(−1)≠0, the sole possible quadratic deck-reflection centre does not make U even, proving birationality as in the audited boundary proof.

Write

    H(u)=(u−3)(2187u²−13962u+22139).

Its three roots are simple and avoid the three critical U-values. The two full discriminants are

    disc_t(U−u)=−(256/27)(3u−13)(3u−8)(3u+19),
    disc_v F=−(1048576/847288609443)
      (3u−13)³(3u−8)³(3u+19)⁵ H(u)².

At each root of H, all parameter branches are analytic in u, and discriminant order two gives exactly one transverse double image. The cusp exponents three, three and five exhaust the remaining critical projections. Thus there are three nodes, not additional folds or unrecorded singularities. The cleared infinity identity Z−X³=8z⁷+O(z⁸), X=z²+O(z³), proves the unique infinity branch has type (2,7). Its delta invariant three agrees with the rational-sextic ledger 10=2+1+1+3+3.

The inherited moving-tube theorem and marked-cluster method are used only as general verification mechanisms. The polynomial, all six paths, every centre and radius, and every local-pair assignment are newly checked. No D4 relation is transported across a cusp degeneration.

## 2. Complete rational path universe and based gauge

All paths start at u*=4+3i and have the six edges through

    u*, c+r, c+ir, c−r, c−ir, c+r, u*.

The exact declared data are:

| Name | c | r | Chronological full word |
|---|---|---|---|
| cusp_plus | 13/3 | 1/32 | −1,−2,−3,2,2,2,3,2,1 |
| cusp_five | −19/3 | 1/32 | 3,2,1,3,3,3,3,3,−1,−2,−3 |
| cusp_two | 8/3 | 1/32 | 3,2,1,1,1,−2,−3 |
| node0 | 3077421/1048576 | 1/64 | −1,3,−2,1,2,2,−1,2,−3,1 |
| node1 | 3 | 1/64 | −1,3,−2,1,1,2,−3,1 |
| node2 | 904195/262144 | 1/32 | −1,3,−2,3,3,2,−3,1 |

The node centres are rational isolating centres, not assumed exact critical values. I independently parsed every decompressed witness. Each row has four rational complex centres; the u and centre-row counts agree; all six declared edges are traversed in exact strict subdivision order; every path closes at the stated basepoint; all six initial labelled centre tuples coincide; and each final tuple is the same unordered tuple. All 1,517 segments are included. No production/scout mode was used in the audit.

The common initial tiny isolating disks determine one contraction from the actual base configuration to the rational centre configuration. The final unordered equality and unique-root isolation make the endpoint contraction the same one with the correct label permutation. Thus the six resulting loops and the open stems have one coherent based fibre gauge. They are not six separately conjugated braid measurements.

## 3. Moving Rouché proof, crossings and witness validity

For a rational complex number z, let A(z)=|Re z|+|Im z| and B(z)=max(|Re z|,|Im z|). These are respectively an upper and lower bound for its Euclidean modulus. On a segment, write u(h)=u0+hΔu and zi(h)=zi0+hΔzi for 0≤h≤1. The source expands the complete literal polynomial as

    F(u0+hΔu,zi0+hΔzi+w)=Σ Cjk h^j w^k.

It retains all powers of h, including nonlinear constant-in-w residuals. If

    B(C01)R > Σ_(j,k)≠(0,1) A(Cjk)R^k,

Rouché against C01w gives exactly one root in the corresponding moving disk for every h. The strict inequality also proves its dominant coefficient is nonzero. Monicity and the four disjoint one-root disks account for all roots of F.

For two affine centres, the exact minimum on [0,1] of max(|Re|,|Im|) occurs among the endpoints, the coordinate zeros, and the Re=±Im breakpoints. This is the complete breakpoint list of a piecewise-linear maximum of four affine functions. The source uses that minimum, not just endpoint distances. Its radii are positive fractions of the pair minima and their sums are strictly smaller than every separation bound.

Consecutive segment disks are concentric at the shared centre value. Both unique-root disks contain the same root, so labels glue even if their radii differ. Contracting each actual root to its affine centre inside these disjoint convex disks produces a configuration isotopy. These contractions agree at the shared vertices and at the base gauge, proving that the polygonal centre braid is the actual braid.

The exact projection is multiplication by 1+i/4 followed by real part. The crossing verifier checks every endpoint pair for distinct projected coordinates, solves the exact rational crossing time for every sign-changing difference, excludes simultaneous crossings, checks adjacency in the current order, and checks nonzero imaginary separation. It updates the order after each crossing and cancels only adjacent inverse letters. I independently reconstructed every full word and incoming stem from the raw rational rows, with no call to the producer's crossing routine.

The convention is chronological representation transport

    H+(a,b)=(aba⁻¹,a),
    H−(a,b)=(b,b⁻¹ab),

where a positive geometric half twist sends the left strand below the right one. This is the inverse-Artin geometric pullback convention for the chosen below-stem meridians, already paid by the inherited method audit. It is used coherently by both crossing extraction and all free-word actions here. The positive core exponents 3,5,3,2,2,2 agree with the counterclockwise local cusp/node models. No push/pull reversal is introduced at the marking step.

The source's method controls are meaningful: an interior affine-distance minimum is retained; a swapped-endpoint collision is rejected; the nonlinear residual in w−h²+h³ is not discarded; and an arbitrarily long common translation of four simple roots cancels exactly and remains certifiable. A separate SymPy substitution checks its affine coefficient routine against the original monomials. The full default replays verify every moving segment rather than trusting a stored word.

## 4. Independent verification of all six whole complex disks

This is the load-bearing local marking sidecar. Over each entire complex disk |u−c|≤r, the source supplies three rational affine holomorphic centres zi(u)=zi0+Li(u−c), one for a pair and two for spectators. For degree m=2 or 1 respectively, it checks

    B(D0m)R^m > Σ_(j,k)≠(0,m) A(Djk)r^jR^k

in the complete expansion F(c+h,zi0+Lih+w). Unlike a check just on a real segment or the diamond boundary, this bounds every complex h with |h|≤r. It gives exactly two or one roots in the entire moving cluster disk. The uniform separation

    B(zi0−zj0)−A(Li−Lj)r > Ri+Rj

keeps the three clusters disjoint throughout that complex disk, including its singular centre.

I wrote a separate exact Fraction verifier which reads only the literal coefficient and centre tables as data. It expands each original monomial directly by the multinomial formula:

    (c+h)^i(v+Lh+w)^l,

choosing a powers of h in the first factor and b powers of h and k powers of w in the second. Its coefficient contribution is

    coefficient · binom(i,a) binom(l,k) binom(l−k,b)
      · c^(i−a) v^(l−k−b) L^b

at h^(a+b)w^k. Thus it does not reuse the producer's base-Taylor or moving-Taylor implementation. This independently passed **all 18 whole-disk dominance comparisons and all 18 cluster separations**.

The pair radii are respectively 1/64,1/64,1/32,5/32,7/64,15/64, and every spectator radius is 1/64. In particular node0 and node1 use their final parameter radius 1/64, not an earlier larger proposal. Pair selection is certified by the whole disk; no nearest-pair heuristic is inherited.

The exact Taylor coefficients of the residual cubic H about each real rational c are

    H(c), 6561c²−41046c+64025,
    6561c−20523, 2187.

At each node centre the linear term dominates on both radius r and radius r/2. Hence there is exactly one node value in the full disk and in the smaller disk contained in the diamond. At each cusp centre a constant comparison excludes all H-roots from the full disk. Other cusp values are excluded by exact distances. I independently checked these nine root-count comparisons and, additionally, all fifteen pairwise parameter-disk separations. Therefore the six diamonds enclose six different, correctly typed singular values.

At the incoming endpoint c+r I also independently re-expanded F at each of the four stored individual centres and verified **24 tiny one-root isolations and 24 unique cluster assignments**. The counts are 2,1,1 in every case. The actual projection-order gaps exceed (5/4) times the tiny-radius sums, and both spectators lie outside the projected pair-disk strip by the same conservative norm bound. These are **30 further projection guards**. Since |1+i/4|<5/4, the estimates pay actual root positions, not merely centre positions.

The pair is therefore adjacent in the endpoint meridian basis and has a distinguished local disk and common access path with no intervening spectator. The independent assignments by original stored strand label are:

| Name | Cluster labels of the four stored strands | Actual pair positions, one based | Actual incoming stem |
|---|---|---|---|
| cusp_plus | 1,2,0,0 | 2,3 | −1,−2,−3 |
| cusp_five | 1,0,0,2 | 3,4 | 3,2,1 |
| cusp_two | 0,1,2,0 | 1,2 | 3,2,−1 |
| node0 | 0,0,1,2 | 2,3 | −1,3,−2,1 |
| node1 | 1,0,2,0 | 1,2 | −1,3,−2 |
| node2 | 0,1,0,2 | 3,4 | −1,3,−2,−3 |

Cluster label zero denotes the pair. A formal decomposition of a full braid word would not have supplied this actual information on its own.

## 5. Exact free-basis algebra and the H4 presentation

Let a,b,c,d be the original four fibre meridians. Applying the free Hurwitz automorphism P0=[−1,3,−2] gives

    x=b,  y=cdc⁻¹,
    z=y⁻¹b⁻¹aby,  w=c.

Its inverse is exactly

    a=xyz y⁻¹x⁻¹,
    b=x,  c=w,  d=w⁻¹yw.

I independently followed the three Hurwitz moves and substituted the inverse. All four new generators are conjugates of positive original meridians. This is a free **fibre** basis; the actual complement group is not asserted to be free. No claim that P0 is a loop of the u-base is needed.

The six full words decompose literally as prefix, positive core, inverse prefix. In the new basis their core pairs are

    (w,w⁻¹yw), (x,w),
    ((xy)z(xy)⁻¹,(xy)y(xy)⁻¹),
    (x,z), (x,y), (z,w).

Simultaneous conjugation by w on the first pair and by (xy)⁻¹ on the third gives

    (w,y), (x,w), (z,y), (x,z), (x,y), (z,w).

These identities hold in the free group before imposing any node relation. The source checks both inverse free-basis substitutions and every individual pair identity.

A Hurwitz square fixes its pair precisely when the letters commute. For odd m, product invariance and the second coordinate of H^m give the exact iff: with g=(ab)^((m−1)/2), the second coordinate is gag⁻¹; its equality to b is ga=bg, the length-m alternating relation. The first coordinate then follows from the fixed product. This applies to both m=3 and m=5 without replacing one by the other.

Consequently the six fixed-word relators become exactly

    yzy=zyz, ywy=wyw, xwxwx=wxwxw,
    [x,y]=[x,z]=[z,w]=1.

All transformations are reversible, so the group defined by these six relations is precisely Artin H4, with chain z–3–y–3–w–5–x. No involution relation has been used.

## 6. Actual inside-pair transport and retained subsets

Comparing actual incoming stems with the formal prefix table, only cusp_two and node2 differ. Each actual stem adds one inverse Hurwitz move **inside its certified colliding pair**. The full-disk and endpoint-strip certificates make this a valid local distinguished-basis change, rather than an assumed loop in the projection base.

For a pair (σ,τ) and its actual retained subsets A⊆Fix(σ), B⊆Fix(τ), that inverse move transports

    (σ,τ)→(τ,τ⁻¹στ),
    (A,B)→(B,τ⁻¹A).

The set B is pointwise τ-fixed, so applying τ proves

    |B∩τ⁻¹A|=|A∩B|.

Individual sizes are retained, hence deleted-overlap cardinality is retained as well. The joint subgroup is unchanged. The subsequent simultaneous conjugations preserve all these counts, full joint fixed sets and moved-support intersections. This proves the announced marking of the three actual cusp neighborhoods and three actual node neighborhoods.

One must apply any local cusp re-access/injection theorem to the directly accessed actual pair, then transport its proven count. No new re-access equation is deduced after rebasing from an arbitrary abstract pair of subsets. Nor are actual retained subsets replaced by full fixed sets or identified between different singularities. The exponent-five cusp already has the direct actual pair (x,w), with no inside-pair correction.

The finite hostile (σ,τ)=((123),(345)) from the audited boundary note still forbids importing the ordinary half-support estimate at this cusp. This supplier makes no such inference and leaves later representation constraints to separate proofs.

## 7. Actual complement generation and the precise quotient direction

Monicity gives a four-punctured vertical fibre over the regular projection base. On a compact disk containing the six paths and all critical values, a sufficiently large constant v* lies outside every fibre root. It supplies an outside-root section which extends across all six critical u-values. Thus each chosen base lasso, when lifted by that section, is null in the actual complement.

The fibre bundle and van Kampen consequently impose each certified **fixed tuple** relation, not merely simultaneous-conjugacy invariance. The common centre contraction and outside-root basepoint pay the common marking. This argument only needs each certified lasso to be an actual loop in the regular base; the tubes prove that condition all along every stem and boundary.

Generation is also paid independently of any claim about these particular stems. The regular-base bundle group is generated by the fibre group and lifts of a base free basis. Restoring critical fibres is surjective on fundamental groups, and the outside section kills those base generators. One may choose any base free basis inside the same large disk. Hence the four positive fibre meridians surject onto the actual affine complement group. It is unnecessary to prove that the six chosen geometric lassos themselves are such a free basis.

It follows that the actual complement is a quotient of the six-relator Artin H4 group. Additional relations are allowed. A finite representation of this quotient supplier would not by itself realize a Keller source. Conversely later exclusions may use these actual marked neighborhoods only with their own proved retention and geometric hypotheses. The current certificate does not supply or assume such a covering.

## 8. Independent execution and frozen artifact pins

Executed independently from `/tmp/math-wt-planar-jacobian-sep06`:

    python3 04-computation/planar_jc48_sep08_mixed_cusp_braid.py
    python3 -O 04-computation/planar_jc48_sep08_mixed_cusp_braid.py

Both default exact runs exited successfully and equal the frozen output byte for byte: **17,100 always-active gates, 1,517 moving segments, 5,138 output bytes**. Neither uses a floating-point root solver. Production and scout modes were not run. I read the entire implementation, including geometry, free algebra, method controls, edge-universe checks, tube bounds, crossing extraction, endpoint isolation and whole-disk marking.

The separate direct-monomial/edge/crossing reconstruction described above was run from a temporary independent verifier. It imports no mathematical implementation from the producer. Its checks covered every raw path row, all six common fibres, all 1,517 segments for crossing extraction, all six full/stem words, all 18 complex-disk expansions, their separations and every endpoint assignment. The persistent producer and six literal witnesses below reproduce the complete accepted certificate; the independent checks are additional verification, not a replacement certificate universe.

| Artifact | Bytes | SHA256 |
|---|---:|---|
| Primary `.md`, pre-promotion | 19,184 | `2ca9d2644531f38559187ef1d80e978a2416137328c825b03448f7d7472a5598` |
| Source `.py` | 21,702 | `b6f884b7bf5871bbf3ba1d12bc94a1f1344c0add07f9832770a7d53ad536a0ff` |
| Frozen `.out` | 5,138 | `94fc466fea4d6bcad8d16c11f0e38a5191addb407533668cb8c9d5684fa883a8` |

All six witness names have prefix `planar_jc48_sep08_mixed_cusp_braid_` and suffix `_certificate.json.gz` in this results directory. Both raw and compressed hashes were recomputed independently:

| Name | Segments | Gzip bytes / SHA256 | Raw bytes / SHA256 |
|---|---:|---|---|
| cusp_plus | 144 | 9,828 / `437f88f34f79d41cac1058db65cd33e8fd7299c00e3f2aba862fab7c50be503d` | 37,461 / `130735471f691aec805ff5fe98b39ed751cc2a911d2a480f16979cf48eae0fc3` |
| cusp_five | 891 | 62,110 / `8b9a6e10554dd20fe382df44eeaf16df4af939312db9fbdcefceb341eb8baef4` | 241,772 / `9fc963f7b0b290ef6315e74d94d03b3f3a0591b7e54d83a0c4b8ca54e528e107` |
| cusp_two | 151 | 9,625 / `bcff0e7ec660bab00365806ebcff35c7d364451c5b9589787fce42d5c1e40c84` | 38,932 / `606536beaeb6301ad214ad0b63b1958c90b948fbd9f604582af67ec94233ad95` |
| node0 | 115 | 8,091 / `ed262265199c118e9fdee078fd8c026377bc82260df43bfef8e62fe99e77c278` | 30,614 / `8cb17de07488147a1a60b6a88f8a2107d3eabcb2e8cd9a427b31da2dec83b018` |
| node1 | 116 | 7,905 / `a1d90904c556445aaeeedf9719552c51875574ed2559c9f6818c72a4dc01381b` | 30,019 / `3e1f0ae0f5a261c964453b522bf44c516afa9e7028b7480f937b7ccf9976d674` |
| node2 | 100 | 6,852 / `71764a003a8e159f6ca395baa1616e91fade4ea85fde54a9a8fa10aca23e1a3a` | 26,463 / `2dab02efcc3b318e85e393999410528eddbf0a061e3aef29585a9400fe52836c` |

**Final acceptance: PASS.** The literal geometry, every moving witness, whole complex-disk local marking, free H4 change of basis, actual inside-pair count transport, positive meridian generation and quotient direction are all accepted. There is no remaining correction and no implied ordinary-cusp estimate at exponent five.
