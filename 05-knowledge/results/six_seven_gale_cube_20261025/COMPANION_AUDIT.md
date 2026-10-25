# Independent geometry audit of the companion character and quotient packets

**Verdict: PASS in the stated scopes.** Audit by the geometry lane, 2026-10-25. No producer implementation is imported by the independent checkers. This audit does not promote numerical coincidence into an isomorphism, or finite checks into an all-order theorem.

Reviewed artifacts and SHA-256 (with the character clarification lineage below):

- `CHARACTER_K3_DECOMPOSITION.md`: `334ba589356c261cb9b4865723fbfe3c83350df28de1584f5de30eb98185cc55`.
- `COMBINATORIAL_BRIDGES.md`: `2089f482bffd76cbd5fc0c34ce8132748842817899db35f85e9dd9ed58e1d659`.
- `BUDGET_AND_GALE.md`: `c2039065eaccdba59c53d5655f7302b56d2fb8fa85cc254755128dc050e5fe7f`.

## 1. Character decomposition and the actual two-form pullbacks

The projective sign quotient has exactly the even-weight characters. The displayed local monomials are a complete free basis, including the cases in which the chart index belongs to the character set; their transition factors are (L_j/L_k)^(|S|/2), giving O(-|S|/2). The multiplication formula retains the required intersection radicands and has the correct denominator exponent.

Finite-flat duality changes O(-r) to O(r−3). The Stacks primary source was opened independently: Section49.15, Lemma49.15.2 and its proof explicitly give f_*f^!O=Hom(f_*O,O). Tensoring with the plane dualizing bundle gives the formula used in the packet. Because every sign character is its own inverse, the dualization does not change the character label. Crepancy then transports the canonical forms to the minimal resolution.

The determinant twist in adjunction is essential and is correctly included. A residue numerator of parity A has character full-set symmetric-difference A; in canonical degree n−6 this character is even. In particular, at n=7 flipping coordinate i fixes exactly the complementary character indexed by its own canonical coordinate section. The naive singleton-character answer would fail to descend through the projective simultaneous sign and would give the wrong fixed-space dimension. The independent checker compares the eigensheaf sum with the complete-intersection Hilbert coefficient and a third count using square-free root monomials for every n=4,...,16.

For each character quotient w²=product L_i, the claim is an equality of spaces of differential forms after pullback, not an equality of surfaces. A common resolution gives a generically finite separable map between smooth surfaces. Pullback of a nonzero rational two-form is nonzero, so it is injective on regular forms. The source quotient's only possible canonical character is its nontrivial double-cover character, and its dimension equals the corresponding eigenspace computed above. Thus containment plus equal dimension proves the asserted identification. The need for a common resolution is explicitly preserved in the packet.

## 2. Triple-branch D4 resolution and H¹

At a triple, w²=xy(x+y). In the blowup chart x=s,y=st,w=sv, normalization removes s² and gives v²=s t(1+t). Thus the exceptional line belongs to the branch and the three distinct tangent directions produce three ordinary double crossings on it, including the direction at infinity in the other chart. These resolve as A1 points; together with the central exceptional curve they are the D4 configuration. The branch half-class changes from rH to rH−E, while the base canonical class changes from −3H to −3H+E. The corrections cancel exactly.

For q, adding a disjoint exceptional divisor E to O(−rH) has quotient O_E(−1), so H⁰ and H¹ of that quotient vanish. Repeating at the distinct triple points preserves the vanishing of H¹ inherited from the plane. A1 resolutions are rational and do not add H¹. At r=3 the resulting canonical bundle is trivial and q=0, hence the smooth surface is K3. The no-fourfold and distinct-line hypotheses are retained; no D4 assertion is used beyond them.

The offending all-three-sign flip on an A1 cone is the image of scalar i on C² under quadratic Veronese coordinates. It fixes the exceptional conic pointwise, giving the stated 1/4(1,1) boundary after quotient. The checker verifies the local substitution and sign action independently.

## 3. Free deck quotients: exact iff and 2→8 boundary

The inertia criterion is both necessary and sufficient on the resolution. Sufficiency follows by mapping a hypothetical fixed point down. For necessity at a node, a nonidentity stabilizer fixes that node, hence preserves its exceptional P¹; every finite-order automorphism of P¹ over C has a fixed point. This covers even the case of trivial action on the entire exceptional curve.

At six roots, only classes represented by triples can be free. A triple and its complement represent the same sign class, so either concurrent triple forbids that one partition. Distinct free classes have odd parity, hence their nonzero sum has an even representative of size two and has fixed points. Thus no free subgroup larger than two exists. The quotient of K3 by a free involution is smooth, has χ=1 and q=0; pullback followed by norm shows its canonical bundle is two-torsion, and χ excludes triviality. The Enriques conclusion is correct.

At seven roots, all branch generators must map to seven distinct nonzero quotient elements, so the quotient rank is at least three. At order-eight kernel, its image labels are exactly the seven nonzero F₂³ vectors. The only additional condition is that every actual concurrent triple be independent; equivalently it avoids the seven Fano triples. The classification is of maximal-degree free deck quotients, not arbitrary automorphisms. An étale free quotient preserves minimality/general type and divides K² and χ by eight. Pullback injects H¹, so q=p_g=0 and K²=2. The numerical Campedelli claim is consequently justified.

An independent enumeration recovers the thirty Fano systems, hourglass16/12/2 bad-triple census, sixteen free kernels, all120 complementary pairs and exactly32 both-free pairs. It checks each direct inertia intersection, not just the Fano line nomenclature. The two actual triple conditions for the outer-deleted six-root K3 leave eight of the ten partitions.

## 4. Normalized fiber product is genuinely the original root cover

For complementary kernels H,H', the fixed fields satisfy intersection K and compositum L by Galois theory. They are Galois of degree8 over K and linearly disjoint, so their tensor product is the field L, of degree64. This field-level statement is sufficient only generically.

The packet correctly supplies the missing scheme argument. Each quotient algebra is an invariant direct summand of the locally free root algebra, hence is locally free over the smooth base; their tensor is locally free too. It therefore injects into its generic field, making the whole fiber product integral. The root algebra is finite over that tensor algebra (it is already finite over the smaller base ring), normal, and has the same fraction field. It is thus precisely the normalization. There are no extra vertical components or hidden generic components.

The warning about raw ramification is also correct: after an étale separation near a generic common branch line, the local product has equation u²=v². The two analytic/étale branches need normalization. Their presence does not contradict global integrality, just as an irreducible nodal curve may have a reducible completed local ring. No unnormalized-isomorphism or rational-point lifting assertion follows.

## 5. Root canonical budget and Gale identity

BG1's linear-general-position hypothesis guarantees at most d vanishing hyperplanes, with independent local equations. Adjoining their r-th roots gives a regular local ring; unit roots are étale. The monic equations are a regular sequence, their degree is r^(n−d−1), and adjunction gives ((r−1)n−r(d+1))H. The general-position and square-root-specific Gale balance scopes are correctly separated. The fourth-root quartic K3 is a valid hostile control.

For BG2, multiplying M by the six columns of C gives exactly (−2A_5,−uA_2,uA_1,−vA_4,vA_3,2A_0), with determinant2uv. Their scalar product is −4u²v². On the dual parameter plane, set the new covector to M^(-t) times the old covector; the discriminant product becomes that scalar times the original product. Dividing the double-cover coordinate by2uv leaves the minus-one twist. The packet correctly refuses to identify the degree-eight K3 intersection with a degree-two double plane, or import a smooth-net equivalence into the special nodal net.

## Clarification amendment

The character author subsequently clarified that §7 concerns the abstract exterior power of the seven-dimensional vector space of canonical sections, not the geometric wedge product of differential forms on a surface. This clarification is correct and does not change any representation or checker result. The amended source hash is `5570d2e294649f9c1deed7440a408dc64eeb190a23d8f52324515d4d067f69fd`; PASS applies to this amendment as well as the preceding hash.

The combinatorics author also made a bounded precision edit in §7: “least extra datum” with full xor values became “a sufficient decision sidecar” with zero/nonzero xor status. I reopened the amended file; this weakens an unneeded minimality claim and preserves the proof. The current combinatorics hash is `c92376d64e55a09b1dffbba3136fc7e5830ed8ff54ae7196b4bb0f87852d3966`; PASS extends to this amendment.

## Reproduction and source receipt

Run `python3 audit_character.py` and `python3 -O audit_character.py`; both produce the same `character_audit_results*.json` (39 explicit gates). Run `python3 audit_free_quotients.py` and repeat under `-O`; their JSON outputs also agree byte-for-byte. These are independent finite controls; the geometric arguments above are the proof audit.

Primary duality reference independently opened: https://stacks.math.columbia.edu/tag/0DWM . The generalized Fermat definition and genus formula used in the separate geometry packet were independently checked in Hidalgo's author manuscript, https://arxiv.org/html/1710.01349 . Existing Session59 theorem T1 and its independent audit supply the root surface/node foundation.
