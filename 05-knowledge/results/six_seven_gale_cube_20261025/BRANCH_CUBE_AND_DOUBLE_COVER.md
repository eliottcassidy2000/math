# One hourglass, two K3 branch languages

Status: **PROVED, independently audited**; see [audit](AUDIT_BRANCH_GEOMETRY.md) and [integration audit](AUDIT.md). Characteristic zero throughout; geometric statements use an algebraically closed field. Arithmetic assertions below retain their stated ground fields. No claim of historical novelty is made for abelian-cover or generalized-Fermat-curve mechanisms.

## Inheritance and live concepts

The closest proved mechanism is Theorem T1 in `05-knowledge/results/session59_f1_geometry_and_information_20261025/toric-shared-roots/TORIFICATION_OBSTRUCTION.md`, independently approved in the sibling `positive-audit/TORIFICATION_AUDIT.md`: the n-root surface is a complete intersection of n−3 quadrics, has only A1 nodes when no four branch lines concur, and its resolution has K=(n−6)H and H²=2^(n−3). The six-root resolution is K3; the seven-root resolution has p_g=7, K²=16 and e=80.

The canonical hostile example is a fourth concurrent branch line, where the inherited A1 argument ceases to apply. The near miss repaired here is treating the map of minimal resolutions as a finite double cover of the original K3: the branch crossings require blowing up the K3 first. The least-used sidecar is the proportionality constant between two old linear forms restricted to the new line: it controls the arithmetic field of the components.

Concept board: labelled line incidence; sign/deck characters; canonical budget; branch normalization; cube incidence; arithmetic descent. The concrete bridge is a change of omitted root, with the same seven-root function field as target. The two branch presentations preserve that function field and the numerical invariants, but lose the other six-root projection unless its omitted coordinate is retained as a sidecar.

## 1. General branch-cube theorem

Let n≥6, and let L_1,...,L_n,L_* be distinct linear forms spanning k³. Assume no four of their lines in P² concur. Let X_n and X_(n+1) be the corresponding square-root surfaces, Y_n and Y_(n+1) their minimal resolutions, and H the pullback hyperplane class on Y_n. Let r be the number of double intersections of the old n lines that lie on the new line ℓ={L_*=0}, and put d=n−r.

**Theorem B1.** Forgetting the new root has degree two. Its branch divisor B on Y_n is in |2H| and avoids the old exceptional curves. It has s=r·2^(n−3) ordinary nodes and no other singularities. Its normalization consists of 2^r smooth connected curves, each of genus

    g=1+2^(d−3)(d−4).

The components are indexed by {±1}^r. Two components meet if and only if their labels differ in exactly one coordinate; every such adjacent pair meets transversely at 2^(d−2) distinct points. Thus the dual multigraph is the r-cube, with 2^(d−2) edges replacing each cube edge. In particular there are no self-nodes.

**Proof.** No old triple point lies on ℓ, since that would create a quadruple point. The r old pairs that meet on ℓ have disjoint index sets, because an old line meets ℓ only once and distinct such intersections cannot coincide. Restricting the n forms to ℓ gives exactly d distinct zeros: r repeated twice and n−2r occurring once. In each repeated pair, the two restricted forms are proportional. Choosing the ratio of their roots supplies one independent binary sign. These give 2^r disjoint components away from the repeated zero positions.

After fixing those signs, normalize by one of the d distinct restricted forms. The other d−1 square classes are independent in k(P¹)*/k(P¹)*²: any selected numerator has an odd valuation at its own distinct zero. Hence each component is a connected degree 2^(d−1) cover of P¹ with group (Z/2)^(d−1), branched at the d distinct points, each inertia of order two. Its normalization is smooth, and Riemann–Hurwitz gives

    2g−2=−2·2^(d−1)+d·2^(d−2)=2^(d−2)(d−4).

This is the classical generalized Fermat curve of type (2,d−1).

At a repeated position, two old lines vanish. The old root surface is smooth there, with local root coordinates u,v. The new line pulls back to αu²+βv²=0, αβ≠0. It has two distinct transverse branches. The corresponding sign is precisely the binary coordinate distinguishing this repeated pair; all other component signs are fixed. Each normalized component has 2^(d−2) points above that position. Thus adjacent components meet at exactly that many points, and nonadjacent components cannot meet. Equivalently the fiber of X_n at that old double intersection has 2^(n−3) reduced points, giving s=r·2^(n−3). No other singularities occur. Finally, the square-root map pulls a base line back to 2H, proving the divisor class. QED.

The geometric component count does not say that every component is defined over the original field. Square roots of the pairwise proportionality constants are the required descent sidecar.

**Theorem B2 (resolved double cover).** Blow up precisely the s nodes of B:

    σ:Z→Y_n,   L'=σ*H−Σ E_j.

The strict transform B' is a disjoint union of the smooth components in B1 and belongs to |2L'|. The finite double cover f:Y_(n+1)→Z branched along B' is the minimal resolution of X_(n+1), and

    K_(Y_(n+1))=f*σ*(K_(Y_n)+H).

The map to Y_n is f followed by σ; it is not finite when r>0.

**Proof.** At a node of B, the new cover has local equation w²=uv after a linear change. Blowing up the base crossing and subtracting twice the exceptional divisor from the branch makes the strict branch smooth and disjoint. The normalized double cover is smooth and resolves w²=uv. The exceptional base P¹ has two branch points, so its preimage is a P¹ of self-intersection −2. At an old node the new radicand is nonzero, so the cover is étale and each old exceptional curve lifts twice. These account for all singularities. The smooth-double-cover identity and the blowup identity give

    K_Z+L'=(σ*K_Y+ΣE_j)+(σ*H−ΣE_j)=σ*(K_Y+H).

The canonical class is nef for n≥6 by the inherited formula, so this resolution has no exceptional (-1)-curves. QED.

For n=6 the four possible branch types are:

| r | d | Components | Genus of each | Branch nodes | Dual multigraph |
|---:|---:|---:|---:|---:|---|
| 0 | 6 | 1 | 17 | 0 | one vertex |
| 1 | 5 | 2 | 5 | 8 | two vertices, eight edges |
| 2 | 4 | 4 | 1 | 16 | square, four edges per side |
| 3 | 3 | 8 | 0 | 24 | cube, two edges per side |

Every row gives e(Y_7)=2(24+8r)−2^r(2−2g)=80 and K²=2H²=16. Also p_g(Y_7)=p_g(Z)+h⁰(Z,K_Z+L')=1+h⁰(Y_6,H)=1+6=7. The last equality follows directly from the six-coordinate complete intersection, and the split pushforward O⊕(L')^(-1) yields the double-cover section formula. These are independent accounting checks, not an inference from numerology.

## 2. Two presentations of the same actual hourglass

Use the source arithmetic exactly:

    c²=E,
    a²=E+x, b²=E−x,
    d²=E+z, e²=E−z,
    h²=E+x+z, k²=E−x−z.

The seven-root surface satisfies a²+b²=d²+e²=h²+k²=2c² and h²=a²+d²−c². Its three branch triples lie at E=0, and its minimal resolution has 24 exceptional (-2)-curves.

### 2a. Delete the outer root k: two genus-five components

The six remaining forms have two triple points, so their complete intersection has eight A1 nodes and K3 resolution. The new line is E−x−z=0. Only one old double intersection lies on it, hence r=1. On the branch locus h²=2c². Over Q(√2), the branch splits as h=+√2c and h=−√2c. Each component is a smooth genus-five hyperplane section of the K3. They meet transversely at eight points.

Each component is isomorphic over Q(√2) to the same curve C, which has a model over Q:

    C ⊂ P⁴[c:a:b:d:e],
    a²+b²=2c²,
    d²+e²=2c²,
    b²+e²=c².

In the chart c≠0, with t=x/E and z/E=1−t, the square classes are t,1+t,1−t,2−t. Thus C→P¹ has degree16 and five branch points {−1,0,1,2,∞}. The normalization calculation in B1 proves smoothness and genus5; it can also be checked from the complete intersection of three quadrics.

The eight intersections of the two components occur at c=h=0. Normalize a=1; then b=±i, e=±1, d=±i, independently. All are Q(i)-points. This illustrates why the component field Q(√2) and intersection-point field Q(i) must not be conflated. The branch has no Q-points: h²=2c² forces h=c=0, and a²+b²=d²+e²=0 then force every coordinate to vanish over Q. It does have real points with c≠0.

Blowing up these eight crossings makes the two genus-five branch curves disjoint. Their double cover is the seven-root resolution; the sixteen lifts of old nodes and the eight new exceptional curves account for its 24 nodes.

### 2b. Delete the center root c: eight rational conics on a smooth K3

The six outer forms have no triple points, so their three-quadric intersection is already a smooth K3. Adding E=0 gives r=3. Over Q(i), its branch divisor consists of eight conics C_ε, indexed by ε=(ε₁,ε₂,ε₃)∈{±1}³:

    b=ε₁ i a,   e=ε₂ i d,   k=ε₃ i h,
    h²=a²+d².

The conic is smooth and rational, with for example [a:d:h]=[1:0:1] a rational point over Q(i). For indices differing in one coordinate, setting the corresponding a, d, or h to zero gives exactly two transverse intersection points. No other pairs meet. There are 12 cube edges and two points on each, hence 24 crossings. Every crossing is defined over Q(i). Complex conjugation changes all three signs, acting antipodally on the cube. There are no real branch points because a²+b²=d²+e²=h²+k²=0 forces every real coordinate to vanish.

The eight conics have intersection matrix

    M=−2I+2A(Q₃).

Its eigenvalues are 4,0,−4,−8 with multiplicities 1,3,3,1; its rank is5. The three zero directions are actual divisor relations: for each coordinate j, each of the two half-sums Σ_(ε_j=±1) C_ε is a hyperplane section H. For example b−ia=0 and b+ia=0 give the two half-sums for j=1. Therefore Σ C_ε=2H, H·C_ε=2, and Σ ε_j C_ε=0. This supplies a geometric explanation of the matrix kernel, rather than treating an intersection-matrix nullspace as an unproved Picard relation.

Blow up all24 crossings, then take the double cover branched over the eight disjoint strict conics. This gives the same seven-root resolution as §2a. Its finite double-cover base is a24-point blowup of a smooth K3, whereas §2a uses an8-point blowup of a nodal K3's resolution.

Thus “two genus-five curves” and “a doubled-edge cube of eight conics” are exact alternate branch languages of one seven-root surface. The operation is changing which root is forgotten. The bases are now known to be nonisomorphic over C: their Picard ranks are 20 and 19, respectively, by the concurrent exact-rank result linked below.

## 3. The dual conic and the A₂ hexagon

The dual coefficient points of the six outer branch lines are

    [1:±1:0], [1:0:±1], [1:1:1], [1:−1:−1].

They all lie on the nonsingular rational conic

    Q(U,V,W)=U²−V²−W²+VW=0.

Its symmetric matrix has determinant3/4. Hence the six original lines are tangent to the corresponding smooth dual conic. In U=1 their (V,W)-coordinates are precisely the six minimal nonzero integer vectors of V²+W²−VW, the A₂ hexagon. The dihedral symmetry of that binary form explains the six labels geometrically; it is not just a count of six.

A parametrization is

    [U:V:W]=[t²−t s+s²:t²−s²:t²−2ts].

The six marked values of t/s are {0,1,2,−1,1/2,∞}. The Möbius transformation t↦(2t−1)/(t+1) cycles them as

    ∞→2→1→1/2→0→−1→∞.

These identities are rational and are checked exactly. The binary sextic

    F(t,s)=ts(t+s)(t−s)(t−2s)(2t−s)

satisfies F(2t−s,t+s)=−27F(t,s). The coefficient −27 is retained: a projective symmetry of the branch set need not lift over Q to an automorphism of y²=F.

The outer-deleted six-line arrangement in §2a has triple intersections and cannot consist of six tangents to a smooth conic: through a fixed point there are at most two tangents to a smooth conic. This is a decisive distinction between the two anchor languages.

**RECOVERED concurrent closure.** Session60 now identifies the central six-root K3 itself with a Kummer surface by an explicit diagonal/Gale coordinate map; see [KUMMER_GALE_BRIDGE](../session60_six_seven_geometry_20261025/deck-k3/KUMMER_GALE_BRIDGE.md). The marked genus-two curve becomes v²=u⁶−1 over C, and its Jacobian is isogenous to two j=0 elliptic curves. The rational Shioda equation retains the factor +27; see [the equations](../session60_six_seven_geometry_20261025/hourglass-symmetry/RESULT.md). This +27 belongs to a different map from the −27 Möbius multiplier computed above. [The exact-rank corollary](../session60_six_seven_geometry_20261025/outer-picard/EXACT_RANK_COROLLARY.md) gives central rank20, outer rank19, and hourglass rank58; the seven deletions form exactly two complex K3 isomorphism classes, in a 1+6 split. These are inherited advances, not conclusions of our branch-count computation. Arithmetic twists and integral-lattice identification still require their own maps.

**Overlap and extension.** The six-to-seven cube, component genera and true hyperplane relations were independently obtained in [DELETION_AND_CUBE, D1–D4](../six_seven_geometry_20261025/root_deletion/DELETION_AND_CUBE.md). That packet also supplies the general Walsh spectrum for the six-root deletion branches. The additions here are B1/B2 for arbitrary old root count n≥6, the exact component/intersection fields, and the [genus-five elliptic factorization](GENUS_FIVE_ELLIPTIC_FACTORS.md).

## 4. Scope, literature, and reproduction

Primary-source imports are classical machinery only. Hidalgo, *Holomorphic differentials of Generalized Fermat curves*, describes type(k,n), its genus, and type(2,4) Humbert curves: https://arxiv.org/html/1710.01349 . The genus formula used here is proved above by Riemann–Hurwitz. Pardini's primary lecture notes, *Abelian covers and complete moduli of varieties of general type*, state eigensheaf and Hurwitz formulas: https://pagine.dm.unipi.it/pardini/files/UGA10B.pdf . Auel–Singer, *The algebraic geometry of 3-by-3 magic squares of squares*, Section4.2, treats the magic octic K3 projections and associated double covers: https://arxiv.org/html/2609.09351v1 . Its discriminant double cover must not be silently identified with a coordinate-character quotient. Existing magic-surface geometry is not claimed as a discovery here.

The exact checker enumerates all252 transitions obtained by selecting six of the nine labelled magic lines and adjoining one omitted label. It checks incidences and node accounting independently from the symbolic branch formulas. It includes the actual center and outer deletions, conic parametrization, genus-five branch equations, cube spectrum and true divisor-kernel relations, plus fourfold-concurrency and coefficient-twist hostile controls. This finite universe supports the applications; the general theorem is proved above.

From repository root run `python3 -B 04-computation/six_seven_gale_cube_20261025/check_geometry.py` and repeat with `-O`; compare with `geometry.json` in this result directory. The [packet replay](../../../04-computation/six_seven_gale_cube_20261025/verify_packet.py) performs both checks. This packet reserves no numeric theorem IDs.
