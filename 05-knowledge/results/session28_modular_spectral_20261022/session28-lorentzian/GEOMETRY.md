# Score geometry at the Leech cusp: reflection rigidity, metric, and extension loss

Status: elementary PROVED statements with FINITE-EXACT controls. This note extends Session27's `LORENTZIAN_SCORE_BRIDGE.md` and independently checks the incoming pre-audit THM-4752 draft (`01-canon/theorems/THM-4752-tournament-scores-in-ii-25-1-backward-length-and-the-a24-deep-hole.md`, incoming commit d9193e48bd3). It does not assert priority over the concurrent draft or the classical ingredients. No theorem ID is reserved here.

## 1. Inheritance and connection contract

- Closest mechanism: for a 25-vertex tournament with outdegree vector s, v=(s;70) belongs to II_(25,1) and v²=-2c3(T). The transitive score lift w=(0,...,24;70) is a Leech cusp.
- Canonical hostile: Session27 has identical 25-vertex scores with Hamiltonian counts657 and643. No function of v alone can reconstruct H.
- Corrected near miss: a cusp pairing is not a finite hyperbolic distance to a boundary cusp. Cyclic v does not define v-perp/Zv.
- Sidecar: the chosen extension from Boolean edge orientations to fractional edge probabilities matters for a Hessian.
- Map: tournaments → score lifts → hyperbolic points for cyclic tournaments. It preserves cyclic-triangle count and score differences, loses higher orientation/cycle information, and needs marked response observables for faithful reconstruction.

All lattice isometries below are **linear origin-fixing** unless explicitly called affine. This qualification is essential.

## 2. The full score-set stabilizer in the Lorentzian lattice

Let S be the finite set of all score lifts v(T) on25 labelled vertices. Then

    Stab_(O(II_(25,1)))(S) = S_25,

acting by permutations of the25 spatial coordinates. In particular, the score converse s↦24·1−s has no linear integral Lorentzian-isometry extension in this fixed embedding.

### Proof

The convex hull of score vectors is the permutahedron P=conv{permutations of(0,...,24)}. One direct proof avoids importing realizability: scores are sums, over unordered pairs{i,j}, of either e_i or e_j. The convex hull is consequently the Minkowski sum of the segments[e_i,e_j]. A generic linear functional chooses the larger of its two coefficients on every edge, which is a transitive tournament. Hence all its vertices are exactly the transitive score vectors.

Write

    c=(12^25;70),  z=(14^25;60),
    U={(a;0): sum a_i=0}.

Then c²=-1300, z²=1300, c·z=0, and the ambient real space is the orthogonal sum U ⊕ Rc ⊕ Rz. The score polytope lies in c+U and has centroid c. A linear Lorentzian isometry preserving S preserves P and c, and restricts to a Euclidean symmetry of the centered permutahedron on U.

For completeness that symmetry group is S_25×{±1}. Its facets indexed by nonempty proper subsets I have |I|!(25−|I|)! vertices. The maximal cardinality is24!, attained exactly by the50 singleton/complement facets. Their outward unit normal directions are ±a_i, with a_i=e_i−1/25·1. Since a_i·a_j=-1/25 for distinct i,j, a Euclidean permutation of these50 normals must have a single common sign after permuting indices. Thus its action on U is ± a coordinate permutation. This also gives a direct elementary proof of the needed classical permutahedron symmetry fact.

Compose by the inverse coordinate permutation. The restriction to U is then ±I and c is fixed. The remaining orthogonal line Rz has just two possible signs. There are therefore precisely four real extensions:

1. I.
2. R_z(v)=v−2(v·z)z/1300 (identity on U, negation on z).
3. C_U(x;t)=(2 mean(x)·1−x;t) (negation on U, fixes c,z).
4. C_z(v)=−v−(v·c)c/650 (negation on U,z, fixes c).

Only the identity preserves II_(25,1). Indeed u=e_0+e_1 lies in that lattice. For R_z(u) the spatial correction is−196/325 in every coordinate; for C_U(u) the spatial correction is4/25; for C_z(u) it is−144/325. None has all coordinates integral or all half-integral. Thus the other three fail even lattice membership of a single vector. This proves the statement.

Consequences:
- The subgroup preserving both S and the particular cusp w is trivial: w has distinct spatial coordinates.
- At the score-polytope level, converse is a perfectly valid Euclidean symmetry. Its failure is specifically the **linear integral Lorentzian extension**, not a failure of converse as tournament operation.
- There is an affine ambient lattice isometry (x;t)↦(24·1−x;t), because translation by(24^25;0) belongs to the lattice. It preserves the bilinear form of differences, not the origin-based norm Q at arbitrary points. Therefore an unqualified 'no Lorentzian extension' would be false.

## 3. Which root reflections remain tournament score operations?

For a norm2 lattice root r, its reflection is v↦v−(v·r)r. If it preserves the entire finite set S, Section2 puts it in S_25. A coordinate permutation is a reflection precisely when it swaps one pair of coordinates. Thus the roots are exactly

    r=±(e_i−e_j), i≠j.

These generate W(A_24)=S_25. General Leech-root reflections leave the score slice.

A concrete hostile is the extra affine-A_24 root from THM-4752,

    r*=(-1/2, (1/2)^23, 3/2; 5/2).

It satisfies r*²=2, r*·w=-1. On the regular score center c, c·r*=-25, and its reflection gives

    c+25r*=(-1/2,(49/2)^23,99/2;265/2).

This is an authentic II_(25,1) vector of norm−1300 but not a tournament score lift: its time coordinate is265/2 rather than70, and its scores are half-integral/out of range. Norm preservation is not score realizability.

## 4. Actual hyperbolic geometry: c3 is a radial coordinate

For cyclic T, put a=c3(T)>0 and p_T=v(T)/sqrt(2a), so p_T²=-1 on the future sheet. Let p_0=c/sqrt1300. Since c·v(T)=-1300 for every score vector,

    cosh d(p_0,p_T) = sqrt(650/a).

Equivalently,

    a = 650 sech² d(p_0,p_T),
    ||s−12·1||² = 1300 tanh² d(p_0,p_T).

Thus regular tournaments occupy the hyperbolic center; the more score-unequal the tournament, the closer its lift is to the null boundary. This is an exact statement about scores/cyclic triples, not Hamiltonian optimization.

For cyclic T,U, with a=c3(T), b=c3(U),

    cosh d(p_T,p_U)
      = (a+b+||s(T)−s(U)||²/2)/(2 sqrt(ab)).

It follows by expanding(v(T)−v(U))²=||s(T)−s(U)||². Equal scores give identical points regardless of H, including the657/643 hostile.

A null cusp w is at infinite hyperbolic distance. The ranking statistic in THM-4752,

    BL_T(pi)=−v(T)·w_pi,

is a well-defined nonnegative cusp pairing, not that infinite distance. For cyclic T a corresponding horospherical height is

    −p_T·w_pi = BL_T(pi)/sqrt(2c3(T)).

Its logarithm is a Busemann function up to additive normalization and sign convention. Minimizing BL over the finite25! ranking cusps minimizes this normalized height because the denominator is independent of pi. Calling U(T) 'Lorentzian distance to the nearest cusp' is only metaphor unless this normalization/convention is stated.

## 5. An exact Hessian bridge, and why it needs an extension sidecar

Choose one orientation i<j for every complete-graph edge. Let B be the vertex-edge incidence matrix whose edge(i,j) column is e_i−e_j. Give edge(i,j) bias x_ij∈[-1,1], meaning P(i→j)=(1+x_ij)/2, independently. Expected scores are

    s(x)=(n−1)/2·1 + Bx/2.

At n25 their Lorentzian lift has

    Q(s(x);70)=−1300+||Bx||²/4.

For actual Boolean tournaments x∈{±1}^300,

    c3(T)=650−||Bx||²/8.

But for fractional probabilities the expected cyclic-triangle count is instead

    E c3 = binom(n,3)/4 + ||x||²/4 − ||Bx||²/8.

This follows either by expanding the two directed orientations of each triple, or by the two-star covariance terms. At n25 its constant is575, not650. The score-norm extension and the true expected count differ by

    (300−||x||²)/4.

They agree on every tournament, but have different Hessians at x=0. The natural expectation has Hessian

    Hess(E c3) = (2I−BᵀB)/4.

If F_n(x)=E H(T) is the independent-edge expected directed Hamiltonian path count, direct expansion over permutations gives

    Hess F_n(0) = ((n−2)!/2^(n−2)) (2I−BᵀB), n≥3.

Indeed a pair of adjacent path edges contributes the signed multiplicity2(n−2)! divided by2^(n−1); disjoint edges cancel, and diagonal entries vanish by multilinearity. Therefore

    Hess F_n(0) = 4 (n−2)!/2^(n−2) · Hess(E c3).

This is an exact second-order bridge, not equality of the full observables. B Bᵀ=nI−J, so BᵀB has eigenvalue n on the cut/score subspace im Bᵀ (dimension n−1) and0 on the circulation subspace ker B (dimension(n−1)(n−2)/2). For n≥3:

- both expected H and expected c3 have negative Hessian on score directions;
- both have positive Hessian on circulation directions;
- the score norm sees only Bx and is completely blind to circulation directions.

The same Boolean invariant can have different off-cube Hessians: adding a multiple of Σ(x_e²−1) changes its Hessian but vanishes on every tournament. Consequently a Hessian analogy must specify the extension, not only the discrete values. This is a concrete mechanism for why score/Lorentzian geometry cannot recover the cycle representation simply by differentiation.

The closest prior historical artifact located by the parallel recovery lane is `05-knowledge/results/ramanujan_expander_tournaments.out`, lines230–255, with n5 Walsh data ±0.75 and a rough Laplacian description. The general signed-incidence formula is independently derived in this session; no novelty claim over all literature is made.

## 6. Incoming A24 cusp: independent scope check

The incoming THM-4752 formulas for nu=(1/2^25;5/2), nu²=0, nu·v=-25, and centered class s−12·1 check directly. The25 listed Leech roots sum to nu and have affine A24 Gram matrix. Its Niemeier quotient construction and root-system inference use classical Conway/Niemeier theory; this note does not recertify the full external classification.

The induced centered-score class lies in A24 (glue class0) inside N(A24), so converse is negation on that quotient. Section2 explains why this quotient symmetry need not lift to a global integral origin-fixing symmetry of the original fixed-time score embedding.

Two wording repairs are appropriate for the draft:
1. Replace literal 'distance to nearest Leech cusp' by 'minimal pairing with the ranking Leech cusps,' or give the horospherical normalization above.
2. In the cannonball uniqueness clause say 'only nontrivial tournament order n>1 compatible with this even-unimodular score lift'; n=1 is also a trivial compatible order. This does not challenge the substantive n25 construction.

## 7. Sources and audit boundaries

Primary classical source: Borcherds–Conway–Queen–Sloane, *A Monster Lie Algebra?*, Advances in Mathematics53(1984),75–79, especially printedp76 for the coordinate lattice and Leech roots. https://neilsloane.com/doc/Me99.pdf . It also discusses the deep-hole/Niemeier relation onp77. Conway's1983 paper is cited there. We use only this limited dictionary, not claims about arithmetic scattering operators.

The source *The Symmetry Group of the Permutahedron*, College Mathematics Journal42(2),135, DOI10.4169/college.math.j.42.2.135, was found but the full text was inaccessible through the web tool; the needed symmetry fact is proved independently above rather than attributed to an unread proof.

The current SRFP311T1 paper https://srfp311t1.com/conway_cusp_scattering was found and passed to the separate Weil/spectral lane. Its arithmetic-quotient/chamber/toy-model distinction is relevant, but none of the geometric theorems above depends on its scattering claims.

Reproduce: `python3 check_geometry.py > check_geometry.json`. The script uses only Python's standard library and exact fractions. It checks300 root swaps, the three forbidden extension witnesses, the r* hostile, all n3..6 quadratic H coefficients from permutation expansion, exhaustive Boolean cyclic-triple identities through n5 plus a stated stride101 sample at n6, fractional probability identities, and100 deterministic order25 metric identities. These finite tests corroborate the all-order arguments; they are not a complete order25 enumeration.
