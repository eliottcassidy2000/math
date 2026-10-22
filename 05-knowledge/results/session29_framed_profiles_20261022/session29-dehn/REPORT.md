# Session29: Dehn invariants, framed angular valuations, and Gaussian spectral decks

Status: PROVED elementary statements below; CITED classical background; no claim of global novelty or a solution of Hadwiger–Nelson. Scope: the exact relation between scissors data and the Session28 rank-88 angular Gaussian operator. No repository files modified by this lane.

## Executive result

The genuine common mechanism is **an additive tensor observable before quotienting by independent rotations of pieces**. It is not that the rank-88 operator is a Dehn invariant. Ordinary scissors congruence identifies rotations, whereas the rank-88 harmonic tensor transforms covariantly. Taking its full orthogonal coinvariant quotient kills it. Retaining a common frame repairs additivity; taking an operator spectrum repairs rotation invariance but destroys a well-defined addition law.

A constructive bridge survives: a direction-enriched Dehn tensor is translation-invariant, rotation-covariant, and additive under polyhedral dissections. Separately, cyclic Gram moments recover all spectral data of the Gaussian angular operator. These are precise versions of “keep the missing sidecar,” and specify exactly which sidecar is needed.

## Inheritance and correction discipline

Read Session28 `session28-rank88/REPORT.md` sections7a–d and `session28-cm-audit/scalar-followon/ADDENDUM.md`. Its fixed operator, not the earlier coarse-grid scalar conjecture, is inherited. At t=1 the scalar ratio r=a/b has a strict local maximum; at least two reciprocal off-center global minimizers exist. Exactly two, their locations, and actual angular eigenvalues are not established.

Recovered historical threads:
- `07-reflections/dehn-scissors-side-channel-addendum-s637.md`: sound missing-coordinate motivation; regular tetrahedron versus box is a useful hostile.
- `07-reflections/equinumerosity-equidecomposability-fiber-bridge-s617.md`: requires a predicate-preserving fiber map rather than equal counts. Sound methodological constraint, not a geometric theorem identifying tournaments with polyhedra.
- `07-reflections/equidecomposability-the-scissors-congruence-of-the-cross-problem-bridge-s599.md`: strong-component reordering defines a possible combinatorial equivalence relation, not an established map to classical scissors congruence. The assertion that a 3-cycle skew matrix has eigenvalues ±2π/3 is false (its nonzero skew eigenvalues are ±i√3); an angle of a root of unity is a different observable.
- `07-reflections/the-hadwiger-needle-niven-equals-dehn-triviality-equals-the-HN-lattice-escape-s685.md`: statements identifying nonzero Dehn data with chromatic forcing cannot be inherited. Classical three-dimensional Dehn data is not even defined for an abstract planar graph without an additional carrier map. An irrational-angle two-edge path has chromatic number2, so irrationality does not imply forcing beyond3. In a fixed lattice, irrational-angle pairs can occur: vectors (1,2),(2,1) in the square lattice have equal norm and angle cosine4/5, an irrational multiple ofπ. Thus “irrational angle escapes every lattice” is false if it means the finite configuration; infinite invariance under the associated rotation is a different condition.
- `01-canon/theorems/THM-418-lattice-chromatic-dichotomy-of-the-plane.md`: the proved square/triangular norm-graph colorings do not make the subsequent Dehn-equivalence prose a theorem. `THM-419-2d-lattice-chromatic-cap.md` explicitly has mixed proved/computational/global scope and is not needed here.

Closest mechanism: Session28 harmonic-theta factorization. Hostile: independently rotate additive pieces. Corrected near miss: retain the frame rather than demand congruence invariance. Least-used sidecar: joint cyclic Gram data, not just scalar theta or marginal norms.

## Concept board and two hypothesis/test rounds

1. Scissors group: dissection plus Euclidean congruence; retained volume and classical Dehn data; orientation is discarded.
2. Harmonic quartic tensor: linear in a measure; rotation-covariant; orientation is retained.
3. Operator spectrum: rotation-invariant; relative eigenspace position is discarded.
4. Gaussian temperature pencil: a(t)I+b(t)A; temperature alters scalars, not eigenspaces.
5. Incidence decks: joint Gram-cycle data can recover spectral moments; radial marginals cannot.

Round1 hypothesis: harmonic fourth moments descend to ordinary scissors classes. Test: explicitly scissors-congruent equal-volume boxes and rotated atomic frames. Result: REFUTED, with exact witnesses below. Strongest survivor: a framed additive valuation.

Round2 hypothesis: passing to eigenvalues restores both invariance and additive gluing. Test: equal-spectrum matrices A and −A obtained by rotating a 2-design. Result: REFUTED; aligned sum has nonzero spectrum, rotated sum is zero. Repair: retain relative-frame mixed moments, equivalently cross-Gram information. Positive continuation: exact Gram-cycle trace formula and a direction-enriched Dehn valuation.

## 1. The orthogonal quotient obstruction (PROVED)

Let H4(R^n), n≥2, denote real homogeneous harmonic quartics, or the equivalent trace-free symmetric fourth-tensor representation. The only O(n)-invariant element is zero: an invariant quartic must be c||x||^4, whose Laplacian is 4(n+2)c||x||².

Choose an O(n)-invariant inner product. The orthogonal complement of

    W = span_R{g·h−h : g∈O(n), h∈H4}

is exactly the invariant subspace. Therefore W=H4 and the O(n)-coinvariant vector space H4/W is zero. The same holds for the rational span of these relations: because h is arbitrary, multiplication of g·h−h by any real a is g·(ah)−ah, another relation. Thus the Q-span is already R-stable and equals W. This justifies the corresponding statement for the tensor-over-Q target in section3. Any real-linear rotation-invariant map from H4 to a trivial O(n)-module is zero. Likewise any equivariant linear map from a trivial module into H4 is zero.

For a spherical 2-design shell X of common norm²m, its isotropic fourth moment subtracted from the fourth tensor is in H4. The operator A on Sym_0 obtained by contracting this tensor transforms by

    A_(gX) = rho(g) A_X rho(g)^−1,
    rho(g)H = g H g^T.

Hence A cannot descend unchanged to an equivalence relation allowing independent rotations of pieces, unless its harmonic part is zero. This rules out the proposed **linear, equivariant, additive** identification with ordinary Dehn data, whose geometric rotation action is trivial. It does not rule out arbitrary nonlinear encodings or maps with additional markings.

## 2. Explicit scissors and spectrum hostiles (FINITE-EXACT + PROVED)

Let p(x,y,z)=x^4−6x²y²+y^4, a harmonic polynomial in R³. For a centered box of side lengths a,b,c and volume1,

    integral_box p = a^4/80 + b^4/80 − a²b²/24.

The unit cube and the 2×(1/2)×1 box are scissors congruent by cutting the cube into two y-halves and translating them end-to-end in x. Both have volume1 and zero classical Dehn invariant, but their p-integrals are respectively

    −1/60 and 611/3840.

Thus even recentering the global bodies does not make harmonic moments a scissors invariant. The integral is additive with a fixed origin; recentering each piece separately is not an additive operation.

For a purely angular hostile, take X={±e1,±e2}. In the orthonormal basis diag(1,−1)/√2 and [[0,1],[1,0]]/√2 of Sym_0(R²), its angular operator is A=diag(1,−1). Rotating X byπ/4 gives B=−A. X and its rotation are both antipodal 2-designs with identical radial data and identical angular spectrum {1,−1}. Nevertheless

    tr((A+A)²)=8, whereas tr((A+B)²)=0.

Consequently addition on conjugacy classes/spectral multisets cannot reproduce addition of actual operators. Relative eigenspace position is indispensable. This is a universal gluing obstruction, not a counterexample involving a rank-88 lattice.

## 3. A repaired direction-enriched Dehn tensor (PROVED)

For a convex three-dimensional Euclidean polyhedron P, define

    D_dir,4(P) = sum_edges e [ length(e) h4(u_e) ] tensor_Q [theta_e]

in H4(R³) tensor_Q (R/πQ), where u_e is either unit direction along the edge, h4(u) is the harmonic projection of u^⊗4, and theta_e is the interior dihedral angle. The even degree makes the choice of sign of u irrelevant. The first factor is a real vector space viewed over Q for the tensor product. Equivalently define the full Sym^4-valued expression before harmonic projection. This is a construction proved here; no novelty relative to general tensor/Hadwiger valuation literature is claimed.

Properties:
- Translation invariance: edge lengths, directions and dihedral angles do not change.
- Rotation covariance: D_dir,4(gP)=g·D_dir,4(P).
- Dissection additivity in a common frame: first use finite dissections into convex polyhedra, and refine all edges to a common edge complex. For a nonconvex polyhedral body, define the value by any finite convex dissection; the following cancellation proof and a common refinement prove independence of this choice. Equivalently use actual interior sectors, including reflex sectors, at each edge. Along a subdivided original edge, lengths and incident dihedral angles sum to the original contribution. Along a newly created boundary seam, angles sum toπ; along an interior edge, to2π. All incident sectors share the same edge direction, so their tensor coefficients agree and their angle sums vanish moduloπQ. Further splitting lengths is additive. This proves the identity after regrouping all subedges.
- It is an invariant under dissections and reassembly using translations alone. Independent rotations are not allowed without transporting tensor values back to their original/common frame.

The ordinary Dehn invariant instead omits h4(u), yielding sum length(e) tensor_Q [theta_e]. Thus the two constructions have the same local cancellation mechanism but different group actions. The harmonic direction component is killed by the real-linear full-rotation coinvariant quotient; the scalar length-angle coordinate survives. Do not confuse the full Sym^4 invariant trace component with its harmonic part.

This construction is genuinely nonzero: take a regular tetrahedron with vertex coordinates (1,1,1),(1,−1,−1),(−1,1,−1),(−1,−1,1). Its six unoriented unit edge directions satisfy sum p(u)=−1 for p=x^4−6x²y²+y^4. All edges have common length l and common angle theta=arccos(1/3), with theta/π irrational. Evaluation of D_dir,4 against this harmonic p therefore gives −l tensor_Q [theta] ≠0. A nonzero simple tensor over Q is nonzero. Boxes give zero because all angles areπ/2.

There is no proposed map here from lattice shells in R88 to polyhedral edge-angle data in R³. The verified common structure is covariance plus additive local cancellation, not an identification of the objects.

## 4. Cyclic Gram decks recover the angular spectrum (PROVED)

Let X be a finite spherical 2-design of N points in R^n of norm²m. Let V=Sym_0(R^n), d=(n−1)(n+2)/2, with trace inner product. Set

    v_x = xx^T − (m/n)I,
    M_X = sum_x |v_x><v_x|,
    c_X = 2N m²/[n(n+2)],
    A_X = M_X − c_X I_V.

The 2-design identity implies tr(A_X)=0 and identifies A_X with the harmonic quartic defect. The centered Gram kernel is

    K(x,y)=<v_x,v_y>=(x·y)²−m²/n.

For every integer k≥1,

    tr(M_X^k) = sum_(x1,...,xk in X) product_(i=1..k) K(x_i,x_(i+1)), x_(k+1)=x1.

Proof: expand the product of rank-one operators and take the trace. Repetitions of vertices are allowed. This is a weighted closed-walk/cyclic Gram-deck statistic, not an induced-subgraph count.

In particular, using the 2-design identity to simplify the quadratic term,

    tr(A_X²) = sum_(x,y) (x·y)^4 − 3N²m^4/[n(n+2)] >= 0.

Equality holds precisely when A_X=0, hence (for an antipodal 2-design) when X is a spherical4-design. The pair fourth-moment excess is exactly the Hilbert–Schmidt angular anisotropy. No scalar theta series determines this excess in general.

For two spherical 2-designs X,Y, of sizes N_X,N_Y and respective norm²m_X,m_Y in the same frame,

    tr(A_X A_Y) = sum_(x,y) (x·y)^4
                  − 3N_X N_Y m_X²m_Y²/[n(n+2)].

This is the missing orientation-sensitive cross-term in tr((A_X+A_Y)²). For equal-radius union, the angular operators add, including multiplicities.

Recover tr(A_X^k) by binomial expansion of M_X−c_XI. Newton identities then recover the characteristic polynomial from k=1,...,d, and therefore every Gaussian Hessian spectrum in the pencil a(t)I+b(t)A_X. For rank88, d=3915, N=168498000, m=8 and c=2723200. This is an exact sufficiency bound, not a claim that enumerating these huge decks is practical, minimal, or needed in a representation-theoretic compression.

Practical next route: exploit Aut(L)-orbits of pairs and tuples or its commutant on Sym_0. Pair inner-product frequencies already give tr(A²). Higher cyclic products require joint data, not merely the unlabelled histogram of pairwise inner products. Even complete operator spectrum does not reconstruct the lattice or its framed angular tensor.

### 4a. A cheap spectral bound from the pair deck (PROVED)

Write V2=tr(A²), with d>1 and tr(A)=0. Then

    −sqrt((d−1)V2/d) <= lambda_min(A) <= −sqrt(V2/[d(d−1)]).

For the first bound, let lambda_min=−a; the other d−1 eigenvalues sum to a, so Cauchy gives V2>=a²+a²/(d−1). For the second, all eigenvalues are at least−a and sum to0; writing them as y_i−a with nonnegative y_i of sum da gives V2=sum y_i²−da² <=d(d−1)a². Both inequalities are sharp for unrestricted trace-zero real symmetric operators. They need not be sharp on shell operators.

Thus a certified scalar lower bound r(t)>=r0 with r0>sqrt((d−1)V2/d) guarantees all-temperature positive Hessians. Conversely, a temperature with r(t)<sqrt(V2/[d(d−1)]) forces a negative mode. Neither criterion is presently evaluated for an actual rank88 shell. The existing positive-semidefinite Gram construction also gives lambda_min>=−c, independently of the pair histogram.

### 4b. Two additional category checks (PROVED, scoped)

Classical Dehn and D_dir,4 are homogeneous of degree1 under scaling the polyhedron; the shell angular operator is homogeneous of degree4 under scaling its vectors. A linear map intertwining these two scaling actions must be zero: L(rD)=rL(D) but scale compatibility requires r^4L(D). This obstructs a naive scale-equivariant linear identification. Restricting to a fixed normalization avoids the scaling mismatch but does not fix the rotation/additivity obstruction.

Also, any continuous additive angle functional f:R→R^k with f(πQ)=0 is zero, since πQ is dense and R^k is Hausdorff. Hence a nonzero Dehn angle coordinate cannot be replaced by a continuous finite-dimensional additive numerical angle statistic. This is only an obstruction to that specific continuous additive angle factor, not a theorem ruling out every nonlinear function of a polyhedron. Smooth harmonic moments and Dehn angle classes live in different regularity categories.

## 5. What transfers to tournament work, and what does not

The map in section4 is actual linear algebra: a point configuration maps to a centered Gram kernel, then to weighted cycle sums, then to an operator characteristic polynomial. It provides a testable pattern for an explicitly defined tournament kernel: closed-walk traces recover that kernel's spectrum, not the tournament itself. To recover tournament isomorphism requires a separate injectivity result; spectral equality alone is insufficient.

For metagraph gluing, the lesson is operational: a scalar or spectrum can be complete for an individual observable yet fail to support composition. Retain the relative transporter/frame/cross-Gram data whenever the construction adds operators. Neither classical Dehn invariants nor the integer24/88 coincidence prove any tournament wiring or Hamiltonian-count theorem.

## Sources and verification boundary

- Rich Schwartz, *The Dehn–Sydler Theorem Explained* (2010), author-hosted exposition of Dehn/Sydler/Jessen: https://www.math.brown.edu/~res/Papers/dehn_sydler.pdf . Used only for the classical three-dimensional volume-plus-Dehn classification and angle tensor convention; our obstruction does not require sufficiency.
- Monika Ludwig and Laura Silverstein, *Tensor valuations on lattice polytopes*, Adv.Math319 (2017),76–110: https://arxiv.org/abs/1704.07177 . Context: tensor-valued valuations with covariance are an established field, not a new vocabulary invented here. Its classification theorem is not used in our proof.
- Daniel Hug and Rolf Schneider, *SO(n) covariant local tensor valuations on polytopes* (2016): https://arxiv.org/abs/1605.00954 . Context: distinguish covariance from invariance and track the precise symmetry group. No classification theorem is imported into the new construction.
- Session28 harmonic-theta derivation and cited rank88 shell hypotheses remain the premises for applying section4 specifically to extremal rank88 lattices.

`check_dehn_angular.py` uses exact rational/algebraic SymPy arithmetic. It checks both hostile examples, harmonicity, trace zero, the pair excess formula, and explicit cyclic sums k=1..5 on a 3D cross-polytope. It does not enumerate rank88 vectors, prove minimum8, compute its actual angular spectrum, or certify global Gaussian stability. The first test invocation needed algebraic simplification in a zero comparison; the mathematical expression was identically zero, and the corrected expanded comparison passes.
