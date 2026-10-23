# Simplicial gauge fields and the exact Fano fourth-order Wilson observable

Status: PROVED elementary finite statements; FINITE-EXACT independent checks; inherited Fano Hamiltonian identity is a separately audited finite identity. No literature-wide novelty claim. This is a finite combinatorial gauge model, not a construction of continuum Yang–Mills theory or a continuum mass-gap result.

## 0. Inheritance, corrections, and concept board

Closest proved mechanisms:
- `01-canon/theorems/THM-1405-star-quotient-is-the-cycle-space-transverse-to-isomorphism.md`: star flips, cut spaces, and binary Wilson loops, there on the tile graph K_n minus a fixed Hamiltonian path. Its carrier is NOT the full K_n used here.
- `01-canon/theorems/THM-2501-switching-fourth-moment-signed-c4-and-gram-energy.md`: symmetric sign fourth moments are C4 Wilson sums and Gram energy. Its explicit skew/symmetric distinction is retained.
- `05-knowledge/results/session33_continuation_states_20261023/FANO_CONTINUATION_SPECTRUM.md`: seven weight-four Walsh characters encode the 171/189 Fano Hamiltonian count. The 128-tournament identity and its prior arithmetic-braids provenance are inherited, not newly discovered.
- `05-knowledge/results/unit_path_ring_kuratowski_pythagorean_yang_mills_20261018.md`, section 5: gauge/flux and reflection-positivity motivations, previously explicitly analogical. No ladder-decay estimate is imported as a mass gap.

Canonical hostile: switching one vertex of the seven-vertex Paley tournament preserves EVERY Wilson loop but changes its Hamiltonian-path count from 189 to 135. Corrected near miss: MISTAKE-664 requires orientation signs in directed path characters. Least-used sidecar: the fixed Fano line partition of the 21 edges. It selects a transverse seven-dimensional slice, not a gauge orbit.

Board: (1) edge cochains and cuts; (2) triangle curvature and tetrahedral closure; (3) Wilson character support; (4) Fano line parity and Hamiltonian response; (5) finite partition functions versus continuum limits. Anchor: an exact fourth-order bridge. Niche: character/surface expansion. Wildcard: nonabelian transported Bianchi identity and its hostile abelianization.

Startup packet read; cached local HEAD was 22d9284f72ed, clean. No Git operations or repository writes were performed. The requested message peek failed because `.machine-id` is absent; no registration was invented. Source search found the two Wilson/switching theorems above; no claim of exhaustive literature priority is made.

## 1. The exact tournament-to-gauge map

Fix a total order on n labelled vertices, n>=2. Set u_ij=+1 if i->j for i<j and -1 otherwise. Regard u as a Z2-valued link field on unordered edges; when traversing backwards set u_ji=u_ij^{-1}=u_ij. This last convention is essential: these are symmetric link signs extracted from a tournament, NOT the skew tournament matrix K, for which K_ji=-K_ij.

A vertex sign g acts by u_ij -> g_i u_ij g_j. On the tournament this reverses precisely the arcs crossing the cut {i:g_i=-1}. Constant g acts trivially; the effective group has 2^(n-1) elements and acts freely. In additive bits x_e in F2, the action is x -> x+delta g.

For an edge subset z define chi_z(u)=product_(e in z) u_e. Gauge transformation multiplies this by product_v g_v^(deg_z(v)). Therefore:

**Theorem 1 (Reynolds / cycle-character projection).** The gauge average R f=2^(-n) sum_g f(g.u) fixes chi_z iff every vertex has even z-degree, and kills it otherwise. These Eulerian characters form a basis for all gauge-invariant functions. The dimension of this function space, equivalently the orbit count, is 2^(E-n+1), E=binom(n,2).

Proof: independent vertex-sign averaging kills exactly odd incidence. Walsh characters are an orthonormal basis on {+-1}^E. The cycle-space dimension is E-n+1. Alternatively Burnside on the full 2^n vertex group has only its two constant elements fixing any edge assignment; (2*2^E)/2^n gives the same count. This is a direct representation-theoretic projection, not an analogy to primes. QED.

A loop Wilson observable is chi_z for the edge-parity support of the loop. Retraversed edges cancel; products of Wilson observables correspond to symmetric difference of edge sets. Degrees mean reduced Walsh degrees, not unreduced word length.

## 2. Curvature is complete; tetrahedral Bianchi is exact

Let f_ijk=u_ij u_jk u_ik for i<j<k. These triangle Wilson loops are gauge invariant. For every four vertices,

f_ijk f_ijl f_ikl f_jkl=1.                                      (1)

Each edge appears twice. Conversely, every assignment of triangle signs satisfying all relations (1) comes from an edge field, uniquely up to vertex switching. To reconstruct choose root 0, put u_0i=1, and u_ij=f_0ij. The relation on {0,i,j,k} supplies f_ijk=f_0ij f_0ik f_0jk. Thus the binom(n-1,2) root triangles are independent coordinates of the gauge quotient. Flatness f=1 implies a pure gauge field.

This uses the full simplex's triangle and tetrahedron incidence, not just graph terminology. For n=3 there is no Bianchi condition, and its single triangle flux is arbitrary. If tetrahedra are absent from an arbitrary complex, the sufficiency assertion must not be imported; topology supplies extra obstructions.

A quadrilateral satisfies the discrete Stokes identity

W_(i j k l)=f_ijk f_ikl.                                       (2)

So four-edge Wilson loops do not introduce a new independent curvature coordinate on the full simplex, but they give genuine fourth-degree edge observables. Any cycle can be triangulated similarly.

**Hamiltonian obstruction.** Let H(T) count directed Hamiltonian paths. For every tournament,

R H = n!/2^(n-1).                                             (3)

For each ordered vertex word, making its n-1 edges all directed forwards imposes n-1 independent equations on the vertex signs along a tree, leaving exactly two of 2^n gauges. Sum over n! words. Thus H is NOT a gauge-invariant observable on all tournaments; its gauge projection is constant. At n=3 a transitive tournament with H=1 can be switched to the directed cycle with H=3. The seven-vertex certificate gives the stronger within-orbit hostile 189->135.

## 3. The Fano rainbow-loop theorem

Use exactly Session33's convention: vertices Z/7, Fano lines L_i={i,i+1,i+3}, and Paley baseline T_0 oriented by differences {1,2,4}. Let u^0 be its upper-order symmetric link field. Reversing every edge on line i introduces z_i=(-1)^(m_i), hence

u_e(m)=u^0_e z_(ell(e)),                                      (4)

where ell(e) is the unique Fano line containing e. Every line has three edges, and the lines partition all 21 edges. The seven line normals in the nonzero-vector point labelling (1,2,4,3,6,7,5) are a=(4,1,3,7,6,5,2). Session33 gives

H(T_m)=693/4+(9/4) sum_(b in F2^3, b!=0) product_(i:b.a_i=1) z_i. (5)

A **rainbow quadrilateral** is an undirected simple four-cycle whose four edges lie on four different fixed Fano lines. Let R be their set. Normalize its Wilson loop by the fixed baseline:

Wbar_C(u)=W_C(u)/W_C(u^0).                                    (6)

This normalization retains the orientation-reference sidecar, avoids a dropped path sign, and is gauge invariant in u. It is not an arbitrary reorientation of the original tournament.

**Theorem 2 (exact fourth-order bridge).** There are 21 rainbow quadrilaterals, with exactly three for each of the seven line-support sets {i:b.a_i=1}, b!=0. On all 128 Fano line orientations,

H(T_m)=693/4+(3/4) sum_(C in R) Wbar_C(u(m)).                  (7)

Proof of the incidence assertion, not just an enumeration: four distinct edge-lines of a quadrilateral have no three concurrent. Indeed any three consecutive edge-lines would otherwise force two of those lines to contain the same two quadrilateral vertices. In the dual Fano plane, a four-point set with no collinear triple is exactly the complement of a line. To see this choose three independent normals a,b,c; the fourth can only be a+b+c, since a+b,a+c,b+c each form a forbidden triple. The complement is {a+b,a+c,b+c}, a line. Therefore the four normals are the affine hyperplane b'.a=1 for a unique nonzero b'. Conversely four lines with no three concurrent have six distinct pairwise intersections. Each of the three unoriented Hamilton cycles of their line-index K4 gives a distinct quadrilateral of intersection points, and every quadrilateral with those four edge-lines occurs this way. This proves both the 7 supports and multiplicity3. By (4), Wbar_C is the product of the four corresponding z_i. Summing and applying (5) gives (7). QED.

The other 84 simple C4s have line-parity weight two: one line occurs twice, and two others once. Two repeated lines cannot support a four-cycle because distinct Fano lines meet only once. Thus the rainbow selection is exactly what removes the lower-degree line characters; summing all C4s would not recover (5).

**What is and is not new here.** The Fano Hamiltonian Fourier identity and the generic C4 switching invariant are prior. The explicit rainbow-loop realization, its incidence proof, the factor3, and the precise transverse gauge extension connect those two prior objects. This is a repository-local synthesis, not a claim of literature-first discovery.

## 4. A constrained gauge model and three distinct symmetries

Let L subset F2^21 be the span of the seven full Fano line edge sets, and B the cut space of K7. Then L intersect B={0}. A cut meets any triangle in zero or two edges, whereas a vector of L meets each designated line triangle in zero or three. Any common vector is therefore empty.

Consequently all128 Fano tournaments occupy DISTINCT switching/gauge orbits. The gauge-saturated family S=u^0+(L+B) has 2^13=8192 fields, with128 orbits of64 points each. The map m -> gauge orbit of u(m) is injective. The 21 rainbow Wilson sum gives a well-defined invariant extension Htilde on S agreeing with H on the chosen slice, but Htilde is not H elsewhere in S. The 189->135 hostile demonstrates this difference.

The set S itself is defined by eight independent gauge-invariant parity constraints: annihilate the 13-dimensional subspace L+B in F2^21. Every vector in (L+B)^perp belongs to the cycle space B^perp and has even parity on each of the seven line edge sets. A basis gives eight Wilson-character equations chi_z(u/u^0)=1. This is a genuine finite constrained Z2 gauge system with an explicitly selected background, not the unconstrained triangular Wilson model below.

Its exactly soluble partition function for the invariant energy -Htilde is

Z_S(beta)=sum_(u in S) exp(beta Htilde(u))
        =64[16 exp(189 beta)+112 exp(171 beta)].                (8)

Equivalently if the action couples to R(u)=sum_R Wbar_C(u), then R=21 on16 orbits and R=-3 on112, so Z_S^R(beta)=64[16 exp(21 beta)+112 exp(-3 beta)]. Normalized Haar/counting conventions divide these expressions by8192.

Do not identify the following groups:
1. Vertex switching B has dimension6 and acts on all21 edge bits; it leaves Wilson loops invariant but not tournament H.
2. Independent Fano line flips L have dimension7, preserve regularity within the selected128 tournaments, and change the Wilson loops.
3. The Session33 Hamming-code translations C=ker(syndrome) have dimension4 inside line-mask space. They preserve the seven displayed quartic characters and H but are NOT vertex-switching gauge transformations. L intersect B=0 proves the distinction sharply. For example a nonzero code translation changes some designated line triangle flux, so it cannot be an edge gauge transformation.

The eight syndrome states forget four line-mask bits. They are sufficient for that Hamiltonian consumer, not for all Wilson loops: the seven designated line triangle ratios individually recover all seven z_i.

## 5. Standard finite triangular Wilson action and its surface expansion

On the full simplex's 2-skeleton let

Z_n(beta)=sum_(u in {+-1}^E) exp(beta sum_(f triangle) W_f(u)).

Write F=binom(n,3), t=tanh(beta). Expanding each factor cosh(beta)(1+t W_f) and summing independent edges gives

Z_n(beta)=2^E cosh(beta)^F sum_(S subset triangles: boundary S=0) t^|S|. (9)

For a cycle c,

E_beta[W_c] = [sum_(S: boundary S=c) t^|S|]
          /[sum_(S: boundary S=0) t^|S|],                      (10)

Boundary is over F2; surfaces may be disconnected or singular. Nothing requires smooth embedded surfaces.

The closed-surface space has dimension binom(n-1,3). A basis consists of boundaries of tetrahedra containing a fixed root: their non-root faces make them independent; coning contracts every closed surface to their span. Thus tetrahedral constraints and the character expansion are two exact aspects of simplicial boundary/coboundary duality.

There is no nonempty closed surface with fewer than four faces; the four-face ones are exactly tetrahedron boundaries. Hence the normalized expansion starts

Z_n/[2^E cosh(beta)^F]=1+binom(n,4)t^4+terms of higher degree.   (11)

For completeness the support argument: select one triangle abc. Its three edges must each occur in another face; no distinct triangle can share two of them. With only three further faces, they must have forms abx, bcy, caz. Closure of ax forces z=x, and closure of bx forces y=x, giving the tetrahedron boundary.

Exact controls:
- n=4: Z_4=64 cosh(beta)^4(1+t^4); E_beta W_C4=2t^2/(1+t^4).
- n=5: normalized Z=1+5t^4+10t^6; the fixed C4 numerator is 2t^2+5t^4+8t^6+t^8.
- n=6: normalized Z=1+15t^4+60t^6+240t^8+392t^10+240t^12+60t^14+15t^16+t^20.

These are three DIFFERENT fourth-order meanings: four edge factors in a Wilson quadrilateral, four independent Fano line-bit factors after restriction, and four face factors in the first closed-surface contribution. Theorem2 makes the first two coincide in one precise carrier. Equation11 has a different grading and must not be conflated with them.

## 6. What survives for nonabelian compact groups

For a compact group G assign U_ij in G with U_ji=U_ij^{-1}, and gauge transform U_ij -> g_i U_ij g_j^{-1}. A based triangle holonomy F_ijk=U_ij U_jk U_ki transforms by conjugation at i, so a character or trace is invariant. A root gauge sets U_0i=1; the remaining binom(n-1,2) links give rooted-gauge coordinates, with residual simultaneous conjugation if quotienting by the FULL vertex gauge group. In particular do not divide cardinalities naively by |G|^n: nonabelian stabilizers vary. For finite G the root-fixed gauge subgroup acts freely and has |G|^(n-1) elements.

The exact transported tetrahedral identity is

F_ijk F_ikl = (U_ij F_jkl U_ji) F_ijl.                        (12)

Both sides expand to U_ij U_jk U_kl U_li. This proves a nonabelian Bianchi relation with transport and order. Dropping transport or commuting factors is wrong in general. The checker exhausts all46656 S3 edge fields on a tetrahedron; the naive untransported identity fails on23328, while (12) always holds. It finds216 flat fields, agreeing with |S3|^3. Flat fields on the full simplex are pure gauge by the same root fixing.

For any finite-dimensional unitary representation rho, a finite simplicial action beta sum_f Re tr rho(F_f) and Haar integration define a finite partition function. Compactness makes the integral finite. These statements require no continuum construction. A nonabelian character expansion can be pursued with representation labels and invariant tensors at shared edges, but this note claims only the explicitly proved Z2 surface formula, not a complete nonabelian spin-foam theorem.

## 7. Source boundary and next decisive work

Public primary source checked on2026-10-23: Wilson, *Confinement of quarks*, Phys. Rev. D10 (1974), DOI https://doi.org/10.1103/PhysRevD.10.2445. Its publisher abstract supports the historical scope: exact lattice gauge invariance, computable strong-coupling limit, and surface expansions; it also explicitly distinguishes that regime from Euclidean invariance. Only the public abstract was used, not the paywalled paper or any PDF.

Clay's official page https://www.claymath.org/millennium/yang-mills-the-maths-gap/ was checked and explicitly labels the continuum existence/mass-gap problem unsolved. We neither import an unreviewed claimed solution nor convert a finite Boolean eigenvalue into a physical mass. No PDF or user-attachment content was read in this lane.

A simplex K_n lacks a specified refining four-dimensional spacetime mesh, lattice spacing, reflection-positive transfer operator with a controlled scaling limit, renormalization estimates, or Osterwalder–Schrader/Wightman reconstruction. Increasing n changes complete-graph coordination, not automatically physical resolution. All such continuum obligations remain open and outside the finite claims here.

Useful next finite target: determine the full unconstrained K7 rainbow-loop partition function and compare it with the constrained8192-field sector. It would measure precisely what the Fano slice excludes; it would not by itself solve any continuum problem. A second target is a chosen sequence of genuine4-dimensional triangulations with local gauge action and an explicit, checked positivity mechanism; no positivity follows merely from calling the building blocks simplices.

## 8. Reproduction

Run `python3 check.py > exact.json` and `python3 -O check.py > optimized.json` from this directory, then compare files. Both outputs agree. The program uses explicit failure exceptions, not optimizable assertions. It performs61275 checks: exhaustive edge fields through n=5 for curvature/root reconstruction and gauge-averaged Hamiltonian count; independent formal polynomial partition identities; n=6 tetrahedral surface enumeration; all128 Fano orientations and all8192 saturated fields; all46656 S3 tetrahedral fields. Fano H is computed directly by subset DP rather than by the claimed Wilson identity. The 21-cycle incidence assertion has the independent geometric proof above.
