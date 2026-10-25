# Shape transitions preserve spectral constraints only with the operator attached

2026-10-25. **Status: PROVED elementary statements below, independently reviewed in
[AUDIT.md](AUDIT.md); CITED background; CONDITIONAL bounded-proof consequence; FINITE-EXACT
controls.** No historical priority is asserted. This packet does not prove RH,
planar JC, Collatz convergence, or P versus NP, and supplies no implication
between their truth values.

## 1. Inheritance and the question changed by the new object

The October 20 [geometric arithmetic session](../geometric_arithmetic_session_20261020/README.md)
and [shape completion packet](../simplex_sphere_cube_completion_20261020/SIMPLEX_SPHERE_CUBE_COMPLETION.md)
already distinguish marked radial charts, probability amplitudes, polarity,
ellipsoid approximation, and moment completion. The new question here is:
**when can an actual linear operator survive the change from a finite shape to
a round one?** This is stronger than replacing a representation by a picture.

Closest proved mechanism: the prior packet's exact norm-section transition and
its complete character-positive matrix. Canonical hostile: its q=17 locally
triangle-feasible PSD object without a global sign law. Corrected near miss:
one-step agreement or pair PSD is not a complete multiplication-compatible
state. Least-used sidecar here: the operator itself and its positive invariant
form, together with its order.

Session 58's [absolute-geometry packet](../session58_hourglass_complement_geometry_20261025/absolute-geometry/ABSOLUTE_GEOMETRY_BRIDGES.md)
already distinguishes fixed-base extension counts, finite Frobenius cycles,
formal absolute zeta, and actual descent. We inherit that firewall rather than
presenting another counting polynomial as an RH construction.

Board: (1) exact linear invariance, (2) finite vertices versus continuous
rotations, (3) Newton support versus coefficients, (4) positive form versus
reciprocal spectrum, (5) proof verification versus proof existence. Anchor:
operator-preserving shape transitions. Niche: coefficient loss at the actual
Keller predicate. Wildcard: a precise common bounded-proof carrier for the
named conjectures.

## 2. Invariant bodies and ellipsoids are equivalent for a fixed linear operator

Let A be an invertible real d-by-d matrix, d>=1. A convex body means a compact
convex subset of R^d with nonempty interior. The following are equivalent:

1. There is a convex body K with A K=K.
2. There is a centered ellipsoid E with A E=E.
3. There is a real positive-definite symmetric matrix G with A^T G A=G.
4. A is diagonalizable over C and all its eigenvalues have absolute value one.

**Proof.** Under (1), volume gives |det A|=1. The centroid mu of K satisfies
A mu=mu. Its normalized covariance

    Sigma = vol(K)^(-1) integral_K (x-mu)(x-mu)^T dx

is positive definite: the square of every nonzero linear functional is positive
on some open subset of K. Change of variables gives

    A Sigma A^T=Sigma.

Thus G=Sigma^(-1) satisfies A^T G A=G. Its unit ellipsoid proves (2), and
(2) and (3) are equivalent by the definition of an ellipsoid. For (3)->(4),
G^(1/2) A G^(-1/2) is orthogonal. Conversely, under (4), a real basis puts A
into blocks +1, -1, and two-dimensional rotations; transporting the Euclidean
form from this basis proves (3). Finally an ellipsoid is a convex body. QED.

This gives an exact way to alternate shapes. After centering K, let E be its
invariant ellipsoid and put

    K_t=(1-t)K+tE,  0<=t<=1.

Linearity gives A K_t=K_t for every t. The operator and its spectral constraint
survive this Minkowski interpolation. The facets, vertex labels, lattice, and
arithmetic predicates need not survive. The ellipsoid usually has more
symmetries than K; the reverse implication for a newly chosen symmetry fails.

For a fixed real q>0 and real invertible T, apply the theorem to A=T/sqrt(q).
Then an exact body relation T K=sqrt(q) K forces every eigenvalue of T to have
absolute value sqrt(q), and forces semisimplicity. A particularly checkable
certificate is

    T^T G T=qG,  G positive definite.

For an eigenvector v over C, this identity gives
|lambda|^2 v*Gv=q v*Gv, proving the spectral conclusion directly.

**Boundary controls.** The unit-spectrum condition without diagonalizability
is insufficient: J=[[1,1],[0,1]] has J^n=[[1,n],[0,1]], so cannot preserve a
bounded full-dimensional body. Conversely

    T=[[0,-2],[1,1]],  G=[[2,1],[1,4]]

has G positive definite and T^TGT=2G; its characteristic polynomial is
z^2-z+2. These are exact finite-dimensional controls, not an identification
of a Frobenius operator attached to Riemann zeta.

## 3. A finite shape imposes the extra condition of finite order

**Proposition.** A preserves a full-dimensional polytope if and only if A has
finite order.

**Proof.** If A preserves a polytope, it permutes its finite vertex set V. Let m
be the order of this permutation. Then A^m fixes every vertex, and the vertices
span R^d, so A^m=I. Conversely, if A^m=I, the convex hull of

    {+/- A^j e_i : 0<=j<m, 1<=i<=d}

is an invariant full-dimensional polytope. QED.

This is the exact finite-versus-round boundary. The rational rotation

    R=[[3/5,-4/5],[4/5,3/5]]

preserves the unit disk. It has infinite order: if it had finite order, its
eigenvalues would be roots of unity, so its trace would be an algebraic
integer. But tr R=6/5 is a rational noninteger. Hence it preserves no
full-dimensional finite-vertex polytope. Circles cannot always be replaced by
triangles or cuboids while keeping the same operator.

The finite families also differ. A3=[[0,-1],[1,-1]] cyclically permutes
e1,e2,-e1-e2 and has order three. A4=[[0,-1],[1,0]] preserves the square and
has order four. No triangle can support A4, since a faithful permutation of
its three vertices cannot have order four. No centered parallelogram can
support A3: its linear symmetry group is conjugate to the signed permutation
group on two coordinates, of order eight. Both operators nevertheless
preserve ellipsoids by section 2.

Thus a simplex has the smallest possible number d+1 of vertices among
full-dimensional polytopes, whereas a ball has infinitely many extreme points
but a short quadratic description. There is no universal complexity order.
In dimension one all three bodies are intervals. Polarity and Boolean order
duality remain the distinct operations already proved in the inherited packet:
a centered regular simplex is similar to its polar; a ball is self-polar;
a cuboid is polar to a crosspolytope. Cubical order complementation does not
make a three-dimensional cube combinatorially polar-self-dual.

## 4. The missing positivity cannot be replaced by duality

The integer matrix and invariant form

    H=[[2,1],[1,1]],  B=[[2,-1],[-1,-2]]

satisfy H^TBH=B, but det B=-5. The eigenvalues of H are
(3+sqrt(5))/2 and (3-sqrt(5))/2, off the unit circle. Their polynomial
z^2-3z+1 is reciprocal. This is a minimal nontrivial hostile to the inference
"self-dual or reciprocal spectral data implies purity."

The strongest survivor is section 2 with a positive-definite form, or its
equivalent bounded invariant body. The missing coordinate is positivity;
the new question is whether a proposed arithmetic construction produces the
specific positive form for the actual operator and the whole required space.
The already-known Weil-positivity approach motivates this question: the
trace-pairing formulation in Connes–Consani–Marcolli is equivalent to RH in
its stated setting. The finite matrix lemma here neither constructs that
pairing nor verifies its positivity. [Primary paper](https://arxiv.org/abs/math/0703392).

## 5. A genuine F1-adjacent language, with a first coefficient hostile

For nonnegative-coefficient polynomials, the Newton-polytope map satisfies

    Newt(fg)=Newt(f)+Newt(g),
    Newt(f+g)=conv(Newt(f) union Newt(g)).

This is an exact semiring map, with the zero polynomial mapped to the empty
set. For the multiplication identity, exposed vertices cannot cancel; for
addition, nonnegative coefficients ensure union of supports. Passing further
to Newt(f)+R^2_{>=0} gives an upward-closed Newton polygon. Such polygons with
convex-hull union and Minkowski sum occur in the reduced square of the
Connes–Consani arithmetic site (Lemma 6.18 and Proposition 6.21). Their
construction connects arithmetic-site points and Frobenius correspondences
to an adelic realization related to completed Riemann zeta. That is a concrete
shared operation language; it is not an identification of planar Keller maps
with their arithmetic site. [Primary paper](https://arxiv.org/abs/1502.05580).

The coefficient loss is visible even before rounding any polygon. Over a
characteristic-zero field, set

    P=x+y^2,
    Q_t=y+x^2+t*x*y^2+y^4,  t!=0.

All these pairs have the same full monomial supports, and therefore the same
Newton polygons. Direct differentiation gives

    Jac(P,Q_t)=1+2(t-2)*y*(x-y^2).

Consequently the pair is Keller exactly when t=2. At t=2 it is visibly a
composition of polynomial shears: Q_2=y+P^2. At t=1 all coefficients remain
positive and the same support semiring is available, but the Jacobian is
nonconstant. This witness refutes only support sufficiency, not planar JC.

First failed implication: multiplicative support compatibility was promoted
to cancellation compatibility in a determinant. Strongest survivor: Newton
polygons control support and exposed degrees. Repaired form: retain the
coefficient relation t=2 for this family, or the actual coefficient equations
for a general family. Future test: before exporting any shape or tropical
obstruction to a live planar-JC row, vary coefficients inside one support
fiber and test the exact determinant predicate.

## 6. A precise common carrier for conjectures is bounded proof search

Fix a formal proof system with a verifier V(phi,pi) whose runtime is polynomial
in the encoded formula and proof lengths. Define

    BoundedProof={(phi,1^B): there exists pi with |pi|<=B and V(phi,pi)=1}.

**Proposition.** BoundedProof belongs to NP. **CONDITIONAL:** if P=NP, a proof
of length at most B can be found, or its absence within that bound determined,
in time polynomial in |phi|+B.

**Proof.** The witness is pi, bounded by the unary input length; verifying it
has polynomial cost. If P=NP, the prefix-extension language for these proofs
is also in P. First ask whether a witness exists, then choose each successive
bit by an extension query, using a padded fixed-length encoding with an end
marker. There are O(B) queries. QED.

This is the standard NP decision-to-search mechanism applied to proofs; Cook's
official account explicitly states the general search consequence of P=NP.
[Primary source, pages 5–6](https://www.claymath.org/wp-content/uploads/2022/06/pvsnp.pdf).

The actual formula, proof system, verifier, encoding, and bound are sidecars.
Writing B in binary changes the witness-size guarantee. More importantly,
P=NP supplies neither a proof of RH, JC or Collatz in that system nor any bound
on one. A system-unprovable statement has no witness at any B. Even for a
provable statement, polynomial time in B need not be small in the statement's
length. Thus this gives a conditional relation between *search procedures*,
not an equivalence of the mathematical conjectures or their truth values.

## 7. Connection ledger and next decisive tests

| Source -> target | Map and preserved predicate | Loss / necessary sidecar | Decisive hostile |
|---|---|---|---|
| Invariant body -> ellipsoid | Covariance inverse; same A and spectral modulus | Facets, lattice and extra symmetries; keep A | R has an invariant disk but no finite invariant polytope |
| Finite polytope -> spectrum | Vertex permutation; finite-order A | Embedding and vertex incidence; retain representation | Order-four square symmetry cannot act on a triangle |
| Reciprocal form -> purity | No valid inference for indefinite forms | Positive-definite form on the actual space | H preserves B with det B=-5 |
| Polynomial -> Newton polygon | Support semiring operations | Coefficients and determinant cancellations | P,Q_1 and P,Q_2 have identical supports |
| Fixed proof verification -> proof search | Unary-bounded NP witness language | Existence, bound, encoding, system | No witness for an unprovable formula |

Revisiting the board: the finite-order obstruction sharpens the inherited
finite-Frobenius boundary; it explains why a permutation model is too small
for general unit-modulus spectra. The invariant-form theorem connects round
geometry to spectral purity only after its strongest missing hypothesis is
proved. The Keller witness shows exactly why keeping shape and operation is
still insufficient for a coefficient consumer. Bounded proof search gives an
honest cross-problem interface without erasing quantifiers or resource bounds.

These elementary mechanisms do not justify a new universal META-PATTERNS card
on their own: further successful use in distinct live threads is required.

## 8. Exact finite controls

The standard-library script `04-computation/shape_spectral_20261025/verify_shapes.py`
checks 80 exact assertions: two finite vertex actions and their orders;
triangle covariance equivariance; four infinite-order rational rotations;
a positive q=2 similitude; an indefinite reciprocal hostile; Jordan growth;
four same-support polynomial pairs; and 36 independent rational point checks
of the determinant formula. Its explicit exception checks remain active with
Python optimization. It imports no repository mathematics and has no hidden
enumeration filters. The all-order claims have proofs above.

From the repository root:

```bash
python3 -B 04-computation/shape_spectral_20261025/verify_shapes.py
python3 -B -O 04-computation/shape_spectral_20261025/verify_shapes.py
```

`shapes-verification.json` records the deterministic output. This computation
is a control for the displayed exact identities, not a finite test of RH,
Collatz, planar JC, or P versus NP.
