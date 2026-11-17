# What a Berggren transplant preserves, and what a dimension analogy loses

**Status:** elementary proved lemmas and explicit controls; sourced boundary statements; research proposal only for the final bridge. No proof of Chan–Pak or planar Jacobian completeness is claimed. The linear-preserver observation is an elementary application of the Hessian, not a priority claim.

## 1. The genuine quadratic boundary

Let K have characteristic zero, r≥2, d≥3, and

    f(x)=sum_i c_i x_i^d,     each c_i≠0.

**Lemma (diagonal-form linear rigidity).** Every L∈GL_r(K) satisfying f(Lx)=f(x) is monomial: it permutes coordinate axes and rescales them. In particular, if c_i∈{1,−1} and L is integral unimodular, its entries are a signed permutation, so the group is finite.

**Proof.** The Hessian determinant is C product_i x_i^(d−2), C≠0. The chain rule gives

    det(L)^2 det Hess(f)(Lx)=det Hess(f)(x).

Consequently the product of the r row-linear forms of L, each with multiplicity d−2, is a nonzero scalar times product_i x_i^(d−2). Unique factorization makes each row proportional to one coordinate. Invertibility forces a permutation. For an integral unimodular monomial matrix the nonzero entries are nonzero integers whose product has absolute value one; each is ±1. Substitution into f determines which signs and permutations are allowed. □

For d=2 the Hessian determinant is constant, and this argument imposes no hyperplane arrangement. The integral matrix

    A = ((1,−2,2),(2,−1,2),(2,−2,3))

has determinant one and preserves Q=x²+y²−z². It sends (3,4,5) to (5,12,13). Its iterates on the positive null cone have increasing hypotenuse, yielding an infinite group and positive orbit. For example z'=2x−2y+3z>z because x+z>y.

This boundary is **degree two versus higher degree of this exact diagonal invariant**, not ambient dimension two versus three. The same A acting on (x,y,z), with w fixed, preserves x²+y²+w²−z². It sends the positive null vector (x,y,w,z)=(1,2,2,3) to (3,6,2,7). Conversely A does not preserve x³+y³−z³ or x⁴+y⁴−z⁴: (1,0,1) goes to (3,4,5).

[Fermat's theorem](https://annals.math.princeton.edu/1995/141-3/p01) excludes nonzero positive integer solutions of x^d+y^d=z^d for d>2. Here d is an exponent, not the dimension: the affine hypersurface in three variables has dimension two for every such d. Changing the number of summands changes the question; 3³+4³+5³=6³ is an immediate control. The rigidity lemma does NOT exclude nonlinear recursion, piecewise maps, maps between different forms, or maps preserving only selected integer solutions.

**Constructive consequence.** A proposed higher-degree Berggren system must escape at least one hypothesis: exact preservation of one diagonal form, linearity, fixed state space, or unimodular integral action. A lifted state with quadratic invariants can be legitimate, but its projection and integrality constraints become part of the theorem, not disposable coordinates.

## 2. The real/complex Jacobian boundary is not optional

Pinchuk maps have positive, nonconstant real Jacobian determinant and can fail injectivity. This is not a counterexample to the constant-Jacobian polynomial question over C. See [Campbell, Pinchuk maps and function fields](https://arxiv.org/html/1202.2949v4). Here “planar” R² has real dimension two, whereas C² has real dimension four. Even the much simpler real bijection (x,y)↦(x+x³,y) has determinant 1+3x² and no polynomial inverse; its complexification has critical points.

**Lemma (critical-divisor persistence).** Suppose a polynomial map F has nonconstant Jacobian determinant. Neither composition on either side by polynomial automorphisms nor adjoining identity coordinates can make its Jacobian determinant constant.

**Proof.** An automorphism has constant nonzero determinant: the chain rule applied to its polynomial inverse makes its determinant a unit in a polynomial ring. Thus for automorphisms α,β,

    J_(α∘F∘β)=c·(J_F∘β),    c≠0.

A ring automorphism cannot take a nonconstant polynomial to a constant. Identity stabilization leaves J_F unchanged. Over C a nonconstant polynomial has a zero, so its complex critical divisor cannot vanish in this category. □

Restriction to a rational chart may discard that divisor, but it also discards the original global polynomial problem. This is exactly the sort of forgotten condition that a recursive structure classifier must retain.

**Lemma (stabilization is upward, not downward).** F×id_k is injective, or polynomially invertible, or proper precisely when F has the corresponding property. Left/right polynomial equivalence preserves each property. If F is dominant generically finite, its generic degree is also unchanged by stabilization.

**Proof.** Injectivity follows by holding the new coordinates fixed. If G is a polynomial inverse of F×id, its first coordinates restricted to the slice with extra coordinates zero give a polynomial inverse of F; the other direction is explicit. For properness, the inverse image of a compact set lies in F^{-1} of its compact first projection times its compact second projection and is closed; conversely use compact sets K×{0}. Polynomial automorphisms are proper homeomorphisms. Generic degree is the degree of the function-field extension; adjoining independent indeterminates preserves that degree. □

Hence higher-dimensional examples cannot be lowered to C² merely by calling the reduction stable. An actual identity-factor splitting, or a different precise descent theorem, is required. Arbitrarily high degree alone is no topological obstruction: (x+p(y),y) is invertible for every polynomial p.

Current source boundary: [Gao, version 2, November 2026](https://arxiv.org/html/2608.00222v2), reports complex constant-Jacobian counterexamples in dimensions greater than two, while explicitly leaving dimension two open. Its stabilization discussion distinguishes actual identity-factor splitting. Its review also recalls that proper complex Keller maps are automorphisms: a proper local biholomorphism covers simply connected C^n with one sheet. This packet does not reprove or independently audit the reported counterexamples.

## 3. A dimension ledger for trees and Gaussian moments

[Random planar trees and the Jacobian conjecture](https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/jlms.70416) uses planar ordered trees for inverse coefficients with vertex types in {1,…,n}. “Planar tree” does not mean a two-variable polynomial map. Its degree reduction concerns all dimensions simultaneously; added dimensions cannot be omitted from a fixed-dimension inference. The January 2026 paper's historical general open-status wording is not used as current status here.

[Derksen–van den Essen–Zhao, §3](https://arxiv.org/html/1506.05192v1), Proposition 3.2 uses GMC(2n), not GMC(n), for the factorial pairing in n pairs of variables. Proposition 3.1 then gives SIC(n); the main theorem is stated with all dimensions. Thus the cited reduction does not by itself supply GMC(2)⇒JC(2). An old repository sentence asserting a same-index implication needs an additional theorem, not merely an analogy. This is a boundary annotation, not an audit of any claimed GMC(2) result.

The useful mechanism to borrow is to record the paired interface and the dimensional cost of eliminating it. Reusing one Gaussian coordinate as two independent resources, just like replacing two endpoint states by their equal scalar Hamiltonian counts, changes the composition law.

## 4. Repository archaeology and a proposed bridge

Read together, these sources identify a useful distinction:

- `05-knowledge/results/pythagorean_form_conic_graphs_rp_margin_20261018.md`: one quadratic form viewed over Z, R and finite fields; the Berggren action is integral Lorentz structure.
- `05-knowledge/results/jcweb_20261020_uniform_rules_branch_integrality_jacobian_web.md`: quotient folding and mixed-weight branching are not automatically global planar Keller maps; it already records the higher-dimensional/planar split.
- `07-reflections/fc-jc2-tandem-the-isolated-coincidence-and-the-fold-no-fold-duality-opus-S4.md`: closing structured subfamilies is not closing a residual global problem.
- `07-reflections/gaussian-moments-the-explicit-failing-mathieu-space-and-why-gmc2-resists-deathstar-S61e.md`: independence of radial and angular resources matters; its same-index GMC/JC wording is not justified by the cited primary reduction above.

**Research principle, not a new mathematical equivalence:** maintain an *obstruction ledger* under every recursive move: exact invariant, arithmetic lattice, forbidden/critical locus, dimension, and physically admissible interface. The theorem in §1 says why one especially tempting linear transplant is unavailable. The theorem in §2 says why a tempting real-to-complex transplant is unavailable. Neither forbids a more expressive recursion.

For Chan–Pak the promising transferable demand is not “all structures have a grammar.” It is a disjoint inverse-chamber decomposition of **target plus interface**, with strict descent and an exhaustive root bank. Berggren supplies all three. A universal tournament grammar supplies neither disjoint scalar chambers nor a realizable predecessor for an arbitrary prescribed odd target. A new lift should be judged by whether its extra coordinates can be solved at each inverse step, not by whether the unconstrained projection looks dense.
