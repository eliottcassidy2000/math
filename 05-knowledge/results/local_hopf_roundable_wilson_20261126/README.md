# Local Hopf, nested repair, and a roundable Moore action

2026-11-26. **SOURCE CLAIM** for the linked curvature construction; **PROVED** below for the abstract fourth-order test and discrete rounding theorem. **FINITE-EXACT** controls accompany the latter. Neither the degree-57 Moore graph nor Chan–Pak is settled here.

## 1. What the new construction changes

Read against the [pinned Local Hopf proof](https://github.com/vixhal-baraiya/local-hopf-problem/blob/3d0a9efeb36155d9d7a8ea02e6eb950815c28ece/proof/PROOF.md), [tensor derivations](https://github.com/vixhal-baraiya/local-hopf-problem/blob/3d0a9efeb36155d9d7a8ea02e6eb950815c28ece/proof/DERIVATIONS.md), and [archived version](https://doi.org/10.5281/zenodo.24534221). The repository claims more than existence of a positive metric: convergence in smooth topology to the round product. The older [Brendle–Hung preprint](https://arxiv.org/abs/2608.19068) is now at v2; the repo's [August audit ledger](../../reference/CORE-PAPERS-BRENDLE-HUNG-S2XS2-2026-08-24.md) concerns v1. Neither the old notebook audit nor this source intake is an independent replay of the new full metric certificate.

The important new issue is the two scales. For each fixed Cheeger parameter t>0, minimize over normal plane tilts before inspecting the s² and s³ coefficients. As t decreases a tilt stiffness degenerates like t³. The linked derivation reports that keeping the original tensors produces a negative quadratic coefficient, with t³ times that coefficient tending to −16/25 at a fixed regular pair. Testing only until the s-linear term vanishes *on* the flat stratum would miss precisely this failure. The new weights and tensors precondition the tilt; the source reports a quadratic bound (3/10)ρ for 0<t≤1/100, where ρ vanishes on the equatorial torus. On that torus the surviving cubic mean is 1/[40000(1+2t)²] in the numerator normalization. The Poisson repair removes the nonconstant part, not the mean. The source leaves the admissible s₀(t) uncomputed and takes a diagonal sequence; the two bounds just quoted are not by themselves a uniform bound on s₀(t).

## 2. What to keep from past Hopf sessions

* **Repair quotient, August:** [THM-3990](../../../01-canon/theorems/THM-3990-componentwise-harmonic-obstruction-and-repair-quotient.md) already isolates the componentwise means as the real obstruction for Laplacian repair. This explains why computing just the torus mean can be enough, provided the completed zero stratum and admissible conformal extension have been identified. The positive mean is not a pointwise curvature bound before the repair. In an integral graph analogue torsion can obstruct exact correction even when every real mean is zero.
* **Completed strata, August:** the separate [S⁶ source ledger](../../reference/CORE-PAPERS-HOPF-S6-2026-08-24.md) and [two-Hopf synthesis](../../../07-reflections/two-hopf-problems-repair-quotients-factorial-cusps-and-frontier-extension-root-20260824.md) supply the useful warning: a normal form on an open chart is not a proof at a completed seam. The local proof handles the coincident and antipodal pairs by a global flat-plane manifold and dense regular charts with continuous reduced coefficients. This is not a statement about integrability on S⁶, nor does it import that manuscript.
* **Phase and transport, March to October:** [early Cayley calculations](../../../04-computation/cayley_hypercomplex_s116i.py) and [October shared incidence](../session40_marked_ports_square_sum_20261023/hopf-incidence/HOPF_SIERPINSKI_SHARED_INCIDENCE.md) motivate keeping chart-transition data. The ordinary quaternion description S³→S² has circle stabilizer; it is not a principal SU(2) bundle over S². October already distinguishes holonomy on a graph from a Chern class, and actual shared transports from phase marginals. The Moore extension below likewise preserves the individual routes, not just their mean occupancy.
* **Hopf algebra, March to October:** the [deletion–contraction checkpoint](https://github.com/eliottcassidy2000/math/commit/132bb761) and the [primitive Hopf–Lie lift](../primitive_hopf_lie_bijection_20261020/PRIMITIVE_HOPF_LIE_BIJECTION.md) concern a different meaning of Hopf. Ordered joins and subset restrictions supply useful exact local operations; the nonmultiplicativity of normalized Walsh under join is a warning against identifying a convenient linear observer with a compositional invariant. No curvature-to-Hamiltonian identity is inferred.

## 3. A fourth-order gate after two normal minimizations

This extends the minimization method, not the claimed curvature coefficients of the source. Suppose all data are smooth with uniform bounds on the compact tubes under consideration, and write the local family as

    F(s z,s) = s² [P₂(z)+sP₃(z)+s²P₄(z)+O(s³)]

with the expansion also valid for the derivatives needed to minimize. Here P₂ has positive definite Hessian H and unique minimizer a=−H⁻¹∇F₁(0) in the rescaled normal coordinates. As in the source, assume the zeroth order has a quadratic normal lower bound and the first coefficient is zero on the base; then the normal minimizers have original displacement O(s) and are interior. Their rescaled position is

    z(s) = a − s H⁻¹∇P₃(a)+O(s²).

This follows by differentiating stationarity, or by the implicit-function theorem. Substitution gives the next coefficient of the first reduced minimum:

    U₄ = P₄(a) − (1/2) ∇P₃(a)^T H⁻¹∇P₃(a).

Here P₄ includes the quartic jet of F₀, cubic jet of F₁, quadratic jet of F₂, linear jet of F₃ and constant F₄. Evaluating the original fourth coefficient on the unshifted flat planes omits the negative correction.

Now suppose U₂≥0, its zero set Σ is a compact embedded submanifold of the first base, and U₂ has positive definite normal Hessian K there. Assume that the repaired U₃ is identically zero on Σ. Minimize again on the normal fibres of Σ, this time in the family divided by s². Away from Σ compactness gives a positive quadratic minimum. Near Σ the resulting minimum is

    s⁴ [ U₄|Σ − (1/2) (∇⊥U₃)^T K⁻¹∇⊥U₃ ] + O(s⁵).

Uniform strict positivity of the bracket gives positivity for small s>0; negativity anywhere gives a negative test plane after the two shifts. The equality case needs higher order. This is a second, independent negative correction, not another copy of H. For a scalar control with F(w,x,s)= w²/2 + s² x²/2 + s²bw + s³cx + s⁴d, the exact minimum is s⁴(d−b²/2−c²/2). The source has a positive cubic mean and does not need this fourth-order test. The test matters for reflection-symmetric alternative tensors or any proposal whose cubic repair has zero mean on its entire remaining stratum. Before using it on another metric, verify normal coercivity of the completed quadratic stratum, not just sharpness on a regular slice.

**Exact symmetry obstruction for alternate tensors.** Let Σ be the product of equators, fixed by simultaneous equatorial reflection. Any metric keeping that reflection as an isometry makes Σ totally geodesic, so Gauss–Bonnet on the torus gives integral zero for ambient sectional curvature on its tangent planes. Such a family cannot succeed at *any* order, whether or not a quartic jet looks positive on a regular chart. More quantitatively, with II the second fundamental form of the fixed embedded torus for an arbitrary metric, the Gauss equation gives

    2 ∫Σ sec(TΣ) = ∫Σ (|II|² − |tr II|²).

This follows by expanding the two-by-two determinant of each normal shape operator and using ∫Σ KΣ=0. Thus a lower bound sec≥κ>0 requires ∫Σ |II|²≥2κ·area(Σ); a successful correction has to create actual extrinsic bending, not merely redistribute intrinsic torus curvature. This is a different obstruction from the normal-minimum correction (which also checks nearby, tilted planes).

## 4. The discrete test: repair Wilson's relaxed zero by keeping route variance

Take the [previous Wilson chart](../moore57_quadratic_tensors_20261125/WILSON-ACTION-AND-RIGID-MOORE-CHARTS.md) with k branches and d=k−1 sheets, and now allow a doubly stochastic link Pᵢⱼ, Pⱼᵢ=Pᵢⱼᵀ, with prescribed matchings fixed and with any additional forbidden entries held at zero. Let Mᵢⱼ be the d route matrices (one direct link and k−2 two-hop products), Rᵢⱼ their sum, and δ(M)=d−‖M‖²_F. The ordinary soft Wilson action E=Σ‖Rᵢⱼ−J‖² can vanish even on the completely diffuse chart Pᵢⱼ=J/d. It is not a rounding certificate.

**Exact repair theorem.** Use
     
    W(P) = E(P) + Σ_{i<j} Σ_{M in Mᵢⱼ} δ(M)
         = 2 Σ_{i<j} Σ_{distinct unordered routes M,N} ⟨M,N⟩_F .

The second identity follows by expanding the square: each route has total entry sum d, there are d routes, and ‖J‖²=d². In a fixed branch-pair the direct and the two-hop routes depend on pairwise disjoint sets of unordered link variables. Thus W is affine in each individual link with all other links held fixed, despite the apparent repeated-link powers in E. It is invariant under independent sheet relabeling, and equals 6T+4S on integral charts.

For the probabilistic reading, use independent Birkhoff draws of the free links with means Pᵢⱼ. The product on a two-hop route is a permutation with mean given by its matrix product. The routes for each output pair are independent even though *different* output pairs need not be. Their entrywise total variance is exactly the added term, so W is the expected discrete action. At each step a linear assignment on the current link support selects a permutation with W no larger than before; Birkhoff decomposition uses only support entries. Thus it preserves all previously fixed matchings and zeros. Repeating gives an integral completion with action ≤W(P); **W(P)<4 is a sufficient exact Moore certificate**, since any nonzero integral action is at least 4. The relaxed global minimum equals the integral global minimum whenever the specified link faces are nonempty. This is not convexity, and it gives no guarantee that blockwise descent finds a global minimizer.

Writing D=Σ_{i<j}δ(Pᵢⱼ), the added variance satisfies

    (k−1)D ≤ W−E ≤ (2k−3)D .

Indeed stochastic matrices contract the Euclidean norm (Jensen), so δ(PQ)≥max(δ(P),δ(Q))≥(δ(P)+δ(Q))/2. Also

    δ(PQ)=δ(P)+δ(Q)−tr[(I−PᵀP)(I−QQᵀ)] ≤ δ(P)+δ(Q),

since the last two matrices are positive semidefinite. Each link is counted once directly and in 2(k−2) two-hop terms. The completely diffuse zero of E actually has W=binom(k,2)d(d−1), or 4,915,680 for k=57. In that case every one-link assignment is tied by symmetry, so an implementation needs actual integral conditioning, not merely a vanishing projected gradient.

**First variation at a zero.** For a true chart P₀ and segment P(s)=(1−s)P₀+sQ, let L=Σ(d−⟨P₀ᵢⱼ,Qᵢⱼ⟩). Since E has zero first derivative there, and the trace correction in δ(PQ) is quadratic in s at permutations,

    W(P(s)) = 2(2k−3)L s + O(s²).

At k=57 the slope is 222 per unit of moved sheet mass. This recovers the previous exact one-link switch cost and shows why looking only at E's vanishing first derivative misses discrete stiffness. It is a quantitative analogue of the linked Hopf proof's lesson about checking the *reduced* coefficient, not an identification of sectional curvature with combinatorial collisions.

## 5. Next concrete tests

For curvature: if attempting a rate s(t) instead of the two-parameter diagonal sequence, bound H⁻¹, the gradient of the cubic rescaled polynomial, and the Hessian along the completed torus in a common tube; the fourth-order losses specify the required derivative ledger. None of the threshold rates are claimed here.

For Moore: freeze the three actual stars and the marked contact, set unsupported edges to zero, compute linear-assignment costs from pairwise route overlap, and track W and E separately. Any decrease should be translated back into actual occupied slots. A fully fractional or coordinate-stationary outcome is not a witness. Run `node wilson_rounding.js` for the all-uniform control, exact blend expectations, all 318 single-transposition slopes, and support-restricted rounding examples on Petersen and Hoffman–Singleton. The saved output is in `wilson_rounding.out`. Degree 57 remains conditional.
