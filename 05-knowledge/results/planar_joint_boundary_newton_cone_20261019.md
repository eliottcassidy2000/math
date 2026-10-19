# A joint-boundary completion extends polynomial rigidity to a Newton cone

2026-10-19. **PROVED + INDEPENDENTLY REVIEWED using the cited fixed-line theorem:** finite-word polynomial rigidity for the Newton cone below. It strictly enlarges the family of [packet 113](planar_boundary_preserving_hamiltonian_cone_20261019.md). **OPEN:** rational realizations, Hamiltonians outside this cone, unrestricted infinite words, source-chart entry, and planar JC. No maximality or priority claim is made.

## 1. The larger family and the missing completion

Use the original source K[x,t] over a characteristic-zero field, q=1+x²t, and {f,g}=f_xg_t−f_tg_x. Assign monomial weight

    w(x^i t^j)=2j−i.

Let R∈K[x,t] satisfy

    every nonzero monomial x^i t^j in R has i≤2j+2.    (1)

Equivalently its minimum weight is at least −2. Put

    H=c+t²qR,   δ={−,H}.

For arbitrary finitely many such Hamiltonians and scalar times, let C be the finite composition of their source t-adic exponentials.

**Theorem.** If C(x),C(t) are polynomials in the original source coordinates, then C(x)=x and C(t)=t. Different invariants are allowed; words which cancel to the identity are permitted.

The family includes packet 113, since a monomial of A(xt,t) has i=a,j=a+b and hence satisfies (1). It also includes R=x and R=x², neither of which lies in K[xt,t]. The old polynomial-coefficient common ring need not support their scalar exponentials: for R=x², the common Hamiltonian in s=xt,τ=t is s²(s²+τ)/τ. Its components are δs=−s⁴/τ and δτ=−4s³−2sτ, which do not raise τ-order. A different completion retaining both the source and the second boundary is therefore essential.

## 2. A precise joint-boundary series ring

Define J to consist of formal sums

    F=∑_{i,j≥0} c_ij x^i t^j

whose support has **finitely many monomials of weight at most M for every integer M**. Set ν(F) to its minimum supported weight, with ν(0)=+∞.

This support condition has the following consequences.

1. Every nonzero element has weight bounded below. For each fixed t exponent j, there are finitely many x exponents: all those monomials have weight at most 2j. Thus J maps coefficientwise and injectively into K[x][[t]].
2. Multiplication is well-defined. If F,G have lower weight bounds a,b, every pair contributing to product weight at most M belongs to the finite sets of terms of F of weight at most M−b and G of weight at most M−a. Thus products satisfy the same local-finiteness condition, with ν(FG)≥ν(F)+ν(G).
3. J is complete for the weight topology. In a Cauchy sequence, each weight-bounded finite truncation stabilizes. The resulting support remains finite below each bound because it agrees there with a sufficiently late element. Polynomials are dense by finite weight truncation.
4. The embedding

       ι:J→K[x][[t]]

   is continuous. Weight at least N implies 2j−i≥N and i≥0, hence j≥ceil(N/2). So a sequence tending to zero in weight tends to zero in source t-order.
5. Boundary evaluation

       ε:J→K((s)),   x↦−1/s, t↦−s²                 (2)

   is a well-defined continuous homomorphism. A monomial becomes (−1)^(i+j)s^(2j−i). Every coefficient receives finitely many terms, only finitely many negative powers occur, and ord_s ε(F)≥ν(F). The same finite-support argument proves multiplicativity.

The map ε is not claimed to be injective; q lies in its kernel. The map ι is injective, and that is what identifies a proposed polynomial realization with the same joint-boundary series before evaluating it.

The weight topology is a linear multiplicative topology, not an ideal-adic topology with x of nonnegative order. Multiplication by a fixed element shifts the weight bound by its finite lower weight. No ideal-adic interpretation is used.

## 3. All scalar-time operations and finite words converge

Differentiation with respect to x raises weight by one, while differentiation with respect to t lowers it by two. Since q is homogeneous of weight zero, the Hamiltonian H−c=t²qR has minimum weight at least two. Therefore

    ν(δF)≥ν(F)+1                                   (3)

for every F∈J. Polynomial differentiation and multiplication extend continuously to J, so δ is a continuous derivation. Its scalar exponential exists on J at every λ∈K, obeys the additive exponential law, is multiplicative, and has inverse obtained by negating λ. The identities follow from the Leibniz and binomial formulas with only finitely many terms below each fixed weight bound, and then by completeness.

Independently H−c is divisible by t² in the original source ring. Hence

    δx∈tK[x,t],  δt∈t²K[x,t],
    δ(t^nK[x][[t]])⊂t^(n+1)K[x][[t]].               (4)

The source t-adic exponential is consequently defined and fixes t=0 pointwise. Because ι commutes with the polynomial derivation and is continuous by §2, it matches the two exponentials term by term. The same compatibility holds for every finite composition C.

In particular C(x),C(t) are elements of J. If their source images are polynomials X,T, injectivity of ι gives the corresponding equalities already in J. This resolves the polynomial-to-boundary transport even when the earlier τ-adic Laurent pole orders are unbounded.

Every Hamiltonian exponential preserves the original source bracket, so C has formal Jacobian one. If X,T are polynomials, their determinant-one and fixed-line identities are polynomial identities by the source embedding.

## 4. The boundary derivation and its own convergent flow

For a single R, define the Laurent polynomial

    v_R(s)=s⁴R(−1/s,−s²).

Condition (1) implies ord_s v_R≥2 unless v_R=0. On source polynomials the chain rule gives

    εδ = Δ_R ε,   Δ_R=v_R(s)d/ds.                   (5)

For example, on q=0 the derivatives of H give

    ε(δx)=s²R(−1/s,−s²),
    ε(δt)=−2s⁵R(−1/s,−s²),

which are exactly v_R times the s derivatives of −1/s and −s². This proves (5) on the generators; density and continuity extend it to J.

The derivation Δ_R raises s-order at least one. Its exponential is therefore defined on K((s)), and termwise continuity gives

    ε exp(λδ)=exp(λΔ_R) ε.

For a finite word, the corresponding ordered boundary exponentials intertwine in the same way. In particular

    ε(C(q))=0.                                      (6)

This is genuine preservation of the second boundary through a retained, convergent boundary dynamics. It is not a q-adic inference in K[x][[t]], where q is a unit.

## 5. A polynomial realization is forced to be identity

Pass to an algebraic closure. By §3, a polynomial realization F=(X,T) has determinant one and fixes the line t=0 pointwise. Gwoździewicz's classical [*Injectivity on one line*, Theorem 1.1 / Lemma 2.3](https://arxiv.org/html/alg-geom/9305008) makes F a polynomial automorphism.

The automorphism preserves the prime ideal (t), so T=at. The determinant identity gives aX_x=1, and pointwise line fixing then gives

    T=t,  X=x+f(t),  f(0)=0.

Apply ε to C(q)=1+X²T using the faithful joint-boundary realization. Equation (6) becomes

    f(−s²)[2−s f(−s²)]=0.

The second factor has constant term two, hence cannot vanish. The first factor must vanish; substitution t↦−s² is injective on the polynomial ring, so f=0. Identity descends to K. ∎

Only the cited fixed-line theorem, not the general Jacobian conjecture, is used. A common invariant of the different flows is not assumed.

### Convergent polynomial limits of finite words

**Corollary.** Let C_n be any sequence of finite words in the stated family, and suppose C_n(x) and C_n(t) converge in J's weight topology to X and T. If both X,T are source polynomials, then X=x and T=t. The words, their lengths, and their scalar times may vary with n.

**Proof.** Differentiation with respect to x and t is continuous, shifting weight bounds by +1 and −2 respectively. Multiplication is continuous. Therefore the determinant-one identities for C_n pass to X_xT_t−X_tT_x=1 in J. Weight convergence implies source t-adic convergence, so the fixed-line identities pass to X(x,0)=x and T(x,0)=0. Finally, continuity of multiplication and ε passes ε(1+C_n(x)²C_n(t))=0 to ε(1+X²T)=0. The fixed-line theorem and the shear argument in §5 apply unchanged. ∎

This excludes polynomial limits in this explicitly specified closure. It does not assert that an arbitrary infinite word converges, that its partial products converge in J whenever they converge in a weaker topology, or that a nonpolynomial limit has an inverse.

## 6. Explicit new controls and the edge of this proof

For R=x², which lies on the boundary i=2j+2 of the Newton cone,

    H=x²t²(1+x²t),   v_R(s)=s²,
    exp(λΔ_R)(s)=s/(1−λs).

Thus the enlarged theorem genuinely includes nontrivial rational motion on the second boundary, while still excluding a nonidentity polynomial map of the whole source plane.

For R=x, the boundary velocity is −s³. These examples lie outside the previous A(xt,t) family. The checker retains their source and boundary vector fields rather than inferring inclusion from a name or a finite carrier.

At R=x³, outside (1), one instead has ν(H−c)=1 and v_R=−s. The Hamiltonian derivation need not increase weight, and its boundary derivation no longer increases s-order. Its weight-completed and boundary scalar exponentials are therefore not justified by the order-raising arguments here; its source t-adic exponential still exists by (4). This is an exact limit of this proof, not a counterexample to the polynomial rigidity conclusion or a claim that analytic or other formal realizations cannot exist.

The second-boundary assumption still matters. Dropping q allows H=t² and the exact nonidentity polynomial shear (x,t)↦(x+2λt,t), although its source flow converges and fixes the line. No claim is made that the Newton cone is the largest possible rigid class.

## 7. Reproduction and limits

The [standalone checker](../../04-computation/experiments/planar_joint_boundary_newton_cone_20261019.py) tests the Newton weights, source convergence, joint-boundary evaluation, induced derivations, finite partial words, derivative continuity, strict family extensions, and the x³ edge case. The [normal JSON](planar_joint_boundary_newton_cone_20261019.json) and [optimized JSON](planar_joint_boundary_newton_cone_20261019_optimized.json) agree byte for byte on **726 exact checks**; the [run receipt](planar_joint_boundary_newton_cone_20261019_run_receipt.json) records both commands. The [independent review](planar_joint_boundary_newton_cone_20261019_independent_review.md) accepts the completion construction, finite-word theorem, and convergence-conditional closure corollary, and independently reproduces both checker modes. Universal claims rest on the proofs above and the cited classical theorem, not on finite experiments.

The theorem covers finite words and the corollary covers only polynomial coordinate limits when convergence in J is established. Existence of arbitrary infinite compositions, nonidentity rational source outputs, chart entry, and planar JC remain open. The joint-boundary ring is an explicit completion chosen to retain precisely the two interfaces needed for this argument.
