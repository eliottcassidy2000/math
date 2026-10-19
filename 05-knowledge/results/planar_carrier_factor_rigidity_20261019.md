# A two-component rigidity gate closes universal single-carrier polynomial repair

2026-10-19. **PROVED:** the finite-order theorem and the universal single-carrier polynomial-specialization obstruction below. **FINITE-EXACT:** the companion controls. **OPEN:** rational specialization for arbitrary carriers; compositions with different invariants; entry into the source-normal chart; planar JC itself. No priority claim is made.

## 1. The precise inherited gap

The [elliptic bridge](../../05-knowledge/results/jacobian_elliptic_bridge_20261016c.md) ends with a requirement to retain the actual source, punctures, and differential. The [long-session board](../../05-knowledge/results/planar_jc_long_20260906_board.md), under “A sharper termination test,” leaves general source-carrier termination open. The [completed Hamiltonian theorem](../../05-knowledge/results/planar_jc_long_20260906_hamiltonian.md) supplies actual completed flows for every

    S = c + p²D R(p,y),  D=p³−y²,  R∈K[p,y],
    p=t(1+x²t),  y=xtp,  {f,g}=f_x g_t−f_t g_x.

The [nonrational-family theorem](../../05-knowledge/results/planar_jc_long_20260906_nonrational.md) excludes rational times for S=f(p^aD), a≥2, by generic genus. Its general-carrier complement was not covered by that family calculation.

Here the retained object is the **complete zero-fibre divisor of the actual source polynomial S−c**. It has two compulsory components, t=0 and 1+x²t=0. A finite carrier alone still gives no birationality. The new argument uses a hypothetical *actual polynomial Keller realization*, its invariant polynomial, and all irreducible components of that polynomial together. It then proves finite order, rather than assuming an inverse.

## 2. Finite divisor permutation without assuming invertibility

**Lemma 1 (invariant fibre factor permutation).** Let K be an algebraically closed characteristic-zero field, let F:A²→A² be polynomial with nonzero constant Jacobian, and suppose H∘F=H for a nonconstant H∈K[x,t]. Write

    H = u ∏_{i=1}^r h_i^{e_i},

where u≠0 and the h_i are distinct irreducibles. Then there is a permutation σ of {1,…,r} and constants a_i≠0 with

    h_i∘F = a_i h_{σ(i)},   e_i=e_{σ(i)}.

**Proof.** F is dominant and étale, hence quasi-finite. Therefore every h_i∘F is nonconstant; two different such pullbacks are relatively prime: a common curve component would map into the finite intersection V(h_i,h_j), contradicting quasi-finiteness. Each h_i∘F is squarefree. This is also the forward direction of Jędrzejewicz, *A characterization of Keller maps*, Theorem 5.1 and Corollary 4.2 ([primary author text](https://arxiv.org/html/1109.2113)); it can be seen locally because an étale pullback of a generically smooth reduced divisor is generically reduced, and every source divisor maps dominantly to a target divisor.

All irreducible factors of every h_i∘F lie in {h_1,…,h_r}, by H∘F=H. Thus r nonempty, pairwise disjoint factor sets partition a set of r elements. Each set is a singleton. Squarefreeness makes its exponent one. The identity for H now gives the multiplicity equality. ∎

For completeness, the quasi-finiteness used here does not assume global properness. Over C, local invertibility implies that each fibre is discrete, and a discrete affine algebraic fibre is finite. Its algebraic form is the standard étale/quasi-finite implication. The proof over a general characteristic-zero field can be checked after embedding the finitely generated field of all relevant coefficients into C, and then descended. The same argument proves the lemma after algebraic closure even when H is not geometrically irreducible over the original field.

This finite-set argument is different from assuming every irreducible polynomial stays irreducible under an arbitrary Keller map. Irreducibility is proved here only for factors of the specified invariant H, by saturation of its full factor support.

## 3. The marked pair has a finite determinant-one stabilizer

**Theorem 2 (two-component finite-order gate).** Let m≥2, let F=(X,T) be a polynomial self-map of A² over a characteristic-zero field, with det DF=1, and suppose

    H∘F=H,   0≠H∈K[x,t],   t(1+x^m t) divides H.

Let r be the number of distinct geometric irreducible factors of H and let N be the order of the permutation from Lemma 1. Then

    F^N(x,t) = (b x,b^(−1)t),   b^(m−1)=1.

Consequently F^{N(m−1)}=id; in particular F is an automorphism, with inverse a positive iterate. For m=2 one has F^N=id and N divides r!.

**Proof.** Work over an algebraic closure. Both t and q=1+x^m t are irreducible and nonassociate. Apply Lemma 1 to G=F^N. It gives

    G_t=a t,   1+G_x^m G_t=d(1+x^m t),   a,d∈Kbar*.

Set t=0 in the second equation to obtain d=1. Then

    G_x^m=a^(−1)x^m.

Unique factorization gives G_x=b x for a nonzero constant b, and ab^m=1. The determinant-one condition on G gives ab=1. Hence b^(m−1)=1 and a=b^(−1). Taking the (m−1)st power finishes. Equality of polynomial maps descends to K. ∎

No properness theorem, rational inverse, generic-fibre genus, or planar Jacobian conjecture is used in this deduction. The resulting inverse is explicitly justified only after the finite-order identity is proved.

**Sharp boundary controls.**

- The exponent m=1 really is exceptional. With F=(2x,t/2), the polynomial H=xt(1+xt) is fixed and is divisible by t(1+xt), but F has infinite order and determinant one.
- With only t as compulsory component, the shear F=(x+t,t) has infinite order, determinant one, and fixes H=t.
- For m≥3, the root-of-unity allowance is real. Choose b of order m−1. Then F=(bx,b^(−1)t) fixes H=xt(1+x^m t), and has exactly that order.
- The full list of factors of H matters. A Keller map could permute the marked factors with other factors before an iterate returns them. One must not assume the two marked factors are individually fixed from the outset.

## 4. Universal single-carrier polynomial specialization is impossible

**Theorem 3.** Let K have characteristic zero, R∈K[p,y] be nonzero, and

    S=c+p²(p³−y²)R(p,y).

For the completed Hamiltonian operation E_λ=exp(λ{−,S}) supplied by the inherited completed Hamiltonian theorem, no λ≠0 has both E_λ(x) and E_λ(t) in K[x,t]. In particular no nonzero time of any nonconstant universal source-preserving Hamiltonian carrier is a polynomial map of the original source plane. The conclusion includes polynomial Keller maps without presupposing a polynomial inverse.

**Proof, divisor part.** Put q=1+x²t. Direct substitution gives

    p=tq,  y=xt²q,  D=t³q²,
    S−c=t⁵q⁴ R(tq,xt²q).                         (1)

The substitution is injective because its rational inverse is x=yp/D, t=D/p². Hence S−c≠0 and t q divides S−c. If E_λ(x),E_λ(t) were polynomials, they would define a polynomial F preserving S and having determinant one, since the completed Hamiltonian operation preserves the source Poisson bracket. Formal identities between these polynomials are polynomial identities: the source embeds injectively into its t-adic completion. Theorem 2, with m=2 and H=S−c, would make F finite order.

**Proof, infinite-order part in the original source completion.** By (1), write directly

    S−c=B(x)t^n+O(t^(n+1)),   n≥5, B(x)≠0.

On K[x][[t]], the derivation δ={−,S} raises t-valuation by at least n−1. Indeed δx=S_t has valuation n−1, while δt=−S_x has valuation at least n. In particular

    δx=nB(x)t^(n−1)+O(t^n),
    E_λ(x)−x=λnB(x)t^(n−1)+O(t^n).              (2)

The higher terms have valuation at least 2(n−1)≥n. The exponential converges in this source completion, obeys the additive group law, and has inverse E_(−λ). For every integer j>0, replacing λ by jλ in (2) yields a nonzero coefficient. Hence E_λ has infinite order directly on the original source completion.

The completion comparison has two precise, compatible formulations. The inherited theorem completes the **unlocalized** ring B₂=K[x,u,p,y]⊂K[x,t], with u=x²t. Here D^N B₂⊂t^(3N)K[x,t] because D=t³(1+x²t)² and every coefficient in B₂ is a source polynomial. Hence this particular inclusion extends continuously to the source t-adic completion and carries the defining partial sums to their source exponential sums. No such continuity is claimed after arbitrary localization of coefficients, and injectivity on the full D-adic completion is not needed.

Alternatively, compare only the actual series in a common ring. Set s=xt and τ=t; then S−c∈τK[s,τ], δs=τS_τ∈τK[s,τ], and δτ=−τS_s∈τ²K[s,τ]. The derivation raises τ-order at least one, so its exponential on x=s/τ belongs to τ^(−1)K[s][[τ]], and its exponential on t=τ belongs to K[s][[τ]]. The substitution

    σ:K[s][[τ]][τ^(−1)] → K[x]((t)),   s↦xt, τ↦t

is well-defined and injective: s^iτ^j maps to x^i t^(i+j), with distinct exponent pairs remaining distinct and only finitely many terms contributing to each t coefficient because j is bounded below. It intertwines the polynomial derivation powers term by term; their partial sums converge t-adically and match the directly defined source exponential. This is a faithful comparison of these series, not of entire generic-cusp or localized completions.

Although t need not belong to B₂, its induced image is E_λ(D)/E_λ(p)² by t=D/p²; multiplicativity of the source exponential identifies this ratio with its directly defined E_λ(t). The denominator is nonzero because E_λ(p)=t+O(t²). If E_λ(x),E_λ(t) are polynomials, substitution into K[x][[t]] is injective and their iterates coincide with the completed iterates. Hamiltonian bracket preservation makes this polynomial map determinant one, and δS=0 makes it preserve S. The finite order from Theorem 2 contradicts (2). ∎

The same logarithmic chart supplies an independent leading-jet control: p=s²+τ, y=s(s²+τ), x=s/τ, and {f,g}=τ(f_sg_τ−f_τg_s). If S−c=A(s)τ^k+O(τ^(k+1)), then k≥1 and δs=kA(s)τ^k+O(τ^(k+1)). Notice that n and k can differ: a monomial s^iτ^j becomes x^i t^(i+j). Neither valuation may be silently substituted for the other.

### An extension and a precise remaining route

Theorem 2 also excludes any infinite-order determinant-one polynomial realization of a completed operation preserving a common nonconstant H divisible by t(1+x²t), whether or not that operation is presented as one exponential. In particular, compositions exp(λ_i{−,f_i(I)}) with a common I=p^aD, a≥2, commute and reduce to the single Hamiltonian sum; if their combined Hamiltonian is nonconstant they are covered. Exact cancellation can of course give the identity.

The useful remaining cut-and-repair question is now stricter: a composition using **different invariants** must evade a common invariant with this compulsory divisor. This note does not prove that such compositions cannot become polynomial. Nor does it exclude arbitrary rational specializations of a general R: rational maps can introduce poles, so the polynomial factor-partition argument does not apply. The previous genus theorem remains stronger in its own special family because it excludes rationality there.

## 5. Exact computational controls

The companion standalone [script](../../04-computation/experiments/planar_carrier_factor_rigidity_20261019.py) imports no repository code; its [deterministic JSON](planar_carrier_factor_rigidity_20261019.json) retains the exact controls. It checks the source divisor identity, Poisson coordinate change, leading valuations for mixed and cusp-divisible carriers, finite factor-support saturation, the diagonal stabilizer equation, and all three sharp boundary controls. Its finite universe is a check of the stated algebra and combinatorial sublemma, not a search over all Keller maps or a proof of JC.

Run the script with ordinary Python and with `python -O`; no check uses Python's removable `assert`. Both runs agree byte for byte on **448 exact checks** over 18 carrier controls, all nonempty factor-support assignments for r≤4, and marked exponents 2≤m≤9. The deterministic JSON SHA256 is `ca335cdf65b8d5e7ae94f5cf124879f753ffd176b6cd24596a572502322cf02d`. The proofs above supply the infinite quantifiers.

## 6. Scope and publication status

This closes a genuine source-polynomial termination branch that the board left open, without auditing or rejecting the inherited finite-row work. It does not produce a Keller pair, establish chart entry, or settle the planar Jacobian conjecture. The statement is author-proved, computationally checked, and accepted by the [independent review](planar_carrier_factor_rigidity_20261019_independent_review.md), which reproduced both runs and added direct-source controls. External input is limited to standard étale geometry, with the squarefree characterization explicitly sourced above; the finite-support and marked-pair arguments are given in full.
