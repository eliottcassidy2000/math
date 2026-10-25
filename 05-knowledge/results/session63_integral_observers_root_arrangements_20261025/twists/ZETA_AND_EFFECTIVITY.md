# Fixed-twist zeta products, character aliasing, and line-complement effectivity

**Status: PROVED below.** These are applications of classical finite-group Fourier algebra, formal zeta series, and the Hodge–Deligne realization. No historical novelty claim is made. The carrier definitions and inherited surface statements are in `RESULT.md`.

## 1. Fixed-twist zeta is not freshly retwisted mass zeta

For a variety W over F_q, define its zeta series formally by

    Z(W/F_q,t)=exp(Σ_(m≥1) #W(F_(q^m)) t^m/m).

Let a nontrivial finite constant elementary abelian2-group G of order h act on W, let B=W/G, and let W_g/F_q be its fixed family of G-twists. Then

    ∏_(g∈G) Z(W_g/F_q,t)
       = Z(B/F_q,t)^h
          · (Z(W/F_(q²),t²)/Z(B/F_(q²),t²))^(h/2).        (Z1)

For the hourglass h=64, so the exponents are64 and32. This applies to W=V,B=U; W=X,B=P²; and W=Y,B=Bl_T(P²).

**Proof.** Theorem B of `RESULT.md` gives the averaged extension count #B for odd m and #W for even m. Therefore

    (1/h)Σ_g log Z(W_g/F_q,t)
      = log Z(B/F_q,t)
        + (1/2)(log Z(W/F_(q²),t²)−log Z(B/F_(q²),t²)).

Multiply by h and exponentiate in Q[[t]]. No convergence or analytic continuation is used. QED.

Thus the omitted factor is exactly the even-extension nontrivial-character information. Freshly choosing all twists over each F_(q^m) instead gives #B at every m and hence Z(B); that is a different operation.

For the good-characteristic open quotient,

    Z(U/F_q,t)=(1−qt)^6/((1−q²t)(1−t)^12).

This formula follows from its actual point polynomial, not from taking the geometric mean of fixed-twist zeta functions.

## 2. General character norm and extension aliasing

Let G be finite constant abelian, and work with a finite-dimensional characteristic-zero representation decomposed into character spaces H_χ over a splitting field. Let F commute with G. For χ of order d,

    ∏_(g∈G) det(1−χ(g)tF | H_χ)
      = det(1−t^d F^d | H_χ)^(|G|/d).                   (Z2)

Indeed χ(G) is the full group of dth roots of unity, each occurring |G|/d times; the identity ∏_(ζ^d=1)(1−ζT)=1−T^d applies to the commuting operator T=tF. Determinants prove (Z2). The same identity applies degreewise with alternating exponents in rational ℓ-adic compactly supported cohomology, ℓ≠char(F_q), after extending coefficients to split G. The Frobenius convention is the cohomological one corresponding to the chosen twist convention; replacing every χ by χ⁻¹ changes no product.

At extension degree m, the fixed-twist trace profile has the form

    N_g(m)=Σ_χ χ(g)^m A_χ(m).

Thus base extension aliases characters by χ↦χ^m. Its mean retains all χ with χ^m=1, not just χ=1. Equivalently it is the trace on G^m-invariants, as in `RESULT.md`. The averaging projector needs |G^m| invertible in its coefficient ring; an integral Z[G] splitting is not asserted.

For our marked seven-line cover the63 nontrivial global even-subset characters have weights2,4,6 in multiplicities21,35,7. An affine Walsh subset S⊂{1,…,6} corresponds to the even subset S when |S| is even and S∪{0} otherwise. Hence the seven weight-six characters are exactly the six affine weight-five modes and one affine weight-six mode. They are the seven canonical/K3 quotient characters inherited from `six_seven_gale_cube_20261025/CHARACTERS.md`. Their open-character variance block includes boundary contributions: it is not automatically a purely transcendental K3 Frobenius trace.

## 3. Classification of split-torifiable projective line complements

Let D be the union of n distinct projective lines over C, and U=P²\D. For each intersection P of at least two lines, let m_P be its multiplicity and put σ=Σ_P(m_P−1). Then

    E_c(U;u,v)=L²+(1−n)L+1−n+σ,
    where L=uv.

Equivalently, with T=L−1,

    E_c(U)=T²+(3−n)T+(3−2n+σ).                          (E1)

**Proof.** The sum of the line classes n(L+1) counts an intersection P precisely m_P times instead of once, so [D]=n(L+1)−σ in the Grothendieck ring. Subtract from [P²]=L²+L+1 and apply the Hodge–Deligne realization. Expansion gives (E1). QED.

Consequently U admits a finite stratification by split algebraic tori **if and only if**

    n≤2, or n=3 and the three lines are not concurrent.

**Proof.** A torus stratification gives a polynomial in T with nonnegative integer coefficients. For n≥4 the coefficient3−n is negative. For three concurrent lines σ=2, and E_c(U)=T²−1 is also impossible. For zero, one, or two lines the complements are P², A², and A¹×G_m, respectively, and have explicit torus stratifications. Three nonconcurrent lines become the coordinate triangle by a projective change of coordinates, whose complement is G_m². These exhaust the cases. QED.

Thus the hourglass coefficient−4 is literally3−7 and does NOT depend on the number of triple intersections. Triple incidence changes the constant term. If seven lines have only double and triple intersections, σ=21−t₃, so χ_c(U)=10−t₃. The hourglass t₃=3 gives7. Seven lines in general position instead give10, while the corresponding root-resolution geometric genus remains7 by the inherited adjunction theorem. The two occurrences of7 in the hourglass are therefore not a universal equality of invariants.

This classifies finite split-torus stratifications of these complex line complements only. It does not classify all F1 frameworks, every rational surface, or the geometry of each root-cover twist from its averaged count.
