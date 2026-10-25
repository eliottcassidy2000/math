# Twist mass, the retained branch algebra, and power-map specialization

**Status:** PROVED below, except the separately labeled finite checks and inherited inputs. Prepared 2026-10-25 from repository commit `eaa81fd8`; historical novelty is not claimed. This packet does not claim an F1 descent, an integer-square solution, or an RH implication.

## 0. Exact inherited carrier and inheritance pass

Closest proved mechanism: `05-knowledge/results/session59_f1_geometry_and_information_20261025/toric-shared-roots/TORIFICATION_OBSTRUCTION.md`, Theorem T1 and Corollary T2. They construct the simultaneous projective root cover, prove the hourglass has 24 A1 nodes and a resolution with p_g=7 and K²=16, and exclude finite split-torus stratifications on every nonempty complex open of that surface. The cited source is the theorem in this file, not a colliding bare theorem number.

Canonical hostile: the full compact cover, its branch-free open, and its resolution have DIFFERENT twist means. Corrected near miss: freshly retwisting over every extension is not the extension count of one fixed family of twists. Least-used sidecar: the power map on the deck group, together with the stabilizer/inertia at each branch stratum.

Board: root-cover torsor; arithmetic twist profiles; Fourier moments; branch strata; equivariant resolution; quotient versus ordinary motive; retained additive multiplication. Anchor: identify precisely what averages to q²−6q+12. Niche: resolution mass. Wildcard: power-map invariants under field extension.

Let k be a field of characteristic not two. On P² with coordinates [E:x:z], take seven forms

    L=(E,E+x,E−x,E+z,E−z,E+x+z,E−x−z).

Let X be the coordinate-square inverse image of this plane in P⁶, let D be the union of the seven lines, U=P²\D, and V=X×_(P²)U. The projective sign deck is G=(μ₂)^7/μ₂≅(C₂)^6. The cover V→U is a G-torsor of degree64: normalize the nonzero zeroth root and adjoin the six independent square roots L_i/E. The square classes are independent in k(P²), by their odd valuations along the six distinct numerator lines. In every odd characteristic all seven lines are distinct; the displayed ratios are units on U.

A twist is indexed by a six-vector of square classes δ_i, equivalently by seven square classes modulo common scaling. The index set is H¹(F_q,G)≅G for odd q. These are the64 **G-twists**; some may be isomorphic after forgetting the marking/action. Such isomorphisms do not remove their multiplicities from twist mass.

## 1. General finite-field twist mass, with the denominator proved

### Theorem A (Frobenius–orbit mass)

Let W be a quasi-projective variety over F_q with an action of a finite constant abelian group G. Write W_g for its twist whose Frobenius descent operator on geometric points is gF (changing the convention to g⁻¹ has no effect on these averages). Assume geometric finite-group quotients in the formulas are used. Then

    (1/|G|) Σ_(g∈G) #W_g(F_q) = #(W/G)(F_q).                 (A1)

Freeness is NOT required. If W→B is a G-torsor, the quotient is B and exactly one twist has a split fiber over each b∈B(F_q); that twist contributes |G| fiber points. Thus the torsor case also has a pointwise proof with no division ambiguity.

**Proof.** A rational point of W/G is one F-stable geometric G-orbit O. For each x∈O, the number of g satisfying gF(x)=x equals |Stab_G(x)|, since F(x) is in O. Therefore

    Σ_g #{x∈O : gF(x)=x}=|O|·|Stab_G(x)|=|G|.

Sum over quotient rational points and divide by |G|. These finite quotients exist for the quasi-projective cases used here; the argument itself requires only that geometric quotient fibers are the G-orbits. QED.

**Stack scope.** The weighted groupoid cardinality of [W/G](F_q) is the same sum in (A1): torsor classes are g∈G, each automorphism group is G, and the twisted-point groupoid contributes #W_g(F_q)/|G|. This agrees numerically with the coarse quotient count in this constant finite-group setting. It does NOT identify the stack with its coarse space when stabilizers are nontrivial. On V the action is free, so [V/G] really is U. The quotient-stack definition and torsor presentation are classical; see Stacks Project [Tag04UV](https://stacks.math.columbia.edu/tag/04UV). The free-action torsor statement is [Tag07S7](https://stacks.math.columbia.edu/tag/07S7). The numerical argument above is self-contained.

### Theorem B (fixed twists and the power-map correction)

Keep one fixed family {W_g/F_q : g∈G}. For every extension degree m≥1, put G^m={g^m:g∈G}, the image subgroup of the m-th power homomorphism. Then

    (1/|G|) Σ_g #W_g(F_(q^m)) = #(W/G^m)(F_(q^m)).           (B1)

In particular, if gcd(m,exp(G))=1 this is the full quotient count; if exp(G) divides m it is #W(F_(q^m)). For our exponent-two group:

    fixed-twist mean over F_(q^m)
      = #U(F_(q^m)) when m is odd,
      = #V(F_(q^m)) when m is even.                         (B2)

**Proof.** Constancy makes F commute with every g, so the m-th iterate of twisted descent is (gF)^m=g^mF^m. The power homomorphism G→G^m has fibers of equal size. Averaging therefore becomes the average of #Fix(hF^m) over h∈G^m, and Theorem A applied over F_(q^m) proves (B1). For a torsor, W/G^m is its intermediate cover of degree |G/G^m| over B. QED.

**Finite hostile controls.** The independently enumerated untwisted open counts over F9,F25,F49 are64,832,3136; the freshly-retwisted means are37,487,2119. Thus the fixed families descending from F3,F5,F7 have these respective untwisted counts as their quadratic-extension means, not the quotient polynomial. The finite source is `session63-recovery/twist_finite_normal.json`, fields `cases[].twist_open_counts[0]`.

The constant-group and abelian hypotheses are substantive: for a nonconstant group F need not commute with g, and in a nonabelian group the power map need not be a homomorphism with equal fibers. The formula is not silently extended to those settings.

**Equivariant realization.** For rational ℓ-adic compactly supported cohomology, using alternating traces and the corresponding cohomological Frobenius convention, the averaged trace is the trace of F^m on the G^m-invariant part: the averaging operator is the Reynolds idempotent |G^m|⁻¹Σ_(h∈G^m)h. On characters, the operation is χ(g)→χ(g^m), the usual power/Adams operation. This is an identity of equivariant traces, not an identity asserting that an ordinary average of geometric varieties is their quotient.

**Decisive motive hostile.** After algebraic closure all64 quadratic twists of V become isomorphic to V. A universal ordinary Grothendieck-class identity (1/64)Σ[V_δ]=[U] is false. For the fixed twists over F5, base change to F25 makes every twist isomorphic to V, but the independently verified point counts are832 for V and487 for U, contradicting such an identity already over finite fields. Separately, for the characteristic-zero version of the same root construction, all twists over C are isomorphic to V and compactly supported Euler characteristic gives χ_c(V)=64χ_c(U)=448 versus χ_c(U)=7. This is a separate characteristic-zero hostile, not a base change from F_q to C. The average is a fieldwise quotient mass, not the geometric average motive of64 isomorphic complex varieties. No compatible fixed64-family zeta function is determined by the polynomial q²−6q+12.

## 2. Branch incidence and the polynomial actually averaged

Away characteristics2 and3 the arrangement has t₃=3 triple points and t₂=12 double points. In characteristic3 it has t₃=5,t₂=6. One direct derivation is to decone at E=0. The three pairs x=±1,z=±1,x+z=±1 meet in three directions at infinity. A finite triple can occur only at x,z∈{±1}, x+z∈{±1}. In odd characteristic this happens exactly in characteristic3, at (1,1) and (−1,−1). No four lines meet. Pair counting then gives t₂+3t₃=21.

For every field in either odd-characteristic regime, all intersections are rational. Inclusion-exclusion gives

    [U] = L² − 6L + (t₂+2t₃−6)
        = L² − 6L + 15−t₃,                              (C1)

where L=[A¹] in the ordinary Grothendieck ring of varieties. This is an ACTUAL class identity for U, independent of the twist-average argument. In particular

    #U(F_q)=q²−6q+12       if char(F_q)≠2,3,
    #U(F_q)=q²−6q+10       if char(F_q)=3.                (C2)

Combining (A1) with (C2) proves the stated freshly-retwisted64 mean of the branch-free covers V_δ. Since the complement includes E=0, this is the affine (x,z)-plane with its six lines removed, not a count of all projective root solutions.

The exact strata are:

| Base stratum | Base count | Nonzero root coordinates | Compatible twists out of64 | Points of one compatible fiber |
|---|---:|---:|---:|---:|
| U | q²−6q+15−t₃ | 7 | 1 | 64 |
| smooth line pieces | 7(q+1)−2t₂−3t₃ | 6 | 2 | 32 |
| double points | t₂ | 5 | 4 | 16 |
| triple points | t₃ | 4 | 8 | 8 |

All fibers here mean reduced geometric/finitary point fibers; the scheme fibers over branch loci are nonreduced. At a point on r lines the inertia subgroup has order2^r. Compatibility fixes the square classes of the 7−r nonzero roots up to common scaling, leaving2^r twists. Each compatible fiber has2^(6−r) rational points. Their average is1 on every stratum. Therefore

    fresh mean #X_δ(F_q) = q²+q+1,                       (C3)

where X_δ denotes the whole singular projective surface. No free-torsor argument is applied incorrectly to its branch locus.

## 3. Resolution changes the quotient, not just the source

Let Y→X be the minimal resolution, obtained by blowing up its reduced A1 node locus. This construction is G-equivariant and twists compatibly. Let T⊂P² be the reduced scheme of triple points.

### Theorem D (resolved quotient)

For this arrangement over every field of odd characteristic,

    Y/G ≅ Bl_T(P²).                                      (D1)

Consequently, over F_q,

    fresh mean #Y_δ(F_q)
       = q²+(t₃+1)q+1
       = q²+4q+1              in characteristic≠2,3,
       = q²+6q+1              in characteristic3.        (D2)

**Proof of the quotient.** Away T the original X is smooth and its quotient is P². Étale-locally at a triple point the relevant cover is the node

    w²=a u²+b v²,       a,b≠0,

mapping to A² with coordinates x=u²,y=v². The remaining roots are units and account for the étale branches; the inertia group is the independent sign group on u,v,w. On the u-chart of the node blowup set V=v/u,W=w/u. Its coordinate ring is

    k[u,V,W]/(W²−a−bV²).

The induced independent signs on u,V,W have invariant ring k[u²,V²]. The map is x=u²,t=V², exactly the x-chart (x,t=y/x) of Bl_0(A²). It is finite of degree8. The v-chart has the analogous description and the charts glue by the usual blowup transition. They cover: on the exceptional smooth conic, u=v=0 would force w=0; similarly they cover off the node. Quotients commute with these étale local identifications because the group order is invertible, and the descriptions agree with P² away T. This proves (D1).

For (D2), apply Theorem A to Y. Blowing up one rational smooth surface point replaces it by P¹ and adds q points. All t₃ triple points are rational here. More generally the correction is q·#T(F_q), not q times the number of geometric points when some triples are nonrational. QED.

**Independent count proof.** Over a rational triple point the average number of rational nodes among X_δ is1, by the last row of the strata table. A rational ordinary double point of each twist is replaced by a smooth projective conic. Over a finite field every smooth conic has q+1 points, split or not; so resolution adds q per rational node. A nonrational node creates no rational points above it. This reproduces (D2) without identifying the quotient. It would be incorrect over an arbitrary field to call every such exceptional conic a split P¹ without justification.

The three good-characteristic formal values at q=1 are therefore:

    branch-free twist mean: 7;
    singular projective twist mean: 3;
    resolved projective twist mean: 6.

These are Euler characteristics of U, P², and Bl_3(P²), respectively. They are not the Euler characteristics of V, X, Y, which are448,56,80 over C. Indeed the branch strata give χ_c(X)=64·7+32·(−19)+16·12+8·3=56; resolution adds24, agreeing with Noether's 12(1+7)−16=80. Nor is the shared number7 a proof that χ_c(U) and p_g(Y) are the same invariant.

## 4. Full Fourier sidecar, moments, and the first arithmetic information beyond mass

Work over F_q of odd characteristic and put f_i=L_i/E on U, for i=1,…,6. Let χ be the quadratic character extended by χ(0)=0; on U every f_i is nonzero. Write

    s(u)=(χ(f_1(u)),…,χ(f_6(u))) ∈ {±1}⁶,
    n_a=#{u∈U(F_q):s(u)=a}.

A twist sign ε∈{±1}⁶ has

    N_ε=#V_ε(F_q)
       =Σ_(u∈U) ∏_(i=1)^6 (1+ε_iχ(f_i(u)))
       =64 n_ε.                                         (E1)

(The inverse twist convention gives the same ε because signs have order2.) For a subset S⊂{1,…,6}, define

    A_S=Σ_(u∈U(F_q)) χ(∏_(i∈S) f_i(u)),   A_∅=|U|.

Then the Walsh transform gives exact inversion

    N_ε=Σ_S ε^S A_S,
    A_S=(1/64)Σ_ε ε^S N_ε.                               (E2)

For uniform ε, all character moments follow. In particular

    E[N_ε]=|U|,
    Var(N_ε)=Σ_(S≠∅) A_S²
            =64Σ_a n_a²−|U|²,                           (E3)
    Cov(N_ε,N_(εh))=Σ_(S≠∅) h^S A_S².                  (E4)

**Proof.** Each factor in (E1) is2 if ε_i=s_i(u) and0 otherwise. Expanding the product proves the first identity in (E2); multiplying by ε^T and averaging uses orthogonality to prove the second. Orthogonality again proves (E3)–(E4); substituting N_ε=64n_ε gives the occupancy form. QED.

Thus the mean destroys all63 nontrivial character sums. The complete labeled twist profile recovers them. All shift covariances recover their squares but not their signs: this is a phase-retrieval boundary, not a full reconstruction. It is an exact connection to the repository's Walsh/local-coefficient methods, with a different underlying six-bit carrier than the three-bit cube; no identification based merely on the number7 is asserted.

If M=|U| and s=#supp(n), Cauchy–Schwarz gives the useful support bound

    s ≥ 64 M² / (M²+Var(N)),

and N_ε≤64M trivially. The second moment measures collision of square-class profiles, not a rational-point solution for any specified twist. A positive mean proves some twist has points, not that the untwisted Diophantine surface does.

## 5. Torification obstruction: exact strength and limits

The actual class (C1), after complex Hodge–Deligne realization, is

    E_c(U;u,v)=(uv)²−6uv+12
              =(uv−1)²−4(uv−1)+7.

A finite stratification by split algebraic tori would give a polynomial in uv−1 with nonnegative integer coefficients. Uniqueness of polynomial coefficients rules this out. This is an obstruction for the actual arrangement complement U, not a blanket no-F1 theorem. U itself is rational, so rationality is only a necessary condition for finite torification, not sufficient.

There is also a precise finite-field consequence for the average. Suppose ALL64 fixed twists over one good odd F_q admitted finite split-torus stratifications. Their count polynomials P_ε(T) lie in N[T−1]. For every ODD extension degree m, Theorem B gives

    (1/64)Σ_ε P_ε(q^m)=q^(2m)−6q^m+12.

Equality on infinitely many values forces polynomial equality. But the left has nonnegative rational coefficients in T−1 while the right has coefficient−4. Contradiction. This argument excludes all64 simultaneously; the mean ALONE does not prove each particular twist non-torifiable.

For the individual complex root opens, the inherited Session59 Corollary T2 supplies the stronger, independent proof: X is irreducible, its smooth projective resolution has p_g=7>0, and every nonempty open has the same nonrational function field. A finite torus stratification of an irreducible surface has a dense two-dimensional torus stratum, forcing rationality. This uses geometry that averaging discarded.

The resolved quotient's positive polynomial q²+4q+1=(q−1)²+6(q−1)+6 is a hostile control: positive coefficients of a quotient count are not a proof that its source twists admit torifications. In this case the three blowup centers are collinear on E=0. Its complement is A². Keep the strict transform of that line as a whole P¹; each of the three exceptional curves, minus its intersection with that line, is A¹. This gives the actual stratification A² ⊔ P¹ ⊔ 3A¹, and refining each piece gives split tori. Nevertheless Y remains nonrational over C. Positivity alone is not sufficient in general.

## 6. Constructive retained algebra: multiplicative roots plus additive relations

The precise positive carrier is stronger than just its bit group. In a chart E≠0 the root algebra is

    O_U[r_1,…,r_6]/(r_i²−f_i).

It has a square-free monomial basis indexed by subsets S⊂{1,…,6}, and multiplication is

    r_S r_T = r_(S△T) ∏_(i∈S∩T) f_i.                   (F1)

The group labels retain S△T; the factor ∏f_i retains the additive branch geometry. Dropping it keeps the character support but changes the algebra, cover, and arithmetic.

Globally the inherited `six_seven_gale_cube_20261025/CHARACTERS.md`, Theorem1, gives

    π_*O_X=⊕_(S⊂{0,…,6}, |S| even) O_P²(−|S|/2),

with the same multiplication expressed in normalized charts. Finite-flat duality yields canonical summands O(|S|/2−3). Exactly the seven six-element subsets contribute one two-form each. These are complementary six-element characters, NOT singleton characters; adjunction carries the determinant twist. This exact algebra-to-character-to-canonical-form map explains p_g=7 independently of the twist mean's Euler7.

The earlier Session59 `toric-shared-roots/RESULT.md`, §5 gives an honest additive blueprint pushout: the free pointed multiplicative root monoid is equipped with preaddition recording equations such as a²+b²≡c²+c² and h²+c²≡a²+d². Ringification recovers the actual shared-root tensor-product algebra. Its nonzero {0,1} homogeneous-coordinate assignment count is1, not7; it is a different specified point functor. This is a constructive equation carrier, not a claim that its existence proves flat descent, a preferred F1 count, Frobenius lifts, or RH.

The operation contract is therefore:

    root equations → equivariant finite algebra → branch-free torsor
      → full twist/Fourier profile → quotient mass;

retain the multiplication factors, inertia strata, and extension power map if the geometry or arithmetic must be recovered. The quotient mean alone forgets exactly the data needed to distinguish twists and singular/resolved carriers.

## 7. Reproduction and scope

Independent exact finite-field twist histograms and all64 shift covariance checks are in sibling `session63-recovery/check_twist_histograms.py`, with byte-identical PASS outputs `twist_finite_normal.json` and `twist_finite_optimized.json` for q=3,5,7,9,11,13,17,25,49. These are controls, not substitutes for the proofs. This directory's `verify.py` independently checks compact/root/node/resolution counts for p=3,5,7,11,13, every inertia fiber pattern, all diagonal exceptional conics over those fields, and4,080 power-map orbit cases across six finite abelian groups. Run `python3 -B verify.py > normal.json` and `python3 -B -O verify.py > optimized.json`; the saved outputs are byte-identical and PASS. The separate `ROOT_STACK.md` retains branch inertia constructively; `ZETA_AND_EFFECTIVITY.md` proves the fixed-twist zeta correction and the general line-complement classification.

Characteristic2 is excluded: square roots are purely inseparable, ± roots coalesce, the displayed G is not an étale constant64-element group scheme, and line pairs coincide. The64 quadratic-twist argument does not continue there. In characteristic3 the seven lines remain distinct and the A1 analysis remains valid, but incidence and every stated U/resolution count must use t₃=5.
