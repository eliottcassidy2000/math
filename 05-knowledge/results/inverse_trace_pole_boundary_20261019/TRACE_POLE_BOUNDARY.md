# Inverse trace poles, surviving sheets, and an alternating-monodromy gate

**Current global scope, 2026-10-20 (CITED).** [Domrina's 2000 theorem](https://www.mathnet.ru/eng/im273) excludes complex polynomial Keller maps of generic/function-field degree four. The repository's [primary-PDF source record](../../reference/CORE-PAPERS-PLANAR-JC.md) confirms that degree convention. Quartic A4 is therefore only a residual of this local method, not a globally open Keller case. The local pole models and conditional trace gates remain valid; the live global program must address higher-degree families, without claiming a particular first open degree.

19 October 2026. **PROVED, scoped synthesis; independent review PASS.** The arguments below use elementary Newton identities, the standard finite-normalization factorization of an étale polynomial map, and the already identified finite-sheet survival mechanism. No historical-priority claim is made. The planar Jacobian conjecture remains open. In particular, no regularity of the trace moments below is asserted for an arbitrary Keller map.

## Main findings

Let F=(P,Q): A²_C → A²_C be polynomial with det DF=κ≠0. Put A=C[u,v], K=C(u,v), L=C(x,y), u=P, v=Q, and d=[L:K].

1. **Exact boundary detector.** Choose generic affine source coordinates (x,y) in which P is monic in y and x is primitive for L/K. If d≥2, the reduced nonproper-value set of F is exactly the union of the pole divisors of

       T_k=Tr_(L/K)(x^k),       1≤k≤d−1.

   Consequently the squarefree least common multiple of the reduced denominators of these d−1 rational functions is an equation of the nonproper-value set. All these moments are polynomial if and only if F is an automorphism. The properness criterion is classical in content; the finite-sheet argument saves the last moment from the usual d-moment Newton test.

2. **Quantitative version.** For any source polynomial h and target prime divisor D, let β_1,...,β_d be its values on the geometric inverse sheets over a splitting field of the completed transverse DVR, with valuation normalized by v_D(t)=1. If m of these values have negative valuation, then m≤d−1 and

       max_i max(0,−v_D(β_i))
       = max_(1≤k≤m) max(0,−v_D(Tr(h^k))/k).

   With m=0 both sides are zero. One may replace m by d−1. Thus cancellation in low trace moments cannot hide an inverse pole indefinitely once actual global sheet degree is fixed.

3. **Extremal blindness is rigid.** If there are m escaping h-values and Tr(h),...,Tr(h^(m−1)) are regular at D, all m escaping values have one common pole order λ. Their nonzero leading coefficients are exactly a scalar multiple of all m-th roots of unity. The first polar moment is Tr(h^m), of pole order mλ, an integer.

4. **Even-degree alternating gate.** In the primitive monic coordinate above, suppose d≥2 is even and T_1,...,T_(d−2) are polynomial. Then the geometric monodromy lies in A_d. More precisely, every nonproper divisor has exactly one finite sheet, T_(d−1) has a pole of some integer order a>0, and the primitive inverse discriminant has pole order exactly da there. Its other affine divisorial orders are even. Thus the discriminant is a square in C(u,v).

   For d=4 this says that S4 monodromy forces a pole in Tr(x) or Tr(x²). The earlier fixed-sheet/normal-generation gate excludes the regular V4 action; combined with it, the only transitive quartic possibility when both moments are polynomial is A4. A4 is not excluded by this method alone; it is globally excluded by Domrina's generic-degree-four theorem cited above.

These are global inverse-field statements, not assertions about a point Taylor jet or a formal Hamiltonian logarithm.

## 1. Inheritance, closest mechanism, and the new object

The preceding packets establish that polynomial symplectic automorphisms can match every finite jet of rational cotangent germs, and that formal, analytic, rational, and polynomial integrability are different categories. The period-three residue identity gives a genuine global sum law, but a residue retains only one coefficient of a principal part. The present change of object is:

       genuine polynomial Keller map
         → finite normalization and all inverse sheets
         → field traces of powers of a source coordinate
         → affine divisorial pole data.

The closest proved inverse-sheet mechanism is THM-2633, `affine-inertia-normal-generation...` (repository filename `THM-2633-derangement-character-obstruction-and-d4-keller-exclusion.md`): every affine target divisor has an actual finite sheet. THM-2621, `planar-degree-four-inverse-spectral...`, retains the inverse quartic and its second-coordinate companion, and proves exactness of the trace–Liouville form. The present detector uses **pole order**, which that exactness and its zero residues do not eliminate.

The relevant underused geometric sidecar is THM-3951, `affine-plane-boundary-incidence-forest`: the actual upstairs boundary has rational components and a forest incidence constraint. Our traces forget that incidence and the individual sheet labels. They do not replace that sidecar. They provide a finite algebraic way to identify which downstairs divisors require it.

The metabelian principal-part criterion is an inspiration to retain whole polar data, not a theorem being transferred to Keller maps. Its maps are birational rational germs; cotangent lifts of nonlinear scalar polynomials have multiple sheets and affine poles; neither class is a class of polynomial Keller maps.

### Concept board

- Source regularity: genuine A² source versus a punctured rational source.
- Surviving sheet: a local bounded inverse branch at every affine target divisor.
- Global moments: Tr(h^k), with degree d fixed before truncation.
- Principal parts: negative valuations, not only residues.
- Extremal cancellation: a root-of-unity pole packet.
- Monodromy: inertia permutation and discriminant parity.

The wildcard connection to tournament spectra is deliberately limited: both use power sums and Newton identities, but adjacency traces lose support-overlap data (THM-500), while our target predicate is only the existence/order of a coordinate pole. No tournament orientation is naturally induced on inverse sheets, so none is imposed. Likewise a Petersen incidence graph cannot be transplanted into an affine-plane boundary: THM-3951 forces the genuine upstairs incidence multigraph to be a forest. A diagram or a spectral coincidence is not a realization.

## 2. The finite normalization and survival of an affine sheet

Let X be the normalization of Spec A in L. It is finite over Spec A, since A is a finitely generated normal C-algebra. Its coordinate ring is contained in C[x,y]: an element integral over A is integral over C[x,y], and C[x,y] is integrally closed. Zariski's Main Theorem gives an open immersion

       j:A²_source → X,       π:X → A²_target finite,       F=π∘j.

The specialization of the standard theorem used here is cited in Section 10; it is not a formal-series assertion.

**Lemma 2.1 (finite sheet).** Every irreducible affine target divisor D meets the open image of F in a nonempty Zariski-open subset. In particular there is at least one bounded geometric inverse sheet at the generic point of D.

**Proof.** F is étale, hence open. If D={g=0} were disjoint from its image, g(P,Q) would be a polynomial without a complex zero, so it would be a nonzero constant. Dominance makes A→C[x,y] injective and rules this out for nonconstant g. Thus the image meets D; openness makes the intersection dense in D. At a source point over the generic point, étaleness makes the corresponding local inverse unramified and bounded. After passing to the strict henselization of the generic transverse DVR, it is an actual branch with both source coordinates integral. ∎

This is the same mechanism as the finite-sheet lemma in THM-2633, with a proof included here. It uses both source polynomial regularity and étaleness. The rational map q↦q^(−d) on a punctured line has no such finite sheet over zero.

**Lemma 2.2 (boundary is precisely coordinate poles).** The complement X\j(A²) is the union of the prime divisors on which x or y has a pole. Its finite image is the nonproper-value set A_F.

**Proof.** If both x and y belong to O_(X,z), they define a morphism to A² on a neighborhood of z. Composing with j is the identity there, since this holds on the dense source open and X is separated. Thus z belongs to the source open. Conversely x and y are regular on that open. A normal local domain is the intersection of its height-one valuation rings; consequently failure of either coordinate to be regular is witnessed by a prime polar divisor through the point. This proves the first claim, including absence of a hidden codimension-two boundary.

Off the finite image of that boundary, F is the full finite map π. Conversely, if F were proper over a neighborhood meeting π(X\j(A²)), then j on that neighborhood would be proper: a map from a proper scheme over the base to a separated scheme over the base is proper. Its dense open image would be closed, hence all of π^(-1) of that neighborhood, contradicting the omitted boundary point. ∎

Thus a nonproper divisor has between one and d−1 finite geometric sheets. Missing sheets are not counted as finite merely because one arbitrarily selected rational observable remains bounded.

## 3. A local moment lemma with exact pole order

Work over a characteristic-zero discretely valued field K_D whose residue characteristic is zero. Fix an extension of its valuation to an algebraic closure, normalized by v(t)=1. Let β_1,...,β_d be a Galois-stable multiset, and suppose at most m entries have negative valuation. Define s_k=Σ_i β_i^k, assumed to lie in K_D.

**Lemma 3.1.** If λ=max_i max(0,−v(β_i))>0 and r entries attain valuation −λ, then

       λ=max_(1≤k≤r) max(0,−v(s_k)/k).

In particular r≤m and the same maximum can be taken through m.

**Proof.** Every term of s_k has valuation at least −kλ, so v(s_k)≥−kλ. After a finite extension allowing a uniformizer power of valuation λ, write the r most singular terms as

       β_i=t^(−λ)(α_i+terms of positive valuation),       α_i≠0.

This notation denotes the corresponding Puiseux scaling, not a claim that t^λ lies in the original field. All other scaled roots have residue zero. If all power sums Σ_(i=1)^r α_i^k vanished for 1≤k≤r, Newton identities in residue characteristic zero would make all elementary symmetric functions vanish. Their product would be zero, contrary to every α_i≠0. For some k≤r that leading power sum is nonzero. Then v(s_k)=−kλ. ∎

**Corollary 3.2.** If m≥1 entries escape and s_1,...,s_(m−1) are integral, then all m escaping entries have valuation −λ, their leading residues are the m distinct roots of Z^m−c for some c≠0, and v(s_m)=−mλ.

**Proof.** If fewer than m entries attained −λ, Lemma 3.1 would detect a pole before m. Thus all m attain it. Integrality of the first m−1 moments kills their leading power sums. Newton identities give the leading polynomial Z^m−c, with c≠0. It is separable in characteristic zero. Its m-th power sum is mc, proving the last assertion. Since s_m is in the base field, mλ is an integer. ∎

This is sharp at every m: the roots of Z^m−t^(−a), together with one zero root, have s_k=0 for 1≤k<m and s_m=m t^(−a). It is a local finite étale algebra over the punctured trait with one retained finite sheet. It is not, by itself, a polynomial Keller realization.

The vector version is also exact. If the leading pole vectors a_1,...,a_m in a residue vector space satisfy Σ_i ℓ(a_i)^k=0 for every linear form ℓ and every k<m, they are all on one line and form a root-of-unity packet there. Indeed for generic ℓ the ratios ℓ(a_i)/ℓ(a_j) belong to the finite set μ_m. A rational function of ℓ with finite image is constant, so a_i=ζ_i a_j, and the scalar result gives all m roots. This describes the strongest possible multivariate moment cancellation; ordinary collinearity alone does not make the moments vanish.

## 4. Exact trace-pole detection for Keller maps

For a source polynomial h, take β_i to be its d conjugate values over K. At every target prime D at least one actual finite sheet has bounded h by Lemma 2.1. Hence at most d−1 of the β_i escape. Lemma 3.1 gives

       λ_D(h)=max_(1≤k≤d−1) max(0,−ord_D Tr(h^k)/k).      (4.1)

This expression uses global field trace, with conjugates counted with multiplicity even if h is not primitive.

For two source affine coordinates define

       Λ_D(F)=max(λ_D(x),λ_D(y)).                         (4.2)

Then Σ_D Λ_D(F)D is a finite effective Q-divisor supported exactly on A_F. It is invariant under invertible affine changes of source coordinates: the maximum negative valuation of a coordinate vector is unchanged by an invertible constant matrix and by a constant translation. It is not invariant under arbitrary polynomial source changes; a high-degree shear can enlarge a pole. The support A_F is intrinsic under every source automorphism. Target coordinate changes transport the prime divisors and their normalized valuations.

This is a genuine inverse-sheet invariant. It forgets residue directions, exact inertia, boundary intersections, and the companion inverse coordinate. Those are separate sidecars.

## 5. One transverse primitive coordinate suffices in the plane

A generic affine source frame can be chosen so that P is monic in y and x is primitive for L/K. Monicity excludes the finitely many directions on which the top homogeneous form of P vanishes; after a constant rescaling, its leading y coefficient is one. A linear form fails to be primitive only when two of the d generic source points have equal value on it; these are finitely many proper linear conditions on its direction. Both conditions can be satisfied at once. These are the same generic-frame facts used in THM-2621.

If x is integral at a valuation lying above a finite target divisor, the equation

       P(x,y)−u=0

is monic in y with integral coefficients, so y is integral as well. Therefore x detects every escaping source sheet in this gauge.

Let χ_x(Z)=Z^d+a_1 Z^(d−1)+...+a_d be its minimal polynomial. Newton identities give

       k a_k + Σ_(j=1)^k a_(k−j) T_j=0,       a_0=1,     (5.1)

for 1≤k≤d. If T_1,...,T_(d−1) are regular at D, so are a_1,...,a_(d−1). By Lemma 2.1, χ_x has a bounded root α in a valued algebraic extension. The identity

       a_d=−α^d−a_1 α^(d−1)−...−a_(d−1)α

then makes a_d regular too. All roots of the resulting monic integral polynomial are integral. Monicity of P in y bounds the companion coordinate, so there is no boundary above D.

Conversely, where F is finite, the trace of every source polynomial is regular. Equivalently, all conjugate coordinate values are integral and their power sums lie in the base DVR. This proves

       D⊂A_F  iff  some T_k, 1≤k≤d−1, has a pole at D.   (5.2)

Since A=C[u,v] is a UFD, (5.2) proves the denominator compiler in the main statement. If every T_k is polynomial, F is finite by Lemma 2.2. A finite complex Keller map is a proper local biholomorphism, hence a covering of simply connected C². Its connected source makes the degree one. Finite birationality onto the normal target gives a polynomial inverse.

The d=1 case is included separately: Lemma 2.1 leaves no room for an escaping sheet at any divisor. Lemma 2.2 again makes F finite and then an automorphism. The empty moment list in degree one is therefore harmless, but it must not be confused with an unknown d for a map supplied only by a germ.

## 6. The even-degree alternating-monodromy gate

Assume d≥2 and T_1,...,T_(d−2) have no affine poles. At a nonproper divisor D, at least one sheet is finite. If there were at most d−2 escaping sheets, Lemma 3.1 would detect a pole in the assumed regular list. Thus there are exactly d−1 escaping sheets and one finite sheet.

Put m=d−1. Corollary 3.2 says that the escaping x-values have pole order

       λ=a/m,       a=−ord_D T_m ∈ Z_(>0),

and their leading coefficients form m distinct roots of one nonzero number. Differences of two escaping values therefore have valuation −λ. A difference of an escaping and the finite value also has valuation −λ. Consequently

       ord_D Disc(χ_x)
         =2[binom(m,2)+m](−λ)
         =−m(m+1)λ
         =−da.                                           (6.1)

At a divisor outside A_F, the finite map is étale, so after strict henselization every sheet is an unramified bounded branch. Differences of conjugate x-values have integer valuation, and the discriminant, the product of their squares, has even valuation. Coincident x-projections can create such even zeros; they are not branch ramification.

If d is even, (6.1) is even at every nonproper divisor too. Every irreducible polynomial in A therefore has even exponent in the numerator or denominator of Disc(χ_x). Constants in C are squares, so

       Disc(χ_x) ∈ K*².

The standard Vandermonde argument identifies its square class with the sign character of the permutation action on the roots: an automorphism acts on ∏_(i<j)(β_i−β_j) by its sign. Thus the geometric monodromy lies in A_d. This proves the gate without using a separate purity argument.

There is a sharper local permutation statement. Transverse analytic continuation t↦e^(2πi)t multiplies the leading escaping coefficients by e^(−2πi a/m). Distinct leading coefficients identify the sheet permutation uniquely. If g=gcd(a,m), inertia has cycle type

       1 · (m/g)^g.                                     (6.2)

It fixes the finite sheet. Equation (6.2) is a local statement after removing codimension-two special points, where the transverse Puiseux branches are defined; no global ordering of roots is chosen.

### Quartic consequence

For d=4, a boundary invisible to T_1,T_2 has exactly three escaping sheets and one finite sheet. Its inertia is the identity if 3 divides a, and a 3-cycle otherwise. In particular it is never a transposition.

Therefore an S4 Keller map of field degree four must have a polar T_1 or T_2. The equivalent discriminant proof above shows this directly. The earlier THM-2633 removes the other even transitive quartic group V4, because its nonidentity permutations have no fixed sheet and affine inertia normally generates monodromy. With that inherited result, polynomial T_1,T_2 leaves A4 as a residual of this method, while Domrina globally excludes generic degree four.

A4 is an actual boundary of this method: a retained root plus a cubic root-of-unity pole packet has precisely the allowable 3-cycle inertia and its first two moments vanish. This local theorem alone neither produces nor excludes a global A4 Keller map; Domrina's theorem already excludes such a map. The model measures this method's missing geometric input, not a live global quartic possibility.

## 7. Polynomial and rational hostiles

### 7.1 Polynomial source and a surviving sheet do not justify shortening d−1

For every d≥2 consider the everywhere polynomial map

       G_d(q,p)=(u,v)=(p, p q^d−q).

Its inverse coordinate q satisfies

       χ(Z)=Z^d−u^(−1)Z−v u^(−1).

The field degree is d: over C(u), the polynomial map q↦u q^d−q has degree d, so [C(u,q):C(u,v)]=d. At u=0 it has the single finite inverse q=−v, while d−1 branches escape. Its trace moments satisfy

       Tr(q^k)=0                  (1≤k≤d−2),
       Tr(q^(d−1))=(d−1)/u.

Thus the d−1 threshold is sharp among polynomial generically finite maps with a surviving finite sheet. The map is not Keller:

       det DG_d=1−d p q^(d−1).

The finite sheet over u=0 is itself unramified (the determinant there is one), so merely checking that one local branch is regular does not repair the missing global étaleness hypothesis. The critical locus elsewhere is the exact failure. For even d the discriminant has additional odd-order critical-value zeros, so the alternating gate cannot be applied to this hostile.

### 7.2 Exact symplecticity does not repair a rational source pole

Let f_d(q)=q/(1+q^d) and use its cotangent lift

       C_d(q,p)=(Q,P)=(f_d(q), p/f_d'(q)),
       f_d'(q)=[1−(d−1)q^d]/(1+q^d)^2.

It has rational Jacobian one wherever defined. The inverse q-polynomial is

       Z^d−Q^(−1)Z+1,

so again the first d−2 traces vanish and Tr(q^(d−1))=(d−1)/Q. There is a bounded branch q~Q at Q=0 and d−1 escaping branches. The map has genuine affine poles in its displayed coordinates and is not a polynomial Keller map. Its critical/deleted divisors and their monodromy cannot be discarded.

This sharpens the earlier polynomial-scalar cotangent hostiles: here the local pole packet includes an actual finite sheet. It still does not provide a Keller counterexample.

### 7.3 Failure of finite-sheet survival costs one full moment

For the local inverse polynomial Z^d−t^(−1), every sheet escapes. The first d−1 trace moments vanish and the d-th is d/t. The generic Newton d bound cannot be reduced without a bounded sheet. A rational symplectic realization is the cotangent lift of f(q)=q^(−d) on its natural punctured domain; it does not have an affine-plane polynomial source map.

### 7.4 Properness is not invertibility without étaleness

The finite polynomial map (q,p)↦(q^d,p) has all coordinate trace moments polynomial but generic degree d. Its determinant d q^(d−1) vanishes. The last step from finite to invertible depends on the Keller condition.

### 7.5 Residue-zero is weaker than pole-free

The exact rational differential d(t^(−a)) has a pole and zero residue. Hence the exactness/zero-residue theorem for trace–Liouville forms does not imply regularity of any trace potential. The present test keeps the polar coefficients that a residue sum discards.

### 7.6 The remaining triple packet passes the local discrepancy test

A local order-of-vanishing argument does not remove the cubic pole packet. In coordinates (z,w), put

       x=z^(−1),       y=z^4 w,
       t=z^3,          s=−w/3.

Then

       dx∧dy=−z² dz∧dw=dt∧ds.

The boundary has ramification index three and different/canonical order two, exactly as required. The source valuation is genuinely realizable by boundary blowups: start on P¹_x×P¹_y with z=1/x and blow up the point z=y=0 four times in the successive charts y=z w_1, y=z² w_2, y=z³ w_3, y=z⁴ w. Every center is outside the original A² and its boundary incidence remains a tree.

The resulting target formula is rational, t=x^(−3), s=−x^4 y/3. It is a local three-sheet model, with no retained finite sheet. The f_4 cotangent example in Section 7.2 supplies a retained finite sheet together with the same cubic escaping behavior, but still has affine poles. Neither example realizes a global A4 Keller map. They establish a specific stopping rule: local volume preservation, ramification order, and a realizable individual boundary valuation do not by themselves exclude the remaining packet. Global polynomial extension and interactions among boundary components must do additional work.

## 8. Differential escalation of a detected pole

The Keller dual derivations are polynomial:

       δ_u=(Q_y ∂_x−Q_x ∂_y)/κ,
       δ_v=(−P_y ∂_x+P_x ∂_y)/κ.

They extend ∂_u,∂_v on K, and trace commutes with them. This follows from uniqueness of derivation extension to a finite separable field extension, or by differentiating the sum over conjugate embeddings.

Suppose T=Tr(H) has a pole of order a>0 along the reduced divisor h(u,v)=0. Choose a target coordinate w∈{u,v} such that ∂_w h is not divisible by h; one exists in characteristic zero. Differentiating a rational function with an exact h-pole increases its pole order by exactly one, so

       ord_D Tr(δ_w^n H)=−a−n       for every n≥0.       (8.1)

Every δ_w^n H is still a polynomial on the affine source. Therefore one detected low-degree moment gives an infinite, explicit polynomial family with increasing trace poles. If the coefficients of δ_w have degree s, then

       deg(δ_w^n H)≤deg H+n(s−1)

whenever the iterate is nonzero. In particular a locally finite δ_w cannot be transverse to a nonproper divisor: its orbit of H would lie in a fixed finite-dimensional polynomial space, whose traces have a uniform denominator bound, contradicting (8.1).

This is a useful connection to the earlier common-boundary-bound obstruction. It does not establish local finiteness of the dual fields; imposing that without proof would reintroduce a Jacobian-conjecture-strength hypothesis. It also does not confuse local finiteness with finite jets or with divergence-free vector fields.

## 9. What would actually advance the planar proof program

The concrete next target is a polynomial regularity statement for genuine inverse traces in selected higher-degree families, not the existence of a formal logarithm. Quartic models are controls because generic degree four is already excluded globally:

- In a selected higher-degree family, compute T_1,...,T_(d-1) in a primitive monic transverse gauge. Polynomiality of all these moments gives the automorphism criterion; for even d, polynomiality of the first d-2 moments gives the conditional alternating-monodromy gate. No exact first open generic degree is asserted.
- If a pole is found, retain its order, the second inverse coordinate, and the finite normalization boundary. Its derivative tower (8.1) quantifies what a proposed bounded-pole approximation must fail to control.
- Retain the quartic A4 cubic root-of-unity packet as a local hostile: it is compatible with this test but cannot globally realize a polynomial Keller map by Domrina. To strengthen the method for higher degrees, retain the rational-boundary/forest sidecar or another genuine polynomial realization obstruction. Repeating scalar residue sums cannot supply it.

The new synthesis does not prove that these trace functions are polynomial for arbitrary Keller maps. The unconditional denominator compiler is an exact reformulation/measurement of nonproperness with fewer moments, and the even-degree gate is a conditional monodromy exclusion. This is the stopping boundary of the present result.

## 10. Sources, dependencies, and reproduction

### Directly checked mathematical sources

- Stacks Project, Zariski's Main Theorem, Lemma 37.43.3, Tag 05K0: https://stacks.math.columbia.edu/tag/05K0 . Used only for the open-immersion/finite factorization.
- Stacks Project, algebraic Zariski Main theorem, Tag 03GT: https://stacks.math.columbia.edu/tag/03GT . Confirms the relative-normalization formulation.
- Stacks Project, normal-ring height-one intersection, Lemma 15.24.18, Tag 0AVB: https://stacks.math.columbia.edu/tag/0AVB . Applied to a normal ring as a module over itself; the elementary UFD intersection for C[u,v] is also immediate.
- Stacks Project, fibers of étale morphisms, Tag 02GL: https://stacks.math.columbia.edu/tag/02GL . The quasi-finite/unramified local-sheet input.

Newton identities, the product-of-root-differences discriminant identity, and its Vandermonde sign interpretation are proved or derived in the text; no specialized classification theorem is silently imported. The final finite-Keller implication is justified by the elementary covering-space argument in Section 5.

### In-repository dependencies read

- `planar_autonomous_lie_jets_20261019`: polynomial autonomous rigidity and fixed-direction jet completion/boundary bound.
- `metabelian_global_invariants_20261019`: orbit-sum principal-part obstructions.
- `parabolic_cotangent_integrability_20261019`: formal/meromorphic/rational invariant distinctions and finite-jet hostiles.
- `finite_jet_bottleneck_clocks_20261019`: finite-jet clock data do not determine analytic embedding.
- `period_three_global_index_20261019`: global residue sum, with polynomial hypothesis explicit.
- THM-2621, inverse spectral coefficient poles and trace–Liouville exactness; THM-2633, finite-sheet survival and monodromy normal generation; THM-1365, repaired polynomial-deck scope.
- THM-3951, actual upstairs rational-boundary/forest obstruction.
- THM-500, tournament power sums versus support-overlap geometry; THM-4529, Petersen/Paley constructions and scope boundaries, used only for the loss ledger.

Reproduce the exact controls from this directory:

    python3 verify_trace_pole_boundary.py
    python3 -O verify_trace_pole_boundary.py

The controls verify Newton sums by both recurrence and companion-matrix traces, the two displayed map families, discriminant valuations in a finite trinomial universe, root-of-unity cancellation, and differential pole escalation. These finite computations supplement the proofs; they do not establish global Keller realization, a general pole bound, or an infinite monodromy classification by enumeration.
