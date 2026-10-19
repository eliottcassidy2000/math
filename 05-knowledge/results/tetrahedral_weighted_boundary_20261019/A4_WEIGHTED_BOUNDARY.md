# The tetrahedral quartic: a global cubic-pole hostile and a weighted-boundary obstruction

**Status:** proved scoped synthesis, with exact algebraic and finite-group verification. No historical novelty claim. This packet is separate from, and does not modify, the frozen trace-pole packet.

Work over **C**. A polynomial Keller map means a polynomial map A² → A² whose Jacobian determinant is a nonzero constant. Regular vector fields on a localization, a constant rational Jacobian, and polynomiality of a few inverse trace moments are weaker conditions.

## Results first

There is an explicit polynomial map

\[
 F(h,e)=(9h^4e^2-8he,\,-27h^6e^3+36h^3e^2-8e)=(u,v)                 \tag{1}
\]

with the following simultaneous properties:

* Its generic degree is four and its field monodromy is **A4**, acting on the four vertices of a tetrahedron.
* Its Jacobian is **64e**. Its restriction to U = A¹ × Gm is étale, with globally regular commuting dual derivations on C[h,e,e⁻¹].
* The finite normalization X of the target in C(h,e) is the quadratic cone b² = ac. The removed cubic ramification curve E is a smooth rational A¹, and its incidence obstruction is a forest.
* Over the generic point of the sole finite branch curve Δ = u³ − v² = 0 there is exactly one retained finite sheet and three escaping h sheets. The exact inverse equation is
  \[
    27\Delta h^4+64vh-64u=0.                                    \tag{2}
  \]
  Consequently Tr(h) = Tr(h²) = 0, while Tr(h³) = −64v/(9Δ).
* Yet X has Cl(X) = Z/2. This excludes **every** affine-plane Keller chart with this exact field extension C(u,v) ⊂ C(h,e), not only the displayed coordinates. The regular dual derivations cannot preserve any polynomial subring C[r,s] with the same fraction field and containing u,v.

Thus the retained cubic packet, global A4 monodromy, rational unibranch boundary, an unweighted forest, and regular commuting dual derivations coexist in a genuine algebraic example. The missing condition is polynomial regularity on an actual affine plane. In this example it is obstructed by a weighted divisor class, independently of the trace calculation.

A useful general consequence is proved below: **if the finite normalization is finitely dominated by a normal factorial affine variety, a polynomial Keller chart forces degree one.** In particular, a nontrivial Keller normalization cannot be a finite quotient of an affine plane.

There is also an explanatory classification: every connected normal quartic cover with A4 Galois closure and affine branch support contained in a polynomial-coordinate copy of the global standard cusp u³ = v² has this quadratic-cone normalization. This does not claim a new exclusion of a Keller map whose entire nonproper-value set is that cusp: Nguyen Van Chau's 1999 theorem already excludes any such curve homeomorphic to C, in every degree. The classification additionally keeps finite branch support distinct from possible unramified omitted divisors.

## 1. The global polynomial-regularity gate

Let B = C[u,v], K = C(u,v), L/K a finite extension, and X = Spec R the normalization of A² in L. Characteristic zero gives unique extensions δu, δv to L of the two coordinate derivations of K. Their commutator vanishes: it is a derivation of L zero on K, and an algebraic separable extension has no nonzero such derivation.

### Proposition 1: stable polynomial models are exactly Keller models

Suppose A = C[r,s] is a polynomial subring of L, contains B, and has fraction field L. Then

\[
 \delta_u A\subseteq A,\quad\delta_v A\subseteq A
 \quad\Longleftrightarrow\quad
 \det\frac{\partial(u,v)}{\partial(r,s)}\in C^*.                 \tag{3}
\]

Indeed, the matrix whose columns are (δu r, δu s) and (δv r, δv s) is an inverse over A to the Jacobian matrix of (u,v). Thus its determinant is a unit in A. Conversely a constant nonzero Jacobian gives the usual polynomial inverse matrix. This uses **no** local finiteness, completeness, nilpotence, or analytic embedding hypothesis.

Every element of R is integral over B and therefore integral over A. Since A is integrally closed in L, R ⊂ A. When (3) holds, the resulting map Spec A → X is the open immersion supplied by the normal form of Zariski's Main Theorem. Hence an obstruction to an affine-plane open in X is an obstruction to every stable polynomial model in (3).

### Proposition 2: the affine-plane boundary classes are a free basis

Let X be a normal integral affine surface and U ⊂ X a dense open isomorphic to A². If E₁,…,Eₛ are the prime divisors in X \ U, then

\[
 \operatorname{Cl}(X)\simeq\bigoplus_{i=1}^s\mathbb Z[E_i].      \tag{4}
\]

The localization sequence for Weil divisors is

\[
 O(X)^*\longrightarrow O(U)^*\longrightarrow
 \bigoplus_i\mathbb Z E_i\longrightarrow\operatorname{Cl}(X)
 \longrightarrow\operatorname{Cl}(U)\longrightarrow0.
\]

Both unit groups equal C*, and Cl(U) = 0. Therefore the middle divisor map is both injective and surjective. Equivalently, a principal relation supported on the boundary would restrict to a unit on A², hence to a constant, so the relation is zero.

**Codimension-two caveat.** No purity assumption about X \ U is needed for (4): codimension-two components do not occur in the Weil divisor sequence. If s = 0, normal Hartogs gives O(X) = O(U) = C[r,s], since removing a codimension-two set changes no height-one valuation. Both X and U are affine, so the given open immersion is an isomorphism. Thus a proper affine-plane open forces s > 0. In the actual Keller normalization, the earlier boundary-pole lemma also identifies the boundary with coordinate-polar divisors; (4) does not depend on that stronger fact.

In particular, nonzero torsion in Cl(X) forbids **any** affine-plane open. If Cl(X) is a torsion group, the only possible affine-plane open is all of X and Cl(X) must be zero.

This is the classical localization mechanism already proved in corpus THM-3922. Its application here is coordinate-independent; detecting a denominator in one formula would be weaker.

### Corollary 3: finite factorial domination excludes nontrivial Keller models

Suppose a normal factorial affine variety V admits a finite dominant map p: V → X, of field degree n. Then n Cl(X) = 0. Indeed, define divisorial pullback using ramification indices. For any Weil divisor D, factoriality makes p*D principal, say div(f). The valuation identity Σ eᵢfᵢ = n and the field-norm identity give

\[
 nD=p_*p^*D=p_*\operatorname{div}(f)
       =\operatorname{div}(N_{C(V)/C(X)}f).                    \tag{5}
\]

These identities are height-one statements and do not require p to be flat everywhere. Thus (4) forces no boundary for an affine-plane Keller model. Its map to the target is finite and étale. A connected finite étale cover of C² is trivial, so [L:K] = 1 and the map is an automorphism.

In particular the conclusion applies if X is a finite quotient of A², or more generally is finitely dominated by a factorial affine variety. It does **not** assert that an arbitrary A4 normalization has torsion class group or admits such a finite domination.

For the norm-pushforward identity see [Stacks, Tag 02RT](https://stacks.math.columbia.edu/tag/02RT). Zariski's Main Theorem and the normal height-one intersection description are recorded at [Tag 03GT](https://stacks.math.columbia.edu/tag/03GT) and [Tag 0AVB](https://stacks.math.columbia.edu/tag/0AVB).

## 2. An exact tetrahedral reflection construction

Put ω² + ω + 1 = 0, δ = 1 + 2ω, so δ² = −3. On A² with coordinates ξ,η define

\[
 r=\begin{pmatrix}1&0\\0&\omega\end{pmatrix},\qquad
 q=\delta^{-1}\begin{pmatrix}1&2\\1&-1\end{pmatrix}.             \tag{6}
\]

Set

\[
 A=\xi^3\eta-\eta^4,\quad
 B=\xi^4+8\xi\eta^3,\quad
 C=\xi^6-20\xi^3\eta^3-8\eta^6.                                 \tag{7}
\]

Direct substitution gives

\[
 r(A)=\omega A,\quad q(A)=A,\quad r(B)=q(B)=B,
 \quad r(C)=q(C)=C,
 \quad C^2=B^3-64A^3.                                          \tag{8}
\]

The projective actions preserve the four zero lines of A. The action of r fixes the line η = 0 and cycles the other three. The action of q exchanges two disjoint pairs. Hence their projective group is A4. A scalar in the group preserves B and C, so its fourth and sixth powers are one; it is therefore ±I. Since q² = −I, G = ⟨r,q⟩ has order 24, center {±I}, and G/{±I} = A4.

This is the **G4 complex-reflection representation**: det r = ω and det q = 1. It is not the determinant-one binary-tetrahedral representation. Confusing these representations would give the wrong quotient surface. The formulas are classical tetrahedral invariant theory; see Elkies's [explicit covariants and reflection twist](https://abel.math.harvard.edu/~elkies/M263.12/index.html), October 31 section. The identities and the full finite group are independently checked here.

Let H = ⟨r,−I⟩, of order six. Its polynomial invariant ring is

\[
 C[\xi,\eta]^H=C[a,b,c]/(b^2-ac),\qquad
 a=\xi^2,\ b=\xi\eta^3,\ c=\eta^6.                              \tag{9}
\]

To see generation, an invariant monomial has η exponent divisible by three and the sum of its ξ exponent and that exponent even. These monomials are generated by ξ², ξη³, η⁶. Thus X = A²/H is a normal quadratic cone. The invariant map to the target is

\[
 u=a^2+8b,\qquad v=a^3-20ab-8c.                                \tag{10}
\]

Eliminating b,c gives

\[
 f(a)=27a^4-18ua^2-8va-u^2=0,\quad
 b=(u-a^2)/8,\quad c=(7a^3-5ua-2v)/16.                          \tag{11}
\]

These expressions are inverse ring presentations. Since 27 is invertible, R is free of rank four over C[u,v]. In particular it is exactly the finite normalization of the target in its fraction field.

The map A²_(ξ,η) → X has degree six. Consequently C(ξ,η)/C(B,C) has degree 24, and the invariance and group order show C(B,C) = C(ξ,η)^G. The core of H in G is {±I}, because a vertex stabilizer C3 has trivial core in A4. Hence the Galois closure of the quartic field is

\[
 M=C(\xi,\eta)^{\{\pm I\}},\qquad\operatorname{Gal}(M/K)=A4.     \tag{12}
\]

This is an exact monodromy identification, not an inference from a square discriminant. The independent discriminant check is

\[
 \operatorname{Disc}_a f=-80621568\,(u^3-v^2)^2.                 \tag{13}
\]

## 3. The actual missing divisor and the retained sheet

Define the regular function and curves on X

\[
 e=a^3-2ab+c,
 \quad E=\{b=a^2,\ c=a^3\},
 \quad C_0=\{b=c=0\}.                                          \tag{14}
\]

Both E and C₀ are smooth affine lines parametrized by a. The cone is singular only at their common point a = b = c = 0. In the original variables e = (ξ³ − η³)². One has

\[
 \operatorname{div}_X(e)=2E,
 \qquad\operatorname{div}_X(\Delta)=C_0+3E,
 \quad\Delta=u^3-v^2.                                         \tag{15}
\]

For example the group H acts transitively on the three lines ξ³ = η³. At their generic points it acts freely, so e has order two along their quotient E. The remaining factor of Δ = 64η³(ξ³−η³)³ descends to C₀ with multiplicity one: the order-three reflection along η = 0 removes the η³ multiplicity in the quotient. This also follows directly from (11).

On Δ = 0, write u = z² and v = z³. Then

\[
 f(a)=(a-z)(3a+z)^3.                                         \tag{16}
\]

The finite companion is C₀, with a = z. The cubic ramification component is E, with a = −z/3. Both residue maps have degree one. The cubic inertia fixes one tetrahedral vertex and rotates the opposite three vertices. A single such three-cycle normally generates A4; no second branch component is forced by normal generation.

Set U = X \ E = D(e). The cone singular point lies on E. On U put

\[
 h=(a^2-b)/e.                                                 \tag{17}
\]

The exact inverse formulas are

\[
 a=h^2e,\qquad b=h^4e^2-he,\qquad
 c=e(h^3e-1)^2.                                               \tag{18}
\]

Thus **U ≅ A¹_h × Gm,e**, as affine schemes. Substituting in (10) gives (1), with

\[
 JF=64e,\qquad\Delta\circ F=64e^2(h^3e-1).                     \tag{19}
\]

The map U → A² is étale, and X → A² is étale off E. Its affine branch support is exactly the cusp Δ = 0. Its nonproper-value set as a map from U is also that cusp, since X is a finite completion and E maps onto the cusp. The origin is omitted by U. In the polynomial extension (1) on all of A²_(h,e), the line e = 0 is contracted to the origin and the Jacobian vanishes there. Therefore (1) is explicitly **not** Keller.

The divisor E is rational, smooth, and unibranch, and there are no two distinct omitted primes giving a repeated-incidence cycle. A full SNC compactification of U also has an unweighted tree: in P¹_h × P¹_e the boundary is e = 0, h = ∞, e = ∞, with chain graph. Thus even full unweighted boundary-tree information does not distinguish this U from an affine plane; the boundary divisor relations and weights matter.

## 4. The frozen trace packet is realized, but its conclusion is not promoted

Eliminating a using (11),(17), or substituting (1), gives (2). There is also the rational recovery identity

\[
 e=\frac{27h^3\Delta+16v}{16(9h^2u-8)}.                      \tag{20a}
\]

Its denominator is not identically zero (evaluate at h = 0). Hence K(h) = C(h,e), and the extension degree four proves that (2) is its irreducible minimal polynomial up to a nonzero scalar. Its monic form is

\[
 h^4+\frac{64v}{27\Delta}h-\frac{64u}{27\Delta}=0.              \tag{20}
\]

Newton identities yield

\[
 T_1=T_2=0,\qquad T_3=-\frac{64v}{9\Delta},\qquad
 T_4=\frac{256u}{27\Delta}.                                    \tag{21}
\]

At the generic point of E, the numerator a²−b has order one and e has order two; hence ord_E h = −1. By (15), ord_E Δ = 3. Therefore the three escaping roots have target-normalized pole order 1/3 and their leading coefficients form a nonzero scalar multiple of the cube roots of unity. The fourth root tends to u/v = 1/z on C₀. No zero leading coefficient is being included among the escaping three.

This is a **global polynomial non-Keller hostile**, with the same A4 packet that survives the earlier sign test. It does not say that first-two-moment regularity is automatic for a Keller map, nor that every A4 extension has this equation. The purpose is to isolate a remaining necessary datum: the polynomial source ring on which the dual derivations act.

## 5. Polynomial dual derivations: the exact denominator debt

On C[h,e,e⁻¹], the target dual derivations are

\[
\begin{aligned}
 \delta_u h&=\frac{-81h^6e^2+72h^3e-8}{64e},&
 \delta_u e&=\frac{81h^5e^2-54h^2e}{32},\\
 \delta_v h&=-\frac{9h^4}{32}+\frac{h}{8e},&
 \delta_v e&=\frac{9h^3e-2}{16}.
\end{aligned}                                                \tag{22}
\]

They commute, preserve this localized ring, and satisfy δu(u,v) = (1,0), δv(u,v) = (0,1). They do not preserve C[h,e]. More strongly, for constants α,β the coefficient of e⁻¹ in (αδu + βδv)h is (−α + βh)/8. It vanishes identically only for α = β = 0. No invertible affine-linear target coordinate change removes the pole from a dual frame on this displayed affine-plane extension.

That calculation alone would not exclude a different polynomial model. Propositions 1–2 and the class-group calculation below do exclude every model with the same L/K. This is the load-bearing global step.

There is even a rational constant-Jacobian presentation: put p = e and q = he, which are birational coordinates on U. Then

\[
 u=9q^4/p^2-8q,\qquad
 v=-27q^6/p^3+36q^3/p-8p,
 \qquad\det\frac{\partial(u,v)}{\partial(q,p)}=64.             \tag{23}
\]

Its poles at p = 0 remain. Thus rational symplectic normalization of the volume form does not repair polynomiality.

## 6. The weighted obstruction that the forest loses

For R = C[a,b,c]/(b²−ac), localization at a gives R[a⁻¹] = C[a,a⁻¹,b], which is factorial. The unique prime over a = 0 is L₀ = (a,b), and div(a) = 2L₀. Hence Cl(X) is generated by [L₀], of order dividing two.

The class is nonzero. Identify the same cone abstractly with C[s²,st,t²] ⊂ C[s,t]. If L₀ were principal, pullback to A²_(s,t) would make its defining rational function have divisor (s = 0). Such a function is a constant times s, whereas every function from the cone is invariant under (s,t) ↦ (−s,−t). Contradiction. Consequently

\[
 \operatorname{Cl}(X)=\mathbb Z/2.                            \tag{24}
\]

One can see the actual omitted class without changing coordinates:

\[
 \operatorname{div}(b-a^2)=E+L_0,
 \qquad \operatorname{div}(e)=2E.
\]

Thus [E] = −[L₀] is the nonzero two-torsion class. The localization sequence for this actual open U records the precise defect: O(U)* / C* = Z generated by e, and its divisor map Z → Z[E] is multiplication by two, with cokernel Z/2. A forest records neither that multiplication nor the nonconstant unit. If an affine-plane open existed, (4) would force all boundary classes to be independent of infinite order. This contradiction is intrinsic to X and excludes every Keller chart inducing the same quartic extension.

The minimal resolution of the cone has one exceptional P¹ with self-intersection −2. Explicit charts are (a,t) with b = at, c = at², and (c,w) with b = cw, a = cw². On the overlap w = t⁻¹ and c = at²; this is the total space of O_P¹(−2). The exceptional intersection matrix is [−2], with determinant of absolute value two. Its unweighted graph is a single vertex, hence a tree. This is an exact example in which forgetting intersection weights loses the torsion obstruction. The class-group result above was proved algebraically, so it is not inferred from a graph determinant without a divisor-lattice identification.

The original THM-3951 forest obstruction remains correct and useful: repeated upstairs boundary incidence or positive genus can exclude an affine plane. This model passes those coarse tests and is excluded by the stronger class-group condition from THM-3922.

## 7. The cubic resolvent and its genuine standard plane

The character by which G acts on A equals the determinant character in (6). Its kernel Q has order eight and lies in SL₂. Its central involution is −I, and its quotient by that involution is V4. Thus Q ≅ Q8. The exact companion verifies six elements of order four and the unique involution.

The cyclic cubic resolvent normalization is

\[
 Y=\operatorname{Spec} C[A,B,C]/(C^2-B^3+64A^3)
      =A^2/Q_8,\qquad (u,v)=(B,C).                            \tag{25}
\]

The hypersurface is normal, singular only at the origin, and finite of degree three over C[u,v], so the invariant-field identification proves the displayed invariant-ring equality. After rescaling A, its equation is w³ = u³ − v². This is the D4 surface, whereas the quartic normalization X is A1. The classical Q8/D4 identification and resolution are also described in [Mikosz–Weber, Section 7.1](https://link.springer.com/article/10.1007/s10801-020-00938-3); the particular invariant ring needed here was derived above.

A nonidentity element of Q8 has no eigenvalue one: it has finite order and determinant one, so an eigenvalue one would make it the identity. Therefore Q8 acts freely on A² \ {0}. Its quotient is Y_reg, and C² \ {0}, which retracts to S³, is simply connected. Hence

\[
 \pi_1(Y_{\rm reg}^{an})=Q_8.                                 \tag{26}
\]

By the complex Riemann existence theorem, finite étale covers are classified by the finite quotients of this group. See [Milne, Lectures on Étale Cohomology, Theorem 21.3](https://www.jmilne.org/math/CourseNotes/LEC210.pdf). In particular its unique connected V4 quotient has kernel {±I} and is

\[
 Z=(A^2\setminus\{0\})/\{\pm I\}\longrightarrow Y_{\rm reg}.    \tag{27}
\]

This is exactly the V4 cover in the quartic Galois closure, not an unrelated D4 diagram.

Moreover Pic(A² \ {0}) = 0 and its units are constants. Descent of line bundles along this free finite cover therefore gives

\[
 \operatorname{Cl}(Y)=\operatorname{Pic}(Y_{\rm reg})
       =\operatorname{Hom}(Q_8,C^*)\simeq (\mathbb Z/2)^2.       \tag{28}
\]

The C3 conjugation action cycles the three nontrivial characters, so the full two-torsion is the required standard F2[C3] plane. The D4 resolution matrix has determinant four and Smith diagonal (1,1,2,2); its three arms are cycled. This recovers the precise THM-2655/2703 carrier in the **actual** cyclic resolvent of this quartic. It does not identify Cl(Y) with Cl(X); the two class groups are different and were computed separately.

## 8. A scoped global-standard-cusp classification

### Theorem 4

Let π: X′ → A²_C be a connected finite normal quartic cover whose Galois closure has group A4 in its transitive degree-four action. Suppose that, after a polynomial automorphism of the target, the **entire affine branch support** of π is contained in D = {u³ − v² = 0}. Then, in those target coordinates, X′ over A² is isomorphic to the cover (9)–(10). In particular Cl(X′) = Z/2, and X′ admits no dense open isomorphic to A².

**Proof.** Let M′/K be the Galois closure and E′ = (M′)^V4. It is cyclic of degree three and is unramified at affine prime divisors outside D. Write E′ = K(g^{1/3}) by Kummer theory. For every affine prime other than Δ, the valuation of g is divisible by three. Unique factorization of C[u,v] and the fact that nonzero complex constants are cubes reduce g, modulo cubes, to Δ^j. The residue j cannot be zero, since then E′/K would be trivial; equivalently a nontrivial connected finite étale cover of A² cannot exist. The fields for j = 1 and j = 2 coincide. Thus E′ is the cubic field in (25), with the same normalization Y.

At D the inertia of M′/K maps nontrivially to C3. Finite inertia in characteristic zero is cyclic, and the only such subgroups of A4 with nontrivial image in C3 have order three. Therefore it meets V4 trivially. Outside D inertia is trivial. Thus normalization in M′ is unramified in codimension one over Y. Purity on the smooth surface Y_reg makes its restriction a connected finite étale V4 cover. This use of purity is exactly [Stacks, Tag 0BMB](https://stacks.math.columbia.edu/tag/0BMB); it does not apply across the singular origin of Y.

By (26), the only connected V4 quotient is Q8/{±I}. Therefore M′/E′ is the cover (27), uniquely as an unmarked cover. The explicit M in (12) is that cover and is already Galois over K. Hence M′ and M are K-isomorphic. The degree-four subfields are fixed fields of order-three subgroups of A4. All these subgroups are conjugate, so they give isomorphic K-extensions. Uniqueness of normalization now gives X′ ≅ X. Equation (24) and Proposition 2 conclude the proof. ∎

The hypothesis is **global**. A local cusp singularity in a larger branch curve does not identify the global cubic resolvent with A²/Q8, does not determine its fundamental group, and does not force (24). Nor does an arbitrary rational branch curve have this normalization.

### Literature boundary

Nguyen Van Chau, *A Remark on Vitushkin's Covering*, Acta Mathematica Vietnamica 24 (1999), 109–115, Theorem 1.1, proves that a Keller map which is an unbranched covering outside a curve homeomorphic to C is bijective. [Primary paper](https://math.ac.vn/uploads/files/9901109.pdf). The standard cusp is homeomorphic to C through z ↦ (z²,z³), so a Keller map whose entire nonproper-value set is that cusp is already excluded, in every degree. We claim no new frontier closure in that familiar special case.

Theorem 4 is retained for its explicit normalization and weighted-class explanation. Its premise concerns branch support of the **finite normalization**. For a general quasi-finite source open, there can also be omitted unramified divisors, so branch support must not silently be replaced by the full nonproper-value set. The class-group proof excludes all affine-plane opens of this exact X without needing that replacement.

## 9. What remains in the A4 lane

A genuine nontrivial planar Keller normalization must evade all the following:

1. It cannot have nonzero torsion in Cl(X); its actual omitted prime classes must form a nonempty free integral basis.
2. It cannot admit a finite dominant morphism from a normal factorial affine variety.
3. Its two rational target-dual derivations must preserve one common polynomial ring C[r,s] containing the target ring and having the full quartic fraction field.
4. Its cubic resolvent must carry the standard F2[C3] plane in H¹_et(Y_reg, μ₂). Under the Kummer sequence, this plane may lie in unit squareclasses or inject into Cl(Y)[2], as in THM-2655; it need not come from the global D4 quotient (25). A local D4 singularity alone supplies no global classification.

None of these items is proved automatic or impossible for all A4 covers. The present packet supplies an exact hostile to discarding the polynomial ring and the weighted divisor lattice, plus a coordinate-independent exclusion of the concrete reflection-quotient family. It does not assume first-two-trace regularity for an arbitrary Keller map, derive global polynomial realization from a formal cotangent model, or settle the general A4 or planar Jacobian problem.

## 10. Dependency and verification ledger

Inherited corpus inputs: THM-3922 (boundary class basis); THM-3951 (upstairs rationality/incidence forest); THM-2633 (allowed quartic inertia and normal-generation scope); THM-2655 (actual resolvent V4 torsor); THM-2703 (weighted three-arm lattice); THM-2690 (do not infer absence of two-torsion from odd degree alone).

Classical external inputs are linked where used: tetrahedral invariant formulas, norm-pushforward, Zariski Main, normal Hartogs, purity, complex finite-cover comparison, and Chau's stronger nonproper-set theorem. The elementary proofs above reconstruct the particular invariant rings, group representation, class groups, and global cover classification needed here. No historical novelty claim attaches to the classical tetrahedral construction or to these standard tools.

The companion checks polynomial identities, inverse presentations, the quartic resultants and traces, the two derivations and their commutator, all 24 exact matrices over Q(√−3), the induced A4 permutations and subgroup cores, the Q8 quotient, and the A1/D4 Smith determinants. These are checks of the displayed finite algebra. The written valuation, normalization, class-group, and covering arguments are mathematical proofs; a finite test is not substituted for them.
