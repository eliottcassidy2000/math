# The trace lattice hidden by a maximally cancelling inverse-sheet packet

20 October 2026. **PROVED LOCAL ALGEBRA; applicability to a polynomial Keller map is conditional on the stated trace regularity. Independent review pending.** No historical-priority claim, new Keller map, or resolution of JC(2) is asserted. **Literature correction:** planar generic degree four is already excluded by Domrina (2000); the quartic A4 calculations below are local-method controls, not a globally open planar degree. See Section 2.

## 1. Result and exact next target

The previous inverse-trace packet detects the nonproper-value set with power traces through order d−1. Its even-degree extremal case leaves a genuine obstacle: d−1 escaping sheets can cancel every lower trace, with one finite sheet left over. Here that scalar cancellation is converted into a precise, basis-independent defect of a trace lattice.

Let R be a complete discrete valuation ring whose residue field has characteristic zero, let K be its fraction field, and normalize v(t)=1. Let β₀,…,β_m be the roots of a monic separable polynomial

    χ(Z)=Z^(m+1)+c₁Z^m+…+c_m Z+c_(m+1) ∈ K[Z],   m≥2.

Let E=K[Z]/(χ), and let h be the class of Z. The algebra trace Tr_(E/K)(h^k) is the power sum T_k=Σ β_i^k in K. Leading coefficients below are taken after a finite splitting extension and reduction in its residue field. Suppose exactly m roots have negative valuation, one is integral, and T₁,…,T_(m−1) lie in R. Put a=−v(T_m).

Then:

1. a is a positive integer. Every escaping root has pole order a/m, and its leading coefficient belongs to one scalar multiple of the full m-th-root-of-unity packet.
2. Let G=(T_(i+j))_(1≤i,j≤m−1), the trace-pairing matrix on the ordered monomials h,h²,…,h^(m−1). Then

       t^a G ∈ GL_(m−1)(R).

   Every Smith exponent of G is exactly −a, and

       v(det G)=−a(m−1).

   Over the stated complete DVR, the least index length of a full-rank sublattice of the original monomial lattice with integral bilinear trace pairing is exactly ceil(a(m−1)/2).

   Thus maximal scalar cancellation forces a uniformly polar rank-(m−1) trace lattice. It does not make that lattice small or positive.
3. The finite root has a canonical rational first approximation

       r=−c_(m+1)/c_m ∈ R,      v(β₀−r)≥a.

   Newton iteration for χ starting at r is defined over K and satisfies

       v(r_n−β₀) ≥ (2^(n+1)−1)a,     r₀=r.

   These rational approximants select the finite sheet locally, with no claim of rational termination or a global inverse.

For a quartic, m=3. Assuming T₁,T₂ regular, the remaining cubic pole packet has

    G = [[T₂,T₃],[T₃,T₄]],
    v(T₃)=−a,     v(T₄)≥−a,
    v(T₂T₄−T₃²)=−2a.

The exact leading term of the determinant is −T₃². Moreover G is congruent over R to

    t^(−a) [[0,1],[1,0]].

This is a hyperbolic trace form, not a positive ellipsoid. A full-rank sublattice on which this bilinear form becomes integral must have index length at least a, and that bound is sharp. The location of that compulsory index is a concrete new sidecar for the local quartic A4 model; the global planar degree-four case is already excluded.

The useful next question is consequently precise: **can the actual polynomial source and its boundary incidence force zero index for the original trace plane, or obstruct the finite-sheet selector's global gluing?** Merely constructing some repaired integral trace lattice is insufficient. Section 7 gives an exact hostile with such a repair and with poles still present.

## 2. Inheritance and scope

The closest proved mechanism is the packet
`05-knowledge/results/inverse_trace_pole_boundary_20261019/TRACE_POLE_BOUNDARY.md`, especially its local moment lemma, finite-sheet survival, and quartic A4 stopping boundary. Its global input is THM-2633, `01-canon/theorems/THM-2633-derangement-character-obstruction-and-d4-keller-exclusion.md`.

For a genuine complex polynomial Keller map F=(P,Q), choose a primitive affine source coordinate h=x with P monic in the transverse source coordinate. At the generic point of a target prime divisor, the finite normalization supplies at least one actual finite inverse sheet. If the first d−2 trace moments are polynomial and a boundary divisor is present, the earlier lemma forces exactly d−1 escaping sheets and one finite sheet. Thus the theorem applies with m=d−1 to that actual inverse algebra. Completion and, when convenient, residue-field extension do not alter the displayed valuations.

For d=4 the earlier trace/discriminant method by itself leaves A4 after its group gates. That is a limitation of that method, not a globally open planar case. A. V. Domrina, *On four-sheeted polynomial mappings of C². II. The general case*, Izvestiya: Mathematics 64 (2000), 1–33, excludes all four-sheeted polynomial plane maps with nonzero constant Jacobian. The [primary publisher page](https://www.mathnet.ru/eng/im273) explicitly confirms this statement; DOI https://doi.org/10.1070/im2000v064n01ABEH000273 . The repository's `05-knowledge/reference/CORE-PAPERS-PLANAR-JC.md` records an earlier direct full-PDF check of the generic-preimage definition on page 1 and the general theorem, distinguishing the more restrictive 1998 paper. In this session the publisher statement was checked directly; the full PDF endpoint timed out. The stale interpretation of the Oct 19 A4 residual as globally open is not inherited.

The local all-m theorem remains nonvacuous as an algebraic statement. For new planar applications one should use a higher generic degree, for example the even d≥6 maximal-cancellation sectors, without asserting their low traces are automatically regular. This note does not establish that premise for an arbitrary Keller map or assert that every local packet has a polynomial source-plane realization.

The canonical hostile is the punctured local algebra

    K × K[Z]/(Z^m−t^(−a)),

which has a finite root zero and the exact root-of-unity pole packet. The corrected near miss is to replace a complex bilinear trace form by a positive Hermitian Gram form without retaining conjugation. The least-used sidecar retained here is the integral lattice and its index, followed by the finite-sheet selector's overlap data.

Concept board: inverse sheets; Newton traces; polar trace lattice; multiplicative closure; finite-sheet selector; real versus complex positivity.

## 3. Recovering the packet, without positivity

Let λ be the greatest pole order of the roots, and let r be the number attaining it. Scale those roots by a uniformizer power of valuation λ in a finite extension, and call their nonzero leading coefficients α₁,…,α_r.

If Σ α_i^k vanished for k=1,…,r, the characteristic-zero Newton identities would make every elementary symmetric polynomial in the α_i vanish, contradicting their nonzero product. Thus one power sum of order at most r detects the maximal pole.

Since T₁,…,T_(m−1) are integral and there are exactly m escaping roots, r must equal m. Their first m−1 leading power sums vanish. Newton identities then say that their leading polynomial is Z^m−c for c≠0. In particular all leading coefficients are distinct, Σ α_i^m=mc≠0, and

    v(T_m)=−mλ.

Because T_m belongs to K, a=mλ is an integer. This proves part 1. Root-of-unity cancellation, rather than an inequality between positive numbers, is the operative mechanism.

## 4. The entire trace plane has one polar Smith slope

Newton identities for χ first give

    c₁,…,c_(m−1) ∈ R,      v(c_m)=−a.

Indeed the m-th identity is m c_m+T_m plus an integral expression. Division by m is harmless because the residue characteristic is zero.

Let β₀ be the integral root. Substitution into χ yields

    c_(m+1)=−β₀^(m+1)−c₁β₀^m−…−c_(m−1)β₀²−c_mβ₀,

so v(c_(m+1))≥−a.

Now use the higher Newton recurrence. For m≤k≤2m−2, every term involving c_m is multiplied by a trace of order at most m−2, and every term involving c_(m+1) is multiplied by a still lower trace (or occurs in the initial Newton identity with its integer coefficient). Those traces are integral. The other coefficients c_j are integral, and induction on k gives

    v(T_k)≥−a,       m≤k≤2m−2.                  (1)

Consequently t^aG is integral. Its entries with i+j<m vanish in the residue field, and every entry on the reverse diagonal i+j=m equals the same unit t^aT_m. Reversing the columns turns the residue matrix into a triangular matrix with that unit on every diagonal position. Therefore

    det(t^aG) mod t
      = (−1)^((m−1)(m−2)/2) (t^aT_m mod t)^(m−1) ≠ 0.   (2)

This proves t^aG invertible, all Smith exponents −a, and the determinant valuation. It also proves that the displayed monomials are linearly independent: a linear dependence would annihilate their Gram matrix. In a degree-(m+1) inverse field this independence also follows from primitivity.

The stronger Smith statement is not inferred from the determinant alone. It uses (1), which places every entry in t^(−a)R, together with (2), which proves the scaled matrix is invertible.

## 5. Quartic hyperbolicization and the sharp index budget

For m=3, write B=t^aG. It has the form

    B=[[u,c],[c,v]],   u∈tR, c∈R*, v∈R.

The polynomial u+2cz+vz² has value in tR at z=0 and has unit derivative 2c there. Hensel's lemma gives z∈tR with u+2cz+vz²=0. Thus e=(1,z) is a primitive isotropic vector for B. With w=(0,1), the scalar B(e,w)=c+zv is a unit. Replace w by

    w−B(w,w)e/[2B(e,w)]

and rescale it by a unit. The resulting pair e,f is an R-basis with

    B(e,e)=B(f,f)=0,       B(e,f)=1.

Therefore an actual U∈GL₂(R) satisfies

    U^T G U=t^(−a) [[0,1],[1,0]].                  (3)

No square root or algebraic closure of the residue field is needed. Completeness is used only for the Hensel lift; the matrix-entry and Smith statements already hold over a noncomplete DVR.

Let L=R h+R h², and let M⊂L be a full-rank sublattice. Its index length is ℓ=length_R(L/M). If the bilinear trace pairing is integral on M, its Gram determinant has nonnegative valuation. Under a basis matrix for M the determinant changes by the square of its determinant, hence

    0≤−2a+2ℓ,       so ℓ≥a.                       (4)

The bound is attained: in the basis (3), take M=R(t^a e)+Rf. Its Gram matrix is the integral unimodular hyperbolic matrix and its index length is a.

If one insists on a uniform scalar repair M=t^bL, then b≥ceil(a/2), and ℓ=2b. For odd a this uses one extra unit of index compared with the sharp asymmetric repair. More generally, in the hyperbolic basis the diagonal repairs t^b e,t^c f are integral exactly when

    b,c≥0,       b+c≥a.

The integer solutions at minimal budget form the segment b+c=a. The symmetric box choice b=c pays a parity cost when a is odd. This is an exact simplex/box relationship in a valuation-allocation domain, not an identification of Euclidean cubes with their polars.

### 5.1 Sharp all-rank repair over the stated complete DVR

Let r=m−1 and retain the complete DVR R with residue characteristic zero. A full-rank sublattice M of the original monomial lattice has determinant valuation −ar+2ℓ, where ℓ is its index length. Integrality of the pairing therefore forces

    ℓ≥ceil(ar/2).                                  (4a)

This lower bound is sharp for the particular unimodular form B=t^aG. Its first floor(r/2) basis vectors span a totally isotropic subspace after reduction: every corresponding entry has i+j<m. Since the residue form is nondegenerate, this is the maximal possible Witt index. The residue form therefore splits into floor(r/2) hyperbolic planes and at most one unit line. Choose a residue isotropic vector and a vector pairing nontrivially with it, lift the isotropic vector by the same one-variable Hensel argument as above, and split off the resulting unimodular hyperbolic plane over R. Repeating in the orthogonal complement leaves at most one unimodular line. This uses the special isotropic subspace forced by the trace matrix, not an assertion that every form over an arbitrary residue field is split. No algebraic closure or positivity assumption is needed.

For each plane of t^(−a) times that unit form, scale one isotropic generator by t^a. This costs index a per plane. For the remaining unit line, if any, scale its generator by t^ceil(a/2), costing index ceil(a/2). The total is

    a floor(r/2)+(r mod 2)ceil(a/2)=ceil(ar/2).

A uniform scalar repair t^bL instead costs r ceil(a/2). When a is odd, its excess over the optimal budget is floor(r/2). These are comparisons among precisely specified lattice repairs. Neither repair implies that the original inverse coordinates have become integral, as the next hostile demonstrates.

## 6. A rational selector for the finite sheet

Since v(c_m)=−a and v(c_(m+1))≥−a, r=−c_(m+1)/c_m lies in R. Dividing χ by c_m gives

    f(Z)=χ(Z)/c_m=Z−r+A(Z),

where every coefficient of A is in t^aR and A has no constant or linear term. Substituting the finite root β₀ gives

    β₀−r=−A(β₀) ∈ t^aR.

For every z∈R, f'(z)=1+A'(z) is a unit. Newton iteration

    r_(n+1)=r_n−f(r_n)/f'(r_n)

therefore stays in R and is rational in the coefficients of χ. If E=v(z−β₀)>0, Taylor expansion around z shows

    v(z−f(z)/f'(z)−β₀)≥2E+a.                    (5)

All second and higher Taylor coefficients contain the factor t^a. Starting with E≥a gives E_n≥(2^(n+1)−1)a.

This is a coefficient-level selector, not an additional inverse sheet. If χ is irreducible over an original global target field K₀, the integral root generally belongs only to its completed extension K, not to K₀; no finite Newton iterate is claimed to equal it. Over the complete field in the theorem, χ necessarily has this one linear factor. The rational formulas have coefficients in K₀ when χ does, but approximate one chosen boundary sheet and need not agree with choices made at other divisors.

For a centered quartic, c₁=0, the coefficients can be read directly from traces:

    c₂=−T₂/2,  c₃=−T₃/3,  c₄=T₂²/8−T₄/4,
    r=(3T₂²−6T₄)/(8T₃).

Centering h by T₁/4 is legitimate globally only when T₁ is a polynomial target function, so h−T₁(F)/4 remains a source polynomial. The uncentered formula r=−c₄/c₃ requires no such change.

## 7. Hostiles: what an integral Gram repair still forgets

Take the exact local algebra with roots 0 and the three roots of Z³−t^(−a). Then

    T₁=T₂=T₄=0,    T₃=3t^(−a),
    G=3t^(−a) [[0,1],[1,0]].

The repaired lattice generated by t^a h and h² has an integral nondegenerate Gram matrix, but h² still has pole order 2a/3 on every escaping branch. Thus **integral, even unimodular, bilinear Gram data on an arbitrarily selected sublattice do not prove bounded inverse sheets**.

There is a stronger hostile even after keeping the original observable and the full inverse-algebra rank. In the same quartic algebra, take

    L′=span_R{1,h,t^a h²,t^a h³}.

It has rank four and contains both 1 and the original h. Its exact Gram matrix in that basis is

    [[4,0,0,3], [0,0,3,0], [0,3,0,0], [3,0,0,3]],

whose determinant is −27, a unit. Nevertheless h is polar, and L′ is not closed under multiplication: h² requires coefficient t^(−a) on its third basis vector. Thus even full rank, original-coordinate containment, and a unimodular trace pairing together do not replace multiplication closure.

The missing coordinate is multiplication closure. If a finite R-module inside an inverse algebra contains 1 and h and is stable under multiplication by h, then h is integral: the characteristic polynomial of its R-linear multiplication matrix annihilates h. That is a valid replacement, but constructing such a module containing the original h is a substantive new obligation; an unrelated repaired trace lattice does not supply it.

This also explains why a Hopf/Lie or primitive-generator analogy cannot be imported just by terminology. The genuine map here is h^k to the power-sum trace T_k and then to the Newton coefficients. It preserves the characteristic polynomial and valuation data when all needed coefficients and divisions are retained. Over residue characteristic p≤m, division by k is no longer valuation-neutral. For example the roots of Z^p−p^(−1), plus a finite zero root, have T₁,…,T_(p−1)=0 and T_p=1, even though p roots have pole order 1/p. The residue-characteristic-zero hypothesis is essential. Keeping characteristic coefficients repairs that diagnostic; an ordinary low ghost-moment integrality test does not.

## 8. The real-positive and complex-Keller boundary

For real numbers β_i, the matrix (Σ β_i^(r+s)) is a real positive semidefinite Gram matrix. In the quartic packet, however, along a real smooth transverse arc with real coefficients, the leading determinant is −T₃²<0. Hence such a packet cannot have all four inverse coordinate values real near that arc.

There is no contradiction for a complex Keller map. A real cubic pole packet has one real escaping value and one conjugate pair; together with the real bounded root it has two real values, not four. The bilinear trace form is hyperbolic exactly because the conjugate pair is retained without complex conjugation in the pairing.

Replacing Tr(h²) by Σ|β_i|² prevents cancellation, but it changes the map. It depends on complex conjugation and the real-analytic location of sheets, and is not the same holomorphic rational trace. An estimate on that Hermitian energy could prove boundedness, but no such estimate follows here from T₁,T₂ polynomiality or from constant Jacobian. A merely positive real Jacobian is also a different hypothesis from a constant nonzero complex Keller determinant.

## 9. Triangle-circle geometry is a finite quadrature, with a precise failure degree

After choosing a complex embedding and dividing the leading packet by one nonzero scalar, the escaping leading coefficients are the vertices ζ^j, 0≤j<m, of a regular m-gon on the unit circle. Equal trace multiplicities give exactly

    (1/m) Σ_(j=0)^(m−1) (ζ^j)^k = 0    for 1≤k<m,
    (1/m) Σ_(j=0)^(m−1) (ζ^j)^m = 1.

Uniform circle measure has zero k-th holomorphic moment at every positive integer k. Thus the polygon exactly reproduces the first m−1 circle moments and fails at the first resonant degree m. For the local quartic A4 model, the triangle reproduces degrees one and two but its cubic moment detects the pole. This is an exact source-to-target map, not an identification of a filled triangle with a disk or ellipsoid.

The same polygon has Hermitian second moment (1/m)Σ|ζ^j|²=1. Hence its bilinear second moment vanishes while its Hermitian second moment does not. Both statements are true; they answer different questions. Trace multiplicities are equal here. Replacing trace by arbitrary signed or weighted moments changes the equations and requires retaining those weights.

## 10. Connection contract and the next falsifiable lemma

Source: the actual finite inverse algebra of a polynomial Keller map, localized at an affine target divisor.

Map: h ↦ powers h,…,h^(m−1) ↦ their bilinear trace Gram matrix ↦ its Smith exponents and determinant.

Preserved predicate: maximal cancellation of the first m−1 traces forces the uniform polar exponent a on an explicitly specified trace lattice. The finite-sheet coefficient ratio recovers the bounded branch modulo t^a.

Lost data: sheet labels, the companion inverse coordinate, multiplication closure, how different target divisors meet, and gluing of the rational selector. These are not recovered from Smith exponents alone.

Concrete local quartic tests (controls for the higher-degree program):

1. Compute T₁,T₂. If either has a pole, the maximally blind cubic-packet premise does not apply.
2. If both are regular, compute D₂=T₂T₄−T₃². Every surviving nonproper divisor must satisfy v_D(D₂)=2v_D(T₃)<0. A proposed source construction with a different order is excluded immediately.
3. Retain a_D=−v_D(T₃), the finite-selector ratio −c₄/c₃, and the actual polynomial-source boundary incidence. Any claimed regular basis change preserving the original trace plane must leave the polar determinant unchanged.
4. Test whether a proposed global source constraint produces a finite multiplication-stable module containing 1 and x, or instead only a repaired two-dimensional bilinear lattice. The latter is explicitly insufficient by Section 7.

**OPEN next lemma (gluing/closure target).** In a specified higher-degree candidate family, for example even d≥6 with T₁,…,T_(d−2) polynomial, determine whether the local finite-sheet selectors and the companion coordinate can be glued with a common finite pole bound into a multiplication-stable module containing the original source coordinates. A positive proof would close that candidate family's properness gap. Neither a family-wide bound nor such a gluing is asserted here. The prior finite-jet/boundary-control packet explains why allowing an increasing pole bound at each approximation cannot replace it.

The next computation should therefore follow the selector's numerator-zero divisors and boundary intersections, rather than increase a point Taylor cutoff or replace the hyperbolic form by a positive one.

## 11. Reproduction

Run `python3 verify_trace_hyperbolic_packet.py` in this packet directory. The exact controls use Newton sums and companion-matrix traces independently, check the whole Gram determinant in finite packet families, verify the quartic pole order and Smith criteria, check finite-selector residual valuations, retain the integral-Gram hostile, and include the mixed-characteristic failure. Finite computations supplement the uniform proofs above.
