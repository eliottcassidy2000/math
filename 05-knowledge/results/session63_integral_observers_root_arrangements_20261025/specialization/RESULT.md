# Integral specialization: lattices, multiplication, jets, and positive cones

**Status:** PROVED elementary algebra below; inherited applications explicitly identified; finite controls are not proofs. No claim of literature novelty.
**Producer:** build_integral_specialization_framework, 2026-10-25.
**Initial repository source:** session63-repo at eaa81fd80089.
**Scope:** an operation audit, not a claim that tournament cancellation, arithmetic torsors, and incidence are one object.

## 0. Inheritance and conceptual board

The nearest inherited mechanisms are Session62's exterior interpretation of clone counting, its integer word-response quotient, and THM-4531's monomial Fourier proof of twisted Mallows–Sloane. The corrected near miss is “the averaging denominator is divisible by p, so modular invariants grow”: this implication is false even in the principal example. The least-used sidecars are the integral form, stabilizer sign, multiplication constants, and full integer Taylor jet.

| Carrier | Operation | Missing sidecar | Hostile |
|---|---|---|---|
| Integral group lattice | Invariants and coefficient change | H¹, actual lattice | C3 augmentation ideal |
| Exterior invariant algebra | Multiply invariant classes | Integral structure constants | xy=3z |
| Dimension polynomial | Reduce coefficients, test t=-1 | All Hasse jets | Depths 5,21,25 |
| Integer constraint matrix | Change coefficients | Smith factors including exponents | [2] versus [4] over Z/4 |
| Finite torsor Fourier basis | Replace point basis by characters | Lattice index | Walsh rank1 mod2 |
| Variety/cover | Count, take invariants, average forms | Carrier and extension-field quantifier | Fresh twists versus fixed family |
| Polynomial class | Evaluate q=1, seek torification | All q-1 coefficients and positive cone | 7-4(q-1)+(q-1)² |

## 1. Exact invariant base-change obstruction

**Theorem 1 (standard, directly proved).** Let A be finite and M a finite-rank free abelian group with A-action. For every prime p there is a natural exact sequence

    0 → M^A ⊗ F_p → (M ⊗ F_p)^A → H¹(A,M)[p] → 0.

Here [p] means the subgroup killed by p, not the full p-primary component. Thus the dimension increase is exactly dim_Fp H¹(A,M)[p].

**Proof.** An invariant residue class represented by x∈M yields the cocycle c_x(g)=(gx-x)/p. Changing x adds a coboundary. If c_x(g)=gy-y, then x-py is an invariant lift. Conversely p[c]=0 means pc(g)=gx-x for some x, which produces an invariant residue class with connecting class [c]. If invariant x equals py in M, torsion-freeness implies y is invariant, proving injectivity on the left. These observations prove every exactness assertion. □

This is also the initial part of the coefficient long exact sequence for 0→M --p→ M→M/pM→0. Kenneth Brown, *Lectures on the Cohomology of Groups*, §3.3, pp.16–20, gives that sequence and Proposition3.3 (Eckmann–Shapiro). The specific assertion is proved above independently of semisimplicity. [Primary author-hosted source](https://pi.math.cornell.edu/~kbrown/papers/cohomology_hangzhou.pdf).

**Hostile.** Let C3 act on Z² by g=[[0,-1],[1,-1]]. Its integral invariants vanish, but mod3 its invariant space is one-dimensional. Since det(g-I)=3 and 1+g+g²=0, H¹(C3,Z²)=Z/3. This augmentation-ideal example excludes a false extension of the next theorem to all integral representations.

## 2. Complete classification for signed permutation lattices

A signed permutation lattice has basis {e_s:s∈S} with g e_s=α(g,s)e_(gs), α(g,s)∈{±1}, satisfying the action cocycle law. For an orbit representative s let H_s=Stab_A(s), and let ε_s(h)=α(h,s). This is a homomorphism H_s→{±1}. Call an orbit oriented if ε_s is trivial. Let a,b count oriented and nonoriented orbits.

**Theorem 2.**

1. M ≅ ⊕_[s] Ind_(H_s)^A Z_(ε_s) as integral A-lattices.
2. For every commutative unital ring R,

       (M⊗R)^A ≅ R^a ⊕ R[2]^b,  where R[2]={r∈R:2r=0}.

   The displayed coordinate isomorphism requires orbit representative/sign choices.
3. H¹(A,M) ≅ (Z/2)^b.
4. The canonical map M^A⊗R→(M⊗R)^A is injective with cokernel R[2]^b, up to those coordinate identifications. It is an isomorphism for every R exactly when b=0.

**Proof.** An orbit of basis lines is the induced rank-one module in1. An invariant vector on that orbit is determined by one coefficient r. Consistency around stabilizer loops says exactly r=ε_s(h)r for every h. This imposes no equation if ε_s is trivial and precisely 2r=0 otherwise. For an oriented orbit, transporting e_s gives a well-defined basis indexed by cosets, genuinely permuted by A. Its integral invariant generator is the sum of those transported vectors, with coefficients ±1 in the original basis. This proves1,2,4.

For3, an independent elementary argument avoids needing Shapiro: if c is an integral cocycle with c(gh)=c(g)+g c(h), put S=Σ_g c(g). Then hS=S-|A|c(h), so v=-S/|A| satisfies c(h)=hv-v. Hence

    H¹(A,M) ≅ (M_Q/M)^A / image(M_Q^A).

The coefficient-consistency calculation above works for the abelian group Q/Z. Each oriented orbit contributes Q/Z, entirely filled by its rational invariant line. Each nonoriented orbit contributes (Q/Z)[2]=Z/2 and has no rational invariant line. The quotient is (Z/2)^b. □

Thus signed permutation invariants commute with every odd-characteristic field, even if that characteristic divides |A|. If |A| is odd, every stabilizer sign is trivial, and they commute with **every** commutative coefficient ring. This is module-specific; it does not make the invariants functor exact on all A-modules.

## 3. Exterior powers and the 50-vertex cancellation example

For A permuting X, the r-subsets S index the wedge basis of Λ^r Z[X]. Their stabilizer sign is

    ε_S(h)=sign of the permutation h restricted to the elements of S.

Theorem2 applies degreewise. If A has odd order, all exterior powers become genuine permutation lattices after orbitwise sign choices. Consequently the canonical map

    (Λ* Z[X])^A ⊗ R → (Λ* R[X])^A

is an isomorphism of graded algebras for every commutative ring R. It is an algebra map because it respects wedge products, and bijective in every degree. This does not say that orbit sums freely generate an exterior algebra.

Every finite tournament automorphism group has odd order: a nontrivial involution swaps two vertices and reverses their edge, while an even-order group has an involution. Therefore the theorem applies to its induced action on an invariant gap set.

For regular A=C5² on25 gaps, the integer polynomial is

    P(t)=((1+t)^25+24(1+t^5)^5)/25.

The cancellation producer owns the intrinsic50-vertex realization and general cancellation law. Its depths are5 overQ,21 after coefficient reduction mod5, and25 mod2. Yet **all degreewise invariant spaces commute with reduction**, and all their H¹ groups vanish. The usual integer-valued Hilbert polynomial, recording actual dimensions, is identical over every field. The depths21 and25 refer to the further operation of mapping that integer polynomial's coefficients into F5 or F2.

Inherited source: 05-knowledge/results/session62_agent_theme_bridges_20261025/SYNTHESIS.md and its clone/ and audits/holonomy-clone/ packets. No tournament realization is proved by the finite controls here.

## 4. A central small hostile: no lifting defect, but xy=3z

Let C3 cycle e1,e2,e3. Its exterior invariants have integral bases

    degree0: 1
    degree1: x=e1+e2+e3
    degree2: y=e1∧e2+e2∧e3+e3∧e1
    degree3: z=e1∧e2∧e3.

The action cycles the indicated basis terms; these really are bases of the entire invariant groups. In x∧y exactly three terms survive, namely e1∧e2∧e3, e2∧e3∧e1, and e3∧e1∧e2. All equal z. Therefore

    x∧y=3z.

Apart from x∧y=y∧x=3z, all products of positive-degree basis elements vanish, either by x²=0 or by dimension. Four layers are now visibly distinct:

* Invariants commute with every coefficient ring; there is no degreewise H¹ defect.
* Multiplication from degree1⊗degree2 into degree3 has image3Zz; it becomes zero mod3.
* P(t)=1+t+t²+t³ has root depth1 at−1 overQ andF3, but depth3 overF2.
* OverQ the algebra has the factorization Λ_Q(x)⊗Q[y]/(y²), deg y=2, with z=xy/3. The analogous integral free factorization fails: changing primitive bases only changes the product3 by a sign, never to1.

An optional derived control is the invariant complex with differential x∧:

    Z --1→ Z --0→ Z --3→ Z.

It has only top cohomology Z/3. Thus torsion can arise after an additional operation even though every term's invariant base change is perfect. This is a scoped example, not an unproved general theorem about invariant Koszul complexes.

The root-cover multiplication law in the torsor lane is another concrete reason to retain algebra structure rather than dimensions alone; no identification of its algebra with this exterior algebra is claimed.

## 5. THM-4531: same induced-sign mechanism, different signs and integral forms

Precise sources:

* 01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md;
* 05-knowledge/results/procgen_tbij_20261001_natural_bijections.md, §1;
* the updated obstruction in 05-knowledge/results/uniform_equivariance_obstruction_20261020/UNIFORM_EQUIVARIANCE_OBSTRUCTION.md.

Let E=E(K_n), U≤F2^E invariant under S_n, W=F2^E/U, and T_U the W-torsor of tournaments modulo U-reversals. The Fourier functions indexed by X∈U^⊥ form a signed monomial lattice

    L=⊕_(X∈U^⊥) Z φ_X ⊂ P=Z^(T_U).

Its stabilizer sign is the parity of **endpoint reversals** on the edges of X:

    ε_X(g)=(-1)^(#{ {u<v}∈E(X):g(u)>g(v) }).

It is neither the vertex-permutation sign nor the sign of the permutation of selected edges. Swapping the endpoints of one fixed edge has endpoint sign−1 and edge-permutation sign+1; the empty graph has endpoint sign+1 even under an odd vertex permutation. The common mechanism with§3 is induction from a stabilizer sign line, not equality of those signs.

Royle–Praeger–Glasby–Freedman–Devillers define this endpoint sign in§4, Lemma4.1, and prove the tournament/even-graph count. Their broader question of a natural unlabelled bijection is distinct from the repository's precise equivariant/type-bijective formulation. [Primary paper](https://arxiv.org/pdf/2204.01947).

**Corollary5A.** If b counts non-even graph orbits in U^⊥, then H¹(S_n,L)=(Z/2)^b. The mod2 invariant dimension of L increases by b: it counts all graph orbits, not only the sign-even graph orbits.

**Corollary5B (integral Fourier sidecar).** Put d=dim_F2 W and N=2^d. The character matrix F satisfies FᵀF=N I, so

    [P:L]=2^(d·2^(d-1))  for d≥1.

L→P becomes an isomorphism after inverting2, but its mod2 reduction has rank1. For d=0 it is the identity. Proof: character orthogonality is the elementary sum Σ_w(-1)^(χ(w)+ψ(w))=N[χ=ψ]; take determinants, then reduce the ±1 matrix mod2.

The honest point-permutation lattice P has H¹=0. Extra mod2 invariants of L therefore do not create new tournament-quotient orbits. They expose the bad reduction of a rationally valid change of basis. The smallest example, n=2,U=0, is

    Z_trivial⊕Z_sign → Z²,  (a,b)↦(a+b,a-b),

with index2. Source and target have mod2 invariant dimensions2 and1 respectively.

This complements but does not replace the repository's stabilizer-containment/Hall proof against natural maps. An equivariant vector-space isomorphism is not an equivariant set bijection.

## 6. Exact integer jet filtration

**Theorem6.** For P(t)=Σ_i a_i t^i∈Z[t], define

    b_j=Σ_(i≥j)(-1)^(i-j) binom(i,j)a_i,
    P(u-1)=Σ_j b_j u^j,
    Δ_r=gcd(b_0,…,b_(r-1)), with gcd(0,…,0)=0.

The order at t=-1 overQ is the least j with b_j≠0; overF_p it is the least j with p∤b_j, with infinity for the zero polynomial. In particular

    depth_modp≥r iff p divides Δ_r,
    Δ_(r+1) divides Δ_r.

More generally u^r divides P(u-1) modulo p^e iff v_p(b_j)≥e for every j<r, or equivalently p^e divides Δ_r. These statements use the convention every integer divides0.

**Proof.** Expand (u-1)^i. Translation is an integral triangular coefficient-basis change with diagonal1, hence commutes with every coefficient reduction. Divisibility by u^r is exactly vanishing of the first r coefficients. Adding a generator to an integer ideal gives the gcd divisibility chain. □

This scalar theorem needs no group, invariant functor, or torsion module. Ordinary derivatives divided by j! are not a valid modular replacement when p≤j; these integral Hasse jets are.

For the25-gap polynomial, independent Python and BigInt-JavaScript controls give

    b_0,…,b_25 =
    0,0,0,0,0,3000,-30000,150000,-495000,1203000,
    -2280000,3486000,-4392000,4621200,-4093200,
    3063600,-1938000,1032840,-460560,169920,-51000,
    12144,-2208,288,-24,1.

The exact filtration is

| r | Δ_r |
|---|---:|
| 1–5 | 0 |
| 6–13 | 3000=2³·3·5³ |
| 14–17 | 600=2³·3·5² |
| 18–21 | 120=2³·3·5 |
| 22–25 | 24=2³·3 |
| 26 | 1 |

Thus prime2 and3 give depth25; prime5 gives21; all other primes give5. The filtration retains strictly more p-adic precision than the field depths: modulo5,25,125,625 the u-divisibility depths are21,17,13,5 respectively.

For any odd-order action on an n-element set, every subset orbit has odd size. Consequently the number of r-subset orbits is binom(n,r) modulo2, and P(t)≡(1+t)^n mod2. This explains maximal mod2 scalar cancellation without adding invariant vectors.

## 7. Genuine constraint defects: Smith data and observer dependence

For an integer matrix D:Z^n→Z^m with nonzero Smith entries s_1,…,s_r, its kernel over any commutative ring R is

    R^(n-r) ⊕ ⊕_i R[s_i].

Here annihilators are actual solution vectors, unlike scalar jet vanishing. OverF_p only divisibility p|s_i is visible; overZ/p^e the exponents matter. Matrices[2] and[4] have identical ranks over every field but kernels of sizes2 and4 overZ/4.

Session62's named positive-word matrix span has quotient (Z/(q-1))^q; adding the identity leaves (Z/(q-1))^(q-1). Its exceptional coefficient primes are genuine cokernel-torsion primes. Yet the whole-word scalar series has a (q+1)-square Hankel minor [[1,1ᵀ],[1,J-I]] of determinant(-1)^q, so its q+1 dimension survives every characteristic. The observer bank and treatment of the empty word changed: there is no contradiction.

Source: 05-knowledge/results/session62_agent_theme_bridges_20261025/audits/palettes/AUDIT.md, §3 and AppendicesA–B. The incoming integral-wrap theorem is prior work owned by the recovery lane, not rediscovered here. Its R[D] terms belong in this Smith column.

**Parameter guard:** palette q is an integer number-of-colors/dimension parameter. Substituting q=1 into a formula is not automatically specialization of a fixed-size matrix family. Such a family must be separately defined. The geometric counting variable q below has a different type.

## 8. q=1, Euler characteristic, and the torification cone

The twist producer identifies U as P² minus the seven lines

    E, E±x, E±z, E±(x+z).

In characteristic0, or finite characteristics other than2,3, there are three triple points and twelve double points. Inclusion–exclusion gives the actual scissor relation

    [U]=[P²]-7[P¹]+(12+2·3)[point]=L²-6L+12.

Thus #U(F_q)=q²-6q+12 in the stated characteristics, E_c(U;u,v)=(uv)²-6uv+12, and χ_c(U)=7. This is a geometric class calculation, not an inference from point counts. The Hodge realization uses only additivity, multiplicativity, and E_c(A¹)=uv; see Cappell–Maxim–Shaneson, *Hodge genera of algebraic varieties,I*, discussion of the E-polynomial. [Primary research source](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/calimash2.pdf).

Writing z=q-1 gives z²-4z+7. A finite stratification by split tori would force its count/Hodge polynomial to be Σ_d a_d z^d with a_d nonnegative integers. The coefficient−4 rules this out. It is an obstruction to the torus-stratum positive cone, not to positivity of point counts.

The general necessary criterion is López Peña–Lorscheid, *Torified varieties and their geometries over F1*, Proposition1.11 and Remark1.21; here the application is coefficient comparison. [Primary source](https://arxiv.org/pdf/0903.2173).

P_U(1)=χ_c(U)=7 is justified by that class calculation. It does not construct a field of size1 or an F1-model. Evaluation at1 forgets precisely the higher jets, including the negative linear coefficient.

## 9. Arithmetic twist averages are not averages of ordinary motives

The degree64 root-open cover V→U is a torsor under G=(μ2)^7/μ2. Over each finite field of odd characteristic, averaging the **fresh64 twists over that field** gives #U(F_q). The twist producer owns the exact torsor/Fourier count, branch corrections, and characteristic3 degeneration.

Fix instead one64-twist family overF_q and extend it toF_(q^m). For even m all its square classes become squares, so its twists become the same V; the average is #V(F_(q^m)), not generally #U(F_(q^m)). For odd m the square-class labels remain distinct. A fresh-field average is not one fixed polynomial-count variety across all extensions.

OverC all64 twists are isomorphic toV. Their ordinary Grothendieck-class average is[V], not[U]. The finite unramified cover has χ_c(V)=64χ_c(U)=448, whereas χ_c(U)=7. Euler characteristic already refutes (1/64)Σ[V_ε]=[U] in K0(Var_C)⊗Q.

The correct characteristic-zero equivariant operation is the projector e_G=(1/64)Σ_g g. Finite-cover transfer identifies its image in cohomology with cohomology ofU. Averaging group-operator traces is not averaging underlying complex varieties. The projector belongs toQ[G] orZ[1/2][G], not generallyZ[G]; no integral or mod2 comparison is supplied by the point-count mean.

The averaged polynomial alone also does not prove each twist open non-torifiable. The inherited individual-cover argument uses its separate nonrationality/geometric-genus obstruction. No unqualified “orbifold Euler” or universal motivic average is asserted.

## 10. Reusable operation-complete integral observers

The framework is a checklist with proved interfaces:

1. Specify the integral carrier and coefficient ring.
2. Retain its constraint/invariant map before decategorification.
3. Retain lattice index, Smith data, and H¹ where applicable.
4. If an algebraic operation will follow, retain multiplication constants.
5. Separately keep the integer polynomial and full translated jet.
6. Distinguish divisibility from positive-cone membership.
7. For averages, retain equivariance and the extension-field quantifier.

The concrete predictions established above are: the50-vertex cancellation jumps have no invariant-lifting defect; THM-4531's Fourier and point lattices differ integrally at2; C3's multiplication degenerates at3 while its scalar root depth jumps at2; and the number7 atq=1 sees neither the negative torification coefficient−4 nor the individual-cover Euler characteristic448.

These compare different operations with explicit preservation/failure data. They imply no new LRC, coloring, or natural-bijection theorem.

## 11. Verification, sources, and task boundary

Run python3 check.py and python3 -O check.py. The independent checker imports neither producer and uses no bare assertions. It compares orbit-stabilizer orientation counts with independent Gaussian-rank constraints overQ,F2,F3,F5 for cyclic actions, symmetric actions, C3², and graph endpoint-sign actions. It also checks Walsh orthogonality/ranks, the augmentation-ideal defect, all26 scalar jets, and the Smith-exponent hostile. Outputs declare the finite universe; general theorems have proofs above.

The separate independent audit agent confirmed the H¹ proof, arbitrary-ring canonical map claims, C3 multiplication/Koszul example, and the exact gcd filtration. Primary sources are used only for their named standard facts; no novelty assertion is based on unsuccessful literature search.

No repository edits or pushes were made by this lane. Startup ran; targeted inbox inspection was blocked by absent .machine-id, and no identity was invented. An execution-service reset interrupted the first prose write; filesystem readback confirmed absence before a bounded create-only retry. The parent owns final integration, audit, and publication.
