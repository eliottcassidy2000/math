# Seven branch lines, simplex subgroups, and a 16 / 12 / 2 surface census

**Status:** PROVED elementary matroid/code/combinatorial statements; PROVED complex-surface deductions using the explicitly named standard surface facts below; FINITE-EXACT independent enumerations. This packet claims no historical novelty, rational square point, or integer magic-square existence theorem. Counts are of **labelled deck subgroups / covers over a fixed plane**, not of abstract surface isomorphism classes.

## 0. Inheritance and live board

Closest inherited mechanisms: Session58 `hourglasses/THEOREM.md` and `MITRE_AFFINE_BRIDGE.md` prove the seven-cell normal form, mitre incidence, non-Fano obstruction and characteristic-3 affine completion. Session59 `graph-magic-degeneration/branch-addendum.md` proves the nine-row arithmetic-matroid degeneration; `distinctness-addendum.md` separates branch avoidance from value injectivity. The current root-cover lane supplies the seven-root four-quadric model and the six-root K3 threshold; the elementary deck presentation is reproved here.

Canonical hostile: the two outer mitre blocks are disjoint, while any two Fano lines meet. Corrected near miss: seven labels alone do not turn either the branch matroid or seven canonical characters into the Fano plane. Least-used sidecar: a **chosen rank-three quotient of the rank-six deck group**, including its local inertia labels.

| Live concept | Representation / operation | Decisive missing coordinate |
|---|---|---|
| Seven square-value forms | Integral rank-three rows; reduce characteristic | Integer minors, not just seven labels |
| Mitre / AG(2,3) | Delete two affine points | Equal-sum incidence is not characteristic-zero collinearity |
| Projective root deck | Seven generators with one relation | Rank six, not three |
| Hamming / simplex | Orthogonal binary codes; choose a subgroup | Which of the 30 labelled simplex subgroups? |
| Quotient singularities | Local inertia rank at each crossing | Independent versus dependent three labels |
| Distinguished central line | Residual monodromy / resolution | Component incidence, not only K-squared |

Anchor: exact six/seven deletion geometry. Niche: code quotients and fixed loci. Wildcard: the central boundary identifies three different birational types.

## 1. Three genuinely different seven-element objects

Index the forms by

    L0=E, L1=E+x, L2=E-x, L3=E+z, L4=E-z,
    L5=E+x+z, L6=E-x-z.

Their integral normal rows are W={(1,0,0),(1,±1,0),(1,0,±1),(1,1,1),(1,-1,-1)}. Let M_k be their vector matroid over a field k. The 35 absolute triple determinants have histogram

    0:3, 1:18, 2:12, 3:2.

The zero triples are

    T1=012, T2=034, T3=056.                         (1)

The determinant-3 triples are 136 and 245. This table follows by direct 3-by-3 determinants; `check.py` records every determinant and independently verifies every subset rank by Gaussian elimination and minors.

**Theorem 1.** In characteristic zero or at least five, the arrangement has precisely three ordinary triple points (1), twelve double points, and no higher crossing. In characteristic three it has precisely the five triples (1),136,245, six double points, and no higher crossing. In characteristic two its normal rows collapse to four distinct vectors of multiplicities 1,2,2,2; the simplification is U_(3,4).

Proof: the displayed determinant table gives every dependent triple. Every pair is independent outside characteristic two. No fourth label can extend any displayed rank-two triple, as the list contains no two triples sharing two labels. Each triple accounts for three of the 21 pairs, leaving the stated double points. The characteristic-two reduced vectors are (1,0,0),(1,1,0),(1,0,1),(1,1,1), any three independent. This also proves that no other characteristic changes occur. ∎

In characteristic three the seven normal points are exactly the affine-plane restriction AG(2,3) with (1,1,-1),(1,-1,1) deleted. Its five rank-two triples are the classical mitre. In characteristic zero, the mitre's two outer blocks are equal-sum triples, **not** rank-two triples: their vector sums equal 3(1,0,0), becoming zero only in characteristic three. This is the exact bridge, not an identification valid in every characteristic.

The five-block mitre is not Fano: its two outer blocks are disjoint. Nor is the three-spoke generic matroid Fano: it has three rather than seven rank-two triples. Nor does reduction modulo two produce Fano: it produces repeated columns. These are distinct failures of three different attempted maps.

For reference, inclusion-exclusion gives the central arrangement characteristic polynomials

    char 0 or ≥5: (t-1)(t-3)(t-4),
    char 3:      (t-1)(t^2-6t+10),
    char 2:      (t-1)(t^2-3t+3).

These count the affine branch complement over corresponding finite fields; they do not count square-root points or distinct values.

### All seven six-line deletions

| Characteristic | Deleted label | Remaining triple points | Characteristic polynomial |
|---|---|---:|---|
| 0 or ≥5 | 0 | 0 | t^3-6t^2+15t-10 |
| 0 or ≥5 | any outer label | 2 | t^3-6t^2+13t-8 |
| 3 | 0 | 2 | t^3-6t^2+13t-8 |
| 3 | any outer label | 3 | t^3-6t^2+12t-7 |

Deletion of 0 is U_(3,6) in characteristic zero. The other six deletions have two triple points. Thus seven six-root models need not be equally smooth: using the local root-cover calculation, the central deletion has no nodes, whereas every outer deletion has 2·2^(5-3)=8 nodes in its maximal degree-32 root cover. In characteristic three these counts become 8 and 12, respectively. The surface/K3 conclusion requires characteristic different from two and the standard complete-intersection/resolution hypotheses established in the sibling root-cover packet.

## 2. The exact binary bridge: thirty extra choices

Work now over C for covers, and over F2 for deck labels. Let X be the normal projective seven-root surface, with coordinates r_i and [r_i^2]=[L_i]. Its full projective sign deck group is

    G=F2^7 / <1>,       1=(1,1,1,1,1,1,1).

The inertia element g_i is the class of e_i; their only relation is Σg_i=0. The dual group is the even-weight code E7={v:Σv_i=0}, of dimension six. The seven canonical characters inherited from the root-cover calculation are c_i=1+e_i. Their only nonempty binary dependence is the full seven-term sum. Hence their column matroid is U_(6,7), decisively different from rank-three Fano.

Because seven is odd, there is a canonical S7-equivariant identification

    G → E7,    [v] ↦ v+(Σv_i)1.                    (2)

It sends g_i to c_i. Dot product restricts to a nondegenerate alternating form on E7: v·v=0 for even v, and its radical would be <1>∩E7=0.

Choose a bijection a:{0,...,6}→F2^3\{0}. The sum of its seven values is zero, so e_i↦a_i descends to a surjection q_a:G→F2^3. Let H_a be its 3-by-7 column matrix. Then

    C_a=ker H_a                     [7,4,3] Hamming code,
    S_a=rowspan H_a                  [7,3,4] simplex code,
    ker q_a=C_a/<1> ≅ C_a∩E7=S_a.    (3)

Here the last isomorphism is (2), not a conflation of a vector space with its dual. Proof: every nonzero row functional takes value one at exactly four of the seven nonzero vectors. Consequently S_a has weights 0 and 4, and H_a H_a^T=0. Thus S_a⊂C_a∩E7; both have dimension three. The Hamming code's weight-3 words are precisely Fano lines, its weight-4 words their complements, and its remaining words are 0 and 1. This proves all asserted dimensions and weights directly.

The standard Fano/simplex and dual-Fano/Hamming terminology is classical; it is explicitly recorded in Freij-Hollanti–Grezet–Hollanti–Westerbäck, *Cyclic flats of binary matroids*, published version §8.2, p34. [Primary published manuscript](https://acris.aalto.fi/ws/portalfiles/portal/62863545/Freij_Hollanti_Cyclic_Flats.1_s2.0_S0196885821000038_main_1.pdf).

**Theorem 2.** There are exactly 7!/|GL(3,2)|=5040/168=30 such kernels K=ker q_a. Equivalently they are the 30 three-dimensional binary subspaces whose seven nonzero vectors all have weight four. They are Lagrangian subspaces of E7. A kernel determines the complete labelled Fano line system, but no kernel is distinguished by the mere existence of the seven branch labels.

Proof: two surjections have the same kernel iff related by a target linear automorphism; the GL action on bijections is free. Since |GL(3,2)|=7·6·4, the count follows. Conversely, a constant-weight-four 3-space S has a generator matrix whose columns give a bijection onto F2^3\{0}: for each nonzero functional on the three-dimensional row space its weight is four, and binary Fourier inversion forces each nonzero column to occur once and the zero column zero times. Self-orthogonality follows from wt(u+v)=wt(u)+wt(v)-2|supp(u)∩supp(v)|: two distinct nonzero words meet in two positions. Thus S is Lagrangian, and (3) recovers the quotient. ∎

The quotient does not preserve the full deck, root signs, or canonical character space. In particular S_a contains none of the weight-six c_i, so none of the seven root-deletion involutions g_i lies in K. The rank-three quotient therefore does not factor through any of the seven maps obtained by forgetting one root. Likewise its pullback character subspace has no canonical weight-six character, explaining why its geometric genus vanishes.

## 3. Why exactly 16 / 12 / 2

Call an actual crossing T_j **dependent** for a if Σ_(i∈T_j)a_i=0. Let d count such crossings.

**Theorem 3.** Among the 30 kernels, d has distribution

    d=0:16,       d=1:12,       d=3:2;             (4)

there are none with d=2.

Proof without enumeration: the three Fano lines through label 0 define a perfect matching on the other six labels. The actual triples (1) define another perfect matching M. Two perfect matchings on six vertices share 0,1, or 3 edges, never exactly two. Of the 15 matchings, one equals M; six share exactly one edge (choose the common edge in three ways and one of the two alternate matchings on the remaining four); the remaining eight share none. Each matching extends to exactly two labelled Fano planes. To see the last assertion, write the three pairs as (b,b+a0),(c,c+a0),(b+c,b+c+a0). The four off-center lines are the transversal triples of one parity; interchanging the two elements in one pair exchanges the two possible parities. Therefore the counts double to 16,12,2. ∎

As a separate incidence-only check, the characteristic-three mitre has intersection-size distribution 0:12,1:10,2:6,4:2 against the 30 Fano planes. This finite census does **not** transfer the complex-surface classification below into characteristic three.

### Fixed loci and the local quotient

A nonidentity k∈K is a sign word with support size four. Its fixed locus in P6 is the union of the coordinate P3 and complementary P2. A point of X in the P2 would require four branch forms to vanish, impossible. A point in the P3 requires the complementary three forms to vanish, possible precisely when their labels form an actual crossing. Its complement is a Fano line by (3). Thus K acts freely on X iff d=0; then it also acts freely on the minimal resolution.

At an independent crossing, choose local branch parameters u,v,u+v. The maximal local root cover is

    s^2=u, t^2=v, w^2=u+v,

an A1 singularity s^2+t^2-w^2=0. In the rank-three quotient it remains one A1 point. At a dependent crossing the local inertia rank falls to two; there are two points above the crossing. Each is the quotient of A1 by simultaneous negation of s,t,w. Writing A1=C2/{±1}, this negation lifts to (a,b)↦(ia,ib), so the quotient is 1/4(1,1). Its minimal resolution has one rational (-4)-curve. Double crossings are smooth because their distinct nonzero binary labels are independent. Therefore

    Sing(Y_a)=(3-d)A1 + 2d·[1/4(1,1)].            (5)

The rank-three cover construction and the Campedelli criterion are classical: Kharlamov–Kulikov, *Surfaces with DIF≠DEF real structures*, §1.3 Proposition1.3, §1.5 definition and Theorem1.10, treat the universal rank-six seven-line cover and its rank-three quotients; Lemma1.7 distinguishes the two local triple cases. Their theorem makes exactly the d=0 labelings Campedelli arrangements. [Primary manuscript](https://webdoc.sub.gwdg.de/ebook/serien/e/mpi_mathematik/2005/76.pdf). Our local calculation and arrangement-specific census are given explicitly rather than inferred from terminology.

### Realized D12 orbit refinement

The independent symmetry packet proves that the full projective automorphism group of this characteristic-zero arrangement is D12. On the present labels its generators act by the image lists

    R=(0,5,6,2,1,3,4),       S=(0,3,4,1,2,5,6),
    R^6=S^2=1,              SRS=R^(-1).

Permuting sign coordinates acts on the thirty kernels. The five orbits are:

| d | Orbit sizes | Stabilizer orders |
|---|---|---|
| 0 | 4, 12 | 3, 1 |
| 1 | 6, 6 | 2, 2 |
| 3 | 2 | 6 |

The certificate lists every orbit and stabilizer permutation. These realized symmetries give isomorphisms of quotient covers after moving the base. For d=0, §4 proves the converse using the intrinsic bicanonical map, giving exactly two abstract surface classes. For d=1 no such converse is asserted. For d=3 the single orbit identifies the two Enriques minimal models with each other.

## 4. From the distinguished line to three birational types

This section is over **C**. Set f:Y_a→P2 and π:S→Y_a its minimal resolution. The binary simplex code continues to be denoted S_a; the smooth surface is S.

All seven nontrivial characters of F2^3 branch over four lines. The finite-cover eigensheaf decomposition is consequently

    f_*O_Y = O_P2 ⊕ O_P2(-2)^7.

One can obtain it by the squarefree character monomials: a weight-four branch product has a square root transforming as a section of O(2). The cover is normal and finite flat (normal surfaces are Cohen–Macaulay; the base is regular of dimension two), and the eigensheaves are rank-one reflexive, hence invertible, so this codimension-one calculation determines them globally. Thus h1(O_Y)=h2(O_Y)=0. A1 and cyclic quotient singularities are rational, so

    q(S)=p_g(S)=0,       χ(O_S)=1.                 (6)

Hurwitz in codimension one gives 2K_Y∼f*H, with K_Y ample and K_Y^2=2. An A1 resolution is crepant; a (-4)-exception has discrepancy -1/2, by adjunction. Hence

    K_S^2=2-2d.                                  (7)

Equation (7) alone does **not** determine the birational type.

### Residual monodromy and rational central components

Let a=a0. On the normalization of R=f^*(L0)/2, quotient out the generic inertia <a>. Its three local residual monodromies are

    m_j = a_(2j-1)+a_(2j) mod <a>,  j=1,2,3.

They sum to zero, and exactly d of them are zero. Thus their span in F2^3/<a> has rank r=2,1,0 for d=0,1,3. The normalization of R has 4/2^r=1,2,4 components, each of degree 2^r over L0. For r=2 the cover is a Klein-four cover of P1 branched at three points; for r=1 a double cover branched at two points; for r=0 degree one. Riemann–Hurwitz gives genus zero in each case.

The incidence can be recovered without pictures. Put

    H=<a, a1+a2, a3+a4, a5+a6>⊂F2^3.

Central components are indexed by G'/H, where G'=F2^3. At a crossing (0,b,c), singular points are indexed by G'/<a,a_b,a_c>. The two branches of R at the point represented by g lie on components [g] and [g+a_b] in G'/H. This gives:

* d=0: one rational component, with two branches at each of the three A1 points.
* d=1: two rational components; each of the two A1 points joins the two components; each of the two 1/4 points has both branches on one component, one such point per component.
* d=3: four rational components; the six 1/4 points join all six distinct component pairs, the edge incidence of K4.

These statements follow also directly from perfect matchings: when d=1, the nonzero residual monodromies equal the class of the uniquely matched pair. `check.py` checks the coset incidence for all 30 quotients.

Let R_i be the smooth strict transforms of the rational central components. Let A_j denote the (-2)-exceptional curves and B_k the (-4)-exceptional curves. At each singular point the two strict branches meet the exceptional curve transversely at distinct points. Vanishing intersection of pullback divisors with exceptional curves gives

    π*R = ΣR_i + ΣA_j + (1/2)ΣB_k.

Since 2K_Y∼2R and K_S=π*K_Y-(1/2)ΣB_k, we obtain the crucial integral linear equivalence

    2K_S ∼ 2(ΣR_i + ΣA_j).                       (8)

Adjunction now determines every central self-intersection. Alternatively compute K_S·R_i as (degree(R_i→L0))/2 minus half the number of its branches at bad singularities. The resulting table is

| d | Number of R_i | Degree each | R_i^2 | Exceptional incidence |
|---|---:|---:|---:|---|
| 0 | 1 | 4 | -4 | three A_j meet R_1 twice each |
| 1 | 2 | 2 | -2 | two A_j connect R_1,R_2; each R_i meets one B_k twice |
| 3 | 4 | 1 | -1 | six B_k connect all pairs, subdivided K4 |

### Independent blowup-base proof and an explicit elliptic pencil

Blow up P2 at the three triple points, obtaining B with exceptional curves E_1,E_2,E_3 and line class H. The strict branch lines have total class 7H-3ΣE_j. Exceptional E_j is additionally branched precisely when its three binary labels sum nonzero. The resulting branch divisor has simple normal crossings; at every crossing its two inertia labels are distinct nonzero vectors, so the normalized rank-three cover g:S→B is smooth. Locally it is exactly the A1 or (-4)-resolution above. Hurwitz gives the stronger exact identity

    2K_S ∼ g*(H-Σ_(dependent j) E_j).             (9)

For d=1, the divisor H-E_bad is the free pencil of lines through the dependent crossing. A general member meets just the four complementary branch lines. Their four labels are the complement of a Fano line and span F2^3, so its pullback is connected. The degree-eight cover of P1 is branched at four points with inertia order two, hence has genus one by 2g-2=-16+4·4=0. Consequently the composite S→B→P1 is the genus-one fibration and (9) says

    O_S(2K_S) ≅ h*O_P1(1),       P_(2m)(S)=m+1 (m≥0).  (10)

The second equality uses connected fibers and the projection formula. Thus κ=1 directly, with no abundance assumption. The central pencil member pulls back to twice the four-cycle described below; no complete singular-fiber or multiple-fiber classification is asserted. For d=3, H-ΣE_j is the strict transform of the central branch line, and (9) immediately becomes 2K_S∼2ΣR_i. This independently checks the cancellation in (8).

This blowup argument was recovered independently by the audit lane and is included to avoid relying on a singular canonical-divisor computation alone.

### Classification theorem

1. **Sixteen Campedelli quotients (d=0).** The resolution has p_g=q=0, K^2=2 and nef canonical divisor, since it is a crepant resolution of ample K_Y. It is a minimal surface of general type, hence numerical Campedelli; the cited classical seven-line theorem supplies the usual Campedelli identification. Its normal canonical model has three A1 points.

2. **Twelve properly elliptic quotients (d=1).** Put F=R_1+R_2+A_1+A_2. These four (-2)-curves form a transverse four-cycle. Hence F·C=0 for each component C of F, and F·C≥0 for every other irreducible curve. Thus F is nef and numerically nonzero. By (8), K_S is nef, numerically nonzero and has square zero. In particular S has no (-1)-curve. The explicit pencil and even plurigenera (10) prove κ=1. Equivalently, the complex surface classification rules out κ=0 and κ=2 from these same numerical and divisor data. The displayed four-cycle is the reduced I4-shaped divisor of the explicitly proved double central fiber; no complete fiber list or list of all other multiple fibers is asserted.

3. **Two Enriques quotients after contraction (d=3).** The four R_i are disjoint (-1)-curves on S. Contract them to a smooth surface β:S→Z. Since K_S=β*K_Z+ΣR_i and (8) has no A_j, it follows that β*(2K_Z)∼0, hence 2K_Z∼0. Moreover p_g(Z)=q(Z)=0, so K_Z is nontrivial (otherwise p_g=1). Thus Z is an Enriques surface. It is already minimal because its canonical class is numerically trivial. Notice that S itself had K_S^2=-4, and only after four explicit contractions does K_Z^2=0.

### Exactly two Campedelli isomorphism classes, not sixteen

For d=0, the complete bicanonical system is intrinsic and recovers the base map:

    H0(S,2K_S)=H0(P2,O(1)⊗f_*O_Y)=H0(P2,O(1)),
    h0(S,2K_S)=3.

Here π_*O_S=O_Y, and 2K_S is the pullback of O(1). Thus any isomorphism between two of these minimal Campedelli surfaces descends, through their bicanonical maps, to a projective transformation of P2 preserving the seven branch lines. By the independent symmetry theorem this transformation belongs to D12. It also conjugates the degree-eight Galois groups and therefore identifies their line-label maps up to GL(3,2), equivalently their kernels. The converse is realized by root permutations. Consequently the two d=0 D12 orbits are **exactly two isomorphism classes** of these complex Campedelli surfaces.

More precisely, for a kernel K in either orbit there is an exact sequence

    1 → (C2)^3 → Aut(S) → Stab_D12(K) → 1.

All stabilizer elements lift by literal root permutations, so the sequence splits. The orbit of size12 has automorphism group (C2)^3, of order8. The orbit of size4 has stabilizer C3; its order-three permutation fixes the central label and has two three-cycles on the other six. Its nontrivial action on the quotient F2^3 splits as a trivial one-dimensional part plus the irreducible two-dimensional part. Therefore this automorphism group is C2×A4, of order24. This conclusion concerns the two d=0 surfaces only; the elliptic orbits need not distinguish abstract isomorphism classes, and no completeness theorem for their automorphism groups is asserted.

The classification inputs in item2 are the usual complex Enriques–Kodaira classification, explicitly summarized in Ravi Vakil, *Complex Algebraic Surfaces, Class19*, §§3–5: the p_g=q=0 minimal κ=0 case has 2K=0, minimal general-type surfaces have K^2>0, and κ=1 surfaces are elliptic. [Primary lecture notes](https://math.stanford.edu/~vakil/02-245/sclass19A.pdf). For item3, Dolgachev, *A brief introduction to Enriques surfaces*, v2, §3 p6, gives the classical q=p_g=0, 2K=0 definition and minimality. [Author paper](https://arxiv.org/pdf/1412.7744). These are CITED classification inputs; (8), the incidence and the contractions are the explicit application.

## 5. Map / loss / restoration ledger

| Source → target | Actual map and preserved predicate | Destroyed information / necessary sidecar |
|---|---|---|
| Integral seven rows → char3 mitre | Reduce coefficients; determinant-3 triples become dependent | Characteristic-zero rank; retain integer minors |
| Mitre → AG(2,3) | Adjoin the two omitted affine points | No square values are manufactured; retain additive and square predicates |
| Root sign deck → Fano labels | Choose one of 30 rank-three quotient maps | Three deck dimensions, root signs, canonical characters; retain kernel |
| Fano quotient → surface type | Retain labels at the three actual crossings and the central line | d alone gives K^2 but not minimality; retain central component incidence |
| Seven roots → six roots | Forget r_j; quotient by <g_j> | One square-membership condition; retain which root was omitted |
| Matroid → dual | Ordinary vector-matroid duality changes ranks 3↔4 or 6↔1 | Not projective point/line duality and not a surface quotient |

Fano's dual matroid has rank four and is associated to the Hamming code, whereas projective Fano point-line duality stays rank three. The branch matroid has rank three but different circuits; the canonical seven-character matroid has rank six and dual U_(1,7). None of these dual operations by itself maps the hourglass to a Fano plane.

## 6. Reproduction, hostile controls and stopping boundary

Run

    python3 check.py
    python3 -O check.py

The standard-library checker uses explicit exceptions, never removable assertions. It checks all 128 subsets in six characteristics by two independent rank methods; all35 determinants; complete affine complements for p=2,3,5,7,11; all5040 binary label bijections; an independent enumeration of constant-weight-four three-spaces; all30 quotient fixed-locus and residual-coset incidences; the full seven-character dependence universe; and the complete D12 action on all thirty kernels. Normal and optimized output/certificates are byte-identical. The finite script verifies combinatorics and arithmetic of adjunction inputs, not the cited surface-classification theorems.

Hostile controls are built into the comparison: char2 repeated columns; char3 extra triples; canonical rank6 versus Fano rank3; the two d=3 labelings with fixed points rather than a free quotient; and the difference between resolution K^2=-4 and minimal-model K^2=0. The sixteen d=0 kernels are not sixteen nonisomorphic Campedelli surfaces; the orbit/descent argument gives exactly two, as explained in §4.

No execution denial occurred in this lane. A requested arXiv HTML version of the code/matroid paper returned cache-miss; one materially different supported retry used the primary arXiv PDF and the authors' published institutional PDF, both successfully inspected. No missing-page content was guessed. The local numerical/code census and all displayed cover calculations are complete; the arithmetic-square existence problem remains separate and open here.
