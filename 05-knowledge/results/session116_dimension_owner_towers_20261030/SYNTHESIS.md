# Session 116 synthesis: retain the multiplication, the lift, and the phase owner

**Status.** Proof-based synthesis of the four local packets, plus the integral full-sedenion norm-shell corollary below. Finite tests corroborate, rather than establish, the infinite statements. “Complete” always refers to the explicitly specified family. This session does not solve Collatz, the general Lonely Runner Conjecture, or the classification of zero divisors at every higher Cayley–Dickson stage. No literature-priority claim is made.

Sources within this packet: [tower](tower/TOWER_DEFECT.md), [arithmetic](arithmetic/ARITHMETIC_SHADOW.md), [second-band owners](lrc/SECOND_BAND.md), [controller completion](collatz/PHASE_COMPLETION_THEOREM.md), and [lineage and corrections](LINEAGE.md). The tower packet's multiplication convention governs all coordinates here, not the different Baez convention used by the concurrent integral-tower companion.

## 1. A dimension chart with its hypotheses attached

For the **standard real** Cayley–Dickson sequence A_n of dimension 2^n, use

    (a,b)(c,d)=(ac−bar(d)b, da+b bar(c)),
    N(a,b)=N(a)+N(b).

Here N is squared Euclidean norm. The first failures are commutativity at dimension4, associativity at8, and alternativity, multiplicative Euclidean norm, and cancellation at16. The conjugation identities, positive N, quadratic relation, power-associativity, and the existence of x*/N(x) as a two-sided inverse for x≠0 survive every doubling. Thus “inverse exists” is not synonymous with “left multiplication is invertible.” For xy=0 with x,y nonzero,

    [x*/N(x),x,y]=y,

where [u,v,w]=(uv)w−u(vw). This equation identifies the exact reassociation that invalidates cancellation. The embedding x↦(x,0) preserves this obstruction at every higher stage. The old raw claim that power-associativity fails at32 must not be inherited.

Dimension alone is insufficient: the standard octonions O have dimension8 and no zero divisors, while the noncanonical subalgebra W⊂S of the same real dimension has many. Its multiplication tensor is different. Likewise Clifford/Bott periodicity is not Cayley–Dickson doubling, and neither supplies an unmarked arithmetic-periodicity theorem.

## 2. From octonion associators to an exact integer bad-prime set

The actual bridge is a sequence of formulas, not a common dimension label. At the octonion-to-sedenion doubling,

    N((a,b)(c,d))−N(a,b)N(c,d)=2<[d,a,c],b>.

Thus an associator at the lower stage becomes a norm defect in its double. At later stages the four inherited lower-stage defects must be retained; this reduced formula is not valid unchanged at32.

For a sedenion x=(a,b), let A=Im(a), B=Im(b), N=N(a)+N(b), and Δ=N(A)N(B)−<A,B>². The exact Gram spectrum is

    N (multiplicity8), N+2√Δ (multiplicity4), N−2√Δ (multiplicity4),

and consequently

    det L_x=N⁴(N²−4Δ)².

For nonzero real x, singularity is equivalent to a,b being pure imaginary, orthogonal, and of equal norm; its right annihilator has dimension4. The packet proves this with a quaternion block decomposition. These real cone phenomena are classical; the packet does not claim their discovery.

In the standard integral basis the determinant is an integer polynomial. Reduction therefore gives the complete criterion, for **every prime p**, including2,

    L_x mod p singular  ⇔  p divides N(N²−4Δ).

The map preserves the determinant and singularity after reduction, but discards real positivity. For example x=(1,i) is real nonsingular but singular modulo2. Reusing the real cone criterion over finite fields would be false. Nor is this16-dimensional determinant formula claimed for every higher stage.

## 3. New corollary: every positive even norm occurs for integral sedenion zero divisors

Let S(Z) be the standard integer-coordinate sedenions. Then

    {N(x): x∈S(Z), x≠0, xy=0 for some nonzero y∈S(Z)} = {2m:m≥1}.

**Proof.** Necessity follows from the complete real cone theorem: x=(a,b) has pure orthogonal halves with N(a)=N(b)=m, a positive integer. Thus N(x)=2m.

For sufficiency, Lagrange's four-square theorem writes m=p²+q²+r²+s². On four fixed imaginary octonion coordinate axes, put

    a=(p,q,r,s),   b=(−q,p,−s,r),

with the other imaginary coordinates and both real parts zero. Then <a,b>=0 and N(a)=N(b)=m. Hence x=(a,b) is singular. Since its multiplication matrix has integer entries, its nonzero real kernel has a nonzero rational vector (rational and real ranks agree); clearing denominators gives a nonzero integral y with xy=0. This proves the full spectrum, not just a subfamily of norms. □

A constructive annihilator is also available. Write H=span{1,a,b,ab}; choose any nonzero integral c∈H^perp. Then

    y=(m c,b(ca))

is integral and xy=0. Such c exists by rational linear algebra. More explicitly, for an imaginary coordinate vector e, its scaled orthogonal projection is

    c=m²e−m<e,a>a−m<e,b>b−<e,ab>ab.

At least one of the seven imaginary coordinate vectors gives c≠0 because H^perp has dimension4. This uses N(ab)=m² inside the octonions; parentheses in b(ca) are essential.

**Why the W shell theorem is genuinely narrower.** Fix H⊂O=H⊕Hℓ and W={(a,sℓ):a,s∈H}. Its zero-divisor shell has an exact Cartesian-product parameterization by an imaginary-quaternion three-square representation and an independent quaternion four-square representation. With r_k(m) counting ordered signed representations,

    #{x∈W(Z):x zero divisor, N(x)=2m}
       =r₃(m)r₄(m)=8r₃(m)∑_(d|m,4∤d)d.

Jacobi supplies the r₄ formula; Legendre gives positivity exactly when m is not 4^j(8k+7). These are classical external dependencies, not session discoveries. Modern primary references verifying their precise forms are [Doležálek, Generalizing Hurwitz's quaternionic proof, introduction](https://arxiv.org/pdf/2409.16970) for Lagrange/Jacobi and [Danilkin–Chevalier, formalized Three Squares Theorem](https://isa-afp.org/entries/Three_Squares.html) for Legendre.

Accordingly W has no norm-squared14 shell, but full S does: use a=(2,1,1,1), b=(−1,2,−1,1) on four imaginary axes. This contrast is a genuine arithmetic restriction induced by the chosen subalgebra, not a global sedenion obstruction. It also separates three different statements: no zero divisors in canonical O; three-square restrictions in noncanonical W; every positive even norm in full S.

`check_synthesis.py` independently implements recursive multiplication and checks the integral projection and explicit annihilator for m=1,…,64, the specified m=7 witness, the W exclusion at7, embedded zero products at32, the parity basis-product map, and the unit/nonunit separation below. Its652 explicit exception gates pass identically under normal Python and `-O`. These finite controls do not replace Lagrange's theorem or the cone proof.

## 4. Characteristic four hides the real obstruction, but not every structure

For any signed-XOR basis product e_a e_b=(−1)^F(a,b)e_(a+b), every coefficient of the norm defect is 0 or ±4. The ideal generated by all coefficients, equivalently all integer defect values, is exactly (0) for a composing table and (4) otherwise. Sparse two-by-two supports realize a nonzero coefficient.

Thus over a commutative unital coefficient ring R an incoherent table composes this norm exactly when 4·1_R=0. At dimensions≥16 all normalized signed-XOR tables are incoherent over the reals, yet all have the characteristic-four norm identity. Modulo8 a sparse ±4 witness detects the failure. This reduction forgets signed rectangles and associators; it does not restore division, alternativity, or a positive definite norm.

The explicit parity map sends an integral signed-XOR vector to its coefficient vector in B_d=F₂[F₂^d], forgetting signs. With XOR convolution *, integer weights satisfy

    wt(a*b)≡wt(a)wt(b) mod4.

On the odd-weight units this is a multiplicative character, not an additive ring homomorphism. The exact integer unit-weight defect ideal is (0) at d=0,1, (8) at d=2, and (4) at d≥3. Thus unit-weight multiplicativity mod8 fails already at dimension8, while real norm composition remains possible there. The two thresholds cannot be conflated: the parity map loses integer carries.

## 5. Uniform odd-eighth extension: a Fourier family becomes an arithmetic subgroup

Put V=F₂^d, N=2^d, d≥3. Session115's complete family of phases having a Walsh coefficient of magnitude N−2 has the unique form

    f_(v,r,c)(x)=c+r·v+r·x+δ_v(x).

Use the truth table as a coefficient vector in B_d. This is an explicit map retaining the whole labelled phase. The new convolution law is

    f_(v,r,c)*f_(w,s,e)=f_(v+w,r+s,c+e).

Hence this family H_d is a subgroup isomorphic to F₂^(2d+1); its arithmetic character is wt(f)≡(−1)^c mod4. Since |B_d^×|=2^(2^d−1),

    [B_d^×:H_d]=2^(2^d−2d−2).

At d=3 H_d exhausts all units; at d=4 its index is64. This is a uniform classification replacing the exceptional dimension8 coincidence. The familiar readings 1/8,3/8,5/8,7/8 are the d=3 Walsh readings, not universal constants forced at every dimension. Convolution is commutative and is **not** octonion multiplication or the sign-gauge coboundary law. Forgetting (v,r) preserves the arithmetic character but loses named Fourier coordinates.

**A sharper separation at dimension16.** Coefficientwise parity reduction is an actual multiplicative map from the standard integral Cayley–Dickson algebra A_d(Z) to B_d: every signed basis product becomes the unsigned XOR product. Moreover N(x) mod2 equals augmentation of its image. Therefore every real integral sedenion zero divisor reduces to a **nonunit** (its norm is even), whereas every member of H₄, and every unit outside H₄, has odd augmentation. Every integral lift of such a unit has odd norm and is consequently nonsingular as a real sedenion. Thus index64 counts cosets among parity units; it is not a zero-divisor count or a stratification of the real singular cone. The map is multiplicative but forgets precisely the signs needed to recover the real obstruction. No converse from nonunit reduction to real singularity is asserted.

## 6. Every declared second-band core: change the owners, not the harmonic law

**Concurrent-source credit.** The same owner-promotion family and harmonic transfer already occur in incoming main commit `98f961a7`, packet `odd_eighth_owner_families_20261030`. Treat this section and the local owner packet as independently audited recovery/refinement, not an exclusive new discovery. Pinned source: [main at98f961a7](https://github.com/eliottcassidy2000/math/tree/98f961a7).

The proved family is explicitly

    C(q,A)={1,…,q−1} ∪ {q+b:b∈A},
    q≥3, 1∈A⊆{1,…,q−1}.

For units u modulo q, define e_A(u)=u+q·1_A(u), E_A={e_A(u)}, and s_A=min E_A. On the strictly bounded chamber 1/(q+2)<r≤1/q the full safe windows are centered at primitive a/q, with left/right owners e_A(a^−1), e_A(q−a^−1). This complete directional atlas is the map that justifies extension. Nonunit attachments are invisible in this chamber only.

There are exactly 2^(φ(q)−1) atlases; possible minima are all units greater than1 and q+1. The old least-missing-prime shortcut is not invariant: q=7,A={1,2,3} has minimum4.

For the harmonic bank {qm,2qm,…,kqm}, k≤q and 2m>s_A, the unchanged functional law is

    M=min{m/(qm+s_A),1/(k+1)}.

For k=q, preservation of harmonic clearance 1/(q+1) holds exactly when m≥s_A. Every census, weight, and period transfers using the actual E_A. Adding q+u replaces only owner u by q+u and deletes precisely the corresponding phase pairs between the two floor cutoffs. It changes multiplicities and mass even when clearance stays fixed.

At q=7 every core therefore has clearance1/8 exactly for m≥s_A; all four odd-eighth phases occur exactly for m≥3s_A. Their phase-weight law stays fixed while the owner atlas selects survivors. This is the requested family-wide extension, not merely another example. The lower endpoint is excluded; below the chamber even the inherited floor equality can fail.

## 7. Controller entries: retain phase cosets, reduce an infinite search to finitely many rotations

**Concurrent-source credit.** Incoming commit `bdbe68566e8`, packet `entry_cosets_resonant_echoes_20261030`, already contains the global completion/rotation and ultrametric result. This section records independently audited recovery/refinement, not an exclusive new theorem. Pinned source: [main atbdbe68566e8](https://github.com/eliottcassidy2000/math/tree/bdbe68566e8).

For any nonempty positive valuation word w, write w=v_j s_j at its d nonempty prefixes. The controller packet, using the Session114 ordinary periodic-mean theorem, proves

    max_u R(u,w)=max_(1≤j≤d)R(s_j,w)
               =max_(1≤j≤d)R(empty,s_j v_j).

Every finite entry u is simultaneously dominated at every 5-adic layer by one of these completion entries. The proof retains its phase coset: choose a phase with maximal membership depth in the cyclic affine subgroup generated by F_w. Through that depth the completion has exactly the same hitset, and above finite depth the original has none while the completion retains a hit. This produces an attained finite maximum, not an unsupported exchange of limits or a p-adic supremum.

The symmetric phase-depth matrix gives a complete finite optimizer classifier. The earlier empty-entry rule does not extend to two resets: w=(3,1,2) has empty-entry rate9/16 but entry(1,2) rate13/16. For nonprimitive w=z^h, use R(u,w)=hR(u,z) and the primitive root when classifying the entry language.

The analogy with owner promotion is methodological only: both reductions work because they retain the coordinate that controls admissibility (directional speed or phase coset). There is no asserted isomorphism between those systems, and no arithmetic conclusion follows solely from a shared denominator8. These are exact local advances on the repository's scoped questions, not convergence of a Collatz orbit or general LRC closure.

## 8. Incoming sparse-tower classification and typical versus worst-case norms

**Concurrent-source credit.** Incoming main `7f70dcd6012439b7cf75fb43f26850f3807126b9`, packet `cayley_dickson_integral_20261030`, C1–C5, already proves the integral two-support block theorem, all-dimension named-pair and census recurrences, every-odd-prime obstruction beginning at dimension16, and random-sign squared-norm variance asymptotic to `2/N`. These are credited source results, not Session116 novelty. [Pinned source](https://github.com/eliottcassidy2000/math/blob/7f70dcd6012439b7cf75fb43f26850f3807126b9/05-knowledge/results/cayley_dickson_integral_20261030/PROOFS.md).

Our arbitrary-multiplier real classification in dimension16 complements that all-dimension integral **two-support** classification. The independent [C1–C5 audit](audit_sparse_merge/AUDIT.md) implements the incoming Baez convention explicitly rather than copying witnesses from our different tower convention. Every basis pair and norm-defect rectangle through dimension128 agrees with the recurrences; direct modular ranks and 336 composite-modulus block-kernel tests also pass normally and under `python -O`.

There is no contradiction between proliferating sparse cancellation failures and improving typical dense norm behavior: the latter concerns two independent random-sign inputs, not an operator bound uniform over every second input. For non-fixed-length independent coordinates, the source's fourth-moment-free extension concerns the composition defect `||xy||²−||x||²||y||²`, not centered product norm `||XY||²−1`. No arbitrary-dimensional dense classification, convention transport of named witnesses, or C6+ automaton claim is imported here.
