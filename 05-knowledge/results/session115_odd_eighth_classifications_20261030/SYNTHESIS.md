# Odd eighths: a formula atlas and two exact cross-field maps

**Session 115, 30 October 2026. Scoped synthesis of the accompanying proofs.** Equal fractions do not by themselves identify observables, measures, dynamics or open problems. Here two connections are literal maps on named objects; the other connections preserve explicitly stated predicates. No major conjecture closure or literature-priority claim is made.

## 1. What is completely classified

| Carrier | Uniform formula or complete classification | Boundary |
|---|---|---|
| Coherent relative Fano gauges | Sixteen classes; relative GL(3,2) orbit sizes 1,7,7,1. The eight cubic classes have precisely the odd-eighth Walsh profiles. | The base multiplication table is retained. Abstract octonion isomorphism forgets these sectors. |
| Boolean single-point defects in dimension d>=3 | Exactly the phases with some Walsh coefficient of magnitude 2^d−2; probabilities lie in {2^−d,1−2^−d,1/2−2^−d,1/2+2^−d}. | This exhausts all odd-weight phases only when d=3. |
| Affine groups V⋊H, V=F₂ᵈ and H odd | Odd-order proportion is the average over h of 2^(−dim Fix(h)). In dimension 3 there are exactly four complement types. | The action, fixed-space profile and multiplicities matter. Group order alone is insufficient. |
| All finite groups with Sylow D8 | Even-order proportion is (2c+1)/8, with c=1,2,3 involution classes. | This is not a realization theorem for tournament switching automorphism groups. |
| Word-address transports q→p | At depth A, products, real values, positivity and permutation-pullback structure are preserved iff p≡q mod 2^(A−1). First failure for distinct multipliers: v₂(p−q)+2. | Linear unitarity alone does not preserve probabilities or nonlinear statistics. |
| Harmonic-tail completion of C_q={1,...,q−1,q+1} | Complete preservation threshold for the q-tail block's optimum: m>=s(q), the least prime not dividing q; all maximizing points and their weights have exact floor-sum formulas. | The full spectrum formula is restricted to the stated upper chamber; the lower chamber remains open. |
| Constant Collatz controllers a^d | Exact all-entry rates and complete equality-entry classification; sharp clock bound 5/8 only at a=1 with an all-one entry. | Fixed periodic controllers, not arbitrary changing words or convergence. |
| Every one-reset controller 1^L a, a>1, L>=0 | Exact maximal all-layer mean over every entry and complete maximizing-entry classification by an explicit mixed-radix digit rule. | Ordinary Lucas support is a special case. The strict clock bound below 1/2 is separately proved for a=2, not asserted for every reset a. |

Proofs: [norm phases](quadratic/ODD_EIGHTH_PHASES.md), [groups](switching/ODD_EIGHTH_GROUPS.md), [word transport](unitary/WORD_UNITARY_CLASSIFICATION.md), [harmonic completion](lrc/HARMONIC_COMPLETION.md), [constant controllers](collatz/CONSTANT_CONTROLLERS.md), [one reset](collatz/ONE_RESET_CLASSIFICATION.md), [every one-reset controller](collatz/GENERAL_RESET.md).

## 2. Exact bridge: affine element order becomes a Boolean measurement

This connects the group and norm-phase lanes through an actual event function on the same named binary cube. Choose the binary coordinate labelling used in the Fano packet; changing that choice transports the phase and its character labels together.

Let V=F₂ᵈ, let h have odd order ℓ, and put

    N_h=1+h+...+h^(ℓ−1),
    f_h(v)=1[N_h v≠0].

The affine element (v,h) satisfies (v,h)^ℓ=(N_h v,1). Consequently

    f_h(v)=1[(v,h) has even order].                         (1)

For a Boolean phase f define the named character reading

    p_f(u)=Pr_uniform v[f(v)+u·v=1].

Then the group observable is exactly

    E(V⋊H)=(1/|H|) Σ_h p_(f_h)(0).                         (2)

Thus this is not an analogy between two numbers: it is equality of predicates before averaging.

Because ℓ is odd, N_h is a projection onto Fix(h). With r=dim Fix(h), its zero-kernel indicator also has the exact Fourier expansion

    1[N_h v=0] = 2^(−r) Σ_(λ∈Fix(h)^*) (−1)^(λ(N_h v)).   (3)

The coefficient at the zero character is 2^−r, yielding the group formula. Equation(3) exhibits the retained linear constraints behind the count.

### The complete three-dimensional mixture

For the odd-order elements of GL₃(2), only orders 1,3,7 occur.

- h=1: f_h=1−δ₀, the raw cubic phase that is 1 at every nonzero cube point. Its named zero-character reading is 7/8. It belongs to the unique relative Fano gauge class with the7/8 exceptional reading, corresponding to A=B=C=D=1. The normalized representative is s=f_h+x₁+x₂+x₃, so this same reading is p_s(111)=7/8, whereas p_s(0)=3/8. Normalization shifts the named character; it does not preserve that label silently.
- h of order3: Fix(h) has dimension 1, so f_h is a nonzero linear function. Its zero-character reading is 1/2.
- h of order7: Fix(h)=0, so f_h=0 and its reading is 0.

The four complete complement types therefore give

    H=1:          7/8;
    H=C₃:         (7/8+2·1/2)/3 =5/8;
    H=C₇:         (7/8+6·0)/7 =1/8;
    H=C₇⋊C₃:      (7/8+14·1/2+6·0)/21 =3/8.              (4)

All four fractions now have one common meaning: the proportion of even-order affine transformations. They are obtained from named Boolean readings by the explicit map h↦f_h and actual element multiplicities.

### The loss exposed by this bridge

The relative norm-gauge quotient identifies Boolean functions differing by a linear function. Thus zero and every nonzero linear f_h belong to the same quotient class, although their **named zero-character** readings are0 and1/2. Their unordered Walsh profiles agree because linear translation permutes character labels. More generally, replacing f by s=f+ℓ requires replacing the observer u by u+ℓ to preserve p_f(u)=p_s(u+ℓ).

Therefore the unordered norm-gauge profile does not suffice to reconstruct (2). The missing data are the named character and the multiplicity of each phase. There is no claimed isomorphism between a finite group and an octonion algebra, nor an identification of their invariant theories. Only the event-function and Fourier maps above are transferred.

## 3. Exact bridge and obstruction: LRC weights through word transport

The harmonic seven-tail family has primitive tail phases t=j/8, j∈{1,3,5,7}. Its reciprocal-owner weight before the physical scaling factor is

    w₈(j)=88/105+(32/105)χ₂(j),
    χ₂(j)=(−1)^((j²−1)/8)=√2 cos(2πj/8).                 (5)

Here χ₂ has values (+1,−1,−1,+1) on(1,3,5,7), so the weights are 8/7 on 1,7 and 8/15 on 3,5. The actual individual physical weight in the q=7 family is w₈(j)/(7m).

For the complete maximizer census, let N(7,m) be the number of physical maximizing points and

    A(m)=Σ_(e∈{2,3,4,5,6,8}) Σ_(odd j<=floor(m/e)) χ₂(j).

The exact boundary slope is

    B(7,m)=88N(7,m)/(735m)+64A(m)/(735m).                 (6)

The signed tally A is bounded by 6 in absolute value and has period 960. It is the missing quadratic-character correction to unweighted witness counting. Reflection preserves χ₂, while multiplication by 3 swaps its signs.

### A precise non-transfer theorem

Consider the word-address unitary U_(3→5) on the odd-residue algebra modulo 8. Its primitive-frequency permutation is

    1→7, 3→1, 5→3, 7→5.

Let χ_(−2) have values (+1,+1,−1,−1) on those odd residues. Since

    χ₂(j)=(√2/2)(exp(2πij/8)+exp(−2πij/8)),

the permutation gives, on odd j,

    U_(3→5)χ₂(j)
      =(√2/2)(exp(2πi·7j/8)+exp(2πi·5j/8))
      =−iχ_(−2)(j).                                      (7)

Hence

    U_(3→5)w₈(j)=88/105−i(32/105)χ_(−2)(j).               (8)

The positive physical weight has become nonreal. The linear transport therefore cannot be used to carry the LRC boundary-weight inequality, even though it preserves the relevant quadratic norm. This is an explicit obstruction on the actual four phase labels, not merely a warning about unrelated fractions.

The restriction to odd residues is essential. The real cosine in(5) is a representative of that restricted function; it is not the zero-on-even extension of the Dirichlet character on the whole cyclic group.

The classification of the unitary gives the exact repair gate. At depth 3, if p≡q mod 4, the transport is a genuine permutation pullback and preserves positive functions. In general, depth-A probability or positive-weight transport through this particular unitary is licensed exactly by p≡q mod 2^(A−1). Even inside that gate, transporting a weight function is not automatically a transport of the physical LRC dynamics or of its source speed set.

There is a sharper distinction for this particular weight: positive transport requires p≡q mod4, but **U_(q→p)w₈=w₈ holds exactly when p≡q mod8**. Within the positive gate, c=q/p mod8 is1 or5 and χ₂(cj)=χ₂(c)χ₂(j). The case c=5 reverses the character correction and exchanges the two weight levels; it does not fix the weight. Outside the mod4 gate the image is nonreal. All16 odd residue pairs modulo8 are checked exactly; residue1 may be represented by a legal multiplier such as9.

## 4. Exact bridge: maximizing return phases become Lucas support

For every unit-reset controller w=1^L a with a>1 and a≢1 modulo4, let K=floor(L/2), δ=L mod2. A maximizing entry is characterized by a phase j=2c+ε with

    0<=j<=L,
    ε<=δ,
    each base5 digit c_i<=K_i,
    1^j u=w^b for some b>=0.                              (9)

The digit condition is equivalent to

    binom(K,c)≠0 mod5,

because over F₅,

    (1+X)^K=∏_i(1+X^(5^i))^(K_i).                         (10)

This maps an actual controller-entry equality predicate to the support of a Frobenius-factorized polynomial. It preserves the condition of attaining every return-layer cap. The number of admissible phases is

    (δ+1)∏_i(K_i+1).                                      (11)

### All reset letters: the intersection coordinate repairs the classifier

The preceding Lucas form is the unit-reset special case h=v₅(2^a−2)=0. The complete result now covers **every a>1 and L>=0**. Put

    P=3^(L+1), Q=2^(L+a), B=3^(L+1)−2^(L+1),
    t=ord₅(P/Q),
    r=v₅(P^t−Q^t) if t>1,
    r=min(v₅(P−Q),v₅(B)) if t=1,
    T_e=t·5^max(0,e−r).

The retained coordinate is the exact subgroup-intersection modulus g_e:

    t>1, e<=h: g_e=5^max(0,min(r,e)−1);
    t>1, e>h:  g_e=2·5^max(0,e−h−1,min(r,e)−1);
    t=1:       g_e=2·5^max(0,min(r,e)−1,e+r−h−1).

Then the maximum ordinary block mean and complete equality predicate are

    max_u R(u,1^L a)=(r+1/4)/t+Σ_e floor(L/g_e)/T_e,
    1^j u=(1^L a)^b and j mod g_e<=L mod g_e for every e,
    0<=j<=L, b>=0.                                       (12)

Empty entry attains the maximum. Divide the block mean by L+a for the original-clock rate. The distinct place values in 1,g₁,g₂,... give these successive radices:

- h=0 or t=1: 2,5,5,...;
- h>0,t>1,r<=h: r−1 copies of5, then2, then5,5,...;
- h>0,t>1,r>h: h−1 copies of5, then10, then5,5,... .

Writing L and j in the same indicated radix system, equality is exactly j_i<=L_i for every digit; the number of phases is ∏_i(L_i+1). The first radix chain has the ordinary Lucas interpretation above (including the t=1 case), but the other chains are mixed-radix no-borrow predicates, not automatically support identities in a single prime characteristic.

At a=5,L=2 the formerly omitted intersection gives g₁=1, correctly retaining all three first-layer positions. At a=5,L=8 the radix begins10 rather than2, and all nine phases j=0,...,8 maximize. Discarding the intersection coordinate would falsely retain only five. Thus this extension repairs a specific failed transfer by carrying its missing sidecar, rather than forcing the unit-reset Lucas formula onto every reset letter.

The count is not the entry itself: the phase j and completion exponent b remain necessary to reconstruct u. Nor does a nonzero binomial coefficient prove a result about arbitrary Collatz trajectories. The exact cyclic-subgroup intersection is what first produces the residue-occupancy inequalities to which Lucas applies.

## 5. Three different arithmetic meanings of an eighth

The constant-controller classification distinguishes quantities that should not be silently equated.

- For the all-one controller, the ordinary all-layer mean per valuation clock is

      5/8=Σ_(e>=1)1/(2·5^(e−1))=1/2+1/8.

  The1/8 here is the higher-layer tail of a five-adic return sum.
- For a=5, the maximal period mean is5/8 but the original clock is5, so its clock-normalized rate is1/8.
- On the stated profinite exponent parameter space, conditional on a being even, the maximal **block-mean expression** has Haar mean3/8. This uses the exact geometric stall-depth law. It does not call a noninteger parameter a legal controller, and does not establish an ordinary Cesàro average over initial integer values of a.

Thus5/8 and3/8 here are not complementary events in one common probability space. The group formula(4) really does put all four odd eighths in one observable; the Collatz formulas use distinct, explicitly retained normalizations and measures.

## 6. What the synthesis adds, and where it stops

Two concrete transfer consequences are now available:

1. The affine-group order count can be expressed as an average of named Boolean/Walsh measurements, with a demonstrable loss when character labels are discarded.
2. The LRC quadratic-character weight can be passed through the word-address operator explicitly; the first nonmultiplicative shell turns it nonreal and prevents an otherwise tempting positivity argument.

The Lucas connection gives a third exact change of representation, from unit-reset maximizing return phases to polynomial support. The complete one-reset extension retains the intersection modulus and replaces that support by mixed-radix no-borrow digits where necessary. These maps may be reused because their preserved predicates and missing information are stated. None identifies the underlying open conjectures or supplies the missing arbitrary-source, arbitrary-controller, or switching-realization step.

`check_synthesis.py` verifies the two new cross-lane maps with exact arithmetic: cyclotomic arithmetic modulo z⁴+1 for the four odd phases; literal affine orders for all odd-order matrices of GL₃(2); the four complement mixtures; and all16 residue-pair tests distinguishing positive weight transport from weight invariance. Its 1,136 gates pass normally and under `python -O`. The all-parameter formulas and classification claims are proved in the linked source packets; this small checker tests the synthesis itself.

## 7. Concurrent merge: a valid conditional-energy decoder

The incoming commit **f128982d824452a919602687554ffea9ef33155c**, published while this session was running, contributes a complementary positive bridge. Its [conditional-energy theorem](../odd_eighths_uniform_families_20261030/COLLATZ.md) and [phase decoder](../odd_eighths_uniform_families_20261030/PHASE_BRIDGE.md) receive their own credit; they are not new claims of this packet. Its multiplier-transport section in turn credits this session's earlier checkpoint.

For its precisely defined q=3 source-character sum S_M, let

    K_M(c)=3^M E[|S_M+1/2|² | R=c mod8],   M>=3.

The expectation is normalized Haar measure on a source cylinder, not evaluation at the integer c. The incoming theorem proves

    K_M(c)=1/4+(-1)^M(3√2/64)χ₂(c).

Consequently the physical LRC phase weight has the exact decoder

    w₈(c)=88/105+(1024√2/315)(-1)^M(K_M(c)−1/4).

For the harmonic family, averaging this expression over its **actual maximizing points**, retaining their phase multiplicities, and multiplying by N_m/(7m) recovers its terminal mass slope. Uniform averaging over four abstract phases generally gives the wrong answer. Clock parity is also essential.

There is no conflict with the earlier unitary obstruction: this map compares two responses to the same named phase bit, whereas U₍₃→₅₎ changes the character realization and fails positivity at that depth. The incoming scalar responses recover the reflection-even bit; the complete oriented phase and physical point still require their other data. Neither map transports a positive integer trajectory or proves fixed-seed convergence.

The incoming commit is preserved in this session's final ancestry. A separate [independent audit](audit_incoming/AUDIT.md) checks its energy identity and decoder; the prior-source theorem remains separately identifiable.
