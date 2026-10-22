# Session 27: repunits, lattice observables, and reconstruction from one bit

This session combines prior-work recovery, four manuscript audits, and exact structural results. It does not claim the complete Hamiltonian-count spectrum, a minimum-10 lattice in dimension96, LRC(14), Collatz, or Hadwiger–Nelson. Classical and recovered results are attributed; independently derived consequences are not advertised as literature-first discoveries.

## 1. Main structural result: a faithful finite tournament observable

For a labelled k-vertex tournament Q, attach an ordered transitive base with a word of k-bit incidence masks. Let F_Q(w) be the Hamiltonian-path count of that extension. These are marked experiments: both the core labels and the attachment addresses are retained.

For every k>=2, experiments with exactly k-1 base vertices reconstruct every arc of Q. More strongly, the single bit

    b(T)=((H(T)-1)/2) mod2

of each response suffices. This is the second least significant bit of H, not its constant parity bit.

Fix core vertices u,v. Arrange all other vertices along an alternating core/base path R beginning at v and ending at the last base vertex. It has2k-3 cross edges. Fix the last base vertex to point toward u and vary only the orientations of the R edges. Then

    XOR of b(T) over these 2^(2k-3) extensions = indicator(u->v).

The proof cancels all odd cycles missing a varied edge. A surviving odd cycle must traverse the whole alternating path, pass through u, and use precisely the target arc u->v. A separate full-H proof gives a signed sum equal to0 or2. Therefore every raw modulus M>=3 also suffices, without dividing by2 in a residue ring. Raw H modulo2 carries no information about a core of order at least two, by Rédei parity.

A fixed reconstruction battery uses at most binom(k,2)2^(2k-3) oracle calls. This is an information theorem, not a practical tournament-enumeration speed claim: computing the queried H values can be expensive. The whole fixed-length response table is faithful under consistent relabelling; adding reversal/complement of the contexts gives the converse quotient used for merged metagraph nodes. It is not a classification of all metagraph edges from their degree counts.

Sources: [full integer reconstruction](session27-response/faithfulness-proof.md), [fixed-endpoint/modular improvement](session27-response/MODULAR_FAITHFULNESS.md), [direct odd-cycle proof](ODD_CYCLE_PARITY_CHARACTER.md), [independent reconstruction audit](session27-response-audit/RECONSTRUCTION_REVIEW.md), [modular audit](session27-response-audit/MODULAR_REVIEW.md).

## 2. Two ports recover the first response geometry exactly

Let S_u and E_u count Hamiltonian paths starting and ending at u. Let C_uv count ordered spanning covers by two nonempty directed paths whose first path ends at u and whose second starts at v. The correct one-port polynomial is

    F_Q(b)=H(Q)+sum_u L_u b_u - sum_(u<v)(C_uv+C_vu)b_u b_v,
    L_u=E_u-S_u+sum_v C_uv.

The two-path covers are essential: deleting the added vertex need not leave a Hamiltonian path in Q. An insertion-only formula fails already on a two-vertex core.

Define J=C-C^T and g=2(E-S)+div J. Then

    F_Q(b)-F_(Q converse)(b)=g dot b.

One port sees a divergence, not the full circulation. Every vertex-transitive core is one-port converse-blind. The directed triangle supplies a nonzero circulation with zero divergence; two ordered attachments distinguish its two labelled orientations by13 versus11.

Two ports recover more than circulation: exactly1+2k+3binom(k,2) queried entries recover H,S,E and the entire ordered matrix C. Its entries are the derivatives of the Hamiltonian-path polynomial with respect to directed-arc variables, including absent arcs. Reversing an actual arc u->v changes H by C_vu-C_uv. The sign of this current does not identify the arc direction: a transitive-three counterexample is explicit.

Complete two-port response distinguishes every labelled core through k=6 in the finite census. Whether two ports suffice in every order remains OPEN. The all-order theorem in section1 uses k-1 ports instead. A separate shared weighted-automaton realization supplies exact equivalence testing for arbitrary candidate rational response series; it does not replace the stronger structural reconstruction proof.

[Proof, tomography and controls](session27-response/proof.md); [independent audit](session27-response-audit/REVIEW.md).

## 3. An actual E8 root can hide different Hamiltonian counts

The score map has a legitimate lattice target, but it loses information. Start with regular Paley Q on seven vertices. Reversing the directed triangle0->1->3->0 preserves all scores and changes H from189 to171. Adjoin the same eighth vertex to both. The resulting tournaments have the identical labelled score vector

    (4,4,4,3,3,3,3,4),

but H=657 and643. Subtracting7/2 gives exactly the same E8 root, with four coordinates+1/2 and four−1/2, squared norm2. This is stronger than merely equal norms or equal theta coefficients: the lattice point itself is identical.

The first score-to-H failure already occurs at order five. Scores determine cyclic triples; for n<=4, H=1+2c3. At n=5, H=1+2c3+2c5. Our smallest pair has the same labelled scores, c3=4 in both, but c5=2 versus1, hence H=13 versus11. The missing information is identified exactly.

[Construction, full masks and first-failure proof](SCORE_ROOT_COLLISION.md). Both eight-vertex counts were independently checked over all40320 vertex permutations.

## 4. Why31,63,127 are rare, with the measure stated

A tournament decomposes uniquely into a totally ordered sequence of strong components. H multiplies across those components; automorphism orders also multiply. For fixed attainable h, let w(h) be the largest number of nontrivial strong factors in an actually realizable factorization of h.

Both the number N_h(n) of ordinary isomorphism classes and the automorphism-weighted count W_h(n)=sum 1/|Aut(T)| are eventually positive polynomials of degree w(h). Labelled counts equal n!W_h(n). For each ordered strong-component sequence of total order M and length r, its contribution is

    binom(n-M+r,r), or binom(n-M+r,r)/product|Aut(component)|.

The generating functions are finite rational sums with denominators(1-z)^(r+1). Their finiteness uses an inherited strong-tournament lower bound, not a conjecture about the full spectrum.

31 and127 are attained primes, so w=1. The composite63 also has w=1: every proper nontrivial factorization runs into7 or21, which are impossible Hamiltonian counts. Thus it is indecomposable in the relevant spectrum even though it is not prime. Meanwhile33=3*11,65=5*13,129=3*43 admit two actual strong factors. Consequently the corresponding31/33,63/65,127/129 ratios tend to zero as vertex order n grows, under both uniform-class and uniform-labelled comparisons.

For63 there are exactly two strong atoms, both asymmetric of order eight. Hence, for n>=8,

    N_63(n)=W_63(n)=2(n-7),
    labelled count=2(n-7)n!.

This is not a universal repunit-rarity law. The counterexample255=3*5*17 has degree3, while257 has degree1, so their eventual comparison goes the other way. Nor does this padding argument explain the distribution among strong classes alone.

[Measure-aware theorem, precise source dependencies and independent tests](session27-repunits/REPORT.md).

## 5. The prior-repunit catalogue and its coverage

The expanded local discovery pass indexed113808 tracked paths and searched72894 physically present text files. It recorded7290 broad matches;249 files contained the literal term“repunit.” Broad results include other subjects and are not all Hamiltonian-count results.40488 tracked paths were absent locally; two oversized files were skipped.

A complementary five-query remote pass returned415 distinct paths, including174 absent locally. Four queries reached their100-result caps. It is therefore not honest to call this an exhaustive search of every historical Git blob. The complete saved catalogues, line previews, hashes, commits, search expressions, and remote limitations are supplied.

The recovered work separates at least four different carriers:

- Mersenne tournament ORDERS in an early coefficient formula;
- k-nacci/necklace TRACE values;
- actual H-values and their multiplicities;
- Collatz repunit digit, debt and clock families.

Their matching integers do not identify their predicates. The catalogue restores March work predating the October synthesis, the first H=63 realizations, the correction lineage, and the S51 independent audit.

[Curated routes and local catalogue](session27-repunits/REPORT.md); [remote recovery](session27-repunits/REMOTE_SEARCH.md).

## 6. The genuine modular and lattice connections

With q=exp(2πiτ), eta=q^(1/24)product(1-q^m), and Delta=eta^24. The factor1728 in(E4^3-E6^2)/1728 normalizes its first coefficient; it is not derived from the adjacent integer1729. The Eisenstein series E8 is E4^2, whereas the E8 LATTICE theta series is E4.

The character identity

    theta_E8/eta^8 = q^(-1/3)(1+248q+4124q^2+...),

explains248 as240 lattice roots plus8 oscillator states. For a rank24 Niemeier lattice with R roots,

    theta_N/eta^24 = j+R-720.

The Leech character is therefore j-720, not j-744. Ordinary binary Golay Construction A has48 roots and gives the A1^24 Niemeier lattice, not the rootless Leech lattice without a further construction.

II_(25,1) is a rank26 Lorentzian even unimodular lattice, not a26-dimensional positive-definite sphere packing. Primitive null quotients give the24 Niemeier classes. Leech deep holes have23 types corresponding to the rootful classes. In II_(9,1), every primitive null quotient is the unique E8 lattice, and the E8 deep-hole geometry has an explicit16-vertex cross-polytope. These are actual constructions with distinct domains and normalizations.

An incoming concurrent proposal also yields a precise intransitivity map: for a25-vertex tournament, the uncentered score lift(s;70) lies in II_(25,1) and has norm-2c3(T). It is null exactly for transitive tournaments; their score permutations give the classical cannonball vector(0,1,...,24;70) and its Leech quotient. The quotient is not defined in the same way for timelike, cyclic cases. Even the Lorentzian point loses H: padding our eight-vertex hostile gives the same norm-40 point with counts657 and643. [Proof, concurrent-source attribution and boundary](LORENTZIAN_SCORE_BRIDGE.md).

The user's Heegner radicand list, including1, is the correct nine-value list. Discriminants are different labels. E8 self-duality of root data, Euclidean lattice duality, and modular transformations should not be conflated merely because they agree in selected examples.

A separate scoped obstruction is proved: no nonzero rational function of q holomorphic in the disk can be a positive-weight holomorphic level-one modular form. It applies to the native fixed-H rational tournament series above. It does NOT prohibit every transform, enriched invariant, or different modular relationship.

[Full dictionary, exact identities and primary sources](session27-modular/report.md); [E8 hole and null-quotient addendum](session27-modular/PRIOR_1729_AND_E8_HOLES.md).

## 7.1729: recovered recurrence, K3 geometry, and arithmetic

The repo already contained a reversible recurrence with denominator1-82x-82x^2+x^3. Its negative-index triple(-9,12,10) recovers9^3+10^3=12^3+1=1729. Ono–Trebat-Leder's associated elliptic K3 construction has geometric Picard number18; the raw homogeneous equationX^3+Y^3=Z^3+W^3 is instead a cubic surface. These must not be identified.

The seven-mode proof of the cubic recurrence identity is already in Hirschhorn's1996 work. We recovered that attribution and proved a uniform synthesis: for characteristic(x+1)(x^2-tx+1), all homogeneous degree-d observables are annihilated by an explicit degree2d+1 polynomial, over every commutative unital ring. Repeated roots, composite characteristic and nilpotents are covered by an integral polynomial-matrix identity. Generic sharpness concerns the entire observable space, not every individual polynomial.

There is also a separate genuine number-theoretic link:1729=7*13*19 is Carmichael and thus passes the inherited tournament-count congruence. Luca's2026 theorem restricts Carmichael numbers of form2^n p^3+1 with prime p to p=3. It does not claim1729 is the unique such number. The Camion converse fails at21. None of this alone identifies the quotient H(Paley11)/|Aut|=1729 with a CM period.

[Recovered sources](session27-modular/PRIOR_1729_AND_E8_HOLES.md); [uniform recurrence theorem](session27-recurrence/RECIPROCAL_POWER_THEOREM.md); [independent audit](session27-recurrence-audit/AUDIT.md).

## 8. Four papers: verified positives and failed implications

All four attachments were read fully. The strongest verified positive is the exact threshold for their SPECIFIED four polynomial ADE moments: no collision through rank42, and one collision at rank43. Independent enumeration covers7122017 formal multisets. Its short mechanism is a rank-mass polynomial

    z^2(z-1)^5(z+1)(3z^2+4z+2),

plus exactly cancelling D-family defects. The fifth moment differs by158592. The factorization proves the witness; the exhaustive census separately proves minimality and uniqueness. This is not a claim about all lattice theta coefficients without an identification of the moments.

The old96-dimensional Gram file is genuinely integral, even, positive definite and unimodular by exact checks. Its minimum10 remains uncertified, not disproved. The older paper's q^6 coefficient is wrong; the corrected formal candidate has2972108280960000. A formal theta candidate, or a bounded short-vector search, is not a minimum certificate.

The Niemeier manuscript's hyperelliptic separation claim is false. Codogni's published theorem makes all five equal-root-count pairs coincide there in every genus, giving19 classes. An explicit coefficient11520 also refutes the manuscript's proposed vanishing hypothesis. The genus5 JACOBIAN separation of all24 classes survives by Codogni's separate established theorem. Ambient coefficients, Jacobian restriction and hyperelliptic restriction are different observables.

[Complete manuscript review](session27-papers-review/REVIEW.md); [independent controls audit](session27-lattice-controls-audit/AUDIT.md); [independent Niemeier restriction audit](session27-niemeier-audit/independent-restriction-audit.md).

## 9. A corrected, exact gluing-realizability certificate

For split orthogonal V=F_q^(2r), q odd, sample BOTH families of maximal isotropic r-spaces. If W is isotropic of dimension d, its containment probability is

    1 / product_(j=r-d)^(r-1)(q^j+1).

The factor2 at d=r matters. For forbidden projective lines B, exact avoidance is a Möbius sum over their FULL span poset, assigning zero to nonisotropic spans. Positivity of that rational sum is an exact finite realizability criterion. The simpler sufficient bound is|B|<q^(r-1)+1.

Pairwise orthogonality alone is insufficient: in F_3^6, two triples both have orthogonality graph K3 but span dimensions2 and3, producing60 versus61 avoiding maximal spaces out of80. The missing coordinate is linear dependence. This is not an analogy: it is a computed and proved obstruction to replacing the span poset by its pairwise graph.

At q=7,r=24, a genuinely proved forbidden-line bound near10^12 would still suffice after correcting the old paper's symplectic count. The complete actual short-vector forbidden set has not been counted; no rank96 extremality conclusion follows.

[Exact certificate and independent finite controls](session27-orthogonal/ORTHOGONAL_AVOIDANCE.md); [audit](session27-orthogonal-audit/AUDIT.md).

## 10. What remains, and what to pursue next

The useful common question is not“do these numbers match?” but“which response experiments preserve enough information to reconstruct the object or certify the predicate?” This session answers it positively for marked tournament response and orthogonal gluing, and exhibits concrete failures for score roots, ambient-to-restricted theta inference, and unqualified modular analogies.

Next obligations:

1. Decide whether two ports always reconstruct a tournament, or find the first counterexample; full finite-context reconstruction is now proved.
2. Exploit the recovered endpoint/cover matrix for verified single-arc and score-preserving cycle-flip analysis before claiming enumeration speedups.
3. Investigate rarity among STRONG classes as h grows; padding-degree laws do not answer that different question.
4. Compute or certify the actual forbidden-line set before applying the rank96 avoidance criterion.
5. Retain phase, owners, integer clocks and geometric realization when returning to LRC, Collatz or unit-distance coloring.

The Smale phrase was provisionally interpreted as h-cobordism pending clarification. The [topological scope note](TOPOLOGY_SCOPE.md) gives a minimal obstruction to forgetting tournament orientations before applying ordinary topology; it supplies no new h-cobordism theorem.
