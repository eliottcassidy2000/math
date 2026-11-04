# From local tournament freedom to global divisor balance

Session 177, 2026-11-04. **Status:** synthesis of proved scoped results, classical constructions, and explicit open proof targets. Odd perfect-number existence and infinitude of amicable pairs remain open. Historical novelty of the new formulations has not been established.

## 1. The useful connection is a sequence of retained constraints

The number 12 supplies a real arithmetic interface: its unit group separates the four prime classes 1,5,7,11 modulo 12, while reduction modulo 4 controls quadratic reciprocity and the Euler prime. This is distinct from an analytic occurrence of 1/12. In particular, primes greater than 3 are not all of the form 12k±1; the classes 5 and 7 cannot be discarded.

The session's central result is a separation, with constructive content on both sides:

1. **Local induced graph:** any finite tournament can be realized by divisor-sum transfers among selected primes all congruent to 11 modulo 12, with one common positive even exponent.
2. **Euler-aware support:** a proposed complete support must pass a stronger order graph and global connectivity test.
3. **Exact neighborhoods:** even allowed edges cannot generally be selected independently; their order labels impose LCM implications.
4. **Closed valuation network:** every divisor-sum prime and its full multiplicity must fit the proposed support and its exact exponent budgets.
5. **Perfect or amicable equation:** the original divisor-sum identity must finally hold. Graph realizability is not a substitute for this identity.

These levels identify what a graph-theoretic argument can preserve and what it must restore. They do not identify the existence of odd perfect numbers with the existence of any arbitrary tournament.

## 2. What the local freedom theorem actually says

The [mod-12 packet](../mod12/NOTE.md) proves that for any finite labelled tournament T there are distinct primes p_i≡11 modulo 12 and a common positive even exponent a such that

    i -> j in T  iff  p_j divides sigma(p_i^a).

One may prescribe either even exponent class modulo 4 and any common lower bound on the valuations along the displayed arcs. The construction uses CRT and Dirichlet to choose Legendre symbols, reciprocity to orient pairs, and multiplicative orders plus lifting to choose exponents.

The key denominator exception is retained: if p≡1 modulo r, then r divides sigma(p^a) exactly when r divides a+1. An order-one label alone is insufficient.

This theorem concerns only the graph **induced on the selected primes**. It leaves all exterior prime factors, exact valuation equalities, and abundancy uncontrolled. In particular it does not show that every finite tournament is a closed component of a perfect number. The explicit local cycle on 11,83,179 in the companion packet has exterior divisor-sum factors. Conversely, a theorem excluding finite tournament shapes alone cannot supply a universal obstruction to these unrestricted induced transfers.

The 12-class bookkeeping has a separate exact use. In a 3-free hypothetical odd perfect number, a nonspecial prime p≡1 modulo 3 has exponent a≡0 or4 modulo 6, not2. Thus the exact component 11² forces 7⁴·19⁴ through sigma(11²)=7·19, whereas the component 13² is impossible because sigma(13²)=3·61. The prime 13 itself is not excluded: it can occupy another exponent or the Euler role. These are conditional exponent costs, not a proof of nonexistence.

## 3. Which global restrictions survive local universality

Write the classical Euler form as N=q^alpha product p^{a_p}, where q≡alpha≡1 modulo 4 and the other positive exponents are even. The [Euler-aware support packet](../oddperfect/odd_order_support_theorem.md) assigns to a finite prime support S and candidate Euler prime q a maximal graph H_q(S):

- ordinary source p permits target r when ord_r(p) is odd;
- Euler source q permits target r when the order has 2-adic valuation at most one.

Every actual transfer lies in this graph. The packet proves that a perfect-number support requires H_q(S) to be strongly connected, and that all its permitted edges can be realized simultaneously with suitable Euler-parity exponents. The latter is again an induced-graph realization, not perfection.

A cheaper preliminary graph uses odd-order edges at every source. Every prime must reach q in this graph; equivalently it has a unique sink strongly connected component containing q. This condition is weaker than strong connectivity of H_q(S).

The illustrative support {3,5,23,31,43}, with sole possible Euler prime5, passes positive indegree/outdegree and the ordinary abundancy ceiling, but fails the global sink test. This example explains the additional information retained by connectivity. It does not improve general prime-count bounds: a five-prime odd perfect number is already excluded by known results.

The [LCM closure supplement](../oddperfect/lcm_neighborhood_closure.md) then classifies exactly which subsets of the allowed targets are realizable neighborhoods. Label target r by k_r=ord_r(p), except replace order1 by k_r=r. Include the role requirement b=1 for ordinary sources and b=2 for Euler sources. A target subset T is realizable exactly when it contains every allowed r whose label divides lcm(b,{k_t:t in T}).

For example, simultaneously requiring 11 and13 to divide sigma(3^a) forces 4561 to divide it, since the three orders are5,3,15. Either first requirement alone does not force the third. This is a conjunction constraint, invisible in an unlabelled directed graph. The theorem accounts separately for the empty ordinary neighborhood, so it does not introduce an invalid zero exponent.

## 4. Exact cuts turn graph paths into arithmetic obligations

The [network packet](../networks/README.md) proves strong connectivity of the actual transfer graph by a short divisor argument. A proper sink block A would satisfy sigma(A)|2A and sigma(A)>A, forcing sigma(A)=2A. Its nontrivial coprime complement B would then have sigma(B)=B, a contradiction. Strong connectivity alone is insufficient: 120 has a strongly connected transfer graph and sigma(120)=3·120.

For a full prime-power block A of a hypothetical odd perfect number, let O collect outgoing odd-prime factors of sigma(A), and I the factors supplied into A from its complement, with multiplicities. If epsilon records whether A contains the Euler component, then

    O/I = sigma(A)/(2^epsilon A),       gcd(O,I)=1.

The reduced fraction determines both O and I. This is stronger than knowing that some edge crosses each cut. Nevertheless a valid cut ratio is not by itself a completed perfect number: all blocks and exponents must be mutually consistent.

For a partial assignment of **exact** components, the residual

    R = sigma(A)/gcd(sigma(A),2A)

must divide the coprime complement. If A contains the Euler component, that complement is a square, so R forces its exponents rounded upward to even values. The resulting abundancy supplies an exact rejection test. Exactness is indispensable: sigma(3²)=13 but sigma(3⁴)=121, so increasing an exponent can remove an old transfer edge.

The packet's eleven conditional fourth-power clauses demonstrate usable, finite conclusions from these budgets. Examples are

- 17||N implies 3⁴|N or13⁴|N;
- 5||N implies at least one of3⁴,13⁴,61⁴ divides N;
- 53||N implies3⁴|N.

The larger table covers each prime q≤97 with q≡1 modulo4 under q||N. Each clause explicitly assumes absence of the listed fourth powers when freezing forced ordinary exponents at2. Its contradiction therefore proves a disjunction, not an all-exponent exclusion.

## 5. A constructive neighboring problem: two-cycle balance

The [amicable packet](../adjacent/AMICABLE_DEFICIENCY_BRIDGE.md) supplies an exact constructive comparison. Perfect numbers satisfy s(N)=N; amicable pairs satisfy s(N)=M and s(M)=N with N≠M. For the particular ansatz N=apq, M=ar with distinct primes p<q, prime r, and gcd(a,pqr)=1, put

    S=sigma(a), D=2a-S, E=S-a.

The amicable condition is equivalent to D>0 and positive factors d<e satisfying

    de=a²,  d≡e≡-E modulo D,
    p=(d+E)/D, q=(e+E)/D, r=pq+p+q,

together with all prime and coprimality guards. The factorization is Euler's classical construction, presented here with unnormalized deficiency coordinates and a reproducible checker. It is not a new amicable-number theorem.

For p,q>3 the last cofactor must satisfy r≡11 modulo12. This is another genuine appearance of the same residue class, but it comes from r+1=(p+1)(q+1), not from the tournament universality mechanism. No equivalence between the two constructions is asserted.

The examples (220,284), (122265,139815), and (69615,87633) are verified directly. The latter two are odd–odd amicable pairs. Their existence is established, whereas infinitude of amicable pairs remains unproved. Each fixed kernel a has only finitely many divisor pairs; a local residue sieve cannot produce the missing infinitude theorem.

## 6. Specific next proof targets

These are proposals, not results of this session.

1. **Order closure forcing exterior support.** For a specified family of proposed supports, combine a required return path to q with the LCM closure rules and prove that at least one required source neighborhood forces a prime outside the support. A useful theorem must cover the entire stated support family, not just enumerate examples.
2. **Forced valuation cost with shared suppliers.** Combine the 3-free exponent costs with exact cut fractions. Seek an inequality that remains valid when several required paths reuse the same prime or supplier. Counting edges as independent resources would not suffice.
3. **Cover higher-exponent branches.** Extend a selected fourth-power clause by analyzing each of its escape branches with the new exact exponents or a proved symbolic exponent family. The target is a genuinely covering disjunction or a terminating branch argument, not repeated use of the squared-component formula after its exponent changes.
4. **Turn graph preprocessing into a transparent filter.** For an explicitly chosen bounded support universe, compare the ordinary sink test, H_q strong connectivity, LCM neighborhood closure, and exact-component residual rejection. Such a result would measure which constraints do work; it would remain a finite classification rather than an OPN existence proof.
5. **Constructive amicable families.** Seek an unbounded family of deficient kernels with explicit divisor-pair parametrizations and fully retained prime guards. An infinitude conclusion would additionally require a theorem producing infinitely many successful prime triples. An algebraic parametrization or congruence admissibility alone does not supply it.

The most promising shared direction is therefore not “exclude a tournament.” It is “force an arithmetic debt across a cut that cannot be paid after every exponent-dependent change of edges is allowed.” The papers here provide exact interfaces for stating that debt; they do not yet prove a global contradiction.

## 7. Evidence and attribution

Proofs and complete validation scopes belong to the linked sibling packets. Their finite checks corroborate implementations; the stated universal lemmas have separate proofs. The synthesis imports no new record bound and asserts no resolution of an open problem.

Primary references:

- Euler, *De Numeris Amicabilibus* (1750), §27, in Jonathan David Evans's translation: https://arxiv.org/html/2409.08783v2 . Exact classical precedent for the amicable factorization.
- Fletcher, Nielsen, Ochem, *Sieve Methods for Odd Perfect Numbers*: https://mathdept.byu.edu/~pace/OPNSieves_web.pdf . Euler form and established sieve/factor-chain context.
- Bibby, Vyncke, Zelinsky, *On the Third Largest Prime Divisor of an Odd Perfect Number*: https://arxiv.org/abs/1908.09420 . Existing work on mutual sigma divisibility.
- Zhu, Li, *Valuation Networks, Order-Filtered Primitive Excess, and Computer-Assisted Obstructions in the Odd Perfect Number Problem*: https://zenodo.org/records/22005604 . Unreviewed August2026 preprint already advertising strong connectivity and exact cuts. The sibling packets give self-contained proofs of the facts they use and import none of its computational bounds.
- Pomerance, *On Amicable Numbers*, manuscript dated March30,2014: https://math.dartmouth.edu/~carlp/amicable032914.pdf . Sparsity and the open infinitude question.
- Pollack, *Amicable Numbers*, October2026 preprint: https://arxiv.org/html/2610.05661v1 . Stronger sparsity bound; infinitude explicitly remains open.
- MIT notes on [quadratic reciprocity](https://ocw.mit.edu/courses/18-781-theory-of-numbers-spring-2012/resources/mit18_781s12_lec9/) and [Dirichlet's theorem](https://ocw.mit.edu/courses/18-785-number-theory-i-fall-2021/mit18_785f21_full_lec.pdf) identify the classical inputs to the local universality construction.
