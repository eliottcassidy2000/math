# Euler-aware odd-order support sieve

**Session 177 research result.** A proved finite-support reduction and an exact realization theorem for its graph relaxation. This does not settle odd perfect numbers, improve their known prime-count bounds, or establish historical novelty. The main practical gain is eliminating impossible prime supports without selecting exponents or factoring divisor sums.

## 1. Setup

Let S be a finite set of odd primes and choose q in S with q = 1 (mod 4). A hypothetical odd perfect number supported on S has Euler form

N = q^alpha product_{p in S, p != q} p^{a_p},

where alpha = 1 (mod 4), and each a_p is positive and even. For distinct p,r define d_r(p) = ord_r(p), the multiplicative order of p modulo r. Construct the directed graph H_q(S) by the rules

- if p != q, put p -> r exactly when d_r(p) is odd;
- if p = q, put q -> r exactly when v_2(d_r(q)) <= 1.

There are no loops. This graph is computable entirely from S and q. One need not factor any sigma value, nor even factor r-1: writing r-1 = 2^t u with u odd, the two tests are p^u = 1 (mod r) and q^{2u} = 1 (mod r), respectively.

## 2. Main theorem

**Theorem.** If an odd perfect number has support S and Euler prime q, then H_q(S) is strongly connected. Moreover, H_q(S) is the exact maximal graph obtainable on S by choosing exponents of Euler parity: there exist alpha = 1 (mod 4) and positive even a_p such that

r divides sigma(p^{a_p}) if and only if p -> r belongs to H_q(S)

for every distinct p,r in S, taking a_q = alpha.

The realization statement concerns the induced graph on S only. The sigma values may have prime factors outside S, their valuations need not balance, and the resulting N need not be perfect. This limitation is essential.

### Proof: local order classification

Let n=a+1. If p is not 1 modulo r, then

r | sigma(p^a) if and only if d_r(p) | n,

by the geometric-sum identity and the invertibility of p-1 modulo r. If p is 1 modulo r, the correct criterion is instead r | n, because sigma(p^a) is n modulo r. This denominator exception must not be replaced by the vacuous order-one condition.

For an even exponent a, n is odd, so every order dividing n is odd. For an Euler exponent alpha = 1 (mod 4), n = alpha+1 is 2 modulo 4, so any order dividing n has 2-adic valuation at most one. The denominator exception is compatible with both cases since r is odd. Thus every actual sigma edge lies in H_q(S).

### Proof: strong connectivity of the actual graph

Use the actual exponents to form the graph with p -> r when r divides sigma(p^{a_p}). Suppose it is not strongly connected. Its condensation has a nonempty proper sink component T. Set A = product_{p in T} p^{a_p} and B=N/A. Since T has no outgoing edge, every odd prime factor of sigma(A) divides A. Also sigma(A) divides sigma(N)=2N, so its odd prime valuations do not exceed those of A, and its 2-adic valuation is at most one. Hence sigma(A) divides 2A.

But A>1 and sigma(A)>A. The only divisor of 2A larger than A is 2A itself. Thus sigma(A)=2A. Multiplicativity and sigma(N)=2N now give sigma(B)=B, impossible since B>1. Therefore the actual graph is strongly connected. Adding edges preserves strong connectivity, so H_q(S) is strongly connected.

This strong-connectivity ingredient is not claimed as new; see the source discussion below.

### Proof: simultaneous exact realization

For each p != q, let L_p be the least common multiple over all outgoing edges p -> r of

- d_r(p), if d_r(p)>1;
- r, if d_r(p)=1.

Every such requirement is odd, so L_p is odd. If L_p=1, replace it by 3 to ensure a positive exponent. Set a_p=L_p-1. This is positive and even. Every allowed outgoing edge is realized by the local criterion. No disallowed edge can appear, since its order is even.

For q, take L_q to be the least common multiple of 2 and the same edge requirements. Its 2-adic valuation is exactly one, because every allowed order has valuation at most one and the exceptional requirements r are odd. Set alpha=L_q-1, which is 1 modulo 4. Every allowed q edge is realized. Any disallowed edge requires an order divisible by 4 and cannot appear.

Choices at different source vertices are independent. This proves the simultaneous statement. QED.

## 3. A whole-support exclusion not detected by simple degree tests

Take S={3,5,23,31,43}. The sole possible Euler prime is q=5. Its maximal graph has edges

3 -> 23;
5 -> 3,23,31,43;
23 -> 43;
31 -> 3,5,23,43;
43 -> 3.

Every vertex has positive indegree and positive outdegree. In particular, 5 receives a permitted edge from 31, so a test merely asking whether the Euler prime can receive a factor does not reject this support. The usual abundancy ceiling also permits it:

product_{p in S} p/(p-1) = 30659/14784 > 2.

Nevertheless {3,23,43} is a proper sink strongly connected component: 3 -> 23 -> 43 -> 3. It cannot reach 5, so the theorem excludes this support for **every** allowable exponent choice.

The exact graph is realized by exponents

(a_3,a_5,a_23,a_31,a_43)=(10,461,20,1154,2).

These are not exponents of a perfect number. They demonstrate that the obstruction comes from the global cut, rather than impossible local edge requirements. The example has only five primes and is already ruled out by known general prime-count bounds; its role is to isolate what the new formulation detects, not to improve known numerical bounds.

## 4. Weaker Euler-independent sink filter

Let G(S) contain only odd-order edges for every source, including q. In any actual sigma graph, a path to q can be truncated at its first visit to q. All earlier sources are non-Euler, so that entire truncated path lies in G(S). Consequently every p in S must reach q in G(S).

Equivalently, G(S) must have a unique sink strongly connected component, and q must belong to it. This is a useful first pass: compute G(S) once, discard all Euler candidates outside its sink, and only then test H_q(S) for the remaining q. The first pass is necessary, not asserted sufficient for H_q(S) to be strongly connected.

## 5. Reciprocity, mod 12, and the limitation of motifs

Every odd-order edge p -> r implies that p is a square modulo r: p lies in an odd-order subgroup of the cyclic group of order r-1. Thus for two primes p,r both 3 modulo 4, quadratic reciprocity forbids simultaneous odd-order edges in both directions. More precisely, because v_2(r-1)=1, odd order is equivalent to quadratic residuosity at such a target; reciprocity selects exactly one of the two directions. The odd-order graph induced on primes 3 modulo 4 is therefore a tournament.

Primes 11 modulo 12 fall in that same class. This is a genuine structural bridge between the mod-12 work and divisor-sum chains, not a justification for classifying all primes by 12k+-1: primes 5 and 7 modulo 12 remain essential.

There is no claim that tournaments alone forbid an OPN. The parallel Session 177 tournament-realization construction shows how flexible such induced graphs can be. The usable restriction here is **global access to the Euler vertex**, combined with eventual exact valuation balance and exclusion of factors outside the proposed support.

## 6. Exact checks and research direction

Run `python check_odd_order_support.py`. The standard-library checker:

1. verifies simultaneous exact realizability for 648 support/Euler-prime cases, using subsets of the odd primes below 30 of sizes two through five;
2. performs 2,952 independent direct geometric-sum comparisons for exponents 1 through 41;
3. verifies the example's sink, degree tests, Euler parity, rational abundancy ceiling, and realizing exponents.

The finite checks corroborate implementation; the proofs above establish the unbounded statements.

A useful next synthesis is to place this support sieve **before** expensive exponent branching. A second layer should retain the exponent congruences L_p | (a_p+1), then integrate exact incoming valuation capacities. Merely adding edges or treating a realizable induced graph as a closed sigma network is insufficient. The outstanding task is a theorem linking the forced LCM growth to valuation overflow or an unavoidable new prime outside S.

## 7. Literature and novelty boundaries

- Fletcher, Nielsen, and Ochem, *Sieve Methods for Odd Perfect Numbers*, gives the Euler form and develops factor-chain and sieve restrictions: https://mathdept.byu.edu/~pace/OPNSieves_web.pdf . Our result is a support-level preprocessing formulation, not a replacement for those methods.
- Bibby, Vyncke, and Zelinsky, *On the Third Largest Prime Divisor of an Odd Perfect Number*, studies mutual sigma divisibility and sigma_(m,n) pairs: https://arxiv.org/abs/1908.09420 . The present note does not claim that sigma graphs, mutual divisibility, or order-based restrictions originated here.
- Zhu and Li, *Valuation Networks, Order-Filtered Primitive Excess, and Computer-Assisted Obstructions in the Odd Perfect Number Problem*, posted August 18, 2026: https://zenodo.org/records/22005604 . Its record explicitly lists strong connectivity and exact cut identities; it is marked not peer reviewed. We give an independent elementary proof of strong connectivity above and do not rely on its computational claims.

No exhaustive literature-priority search has been completed. The specifically delivered contribution is the proved Euler-aware maximal graph, its simultaneous realization including the order-one exception, the two-stage sink algorithm, and the exact support example. Treat historical novelty as unestablished.
