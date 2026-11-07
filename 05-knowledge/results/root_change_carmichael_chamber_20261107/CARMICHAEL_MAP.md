# Carmichael's totient conjecture: literature, repository work, and an exact mechanism dictionary

kind-pasteur-2026-11-07-S14. Companion to [PROOF.md](PROOF.md) (THM-4904).

This note answers the owner's request to investigate prior work on the Carmichael conjecture in depth. It separates
three layers:

1. what the number-theory literature proves;
2. what this repository already did;
3. which mechanisms transfer to the root-change word grammar, and which do not.

The four objects called "Carmichael" (totient conjecture, Carmichael numbers, Carmichael lambda, and the repo's Camion
pseudoprimes) are kept apart exactly as in Session218's
[provenance note](../session218_smith_repunit_totient_20261107/provenance/PROVENANCE.md). Only the totient conjecture
is discussed here.

## 1. The conjecture and its literature

Write A(m) = #{x : phi(x) = m} and V(x) = #{m <= x : A(m) > 0}. Let V_k(x) count the totients up to x with
multiplicity exactly k.

- **Carmichael's conjecture (OPEN).** A(m) is never 1. Carmichael announced a proof in 1907 (Bull. AMS 13, 241-243).
  The proof was flawed. He later showed that no m < 10^37 has multiplicity 1, and posed the statement as a conjecture.
- **Sierpinski's conjecture (THEOREM).** Every k >= 2 is a multiplicity.
  - Schinzel deduced it from Hypothesis H (1961).
  - Ford deduced it from the prime k-tuples conjecture.
  - Ford and Konyagin proved it unconditionally for even k.
  - Ford proved it for odd k: "The number of solutions of phi(x) = m", Ann. of Math. 150 (1999), 283-311.
- **Lower bounds for a counterexample.**

  | Bound | Source |
  | --- | --- |
  | 10^37 | Carmichael |
  | 10^400 | Klee (1947) |
  | 10^10,000 | Masai and Valette (1982) |
  | 10^10,000,000 | Schlafly and Wagon (Math. Comp. 1994) |
  | 10^(10^10) | Ford (1998), Theorem 6 |

  Ford's source is "The distribution of totients", Ramanujan J. 2 (1998), 67-151, arXiv:1104.3264.
- **Ford's dichotomy (Ford 1998, Theorems 2 and 5).** If A(d) = k, then V_k(x) >> d^(-1-eps) V(x). Also
  limsup V_1(x)/V(x) < 1, and Carmichael's conjecture is equivalent to liminf V_1(x)/V(x) = 0. In Ford's words, the
  set of counterexamples "cannot be a 'thin' subset": either there are none, or a positive fraction of totients are
  counterexamples.
- **Ford's other quantitative results.**
  - Theorem 7: liminf V_1(x)/V(x) <= 10^(-5,000,000,000), improving Pomerance's 1/2 (1991).
  - Theorem 8: for each k, the totients all of whose preimages are divisible by k either never occur or have positive
    proportion.
  - A unique preimage must be divisible by 6. So if no totient has all of its preimages divisible by 6, Carmichael's
    conjecture follows.
- **The Carmichael-Klee lemma (Ford's Lemma 7.2).** Suppose phi(x) = m and A(m) = 1. Suppose also that:
  - d divides x and gcd(d, x/d) = 1;
  - the squarefree kernel of phi(d) divides x;
  - e divides (x/d)/s(x/d), where s(n) is the squarefree kernel of n;
  - P = 1 + e phi(d) is prime.

  Then P^2 divides x. This lemma forces 2^2 3^2 7^2 43^2 | x. Carmichael's case split is
  3^2 || x, which also forces 13^2, against 3^3 | x. Every computational bound iterates the lemma to build an ever
  larger set of primes whose squares divide x. The elementary first step is the doubling twin: phi(x) = phi(2x) for
  odd x, so a unique preimage has 4 | x.
- **Pomerance's criterion (Proc. AMS 43, 1974).** If every prime p with (p - 1) | phi(n) satisfies p^2 | n, then n is
  a counterexample. Pomerance also showed such n are very improbable.
- **Erdos's propagation.** If A(m) = k, then A(m(p - 1)) = k for most primes p. This is the mechanism behind
  "positive proportion".
- **Other functions.** For sigma, A_sigma(3) = 1, so there is no Carmichael conjecture for sigma (Ford, section 8).
- **Recent work.** Shannon, Shiue, He and Saito (Notes Number Theory Discrete Math. 31(3), 2025, 504-534; Saito's
  UNLV thesis, 2025) verify the conjecture for totient values of special shapes: 2^t, 2p^s, 4p and 2pq. No new general
  bound or counterexample appears through 2026-11-07.
- **Nontotients.** The empty fibers are 14, 26, 34, 38, 50, ... (A005277).
  - Every odd number above 1 is a nontotient.
  - 2p is a nontotient exactly when 2p + 1 is composite (p > 3 prime).
  - An elementary infinite family: 2 * 7^k is a nontotient for every k >= 1. A prime p with p - 1 | 2 * 7^k is 2, 3,
    or 2 * 7^j + 1, and the last is divisible by 3. A power of 7 in x would put a factor 3 into phi(x).

  Session218 proves the Mersenne-Wagstaff family and the exact dyadic repair law. As Session218 stresses, empty fibers
  are not counterexamples to Carmichael's conjecture, which concerns singleton fibers.

## 2. Repository work on Carmichael-type fibers (read in full, credited)

| Where | Fiber notion | Exceptional (singleton) set | Status |
| --- | --- | --- | --- |
| HYP-9477 (A) / THM-4884 (opus-S81; scope repaired by Session218) | strong tournament classes with H = h | {3, 5, 11, 13} | FINITE-EXACT below 625 |
| Session200 `carmichael-fibers` | numerical source profiles s with F(s) = n (unbounded J) | {3, 11} | FINITE-EXACT to 50,001; all-order OPEN |
| this session, THM-4904 | root-change words with at most K letters J | K = 0: about 0.5 sqrt(X) values up to X; K = 1: sporadic; K = 2: none above 1 | PROVED structure, FINITE-EXACT counts |
| Session218 | phi itself: dyadic repair, Mersenne-Wagstaff empty fibers, Smith and repunit bridges | (empty fibers) | PROVED |

Session200 introduced the transfer: its O/J twin F(O(s)) = F(J(s)) = G(s) plays the role of phi(x) = phi(2x), and its
"nonbalanced chamber" plays the role of 4 | x. Its two-run chamber values are the Mersenne products
2 + (2^(a+1) - 1)(2^(b+1) - 1), with hostile diagonal 11, 51, 227, 963, 3971. HYP-9477 (A) records strong-class
singletons; Session218 withdrew a bounded-census Ford-type inference there (MISTAKE entry in its provenance folder).
The older `carmichael_*` packets of 2026-10-14..19, HYP-9358 and THM-4746 concern Carmichael numbers and Korselt
criteria, a different question.

## 3. Mechanism dictionary (proof mechanisms, not number-preserving maps)

| Totient side | Root-change word side (THM-4904) |
| --- | --- |
| doubling twin phi(x) = phi(2x), x odd | final-letter twin: wZ and wJ have the same value F(t.w) |
| a unique preimage has 4 \| x | a singleton word has exactly K letters J and its value lies outside F_<=K-1 |
| Carmichael-Klee: P^2 \| x forced for many primes | Session199's adjoint w -> w^dagger pairs the remaining words; a singleton word is dagger-fixed, of the form A c A^dagger, and (Session200's O/J twin) neither begins nor ends with J |
| Ford: unique preimage divisible by 6; Theorem 8 classes | singleton values lie in F_<=K minus F_<=K-1 (all their words have exactly K letters J) |
| Ford: no counterexample, or a positive proportion | singleton values number O_K(sqrt(X) (log X)^floor(K/2)): always density zero. Sharp for K = 0 (about 0.49 sqrt(X) observed) |
| Sierpinski: every k >= 2 occurs | M_K(N) = [N = 1] + Fix_K(N) (mod 2): odd multiplicities need an odd number of chamber words, a density-zero event |
| sigma has singletons; phi conjecturally none | K = 0 and K = 1 grammars have singletons; K = 2 has none through 10^9 (pre-registered P3) |

The dictionary transfers the shape of the argument, not numbers. As Session218 records, no map takes totient gaps to
Hamiltonian-path gaps: halving 14 and 26 gives 7 and 13, but 7 is excluded from the full spectrum while 13 is attained.
Inside the two-root-change grammar both 7 and 13 are missed. The first nine misses of F_<=2 are
7, 13, 21, 41, 43, 55, 67, 77, 83. Doubling them gives 14, 26, 42, 82, 86, 110, 134, 154, 166, of which five are
nontotients (14, 26, 86, 134, 154). This is unremarkable: about 30% of the even numbers in [14, 166] are nontotients,
and for a prime p > 3 the double 2p is a nontotient exactly when 2p + 1 is composite, which holds for most primes. It
is reported as **no fit**.
