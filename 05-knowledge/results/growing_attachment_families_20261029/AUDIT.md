# Proof and computation audit

2026-10-29. Same-agent review, not independent-agent certification.
**PROVED scoped deductions from CITED inputs / FINITE-EXACT controls.**

## Proof gates

1. **P1 rank and realization.** Rank is polynomial squareclass rank after
   removing constants; nonsquare constants translate the answer code.
   Every nonzero basis combination has a nonconstant squarefree
   representative, so the cited Weil hypothesis holds. Products may
   acquire repeated factors, but parity reduction is justified only off
   the complete deleted set. D ignores even-only factors, Z does not.
   Rank 0 is treated separately. The claim covers all finite banks; it
   does not claim every algebraic word occurs without the lower bound.
2. **P2 new roots.** Each earlier actual vertex lies outside Z, so its
   linear factor is coprime to the bank. Distinctness gives exactly t new
   independent factors, not merely an upper bound on rank. Sequential
   counts retain their common actual frame. Multiplication of uniform
   lower/upper extension bounds proves the joint bound; no assumption
   of independence is used. Root-root incidences are held fixed.
3. **Growth budget.** The per-step relative error against q/2^(r+t)
   is bounded by a geometric factor times (D+z+t+1)/sqrt(q), plus deleted
   points/q. Summing over t=h,...,h+k-1 proves the stated sufficient
   budget. Summing equal-main-term full patterns also proves Boolean
   closure without enlarging relative error. Arbitrary identity-based
   predicates are outside that answer space.
4. **P3 universal codes.** Choose r distinct roots and a full-rank matrix;
   distinct linear prime divisors certify independence algebraically.
   The explicit positive lower bound, not a generic rank claim alone,
   proves that every word occurs. Constant factors implement affine
   shifts. Code image and square-product kernel are dual. Characteristic
   two is excluded. Odd extension degree preserves base constant signs;
   Frobenius transports the entire named source together.
5. **P4 quantifiers.** Fixed A precedes X; the bank, code and predicates
   may then vary subject to Q<=log^A X. Summing AP errors costs at most
   phi(Q), absorbed by an arbitrarily large logarithmic saving. This is
   exactly where fixed-modulus PNT was insufficient. Siegel--Walfisz
   constants are ineffective. The spectral pattern estimate is inherited
   with its diagonal and collision correction. Multiplication by 2^m
   is necessary for relative, rather than absolute, pattern error.
6. **P5 specialization.** Direct substitution in the inherited six defect
   polynomials gives the four displayed classes. Exclude 2,3,5,23: all
   three quadratics are then squarefree and pairwise coprime. The constant
   shift has even parity. Source t=0 is removed for noncollinear centers;
   the statement remains a discriminant-character law, without asserting
   eight distinct circles or closing the older full geometric problem.
7. **P6 class normalization.** 81 and 27 are inherited automorphism orders.
   The current script only checks their factor-three orbit normalization.
   It does not recompute their groups, paths or response modules. Colored
   or rooted class normalization uses the appropriate stabilizer.
8. **P7 arithmetic sources.** Every intermediate value is an exact affine
   function of the same integer cylinder parameter. A fixed choice of
   plus or minus gives a nonzero same-sign segment intercept, hence every
   determinant is nonzero. The excluded primes include 2,3 and every root
   collision. Sampling t through one full prime residue interval is an
   explicit ensemble of legal integer trajectories, not a claim about
   temporal distribution. A finite collection has one exceptional-prime
   certificate; no product-independence between different words follows.

## Hostiles and strongest survivors

* Squared factors add no bit but add excluded sources. All (X-a)^2 over
  the whole field leave an empty domain. Survivor: the free-query claim
  is valid for products supported on a retained finite divisor set.
* A root at a bank zero does not necessarily add rank. X and the root0
  incidence are the same sign. Survivor: the coprimality/zero-set guard
  in P2 makes rank additive exactly.
* At q=7 the Apollonius rank 3 family realizes only one of its eight formal
  words. This is compatible with the estimate, whose lower bound is
  nonpositive there. Algebraic code size is not finite realized size.
* In F27 moving x by Frobenius while leaving a=3 fixed changes the arc
  for x=5. Simultaneous transport preserves every arc and every response.
* Characteristic 23 identifies C with A in P5. The factor-parity rank
  must be recomputed there; the generic eight-word conclusion is withheld.
* The 81/27 pair still distinguishes structural response data that this
  static answer carrier never promised to reconstruct.

No established false theorem was discovered or left live, so no new
MISTAKES entry is needed. Failed tool lookups were navigation failures,
not mathematical evidence. Current result packets take precedence over
historical claims about unresolved prefix equidistribution.

## Reproduction and provenance

The verifier uses integer arithmetic and explicit exceptions under both
normal and optimized Python. Irrational envelopes are tested by moving
their rational term to the other side and squaring the remaining
nonnegative quantity. Joint product envelopes use the rational outer
bound ceil(sqrt(q)); exact sequential and literal count dictionaries agree.
F27 arithmetic is F3[u]/(u^3-u+1), an irreducible cubic (no F3 roots).
Euler's criterion is independently compared with the set of literal
squares, and Frobenius is checked against addition, multiplication and
all arcs. No floating-point spectral calculation is used.

Receipt SHA256 (normal and optimized):
`add2b3ce361f01307272be887f3d54a3c201dcc5c681c3c0e73143aa438c4b95`.
The receipt declares every finite universe. Source digests are in
`sources.json`; they are provenance pins rather than a sandboxed replay
manifest. Published inputs are Weil's bound, the classical Paley extension
estimate, Siegel--Walfisz and the already audited quadratic large sieve.
The general methods are standard. The contribution here is their explicit
family interface, budgets, source-preserving word application and recovered
discriminant specialization, with no claim of historical novelty.
