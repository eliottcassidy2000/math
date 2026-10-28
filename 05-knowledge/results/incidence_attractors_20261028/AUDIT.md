# Proof and computation audit

This is same-agent adversarial review with independently implemented
arithmetic checks. It is not an independent-agent or external peer audit.

## Types and quantifiers

- F,U,S are counts on1<d<n; sigma is the weighted sum including1,n.
  F=U+S is a Boolean predicate in our87-state model, not the numeric defect.
- B has seven labelled atoms but only six prime supports. The forbidden pair
  is2/4. Its parity-incidence graph vertices0,2,3,7 are a different set from
  the atom labels. Repeated support cannot be ignored because sigma is not
  completely multiplicative.
- The96 states are legal subset products only; operations outside the seven
  named multipliers are absent. Failed requests return a visible rejection
  and keep the current state. The86/87/96 counts depend on those conventions.
- Complete subset context tables suffice because every successful request
  consumes a fresh prime support. Rejections are visible self-loops; successful
  operations commute in this universe. State minima are distinct from ranks.
- The Fano bank consists of seven distinct primes, all verified. Its code
  covers all128 subsets. The non-graphic proof concerns cycle matroids;
  it does not deny other useful graphs representing finite incidence systems.
- The polynomial basis, normal basis, Singer order and orientation are
  explicitly chosen. The trace-zero orbit(2,4,6) fails the basis test, while
  (3,5,7) passes. Frobenius preserves Boolean acceptance, not sigma magnitudes.
- The Fano automorphism order168 is both elementary and finite-exact; the
  chosen Paley H189/Aut21 counts are finite-exact and inherited classical
  controls. The original81/27 packet remains a separate, pinned control.
- The minimum364 is exhausted below a witness; the context minimum205 is
  exhausted through its witness. Neither requires an assumed asymptotic bound.
- Divisor-height1 is terminal, not assigned a finite value of h_N(1).
  Its Monster provenance is not needed to verify the displayed integer N.
- Berggren's Q=0 and Q=1 shells share branch addresses, not their norm values.
  The local absorber checks actual positive distinct vertices. It supplies
  neither a prescribed Hamiltonian host nor a general graph existence proof.
- Petersen local lifting includes two runs through a cluster and separately
  tests connectivity. K3,3 keeps the three suppressed-edge marks. Graceful
  overwrite connectivity is not the Graceful Tree Conjecture.

## Exact universes and independent paths

`check.py` computes sigma from prime factorization and checks:

1. All128 subsets of B, the96 legal states, the16 formal kernel elements,
   the12 legal square products, every minimal circuit, and every transition
   until partition refinement stabilizes.
2. All integers2..364 for the p²qr minimum and all common legal multipliers
   1..205 for the4/13 square-output separation.
3. All128 Fano-bank subsets; all5040 permutations for each of the Fano and
   Paley automorphism tests and the Paley path count; all8 normal coordinates.
4. All5759 nonunit divisors of the displayed N by repeated integer divisibility.
5. All364 words of length0..5 in the three Berggren matrices, in both shells.
6. All1024 K5 edge subsets for the local Petersen lift, independently compared
   to all binom(15,10)=3003 native candidate2-factors; all rooted cyclic orders
   of K3,3 for its six cycles and two marked lifts.

`audit.py` imports no producer functions. It independently adds every divisor
into its multiples through60060, then checks all96*128 residual entries and
both extended observer partitions. The maximum legal bank product is60060;
noncoprime contexts are rejected before array access. It verifies all128
Fano products by direct divisor-pair enumeration, independently of the parity
matrix. It also directly enumerates all positive norm-one triples with c<=150
and verifies their unique descent to one of the three inherited roots.

Every check uses an explicit exception, not a removable assert. Producer and
audit each run normally and under-O, with equal JSON receipts. `verify.py`
checks receipts, immutable source pins, local links, the sealed file set, and
optionally replays all four runs. No random sample substitutes for a complete
finite universe. The source receipts retain the full list of merged states.

## Literature and novelty boundaries

Primary pages checked: Lovasz's connection-matrix abstract; Courcelle–Weil's
recognizability abstract; Romik's modified-Euclidean/tree-address abstract;
Utrecht's bibliographic record for Beukers–Luca–Oort. The Utrecht PDF endpoint
did not load in the browser tool. We therefore do not claim a new independent
reading of its full proof; the arithmetic parity argument used here is proved
directly, and the prior Session82 attribution is retained.

No universal graph-homomorphism representation theorem is applied to the new
carrier. No Hamming-code, Fano-plane, Myhill–Nerode, or Berggren novelty claim
is made. The finite arithmetic realizations and operation-relative comparisons
are the present deliverables. Unbounded prime realization, unrestricted
tournament composition, Collatz/aliquot termination, and LRC(14) remain outside
the conclusions.
