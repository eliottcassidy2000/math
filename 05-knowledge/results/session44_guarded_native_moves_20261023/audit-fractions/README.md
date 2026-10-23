# Independent audit of observer fractions and ternary continuation

23 October 2026. [Audited producer packet](../fractions-ternary/README.md).
**PASS, scoped.** The producer's finite weighted-incidence
theorem, order4/order5 spectral certificates, reflected-converse theorem and
native Jacobsthal continuation law withstand independent mathematical and
exact computational checks. No all-order native gap classification or
original electrical-guard factorization is inferred.

## Findings first

1. The trace defect really is a nonnegative weighted profile variance.
   Vanishing is equivalent to the full factorization K=J Kbar E, not just
   lumpability or an equality of ranks. Equality of nonzero spectra implies
   equal traces and therefore forces the same condition. Both directions
   and multiplicities are justified.
2. Independent literal arc-set reconstruction gives the same order4 and
   order5 incidence tables. A different characteristic-polynomial route,
   SymPy's exact Berkowitz algorithm, confirms every full, compressed,
   joint-refinement and complement-paired matrix in the frozen receipt.
   Independent finite-field and rational factorization confirms algebraic
   degrees 7, 7 and 5; exact root counts confirm all stated largest-root
   brackets and the strict reversed gap comparison.
3. The all-order converse argument is sound because it uses the actual
   reflected-converse involution. Native chord (i,j+1) becomes
   (n-2-j,n-1-i), exactly the image of compressed edge (i,j) under
   i->n-2-i. Ordinary source classes remain unchanged while the tournament
   class becomes its converse. Equal weighted columns follow only if the
   accepted word set and weights are invariant under that bijection.
4. The Jacobsthal first-step law follows from the exact binary-residue
   counts. The stated halves, reduced unequal pairs, determinant-two
   continued fraction, A-spine depth and ternary-heap address are correct.
   Neither 7/15 nor 8/15 belongs to these unguarded continuation pairs.
5. Two prose clarifications were requested and accepted, with no change
   to the mathematics or computation: cancel the common factor 2 before
   describing the displayed unequal probabilities as reduced; and explain
   that an arbitrary remaining clock r is a continuation suffix inside a
   sufficiently large odd-cell regular host, not necessarily its full
   (even) number of phase differences.

## Independent test universe and methods

The audit imports none of the producer's code. It reads the pinned producer
receipt only as a target for comparison.

- Rebuild every one of the 8 order4 and 64 order5 words directly as sets of
  native arcs. Canonical graph labels are found by literal permuted edge
  sets; canonical tournament labels are found by literal permuted arc sets.
  This avoids the producer's affine bit-permutation implementation.
- Build all eight incidence matrices from these independent labels and
  compute their characteristic polynomials by Berkowitz, rather than the
  producer's Newton trace recurrence.
- Use independent exact symbolic factorization over Q and over the stated
  finite fields, and exact root counting in all isolating intervals.
  No approximate eigenvalue is used to establish a degree or inequality.
- Exhaust all 3-by-2 nonnegative incidence matrices with entries in {0,1,2}
  and positive marginals, merging rows 0 and 1. Check the compression
  equation, trace/variance equality, full factorization, and equality of
  nonzero spectra. This is evidence for the elementary all-size proof,
  not a replacement for it.
- Exhaust every fixed-path word at native orders 2 through 6, checking
  both native and compressed reflected coordinates literally.
- On actual regular hosts with 3, 5 and 7 C3 cells, enumerate **all** 3^t
  phase tuples rather than generating only already-accepted binary words.
  Verify native Hamiltonian validity iff every successive difference is
  0 or 1; then recover every endpoint-conditioned first-step count.
- Independently sum binomial residue classes at every remaining clock
  r=3,...,120 and every endpoint residue. Check reduced pairs, determinant
  J_r, Cassini's determinant-two identity and primitive adjacent ratios.
  Symbolically verify that the inherited A map sends k/(2k+1) to
  (k+1)/(2k+3), proving the address rule for arbitrary k, not just samples.

The independent checker uses installed SymPy 1.14.0 and exact integers and
rationals. It has no optimization-sensitive assert statements. Normal and
optimized outputs byte-match the frozen audit receipt, with **7,992** exact
checks in each recorded run (plus its replay equality check).

    python3 -B check.py --producer <fractions directory>
    python3 -O -B check.py --producer <fractions directory>

`source_pins.json` identifies the exact producer README, checker and receipt
bytes audited. `manifest.json` pins this separate report, independent code,
receipt and source-pin file. The producer's own normal/-O tests remain an
additional, distinct replay path.

## The theorem/evidence boundary retained

The finite audit verifies orders4 and5, not higher-order native spectra.
The symbolic compression identity and reflected-converse proof are all-size
within their hypotheses. The Jacobsthal law is exact for the stated native
phase transfer and its suffix continuations, not for arbitrary electrical
filters. Shared endpoints, the common-phase lift, actual remaining clock
and integer displacement are not reconstructed from a rational-tree address
alone. Positive incidence rank and rational matrix entries do not imply a
rational gap. Nothing here identifies a Farey tree with a native move graph.
