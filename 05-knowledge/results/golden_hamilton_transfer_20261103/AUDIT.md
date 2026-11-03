# Scope and independent computational checks

**Same-agent proof audit, with separate combinatorial and computational
paths. No external or cross-agent review is claimed.**

The all-order proof is (G1)'s path-run bijection, the three-response
reconstruction (G3), and the matrix cone and monotonicity proof (R1).
Finite tests check their implementation; finite sample agreement is not
used as the proof of closure or of the unrestricted spectrum conjecture.

## Types, hypotheses, and equality boundaries

- Vertices are named, arcs are Boolean, and all pairwise arcs are present.
  A block join reverses exactly one otherwise forward cross-edge.
- The source word has a fixed Hamilton spine. K is actual path complement,
  not global converse or a resampled class-level transition.
- The seven-symbol grammar retains all five Q4 charts. Processing reverses
  the source-component order and the physical vertex labels.
- The scalar output is H; the countable three-coordinate quotient is exact
  for future single-bridge H queries. Vertex cost and actual source words
  are separate sidecars when those outputs are required.
- C counts a marked path and at most one other path, with no double-counting
  of the unordered pair. The empty complement contributes H(empty)=1.
- Three is the minimum dimension over Q of linear representations; arbitrary
  encodings of integers or the cardinality of a finite automaton are different
  questions. The quotient is not a complete isomorphism invariant.
- All-order family absence uses nondecreasing H, not an assumption that
  bounded order exhausts bounded H. W increases even on equal-H transitions.
- The extended family is all marked source SCC words with component size<=5.
  Its failure at47 is not a failure of unrestricted constructions.

## Exact universes and independent paths

`verify.py`, normal and optimized Python, gives identical receipts with
44883 always-active checks. It verifies:

1. All692 marked golden words through order10, including the empty word.
   Literal ordered source joins, actual K, and reversed labels agree with
   the independently generated backward-edge descriptors. Subset path DP
   verifies all three coordinates, not merely H.
2. All57 block matrices independently reconstructed from literal graph
   extensions of three linearly independent prefixes. The forest compiler
   and endpoint/path-cover DP do not share their transition mechanism.
3. Every pair of arbitrary labelled tournaments of orders1 through4:
   75^2=5625 one-bridge joins, including blocks without the fixed spine.
   Both (G1) and (G4) agree with actual graph DP.
4. All five Q4 phases distinguished by short contexts; the exact determinant
   -4 proving the linear lower bound; entrywise cone inequalities for all
   matrices; the special singleton cone argument.
5. Every successor of every state in the golden caps37,133,343 and extended
   cap47 certificates. These contain61,206,536,171 states respectively.
   Forward invariance plus the initial state is sufficient for the stated
   negative witnesses, regardless of how a candidate certificate was found.
6. Literal strong targets for the selected37,43,63,133,343 repairs, including
   their full arc rows. Literal permutations independently check three
   selected targets of orders6/7.
7. All99 words with H<=55 and all128 new-vertex cuts on each of the two H55
   words. The inherited positive55->111 source and mask96 are also replayed.
8. Actual81/27 automorphism groups, equal67 cyclic triples, common H462915,
   and different named two-birth responses2390013/1754919. The graph
   reconstruction and exact automorphism enumeration are inherited controls.
9. The rational aggregate series by symbolic matrix inversion, coprime
   numerator/denominator, direct weighted transfer through cost25, and
   complete word sums through cost10. These checks count presentations,
   never silently distinct target isomorphism classes.
10. Square-sum edge counts108/111 at55/56, and no earlier surplus for2<=n<56.

`reach.py` additionally completes golden closure through H100000, recording
the state hash, all selected witnesses, minimum vertex costs, and lane totals.
`extension.py` enumerates all64 strong5 candidates and filters exactly50;
it completes the57-symbol closure through H10001. The mathematical termination
proof applies to both routines. The large state sets are reproducible rather
than stored redundantly; small decisive closed sets are included literally.

## Inherited claims and current limits

The singleton tribonacci recurrence, golden source language, nonstrong-source
strong-complement lemma, prime/7p reduction and missing7/21 theorem are
inherited, with links in RECURRENCE and REACHABILITY. The111-dimensional
fixed-complement calculation is inherited, not recertified here. The square
graph's large-order surplus persistence is inherited; our fresh check is
only its first boundary through56.

The new material is the explicit single-bridge observer proof, the seven
matrices with full Q4 phase, the aggregate weighted series, complete bounded-H
family classification, strong5 extension, and golden55 parent obstruction.
No literature-priority claim is made. Existing general endpoint/path-cover
tools are credited rather than renamed as new theory.

Full Hamilton spectrum completeness, LRC(14), arbitrary-port minimality and
source-preserving inverse construction for every prime/7p remain OPEN.
