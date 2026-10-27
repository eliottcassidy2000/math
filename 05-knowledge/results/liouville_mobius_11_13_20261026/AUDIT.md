# Proof and computation audit

2026-10-26. Reviewed base `4959b4933cba`. The all-order arguments have a
same-agent proof review; the finite claims have separately implemented
computational paths. There was no independent-agent proof audit in this
continuation. Earlier independent audits are cited only for their own
inherited results.

## Proof obligations

- **Arithmetic types:** pointwise mu² differs from convolution squaring;
  the unit is epsilon, not the constant sequence1. A1 proves all inverse
  identities prime by prime, including n=1. Dirichlet series are restricted
  to absolute convergence; ordinary products are formal and locally finite.
- **Contexts and information:** A4 quantifies over every positive integer
  context for its arithmetic classification. Its tournament application
  uses the explicitly attained alphabet11,13 and the proved order-join
  product law. No all-primes attainment or arbitrary-substitution product
  law is assumed. Rank over every field follows from an integral inverse.
- **Necessary and sufficient:** primitive-period inversion is derived from
  the partition by least period, not a division of nonfree orbits. At a
  squarefree period lambda happens to agree with mu; the proper-square
  terms give its exact failure elsewhere. The constant words supply the
  hostile at prime squares.
- **Symmetry induction:** the fixed labelled C3 realization theorem is a
  cited structural input, not inferred from a census. A squarefree parent
  has squarefree child automorphism orders because each divides the product
  formula. Induction handles child p-valuations before the quotient step.
  A p-subgroup failing to preserve the other colouring must have a
  nontrivial orbit, of at least p equally sized children. Below5p they are
  size4; their common ambiguous C3 type forces repeated factors3 in the
  squarefree realization. Both prime quotient orientations have the same
  permutation automorphism group. Linear nodes have no quotient symmetry.
- **Sharpness:** D and D* are rigid and distinct; C_p acts regularly on
  carousel ports. Thus the two orders are p and1, not merely quantities
  with differing p-parts. At15 the uniqueness argument first rules out
  changes below the root, then p>=5, then size4 children. Only the
  exhaustive finite order5 classification is computationally supplied.
- **Orbit classification:** labelled lifts, ordinary classes and hypergraph
  orbits are distinguished. Pair codegree intrinsically recovers outer
  blocks, yielding the full384-element hypergraph group. The displayed
  mass is orbit-stabilizer, not division by one uniform tournament order.
- **Historical recovery:** equal score and equal number of cyclic triples
  is not equality of labelled triple sets. The H11/13 pair is identified
  by full relabelling. The old THM-466 inverse wording is explicitly
  repaired with a minimal-order collision; its forward formula survives.

The scope checks expose no remaining gap in these proofs, conditional on
the cited C3 realization input. They do not certify an all-order metagraph
classifier, an analytic cancellation conjecture, or LRC(14).

## Finite universes and controls

`check.py` uses only the Python standard library and explicit runtime gates
that remain enabled under `python -O`. There are no inherited census filters.

| Universe | Positive and hostile controls |
|---|---|
| Integers1..1000; all ordered pairs1..100 | Pointwise mask, convolution inverses, square lift, summatory identities; overlap annihilation versus complete multiplicativity |
| All subsets of r prime coordinates, r=1..5 | Integer inverse of mu matrix; exact rational rank one for lambda |
| Formal series through degree80 | Direct product and ghost-recursion coefficients; square dilation and its inverse; prime-local ghost formulas |
| Every labelled order4 tournament (64), all24 relabellings | Only ambiguous C3 type has two ordinary classes, each Aut3 |
| Every labelled order5 tournament (1024), all120 relabellings | H/odd-cycle identity, full H census, historical H11/13 isomorphisms; all30 rigid ambiguous labelled C3 types are D/E |
| Four actual order10 joins | Independent Hamilton-path counts121,143,143,169; mu separates the cross product from squares |
| Two explicit order6 tournaments | H9 for both, cycle-collection profiles `(1,2,1)` and `(1,4,0)`; hostile to inversion from H digits |
| Full first order15 C3 fibre | All1024 labelled lifts,384 hypergraph permutations, four disjoint orbits128/384/384/128; all four Hamilton counts and automorphism orders |
| Sharp constructions for p=3,5,7,11,13 | Full expanded adjacency automorphism searches at5p vertices, and literal equality of labelled C3 hypergraphs; orders p/1 |
| All binary words of lengths1..13 | Least period tested directly against every nontrivial rotation; primitive counts and necklaces agree with mu inversion |
| Periods121,169 | Exact huge integers; true numerator divisible by period, lambda replacement has remainder2 |

`audit.py` imports no producer functions. An independent sieve verifies
mu and lambda divisor sums to1000. For each of the four15-vertex rows it
checks the tournament axioms, refines vertex colours by adjacency counts,
then exhausts every within-cell permutation. The homogeneous classes
require7776 candidates each; mixed classes refine to singleton cells.
A backward memoized endpoint recurrence independently reproduces the
four Hamilton-path counts. These paths use no predicted module partition.

## Frozen receipts and source provenance

`SHA256SUMS` covers every packet file except itself and ignored bytecode.
`source-pins.json` records the reviewed revision, inherited-file hashes
with normalized line endings, and before/after hashes for the repaired
THM-466 source. Source pins provide provenance;
`verify.py --check-sources` checks for later drift explicitly. It is
separate from replay so a later linked source correction does not silently
rewrite a receipt.

`verify.py --replay` verifies hashes first, runs both programs in normal
and optimized modes without rewriting committed expectations, and compares
each output byte-for-byte with its frozen receipt. It also checks relative
links in the packet. Numerical agreement supports the specified finite
claims; the all-order conclusions rely on the displayed proofs.
