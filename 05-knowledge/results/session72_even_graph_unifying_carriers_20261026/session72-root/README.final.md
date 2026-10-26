# Session72: unifying theorems for symmetry, observations, and alignment

2026-10-26. **PROVED** scoped results and **FINITE-EXACT** certificates below;
**INHERITED/CITED** results and one new **OPEN** conjecture are distinguished.
This is a proposed unifying research programme, not a claim that the repo's
open problems are equivalent or share an unproved complexity reduction.

## The proposed bigger picture

I propose organizing these advances as **operation-relative representation
theory**: specify the actual state, allowed operations, observation and level;
then determine its exact quotient, response space, realization constraints
and alignment data. The earlier `(X,M,O,ρ)` reflection already supplied the
framework. This session gives four concrete theorem families within it.

[Full synthesis, proofs map and next conjecture](session72-root/UNIFYING_THEOREMS.md).

### 1. Harmonic visibility: counting is a dimension shadow

Tournament/sign-even-graph equinumerosity is an equality of SAME-ORDER
UNLABELLED types. Here “even” means automorphism-sign-even, not even vertex
degrees. The natural linear bridge is a signed Walsh basis of invariant
observables. Explicit inherited type codecs also exist, but their target
operations and symmetry scopes must be retained.

The operation-envelope theorem identifies all invariant source-independent
difference kernels, their invisible part, and their behaviour under quotient
averaging. Keeping the SAME weighted single-arc-reversal operation gives:

- Ordinary tournament types: inherited Laplacian gap4, n≥3.
- Switching-class types, modulo cuts: gap8, n≥4, with a simple quartic mode.
- Score-parity classes, modulo Eulerian edge sets: gap2(n−1), n≥3, with an
  explicit birth–death chain and complete spectrum.

Loops and multiplicities are essential. These are arc-flip quotient operators,
not the literal Hamilton-path-complement metagraph, nor Boolean simple graphs.
Source-dependent guarded operations need additional analysis.

[Operation theorem and inherited correspondence](session72-even/EVEN_GRAPH_OPERATION_ENVELOPE.md)
| [Score-parity theorem](session72-even/SCORE_PARITY_QUOTIENT.md)
| [Independent audit](session72-even-audit/REPORT.md).

### 2. State, response rank, and positive realization are different minima

For a finite state set and q-dependent reset rows plus deterministic actions,
the exact state
signature and the exact terminal linear response need not have equal size.
For a separately allowed involution, the signature is(q(x),q(τx)). Its realized
pairs form a bipartite graph. The terminal rank is the vertex count minus
component count; the missing joint-law directions are exactly its cycle space.
A spanning forest gives an optimal repair. The positive response cone is an
orthant up to linear isomorphism exactly when this graph is a forest.

Native order5 tournaments give50 pair cells, rank23, and27 invisible pair-law
directions. The normalized positive polytope has50 vertices in dimension22.
The full64-state literal-law kernel is41-dimensional, not27. These statements
concern preselected words and TERMINAL observations, not observed trajectories
or conditional tests.

[Reset/action, cycle and positivity proofs](session72-response/REPORT.md)
| [Independent combined audit](session72-unification-audit/REPORT.md).

The full Walsh basis has a different positivity issue: no information is lost,
but graph-indexed moments are not probabilities on graph types. Converse has
a negative diagonal coordinate already at order4 while preserving the correct
transported cone. [Exact matrix and cone](session72-root/WALSH_REALIZABILITY.md).

### 3. Reversible test-return converts state separation into law recovery

On a finite valid state set, suppose every successful guarded test can
restore the original state and every failure is absorbing. Then concatenated
test-return words realize products of guards. Separating guards then give full
terminal linear rank on all valid states. Named Hamilton re-marking followed
by its positional inverse supplies exactly this mechanism.

At order5, raw single-guard rank is44 over Q and43 over F2. Two restored tests
recover all64 coordinates: an explicit determinant−1 certificate works over
every field. At order6, two-test completeness over Q and F2 has an independent
1024-row certificate.

**New sharp-depth conjecture (OPEN):** over Q, the least number d_n of restored
whole-permutation tests needed to linearly span all functions on X_n is
ceil((n−2)/2), for every n≥3. The lower bound is proved by Boolean polynomial
degree; equality is verified through n=6. The first untested prediction is
full rank32768 at three tests when n=7. This is not a claim about elementary
move counts, and the other occurrences of44 are not identified with this rank.

### 4. Coherent alignment is controlled by prime-power layers

For r≥1 cyclic factors of periods d_i, the joint phase classes number
∏d_i/lcm(d_i). Full pairwise relative residues, with compatibility conditions,
recover the joint phase class. A sparse set of pair tests suffices exactly
when every nonempty prime-power layer induces a connected graph.

This explains precisely when a spanning tree is enough and when it is not.
The same phase action appears in necklaces, Frobenius coordinates, finite-field
tensor factors and an odd-period native carousel-tournament family. Arbitrary set-families
can occur as prime layers; no efficient sparse-optimization or complexity
claim is made.

The cyclic restriction is real: Paley7 triples(0,1,3) and(0,1,5) have identical
pair orbits but different triple orbits and different common-neighbour guards.
Their F8 realizations have F2 ranks2 and3. Deterministic coherent phases must
also be distinguished from pairwise marginals of probability laws.

[All-arity phase and sparse-test theorem](session72-cross-domain/THEOREM-CYCLIC-ALIGNMENT.md)
| [Cross-domain mechanisms and boundaries](session72-cross-domain/REPORT.md).

## What this unifies—and what it does not

The common success is finding the right representation AND proving what it
preserves under future operations. The limits are just as useful:

- A complete identifier always transports operations abstractly; apparent
  composition failure signals extra markings, a different state, efficiency,
  or a demanded simpler operation law.
- Affine carries without exact domains do not decide legal arithmetic or
  geometric continuation.
- The recovered factorial-Hankel theorem excludes every fixed finite-dimensional time-homogeneous linear
  carrier for the full characteristic-zero moment tower; characteristic p
  has a different sharp p-dimensional carrier.
- Identification of an unmarked source is not identification of its selected
  ports, and finite injectivity is not an all-order classification proof.

No Collatz, LRC, dissection-optimality, prime-atlas or other prize-problem
closure follows from these comparisons. The new results are exact unifying
subtheorems and a falsifiable next target, not a universal analogy.

## Sources and reproduction

The original counting theorem is
[Royle–Praeger–Glasby–Freedman–Devillers](https://arxiv.org/abs/2204.01947v2),
J. Algebraic Combin.57(2023),515–524. The lane reports credit prior repo
THM-4531/4721/4724/4725, THM-584/587/588, THM-840, the existing codecs, and
Sessions69–71 rather than presenting their results as new discoveries.

Run `python3 -B verify_packet.py --replay` in the packet. Curated native files,
exact source pins, independent normal/optimized receipts, and hostile controls
are frozen by the packet manifests. Numeric universes and coefficient fields
are stated per result; an empty CI check list is not described as a CI pass.
