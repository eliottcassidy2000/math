# Independent review: conditional indexed observer compiler

Verdict: **PASS** on the stated mixed-depth ternary automorphism-orbit
carrier, after a selector-scope clarification made by the author. Reviewed
`../compiler/README.md` sections 2–6, their surrounding carrier/guard scope,
and `../compiler/check.py`. The inherited full hierarchical-group
classification is a supplied premise, not a new result of this audit.

## Mathematical checks

- An indexed owner answer fixes exactly its input-prefix parameter
  coordinates. Shared coordinates must agree, and the native action has no
  further constraints. Freeness on the full ordered spanning path makes the
  completion count `3^(D-|F|)` exact. The group law is not being mistaken for
  coordinatewise addition.
- For a target observation tuple, distinct assignments of newly exposed
  coordinates give distinct answer tuples: every coordinate occurs in an
  actual requested digit and is recovered by subtracting the source digit.
  The fibre multiplicity and coverage-function dependence are exact.
  A forced deeper coordinate forces its ancestors because threshold
  components are nested and the input prefixes agree. Thus the single-owner
  prefix-cylinder statement follows in this owner-pin model.
- A previously untouched required bottom slot retains a rotation that
  changes the requested owner and fixes all observed source owners. For
  nonzero answers the correct action is `g0 k`, as stated. Every required
  slot must therefore be touched by an actual indexed flag. The minimum
  cover is `|R|-nu(J[R])`, treating flags to nonrequired slots as singleton
  covers. The source positions exposed by a query do not depend on answer
  values, so adaptive queries do not improve this optimum within the pure
  orbit model.
- If one input-prefix coordinate changes, `g_(a') g_a^-1` really changes
  just the corresponding digit on the common **target** prefix subtree.
  Later target digits stay unchanged because both frames evaluate those
  digits from the same original input address and unchanged parameters.
  Nested components make the earlier prefix image common across the
  affected named cells. This validates the actual native realization of
  the Hamming coordinate edges, not merely their count.
- The compatible-pin move graph is exactly `H(r,3)`. Its random-scan
  heat-bath spectrum is `1-j/r` with multiplicity `binomial(r,j)2^j` for
  `r>0`; the dimension-zero case is separately handled. These automorphisms
  preserve the ordered arc bits and backward free word, so the result
  concerns owner lifts in one target isomorphism class, not arbitrary
  changes between target classes or elementary metagraph moves.
- If an authenticated guard factors through a query footprint exposing
  `q` free coordinates, the induced graph is exactly
  `H(q,3)[A] Cartesian-product H(r-q,3)`. Count, components and component
  diameters follow directly, including the empty and zero-dimensional
  boundaries. Additional invisible coordinates cannot reconnect distinct
  visible components. The necessary guard factorization is expressly a
  separate obligation; it is not assumed for electrical or mod-eight data.

## Selector clarification and corrected scope

The first wording could have read the supplied reference path `P` as part
of the selector input while simultaneously claiming residual symmetry could
not select a path. A retained reference path or parameter origin already
provides a gauge: with no pins, choosing `P` is a completion rule.

The author clarified section 6 before final replay. The no-equivariant-
selector claim takes the **unpointed candidate fibre plus observations** as
input; the reference path is only an auxiliary chart. In this scope the
target-owner stabilizer acts freely and transitively on the fibre, and no
equivariant choice exists when that stabilizer is nontrivial. Chart-dependent
completion algorithms using an additionally retained `P` are explicitly
allowed. The abstract `H -> H/K` argument is correct.

## Reproduction and finite universe

Final normal and optimized runs both passed **117,691** explicit checks:

```text
python -B 05-knowledge/results/hamiltonian_observer_restrictions_20261023/compiler/check.py
python -O -B 05-knowledge/results/hamiltonian_observer_restrictions_20261023/compiler/check.py
```

The controls retain actual source arcs and path owners. They include all
`9!` permutations of the regular-nine carrier, all `81^3` internal-chart
triples for the named mixed order-27 source, all 32,768 required-slot subsets
of the order-63 hostile, 30,618 actual directed coordinate moves, and all
512 indexed owner-clock subsets on regular nine. The final footprint lane
also exhausts all 512 subsets of one actual owner's nine images. Frozen
receipt equality adds one check after the receipt's internal subtotal;
the different displayed subtotal is intentional, not a lost test.

The positive footprint example has 18 frames in one diameter-three
component; the disconnected example has three nine-frame components of
diameter two. The guard hostile is retained rather than silently assuming
that all consistent pins imply connectivity under arbitrary guards.

No claim that this orbit exhausts all Hamiltonian paths, all tournament
isomorphism classes, or all native metagraph wiring was found in the final
scope. No blocker remains.
