# Validation receipt: pure-blue atom classification

Status: PASS, 2026-10-20.

## Mathematical result

Every normalized odd circulant tournament of order at least seven has
an explicitly constructed directed Hamiltonian path whose positional
reversal is not an anti-automorphism. The three prefix rerouting
families and unit reindexings by 2 or -2 prove this uniformly.
No primeness, size cutoff, cycle-count estimate, or large census enters
the proof.

With the independently audited strong-pure-blue odd-circulant reduction,
the nontrivial strong pure-blue atoms are exactly C3 and the regular C5.
The singleton is included only when the atom convention permits the
trivial strong component.

## Producer check

Command:

    python3 verify_prefix_classification.py --max-exhaustive 31 --large-max 501

Result:

- All 32,760 normalized skew connection sets at odd orders 9 through 31
  have valid constructed counterexample paths.
- All 1,235 structured and deterministic random examples at every odd
  order 9 through 501 pass.
- Every case checks the actual failed reversal arc, not just an
  indirect statistic.
- The search-free failed-arc rule is independently checked against the
  entire path word.
- Permutation orders, disjoint supports, and the conjugation identity
  are independently evaluated as full arrays.
- The four order-seven certificates are checked against all four
  normalized connection sets.
- All Hamiltonian paths of C3 and both normalized C5 presentations are
  enumerated and have anti-automorphic reversal.
- At most two unit reindexings occur; the two-step boundary is order nine.

Machine-readable output: `verification.json`.

## Separate independent audits

The parent independently implemented the witness construction without
producer imports and checked every normalized skew set at odd orders
9 through 23. Its symbolic review checked all parameter inequalities.

A separate adversarial auditor independently implemented the proof
families in `independent_prefix_obstruction.py`, with no producer imports.
It checked:

- Every normalized skew connection set at odd orders 9 through 31:
  32,760 cases.
- All four order-seven certificates.
- Every Hamiltonian path in C3 and both normalized C5 presentations.
- All three residual patterns A/B/C for every m from 4 through 200:
  591 cases.
- The actual failed anti-automorphism arcs, permutation orders, support
  identity, unit changes of coordinates, and all boundary cases.

The auditor separately reviewed the all-order proof and the direct
violated-arc table, and approved the mathematical classification.

## Failure boundaries preserved

- Circulancy alone is not purity.
- C5 is pure-blue but nonrigid; the proof never identifies these notions.
- The prefix lemma requires m>=4; order seven is handled by four
  completely explicit counterexample words, not by extrapolation.
- The unit changes of coordinates use that 2 is coprime to every odd n.
- The even-order exclusion remains a dependency of the companion
  foundations packet, including its cited crossing input.
- The result is not a classification of all regular tournaments, nor
  of all vertex-transitive tournaments without the purity assumption.

## Counting and branch-geometry controls

The separately audited corollary note
`PURE_BLUE_COUNTS_AND_BRANCH_GEOMETRY.md` uses the SCC center law and
the all-order rigid alphabet in addition to the atom classification.
It gives the rational counting series, mass-one and mass-three counts,
the exact pure-blue-to-mixed blue-line budget for n>=4, and the proof
that every trivalent pure-blue node has three distinct mixed neighbours.
The order-three pure-blue-to-pure-blue edge is retained explicitly.

Command:

    python3 verify_branch_corollaries.py

This passes:

- Coefficient identities and nonnegative line budgets through order100.
- Actual complement construction for 188 C5-centered normal forms,
  with 564 complement partners, through half-word weight12.
- The exact central-score-shell identity in every nontrivial outer
  stack, and independently enumerated center invariants
  (H,|Aut|)=(11,1),(15,3),(13,1).
- Every suffix-map bijection at m=0,...,12.
- Direct staircase-reflection orbit counts at n=3,...,50.

Output: `branch_corollaries_verification.json`. These are controls of
the symbolic proofs, not finite substitutions for their quantifiers.
