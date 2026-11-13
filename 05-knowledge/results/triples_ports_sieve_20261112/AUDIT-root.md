# Independent audit of the dihedral and triple lanes

Status: **VERIFIED** proof review and independent finite controls by the
arithmetic lane. No unresolved defect was found in the elementary theorems
after the residual-demand clarification below. This audit does not promote
the separate cited cluster-algebra or root-lattice literature to new proofs.

## Dihedral boundary states

The proof in [dihedral/PROOF.md](dihedral/PROOF.md) has the required
quantifiers: finite connected support, at least one observed boundary
vertex, signed-affine bijective transport, and ordered boundary labels.
Tree transport determines every interior value from the first boundary
value. A nontrivial translation has no fixed point; distinct reflections
have incompatible fixed equations. This leaves exactly the listed full,
singleton, antipodal-pair, or empty relations. The exceptional modulus 2
correctly identifies the two signs as permutations.

The state counts do not identify boundary relabellings. Distinct proper
relations are distinguished by a realizable pinning context: a singleton
for odd modulus, an antipodal pair for even modulus, and a full translation
orbit at modulus 2. This supports the stated contextual minimality, with
its explicit arbitrary-compatible-context scope. It does not claim that
every state is reached from one fixed initial gasket cell.

The square-label realization is valid in the abstract constraint category:
0 and 1 are square residues, and `R_1 R_0` is translation by 1. Subdivision
realizes words without asserting an injective lift to integer square-sum
vertices. The even-modulus circle counterexample and the loss of sign at
modulus 2 are necessary boundaries. The CRT actions at 30 and 42 have the
claimed index 2 because the odd coordinates retain a common sign; a free
2-coordinate and free 3-coordinate with a pinned 5- or 7-coordinate indeed
cannot be one global nonempty fixed set.

## Theta matrices and gasket scales

The proof in [triples/PROOF.md](triples/PROOF.md) correctly grounds one theta
hub and subtracts one unit at the other. Direct edges and parallel length-1
arms contribute to the grounded diagonal, so the root remains diagonal 2
after subtraction. The path blocks are positive definite, and their root
inverse entries give the stated Schur complement. The determinant, complete
inertia split and equality vectors follow by block congruence. The sorted
positive/zero/negative triple classification includes the empty-arm cases.

For the gasket, the proposed midpoint value satisfies every harmonic
equation and gives the displayed boundary current. The eliminated interior
Laplacian is positive definite. Iteration of the effective boundary form
therefore gives the conductance and resistance factors. The mass-resistance
factor is explicitly an elementary network scale, not a spectral-limit
claim.

One underspecified sentence was repaired during review: the theta parity
hostile requires **internal residual demand 1**, in addition to hub demands
`(2,2)`. With those demands it has two degree completions for arm triples
`(2,3,4)`, `(2,3,6)`, and `(2,3,8)`. Ordinary degree-2 interiors would be a
different problem. The current proof now states the correct demands and
retains all residual demands in its combined state.

The square-sum cross-reference now exists and explicitly supplies both
marked theta kernels and the Q25/Q26 outside-partition witness. This audit
checks consistency of that reference; its full first-25 graph computation
is audited and reproduced in its own lane.

## Independent controls and scope

Both root checkers were rerun normally and under `python -O -B`; all four
runs passed. Their reported exact finite universes were:

- Literal dihedral triangle relations for moduli 2 through 11; all 125
  joins at modulus 2 and all 912,673 joins at modulus 4. The four additional
  5,000-join samples remain explicitly sampled, not exhaustive.
- All 364 sorted arm triples in `1<=p<=q<=r<=12`, with exact determinant,
  inertia and null-vector checks; literal gasket elimination for 2 through
  9 corners.

The additional independent [audit control](arithmetic/audit_root.py) was
also run normally and with `-O`, both passing. Its
[receipt](arithmetic/audit_root.json) records:

- All 364 literal theta multigraph Laplacians, grounded and root-subtracted,
  exactly equal the separate star-matrix construction.
- Complete edge-subset counts for the three corrected parity hostiles.
- All 21 combinations of moduli 2 through 8 and boundary counts 1, 2, 3,
  enumerating two-holonomy fixed sets and all transport tuples independently
  of the normal-form state generator.
- Explicit group closure under the two square-labelled reflection
  generators for moduli 2 through 8.

Reproduce the extra controls from the repository root with
`python -B 05-knowledge/results/triples_ports_sieve_20261112/arithmetic/audit_root.py`
and the same command with `-O` before `-B`. All checks use explicit
exceptions, so optimization cannot disable them. These finite controls
corroborate the proofs; the unbounded statements depend on the proofs and
their stated categories, not on extrapolation from these ranges.
