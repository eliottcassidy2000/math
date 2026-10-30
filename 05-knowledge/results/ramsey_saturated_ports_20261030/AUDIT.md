# Audit of the port assumptions

**REFUTED local implication, repaired before promotion.** The discovery
model initially called all four R-by-G row/column margins equal to four
necessary. P1 only forces red row sums and green column sums. It forces
the other two averages, not their individual values.

The first failed inference used regularity of an induced neighborhood
at an endpoint that did not belong to that neighborhood. On a red edge
uv, `|G(u) intersect G(v)|` is a variable t, not automatically four.
The correct table is in [P1](PROOFS.md).

`unequal_ports.json` is a27-vertex witness satisfying the repaired local
axioms and avoiding red/green triangles and blue K5. Its green row sums
are

    3,4,4,4,4,4,5,4,3,4,4,5,4.

They sum to52. Red column sums are nonuniform as well. This refutes
the implication from the **local** constraints to equitability; it is
not a56-vertex counterexample and does not refute a possible additional
global theorem forcing more uniformity.

`balanced_ports.json` imposes the stronger optional equitability ansatz.
It is also satisfiable. Therefore neither tested two-port model alone
produces an upper-bound certificate. The first discovery run used a
proof-output option and returned a nonzero process status after writing
its model. Re-running without that unused option succeeded; the stored
models were then independently checked. No solver status from the
failed process is used as a proof.

Recreate the two experiments with:

```text
python -X utf8 -B 05-knowledge/results/ramsey_saturated_ports_20261030/port_probe.py --equitable --prefix balanced_ports
python -X utf8 -B 05-knowledge/results/ramsey_saturated_ports_20261030/port_probe.py --force-defect --prefix unequal_ports
```

The first CNF has4251 variables and17966 clauses; the second has2409
variables and14289 clauses. SAT models can vary by solver version.
The native verifier checks the semantic coloring conditions directly,
including every triangle and every five-element subset. The independent
audit instead enumerates maximal cliques with NetworkX. Both witnesses
have color clique numbers `(2,2,4)`.

## Proof boundaries

- A52-element local automorphism group is not the global symmetry group.
  P4 supplies an injective restriction, and P5 excludes13-torsion using
  the actual fixed vertices. The weaker preliminary896 bound is subsumed
  by the single-port224 bound.
- P5 applies to color-preserving automorphisms. Allowing the red-green
  swap changes the bound to448. It does not turn the object into a tournament.
- Checking the members of a graph bank does not prove completeness.
  Completeness is a cited dependency; membership, isomorphisms and orbit
  counts have exact independent checks.
- The uniform formula in P2 assumes the specified constant degrees and
  local regularity at every vertex. It is not inferred from a single
  root's27-vertex witness.
- The P5 residue condition fails on the GF(16) positive control. Its
  order-five symmetry therefore does not contradict the general lemma.
- The81/27 control still has the same67 labelled cyclic triples and
  different full groups. Its67 is unrelated to the67 rooted Ramsey types;
  no map identifying those two sets has been supplied.
- Literature search and derivation do not establish novelty. These are
  recorded proofs of scoped implications, without a priority claim.

## Validity gate

Vertices are the original colored-graph vertices; the pairwise observable
is one of three named edge colors. There are no ties and no tournament
orientation is introduced. The preserved target is avoidance of specified
monochromatic cliques. Passing to counts loses joint adjacency; passing
to unrooted types loses the chosen ports. The saved witnesses, pair-orbit
inventory, and rooted recurrence explicitly retain those coordinates.

No random search is used to certify completeness. Exact universes are
the supplied one/seven/one graph banks, both fixed27-vertex matrices,
the GF(16) coloring, the inherited12-vertex pair, and all labelled growth
cases through five vertices. No exhaustive56- or90-vertex search was run.
