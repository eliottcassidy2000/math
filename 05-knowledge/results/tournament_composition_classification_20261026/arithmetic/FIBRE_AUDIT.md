# Independent audit of the cyclic triple fibre results

**Status: VERIFIED.** Arithmetic-agent review of
[fibre/RESULT.md](../fibre/RESULT.md), with the reviewed file SHA-256 recorded
in the two `fibre-audit-*.json` receipts. No mathematical defect found in
F1-F5. One type clarification was requested: `Aut(Q)/Stab(0)` is a coset
space; it is not generally a quotient group.

The primary [Boussaïri et al. PDF](https://arxiv.org/pdf/1805.04607) was
independently reopened. Theorem 10 states equality of the strong modules;
Theorem 24 gives the two labelled realizations of a prime tournament's
C3 structure; Proposition 51 gives the bijection with independent choices
at the decomposition nodes; Corollary 57 counts those choices. These are
the needed inputs, rather than a claim that all ordinary modules agree.

The sharp `4p` proof is valid. The size-preserving ambient quotient group
preserves the set of large child ports. A p-group acting on fewer than p
such ports fixes each one. It therefore preserves every possible refined
coloring on those ports and has the same p-part in either child-color
stabilizer. The remaining children have size at most three and their
ordinary isomorphism colors cannot vary inside their C3 fibres. Induction
handles the kernel. At `4p`, changing one of p nonisomorphic A/B colors
reduces the transitive quotient action to a point stabilizer and loses
exactly one factor p; there is no assumed classification of Aut(Q).

The uniqueness argument at twelve is valid. Any variable prime must be
3, and all proper child orders are stable by the smaller-order theorem.
Variation must therefore occur at a node with three size-four children;
the size budget makes this node the root and its quotient C3. To change
the color stabilizer, one realization must have all three colors equal.
The four-vertex classification leaves exactly the common A/B fibre, so
there is one such cyclic-triple hypergraph up to isomorphism.

The ordered path-cover proof is also valid. Run reversal at a prime node
and run reordering at a linear node preserve each covering path's vertex
set. For a move inside a proper module, this stronger preservation applies
to each collected run, so the inverse recovers the same positions and the
outside entering/exiting arcs remain valid by uniformity. Repeated visits
to a child do not collapse adjacent runs or change their external slots.

The independent checker `audit_fibre.py` uses the arithmetic agent's
separately written automorphism code and a direct subset Hamilton-path
recurrence. It checks the four stated twelve-vertex representatives on
identical labels. All are strong, all have the same 67 labelled cyclic
triples, the orders are `81,27,27,81`, the Hamilton count is `462915` in
each case, and the path-orbit counts are `5715,17145,17145,5715`.
The pair-codegree histogram is `{0:9, 1:9, 4:48}`, confirming the intrinsic
recovery of the three four-vertex blocks used for the hypergraph group.
These checks do not substitute for the all-order proofs above.

Reproduce from the repository root:

```text
python3 -B 05-knowledge/results/tournament_composition_classification_20261026/arithmetic/audit_fibre.py --out 05-knowledge/results/tournament_composition_classification_20261026/arithmetic/fibre-audit-normal.json
python3 -B -O 05-knowledge/results/tournament_composition_classification_20261026/arithmetic/audit_fibre.py --out 05-knowledge/results/tournament_composition_classification_20261026/arithmetic/fibre-audit-optimized.json
```

Normal and optimized outputs agree byte-for-byte. Checks remain active
under optimization, and all source/output files use LF line endings.
