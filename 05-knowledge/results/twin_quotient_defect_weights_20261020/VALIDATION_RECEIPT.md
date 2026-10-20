# Validation receipt

**Status: FROZEN_REVIEWED.** 20 October 2026.

Both proof notes have received fresh independent review. In particular, the
sharp <=11 bound was checked for arbitrary small-class automorphism groups;
the order12 proof exhausts the cases by number and size of big twin classes;
and the component-product argument is used only for nonisomorphic connected
nontrivial components. The identical-component example 2K_(3,3) has its own
direct packing proof.

| Path | Exact scope | Runtime checks | Result |
|---|---|---:|---|
| `verify_twin_quotient.py` | 8,468 weighted quotient cases; 142 explicit contraction tests; balanced c=1..60; split/join r=1..20 | 25,986 | PASS |
| `independent_packing_audit.py` | All 208 unlabelled graphs on 1..6 vertices via full symmetric groups; 534 actual packing-stabilizer sign tests; all 32 labelled false-twin-free four-vertex quotients with weight3 | 994 | PASS |
| `verify_sharp_defects.py` | All 16,798 labelled false-twin-free quotients with 1..4 vertices, positive integer weights totaling <=12; opposite odd-sign families r=1..40 | 77,005 | PASS |

All three paths were run both normally and with `python3 -O`; each pair of
outputs is byte-identical. Validation uses explicit exceptions, not removable
assert statements. No external package, graph database, or network is needed.

The literal packing program imports neither quotient verifier. It constructs
the marked triples as actual vertex subsets, finds the actual stabilizers from
all vertex permutations, and computes the contracted graph sign explicitly.
The general sharp theorem is a proof, not an inference from either restricted
finite universe. The q<=4 weighted scan has 16 labelled first-boundary
presentations: one K_(6,6), twelve P4[3], and three 2K_(3,3), representing the
three unlabelled types explained in the theorem.

## Scope boundaries

This packet does not furnish a graph-structural involution on all decorated
defects. It does not refute existence of such an involution, and it does not
reopen the already-proved global p-group equality. It refutes the specific
unweighted all-order extrapolation and its even/odd sign partition, preserves
the correct finite observations, and identifies the missing packing-orbit
and stabilizer data.

`SHA256SUMS.txt` pins the deliverable files.
