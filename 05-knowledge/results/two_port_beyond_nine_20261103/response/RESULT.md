# Full joint endpoints and three paths, with actual query checks

**PROVED evaluation formula:** [algebra A4](../algebra/PROOF.md).
**FINITE-EXACT independent controls:** eight named tournaments, including six
at order10. These controls are not an order-ten exhaustive response census.

The saved `tensors.json` retains every named entry of S,E,M,C,U,V,D, rather
than just their marginals or their supports. Here M is the joint start/end
Hamilton matrix, C is the ordered two-path matrix, U and V retain three joint
endpoints of ordered two-path covers, and

    D[i,j,k,l] = # ordered spanning three-path covers(P,Q,R)
                 with end(P)=i, start(Q)=j, end(Q)=k, start(R)=l.

The middle singleton j=k is included. V is the tensor called W in A4's
prose; the notation change is only a field name. Arrays are flattened in
the displayed index order. Full support is implicit in the spanning condition.

Two algorithms compute every entry: literal vertex permutations with zero,
one or two cuts; and independent squarefree subset multiplication of the
simple-path matrix. For a literal permutation, all backward consecutive
adjacencies must be cut. This is a direct, complete enumeration of ordered
path slots and supplies an independent check of disjointness and normalization.

The tensor formula is then evaluated on every mask pair of weights<=2,
giving the complete bidegree-(2,2) battery. Each value is compared with a
fresh Hamilton-path subset DP on the actual tournament with two added ports.
Nine additional high-weight queries per control compare direct counts,
the tensor formula, and Boolean Moebius interpolation from that battery.
Mask bits here mean old vertex->port, and the fixed port arc is x->y.

| Control | H | Ordered two-path covers | Ordered three-path covers |
|---|---:|---:|---:|
| Cyclic5 and its converse, each | 15 | 90 | 210 |
| Order10 bank index0 | 3 | 1530 | 67572 |
| Order10 bank index100000 | 10019 | 166594 | 1246428 |
| Order10 bank index4866528 | 6991 | 126986 | 1019052 |
| Order10 bank index9733055 | 1971 | 63378 | 674748 |

Two additional order-ten controls are transitive10 and a universal source
above C3[C3,C3,C3]. Their sidecar automorphism orders are1 and81. At every
order-ten control, the census implementation's S,E,C entries match literal
counts, and its refined canonicalization has the same symmetry order and
extra-symmetry verdict as enumeration of all10! labellings. The brute audit
only rejects a labelling early when its full vertex-colour prefix is already
lexicographically worse; no potential minimum is dropped. The n5 output's
`sidecar_aut:0` is an unused sentinel; that canonical audit is performed only
on the six n10 controls.

Each n5 control checks256 battery entries plus9 additional queries; each
n10 control checks3136 plus9, for19,400 actual two-port checks altogether.
The cyclic-five pair has different full batteries despite satisfying the
earlier robust changed-arc conditions. The saved normal and optimized
artifacts agree exactly. Both wrapper modes compile C++ with-O3; this is not
an audit of different C++ optimization levels.

The tensors are sufficient to evaluate every two-port query. The converse
is deliberately not asserted: two probes expose combinations of these
tensor entries, and equality of the response polynomial need not imply
entrywise equality of the entire tensor. A proposed universal proof must
justify that additional inference rather than assume it.

Reproduce from the repository root:

```sh
python3 05-knowledge/results/two_port_beyond_nine_20261103/response/run.py
```

This prints the complete artifact, including all tensors and query tables;
use the packet verifier for a compact frozen comparison. Source and bank
hashes are included. All exact quantities are bounded by the number of
permutations and cuts at these fixed orders; unsigned64-bit storage suffices.
