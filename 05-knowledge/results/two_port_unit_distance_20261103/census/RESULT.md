# Two-port reconstruction is exact through order nine

**FINITE-EXACT: every labelled tournament on at most nine vertices is
determined by (H,S,E,C), hence by its two-port response.** The universe is
the complete inherited class bank, not a sample of labelled tournaments.
No unbounded two-port reconstruction claim follows.

Here S and E count paths starting and ending at each named vertex, and
C_uv counts ordered spanning two-path covers whose first path ends at u
and second path starts at v. H=sum S=sum E. Session27's actual two-port
queries recover this entire sidecar. We keep its ordered C entries; the
symmetric kernel or skew-current signs alone are not substituted for it.

## Why a class census proves a labelled statement

For every inherited tournament representative, form the complete directed
weighted graph with vertex colours (S_u,E_u) and edge labels C_uv. Compute
an exact canonical representative of this sidecar, retaining **every**
minimizing labelling. Two conditions are checked:

1. Distinct tournament isomorphism classes have distinct sidecar classes.
2. All sidecar-minimizing labellings produce the same labelled tournament
   code. Equivalently, every automorphism of the sidecar preserves the
   original tournament.

The second check is essential: unlabelled class separation alone would not
rule out two different namings of one tournament with identical response.
Together the conditions prove the stated labelled injectivity, by first
aligning sidecars and then using the automorphism condition. Sidecar
automorphism orders also agree with the independently inherited tournament
automorphism orders in every case.

| Order | Tournament classes | Sidecar classes | Extra sidecar symmetries |
|---:|---:|---:|---:|
|1|1|1|0|
|2|1|1|0|
|3|2|2|0|
|4|4|4|0|
|5|12|12|0|
|6|56|56|0|
|7|456|456|0|
|8|6880|6880|0|
|9|191536|191536|0|

Thus a counterexample to universal two-port reconstruction must now have
at least ten vertices. A collision of the smaller (H,S,E,C) sidecar beyond
this range would still need testing against the complete bidegree-(2,2)
response: ordered three-path-cover information may distinguish it.

## Exact computation and independent controls

The source banks and their full-enumeration provenance are
[the October26 enumeration packet](../../hn_metagraph_information_20261026/enumeration/RESULT.md).
`representatives.json` supplies orders through8; `representatives9.json.gz`
supplies orders through9. This session does not re-claim their enumeration
as new. Their byte hashes, coding convention and counts are in the receipts.

The producer computes endpoint subset tables, then

    C_uv = sum_(empty!=A!=V) end_A(u) start_(V\A)(v).

For canonicalization, equitable refinement uses the actual vertex colours
and the multiset of ordered incoming/outgoing C labels in every cell.
Refinement is label-equivariant. Every permutation inside each remaining
cell is explicitly visited; there is no heuristic search cutoff. All
minimum leaves are retained, recording minimum and maximum induced
tournament codes. At order9 this uses1,362,537 leaves, with worst case9!;
equal min/max codes certify the automorphism condition without assuming
the inherited automorphism count is correct.

`check.py` is the Python implementation through8. `exact.cpp` and its
temporary-build wrapper `run_cpp.py` give the full order9 run. Their entire
order-by-order results through8 agree. They implement the same mathematics
in different languages, not independent mathematical algorithms.

`audit.py` separately counts by literal vertex permutations and admissible
cuts, then canonicalizes by **all n! labellings without refinement** for
every class through6. Both data and automorphism tests agree. Six specified
larger classes through9 also have independent literal path/cover counts.
This is a separate algorithmic control; the order9 completeness conclusion
still uses the inherited complete bank and the exhaustive canonicalization
proved correct above.

All path and cover counts are bounded by8*9! in this universe, so unsigned
64-bit arithmetic is exact. Adjacency codes use at most36 bits. No floating
point, random hash comparison, skipped symmetry case, or quotient merging
is used. Standard-library Python plus an available C++17 compiler suffice;
the wrapper installs nothing and builds in a temporary directory.

Reproduce from this directory:

```sh
python3 check.py
python3 run_cpp.py
python3 audit.py
```

Their normal/optimized Python receipts agree. The C++ compilation uses-O3
in both wrapper modes; those are not falsely described as two different
C++ optimization audits. The proof does not depend on performance timings.
