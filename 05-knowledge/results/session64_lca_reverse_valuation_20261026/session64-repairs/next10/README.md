# Complete order-ten closed endpoint-four replacement

2026-10-26. **FINITE-EXACT:** every closed blue marked Hamiltonian path
at order10 with endpoint outdegrees `(5,4)` has three directed black
K moves to its exact labelled target. This packet extends the separately
sealed order11 and n>=12 range. It does not claim full black-core
connectivity or cover open marked paths.

## Exact source partition and the nine exceptions

The complete universe imposes the path `P=(0,...,9)`, its closing arc,
its actual positional anti-reversal, `N-(0)={9} union I`, with I any
three-element subset of `{2,...,8}`, and both bits `T(1,8),T(2,8)`.
No strongness, clone, module, shared-port or isomorphism filter occurs.

All 140 source cells are nonempty. The inherited base bank covers 130:
122 are uniformly three-black and eight admit a blue word. The other 10
closed-endpair cells split into 21 nonempty disjoint source leaves:
20 are uniformly three-black and one admits a blue word. Every inherited
word in all 151 leaves is universally directed and has the exact target.
Each of the nine feasible blue loci has dimension 0, and the nine sources
are distinct.

| I | Middle | x,c chord | Old branch | Blue step |
|---|---|---|---|---:|
|{2,3,4}|open|x->c|open_release|1|
|{2,3,6}|open|c->x|open_release|1|
|{2,3,6}|closed|c->x|closed_zero|1|
|{2,5,6}|closed|c->x|closed_zero|3|
|{3,4,6}|closed|x->c|closed_one|3|
|{3,4,8}|open|x->c|open_release|1|
|{3,6,8}|open|c->x|open_release|1|
|{3,7,8}|closed|x->c|h3_forward,1->3|3|
|{4,5,6}|closed|c->x|closed_zero|3|

The old ledger applies off each blue locus. On it, `found.json` gives
an exact three-black replacement. `full_check.py` reconstructs the
complete raw affine systems, checks exact equality of their solution
union with the nine rows, and replays every new word on literal vertices.
The inherited endpair splits are rechecked at10, not assumed valid below
their published uniform thresholds. Their disjointness and exact
cardinality sum are checked inside each original source cell.

## Exact clone transport is substantive

Every state records its complete rows, all labelled endpoint ports,
all flipped pairs, an anti-reversal mismatch proving blackness, global
twin pairs, ordered maximal transitive modules and their full intersection
matrix through each move. Of the 36 states in these repairs, 17 have
nontrivial clones; seven states have five distinct size-two clone blocks.
There is no contraction, isomorphism quotient, guessed merger or lifted
word. Internal clone arcs remain explicit in every K move.

## Reproduction

```
python full_check.py > full-normal.json
python -O full_check.py > full-optimized.json
cmp full-normal.json full-optimized.json
```

Both receipts are PASS and byte-identical. The sealed parent Gaussian
engine and old-bank word definitions are imported through a checked
SHA256 of `../check.py`. All new finite partitions and literal replay
are local. Only the files in MANIFEST.json belong to this proof packet;
exploratory diagnostics and search scripts are not dependencies.
