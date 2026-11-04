# Three complete named ternary-return rounds

**FINITE-EXACT.** Starting from the actual named H170807 source, both
`3H+2` and `3H+4` can be realized for three consecutive rounds. This is a
complete binary tree with14 actual strong endpoint tournaments and14
explicit native programs. Every physical old name survives along every
branch. The first two endpoints are exactly the already fixed21-name
sources in [named_targets.json](../inverse/named_targets.json).

| Parent H | Child3H+2 (order) | Child3H+4 (order) |
|---:|---:|---:|
| 170807 | 512423 (21) | 512425 (21) |
| 512423 | 1537271 (22) | 1537273 (21) |
| 512425 | 1537277 (22) | 1537279 (22) |
| 1537271 | 4611815 (23) | 4611817 (23) |
| 1537273 | 4611821 (23) | 4611823 (23) |
| 1537277 | 4611833 (23) | 4611835 (24) |
| 1537279 | 4611839 (23) | 4611841 (23) |

This establishes actual repeated returns at this finite source tree. It
does not prove arbitrary-depth closure, prime/7p coverage, or the complete
spectrum. The components change between rounds; equal field coordinates
or response profiles are never used to identify two source graphs.

## Construction and preserved interface

For rounds2 and3 the inverse bank consists of every labelled strong seed
through order5, reduced by its eight retained-reader fields, together with
the35 inherited strong six-vertex representatives in
[low_sources.txt](../../metagraph_primegame_20261102/low_sources.txt).
Each eight-field key retains one literal graph and marked vertex. There
are219 keys. The finite inverse searches perform4,583,761 arithmetic
candidate checks across the12 later targets; they do not enumerate the
21–24 vertex tournament universe.

The reader appends sinks or vertices beating only the selected mark;
its exact eight-field transfer is in [reader.py](../inverse/reader.py).
A selected word ends at a strong component B. The outer graph is the
literal zero-ready wrapper with a transitive prefix and one final root.
Every selected witness has prefix length at least2 and total order at
least its actual parent order.

The physical root is39, the component mark31, and the spare port6. The
first two prefix vertices are50,51 in that order. Remaining old names
are assigned deterministically to the unused prefix, seed and tail slots;
any additional permanent names are fresh integers below1000. The full
relation, seed graph/root, word, eight fields, prefix, old-name identity
map, and all changed old arcs are saved in
[named_tree.json](named_tree.json).

Names persist, but old arcs and joint observations may change. These are
source-changing excursions, not induced extensions or H-neutral paths.
The two endpoint joint17 observers are retained separately for the named
pairs(39,6) and(31,6), with their ordered cross-bin attachments. Integer
golden coordinates and the inherited F49/F343 charts are recorded for39,
31 and6 on every endpoint. No field chart is substituted for an actual
source or a legal attachment.

## Native programs and independent verification

For each changed old arc, the inherited one-helper compiler supplies two
legal native complements followed by deleting only that helper. The
program explicitly stores the helper's out-mask and both full Hamilton
words. Every native guard and resulting old arc is replayed, and the
completed endpoint is compared with the exact named target relation.
A bounded search tries both orientations before the inherited complete
search; finding a path never bypasses its literal legality check.

[native_tree.json](native_tree.json) contains all14 complete programs,
with944 native complements in total and peak order25. It also records
permanent births, all temporary helper names, and the hash of every
completed arc-gadget endpoint. Every old name survives. Intermediate
strongness, H, port stars and field phases are not asserted invariant.

[verify.json](verify.json) independently recounts each actual graph using
the transitive-gap engine, verifies its full named wrapper relation and
strongness, checks component eight/nine fields against the seed-reader
replay, and matches full-graph H, P2 and root fields with the wrapper lift.
It then computes both named joint17 vectors and all three retained-port
charts. The first-round pair data agree exactly with the independent
[interface receipt](../interface/phase.json).

The mathematical return ingredient is now narrower: obtain a reusable
source bank or a source rewrite theorem that supplies such named reader
witnesses at arbitrary depth, with the joint interfaces needed by the
next step. The finite tree supplies matched source objects for three
rounds, but no uniform bank-closure theorem. Bare ternary branches still
omit a residue class, and prime/7p predecessors need not remain in the
same lane.

## Reproduction

Run in this directory:

```
python3 -B build.py
python3 -B native_tree.py
python3 -B verify.py
```

The native script replays cached programs before accepting them; deleting
its JSON reruns the constructive compiler. Source hashes are pinned in
all three receipts. The same scripts also run with `python3 -B -O` and
use explicit validity gates rather than Python assertions. All three receipts
are byte-identical under normal and optimized replay; see
[replay.json](replay.json). The computation stops at depth3; there is no
fourth-round extrapolation.
