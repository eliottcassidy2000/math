# Independent audit: count-free local tournament enumeration

Status: **PASS**. Mathematical conclusions are all-order corollaries of the
inherited relocation-connectivity proof. New finite checks cover the explicitly
listed ranges below; no practical all-order runtime claim or novelty claim is
made. The audited production source is `enumerate.cpp`, SHA-256
`21a1c77169aa66e45b7d868e1e23ff94fc7ceb7bb9bea960ca1e758a941138c3`.

## 1. Inherited mathematical dependency

Reviewed `THEOREM.txt` and `LOCAL_FOREST_ADDENDUM.txt` in repository directory
`05-knowledge/results/hamilton_path_relocation_connectivity_20261020/` at source
commit `b4234f728c0516b5fc521c001d69e26f3d2508d5`.

The proof is valid for every finite nonempty tournament. Its potentially delicate
steps are sound:

- Common deletion is directed **only for distinct legal relocation neighbours**;
  this is not an assertion that arbitrary vertex deletion preserves a path.
- In edge lifting, a shared insertion gap accepts the two inserted vertices in
  their mutual arc order. Both deletion identities hold, including endpoint gaps.
- Front promotion in a strong tournament has a legal last-in-neighbour insertion.
  A vertex eventually absent from the front has a stabilized nonincreasing
  position. Dominating a recurrent front vertex would force that position lower.
  The recurrent front set therefore dominates its complement and, by strong
  connectivity, is the entire vertex set.
- The initial SCC is a contiguous first block of every Hamilton word. Promotion
  within it preserves all outgoing joins. The chosen vertex can be made first,
  connecting every word to the connected removable-vertex lift subgraph.
- Automorphism quotienting identifies precisely equal normalized path tilings.
  Collapsed relocation loops cause no problem: C3 is the smallest control.
- Reversal plus converse permutes the free coordinates without flipping bits.
  Adding those reflection edges merges exactly converse-paired ordinary classes.

No external reconfiguration result or supplied class count enters this proof.

## 2. Count-free recognition: necessity and sufficiency

For a partition of the **whole** fixed-path cube, requiring each block to be
connected by permitted auxiliary moves and requiring both endpoints of **every**
permitted edge to have the same block label is necessary and sufficient for the
component partition. Connectedness puts a block inside one component; full edge
closure prevents the component from leaving that block. Reflection must be
included in both conditions for converse-merged classes.

The equivalent expanded half-edge criterion still requires a bijection onto the
whole cube, exact free-bit complement pairing, and exact reflection-fixed colour.
It no longer needs a supplied value of K_n. Neither the node count nor a census
is hidden in the closure check. Relocation/reflection edges are transport
witnesses, not the literal complement lines.

All **4140** partitions of the order-four cube were tested, separately for
ordinary and merged recognition. In each case exactly the true partition passed
both conditions. Ordinary: 34 partitions pass connectedness, 15 pass closure,
and only one passes both. Merged: 68 pass connectedness, 5 pass closure, and only
one passes both. Splitting into singletons is connected but not closed. Merging
the whole cube is closed but not connected. See `recognition_hostiles.json`.

It is sufficient to enumerate rightward relocations globally and compare their
two endpoint labels: this is undirected edge closure, not one-sided forward
reachability. The inverse of every leftward move is among those rightward moves.

## 3. Independent finite verification

`independent_check.py` was written before inspecting the production map code.
It decodes each free-bit mask into a full adjacency matrix. It constructs moved
words explicitly, checks every consecutive arc, and encodes the complete new
relation. No imported class atlas, automorphism routine, Burnside computation,
or A000568 value is used to generate its partitions.

- Orders **1 through 7**: exhaustive cube enumeration, every rightward legality
  predicate versus full-path validity, ordinary and merged components, reflection
  involution/fixed count, complement colour consistency, and fibre-size census.
- Orders **1 through 6**: every vertex permutation is inspected. The ordinary and
  converse-merged partitions agree **pointwise**, not merely in cardinality,
  with the exhaustive Hamilton-word normalization oracle. Every possible legal
  leftward relocation is checked against its inverse rightward move; adjacent
  transpositions are rejected.
- Orders **1 through 7**: production ordinary representatives, every tc value,
  merged identity, self-converse flag, and the complete coloured endpoint-line
  tables agree with independent output (node names normalized before comparison).
- Orders **1 through 8**: `production_direct_audit.cpp` instruments the actual
  production LUT code. For every mask and every proposed move it checks the full
  consecutive-arc path predicate; for every legal move it compares the actual
  LUT destination with full pairwise re-encoding. Every reflection output is
  independently compared with explicit reversal/converse. Its complete census
  output equals production output.
- Order **8**: the gzip ledger has 438374 coloured endpoint records representing
  1048576 literal lines. Every merged-node degree, counting a loop twice, equals
  its cube-fibre mass. Exact totals: 2097152 half-edges, 6880 ordinary classes,
  3528 merged classes, 176 self-converse ordinary classes, 2048 blue lines,
  1046528 black lines, and 622 loop lines.

Ordinary counts at n=1,...,7: **1, 1, 2, 4, 12, 56, 456**.
Merged counts: **1, 1, 2, 3, 10, 34, 272**.
Self-converse counts: **1, 1, 2, 2, 8, 12, 88**.

All Python correctness guards use explicit runtime checks; `python3 -O` does not
remove them. Optimized negative controls corrupt one node mass and one line
multiplicity in separate temporary copies; both are rejected. Production C++
uses exception-producing `need`, not removable assertions.

## 4. Binary relation-space interface and smallest hostile

For any finite graph and any field, the edge-difference span has codimension
its number of components. Apply this to the relocation graph to obtain A000568
without a pre-supplied count. Over F2 the generators are e_x+e_y. The kernel of
edge-difference **evaluation on functions** is the algebra of functions constant
on components, with pointwise multiplication.

Keep two constructions distinct: the quotient by the edge-difference span is a
vector-space quotient; that span is not generally a pointwise ideal. Indeed,
e_x times (e_x+e_y) is e_x, which need not lie in that span. The invariant-function
algebra is instead its annihilator/kernel of edge differences.

The F2 graph Laplacian is **not** a replacement for the incidence constraints.
The smallest actual relocation hostile is n=4. In bit order
((0,2),(0,3),(1,3)), its nontrivial ordinary component is
{0,1,2,4,5}, with edges
{(0,1),(0,4),(1,2),(1,5),(2,4),(4,5)}.
It is K_(2,3), with shores {1,4} and {0,2,5}. Incidence rank is four; Laplacian
rank over F2 is two. The indicator of {0,2} is Laplacian-harmonic but is not
constant on that component. The whole cube therefore has four ordinary classes
but F2-Laplacian nullity six. At n=5 the figures are 12 versus 25; at n=6 they
are 56 versus 158. Orders n<=3 have no discrepancy. These are simple undirected
relocation graphs, omitting collapsed loops and duplicate edges.

## 5. Corrections and scope

The original draft's direct-encoding O(n^3 2^m) bound was corrected before
acceptance: there are O(n^2) possible moves and O(n^2) output chord bits per
straightforward re-encoding, giving conservative O(n^4 2^m), up to union-find
overhead. The fixed n<=8 LUT implementation is valid and deliberately bounded.
No graph-isomorphism call enters the production algorithm; the independent
small-order permutation oracle is an **audit**, not a generation dependency.

No flaw was found in the inherited proof, the corrected count-free recognition
statement, the bit mapping, or the reported finite class/line data. The targeted
agent-message command was unavailable because the sparse checkout lacks
`.machine-id`; no registry was created or altered. This did not block any audit.

## Reproduction

From this audit directory, set `OUTPUT` to the production `output/` directory:

```sh
python3 -O independent_check.py --compare-dir "$OUTPUT" > independent_check.out
python3 -O recognition_hostiles.py
python3 -O binary_laplacian_hostile.py
c++ -O2 -std=c++17 -Wall -Wextra -pedantic production_direct_audit.cpp -o /tmp/session67-direct-audit
/tmp/session67-direct-audit 8 > production_direct_audit.json
python3 -O compare_outputs.py "$OUTPUT"
python3 -O optimization_negative_controls.py "$OUTPUT"
```

The instrumented C++ source is self-contained. Production must supply all node
tables through eight, line tables through seven, `census.json`, and
`lines_8.tsv.gz`. The independent Python oracle defaults to n<=6 while cube
checks default to n<=7; these bounds are deliberate and recorded.
