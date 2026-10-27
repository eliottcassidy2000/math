# Order nine: a complete blue quotient from a certified black tree

**FINITE-EXACT.** Every tournament class at order nine is included. Every
blue line is counted. A spanning tree certifies connectivity of the entire
black-active locus. The complete black-line multigraph is **not** enumerated.
No all-order connectivity or no-nonleaf-bridge theorem is claimed.

## 1. What the larger test establishes

At tournament order nine, the converse-merged complement-line metagraph has:

| Object | Exact result |
|---|---:|
| Ordinary tournament classes | 191,536 |
| Self-converse classes | 2,752 |
| Converse-merged nodes | 97,144 |
| Blue tilings / blue lines | 65,536 / 32,768 |
| Pure-blue nodes | 6: five leaves and one trivalent node |
| Black-active nodes | 97,138 |
| Black-active connected components | 1 |
| Black components including pure-blue isolates | 7 |
| Full metagraph connected components | 1 |
| Full metagraph bridges | 5, all incident to leaves |
| Non-leaf bridges of the full metagraph | 0 |

Thus no genuine non-leaf bridge or black-active disconnection occurs at
order nine. This extends the inherited complete checks through order eight.
The result is a finite certificate rather than evidence inferred from a
sample of black lines.

After contracting the black components, the blue quotient has seven
vertices. Its hub has 32,760 blue loops. Five pure-blue leaves join the hub
by one line each; the remaining pure-blue node joins it by three parallel
lines. The trivalent node has three distinct mixed neighbours before
contraction, independently agreeing with the all-order branch theorem.

## 2. Inheritance and the information retained

The previous
[bridge criterion](../boundary_support_20261020_metagraph.md) proves that
every full-metagraph bridge is blue and is a bridge precisely after black
components are contracted. This uses the sigma-twin theorem from
[THM-4731, self-lines, twins and depth](../../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md):
each individual black line has a distinct parallel twin.

The computational inheritance is the complete order-eight class table and
the score-block canonicalizer in
`04-computation/experiments/tiling_line_metagraph_census_20261020.cpp`.
The old full order-eight line TSV was not committed, so reading that
missing file is not presented as a replay.

The incoming
[pure-blue classification and branch geometry](../pure_blue_atom_classification_20261020/PURE_BLUE_COUNTS_AND_BRANCH_GEOMETRY.md)
predict five leaves and one trivalent pure-blue node at order nine. Here
these counts are checked independently by comparing each self-converse
class's complete blue mass with its Hamiltonian-path count divided by
its automorphism order.

The incoming
[distance-moment packet](../merged_metagraph_distance_moment_constraints_20261020/PROOF.md)
adds necessary joint Fourier constraints and an explicit rewiring hostile.
It does not make a degree table or a blue graph determine the black
components. The 19 non-leaf blue-only bridges below show why those missing
connections still matter.

Anchor: the complete blue quotient. Niche: a compact connectivity certificate
instead of a complete line census. Wildcard: canonical Hamiltonian-path
choices and their triangle current. Closest mechanism: actual path-level
complement operations. Hostile: confusing a blue-only bridge with a bridge
after black paths are allowed. Sidecar: the Hamiltonian word and explicit
partner isomorphism for each retained black edge.

## 3. A complete class universe without the 268-million-tiling cube

The producer first validates all 6,880 inherited order-eight representatives:
its independently implemented refinement canonicalizer reproduces every
automorphism order and the full labelled mass `2^28`.

It then adjoins a ninth vertex in all 256 ways to each representative.
This is a complete class generator: deleting any vertex of any order-nine
tournament gives an order-eight class, and relabelling that deletion to the
chosen representative turns its incident arcs into one of those 256 words.
The resulting canonical representatives number 191,536.

The independent verifier uses the inherited, different algorithm:
enumerate every permutation within score blocks and take the least arc
word. It confirms that all 191,536 representatives are pairwise
nonisomorphic and reproduces every automorphism order. Their labelled
orbit masses sum to

```text
sum_T 9!/|Aut(T)| = 68,719,476,736 = 2^36.
```

Distinctness plus that mass independently confirms completeness, without
using a known class count as the sole completeness test. Explicit
permutation witnesses authenticate every converse pairing; 2,752 classes
are self-converse.

Canonical words use the upper-triangular arc bits in lexicographic pair
order, with the first bit most significant and `1` meaning `i -> j`.
They are implementation addresses, not unexplained historical class IDs.

## 4. Every blue line and an authenticated black spanning tree

Reflection has 16 coordinate orbits on the 28 free path-arc bits at order
nine. Enumerating their `2^16` assignments gives every blue tiling.
Complement pairs them into 32,768 lines, occupying 31,832 unordered class
pairs. Both canonicalizers independently reconstruct the complete matrix
and every self-converse class's blue mass.

For each merged black-active class, choose its smaller canonical
representative and enumerate Hamiltonian paths lexicographically in those
canonical vertex labels. Every nonsymmetric path gives an actual black
line: reverse all arcs except its eight path arcs, canonicalize the
partner, and merge converse classes.

The first path at each of the 97,138 active nodes leaves 237 components.
Considering up to the first four paths where components still need to
join produces a spanning tree with 97,137 edges. In total the producer
examines 97,719 Hamiltonian paths, of which 97,601 are nonsymmetric.

Each recorded tree edge contains:

1. the source canonical tournament word;
2. the target canonical tournament word;
3. its nine-vertex Hamiltonian order;
4. a permutation explicitly identifying the constructed complement
   partner with that target.

Orders and isomorphisms are packed in four bits per vertex, with the first
entry in the lowest nibble. The independent audit checks every path arc,
checks that the normalized tiling is nonsymmetric, applies the path
complement, verifies the stated isomorphism, and checks that the edges
form a tree on **all** 97,138 active nodes.

This is why the connectivity conclusion is exact even though most black
edges were never generated: a spanning connected subgraph of the full
certified vertex set proves connectivity. No missing black edge can split
that tree. Pure-blue nodes have no black half-edges, so the seven black
components are also determined exactly.

## 5. Nineteen genuine blue-only bridges disappear in the full graph

The complete blue graph on the 2,752 self-converse classes is connected.
It has 24 bridges: the five pure-blue leaf edges and **19 mixed-to-mixed
bridges**. Those 19 have neither endpoint a full-metagraph leaf.

For example, the blue bridge between canonical words

```text
147968 and 17181073669
```

has full fibre masses 5 and 191. It is a non-leaf bridge of the blue
subgraph. The certified black tree connects its endpoints without using
that edge, so it is not a bridge of the full metagraph.

All 19 are listed in `n9_analysis.json`, under
`quotient.blue_only_nonleaf_bridges`. This is an explicit hostile to
classifying full bridges from the blue graph alone. After the correct
black contraction, exactly the five leaf bridges survive.

## 6. What the path choices suggest, and what they do not prove

The second round adds 236 tree edges from 219 source classes:

| Lexicographic Hamiltonian-path rank | Added tree edges |
|---:|---:|
| 2 | 187 |
| 3 | 20 |
| 4 | 29 |

Of those 219 source classes, 115 are self-converse. The full records in
`n9_second_round.json` retain scores, strong-component sizes, automorphism
order, the path, and the triangle-count change

```text
outdegree(first) - outdegree(last) - 1.
```

The choices depend on canonical vertex labels. No intrinsic four-path
theorem or all-order connectivity claim follows from this finite success.
The records are a concrete test bank for such a proposed theorem.

The **certificate tree** has diameter 154 and radius 77; these are not
claimed distances in the full metagraph. Its maximum degree is 49.
Extra black edges can only shorten distances. The tree therefore gives
the finite upper bound 154 on distances inside the black-active locus,
without determining its actual diameter.

## 7. Reproduction and compact certificates

From the repository root:

```text
python3 05-knowledge/results/blue_quotient_20261020b/replay_census.py
python3 05-knowledge/results/blue_quotient_20261020b/replay_census.py --rebuild
```

The first command verifies source hashes, expands the compressed
certificates in a temporary directory, compiles the independent verifier,
rechecks all classes, the complete blue cube and every black-tree edge,
then compares normal and `-O` Python analyses.
The second also regenerates the producer's whole class universe and
certificates and compares them after newline normalization.

Files:

- `census.cpp`: complete augmentation, blue enumeration and black-tree
  production;
- `audit.cpp`: independent score-block canonicalization and explicit
  witness replay;
- `analyse_census.py`: bridge quotient, tree geometry and path-rank audit;
- `n9_classes.tsv.gz`, `n9_blue.tsv.gz`, `n9_black_tree.tsv.gz`: complete
  compact certificates, deterministic gzip with timestamp zero;
- `n9_independent_audit.txt`, `n9_summary.txt`, `n9_analysis.json`:
  exact receipts;
- `census_manifest.json`: owned and inherited source hashes.

Compiled with GCC 15.2.0, `-O3 -std=c++17`. All checks use explicit errors;
no essential Python `assert` is removed by optimization. Text hashes use
canonical LF; archived producer TSVs retain their exact Windows newlines.
Executables and expanded replay files remain in a temporary directory.
