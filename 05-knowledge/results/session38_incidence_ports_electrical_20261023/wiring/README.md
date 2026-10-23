# Complete order-eight complement-line wiring

23 October 2026. **FINITE-EXACT.** All 6,880 ordinary tournament classes and
all 2,097,152 fixed-path presentations are enumerated. This is the complete
converse-merged complement-line multigraph, not a blue-only slice, class
signature census, or sampled transport experiment. All parallel lines,
loops, endpoint ports, and reflection twins are retained.

## Inheritance and boundary

The source mechanism is the fixed-path complement/reflection involution of
THM-4731, `01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md`.
The current source/path/target comparison is THM-4739,
`01-canon/theorems/THM-4739-kirchhoff-and-hodge-laws-of-the-tournament-metagraph-intransitivity-currents-second-order-determination-and-the-consecutive-start-port-ledgers.md`.
The incoming `metagraph_bucket_obstruction_20261023` packet has a complete
order-eight class-signature census but explicitly does not replay the
order-eight transport. This packet supplies that missing replay.

Canonical hostile: an integer transportation table can satisfy all class
margins and black parity while changing the native wiring. Therefore this
consumer does not infer wiring from margins: it constructs every actual
source path, its literal complement, and its actual destination. The
endpoint-capacity checks are subsequent conservation tests.

Corrected near miss: simple support edges, literal line multiplicity, and
paired-black multiplicity are different. A doubled black support bridge is
not a bridge of the literal multigraph. Loops contribute twice to degree
and zero to a Laplacian. Least-used sidecar: the complete mask-to-class
assignment, retained compactly so future signatures can be tested without
rerunning isomorphism classification.

## Exact result

| Quantity | Count |
|---|---:|
| Ordinary classes | 6,880 |
| Self-converse classes | 176 |
| Converse-merged nodes | 3,528 |
| Fixed-path presentations | 2,097,152 |
| Blue presentations | 4,096 |
| Literal complement lines | 1,048,576 |
| Blue lines | 2,048 |
| Literal black lines | 1,046,528 |
| Paired-black lines | 523,264 |
| Fundamental blue-plus-paired-black lines | 525,312 |
| Literal loops | 622 = 4 blue + 618 black |
| Fundamental loops | 313 = 4 blue + 309 paired black |
| Directed merged closed-to-open endpoint pairs | 460,699 |
| Ordinary endpoint-port transport entries | 956,339 |
| Undirected endpoint pairs, including loops | 438,216 |
| Distinct nonloop supports | 437,981 |

The full graph is connected. The blue graph connects all 176 self-converse
nodes, with the remaining 3,352 nodes isolated. The black and paired-black
graphs have one 3,525-node component and three isolated nodes, namely
0, 331, and 3357 in this packet's canonical numbering.

There are eight full-graph **support bridges**. Three carry one blue line;
five carry two black lines. Thus the literal full multigraph has three
bridges, while the fundamental graph has eight. Every full support bridge
is incident to a node with ordinary fibre h=H/|Aut|=1, and all eight merged
h=1 nodes occur this way. There are 13 corresponding ordinary classes.
There is no additional interior support bridge at order eight.

| Node pair | Color | Literal multiplicity | Leaf canonical code |
|---|---|---:|---:|
| (0,3376) | blue | 1 | 0 |
| (1,2388) | black | 2 | 2 |
| (3,2540) | black | 2 | 16 |
| (10,1856) | black | 2 | 256 |
| (40,1140) | black | 2 | 8194 |
| (329,2383) | black | 2 | 524290 |
| (331,2432) | blue | 1 | 524304 |
| (2387,3357) | blue | 1 | 67108866 |

The paired-black graph has six bridges. Its extra bridge (30,1161) is
bypassed by blue wiring in the full graph. The blue graph's three bridges
are precisely the full literal multigraph's bridges. Complete component,
loop, support-cycle-rank, multigraph-cycle-rank, and bridge inventories for
all five graph conventions are recorded in `receipt.json`.

These are finite statements at eight, not an all-order bridge or
connectivity theorem.

## Universe and canonical labels

A tournament's canonical code is a 28-bit, **big-endian** upper-triangle
word: pairs (i,j), i<j, occur in lexicographic order, with bit one for i→j.
Vertices are first ordered by nondecreasing outdegree; among permutations
within equal-score blocks, choose the minimum word. `class_id` is the
zero-based rank of this code among all 6,880 codes. A merged `node_id` is
the rank of min(class_id,converse_class) among distinct such minima.
These are stable canonical labels, not the historical `electric_n7.json`
IDs or the original input TSV row IDs. The original row ID is preserved
as `input_class_id` in `classes.tsv.gz`.

Fixed-path masks use the **descending** Hamiltonian path
7→6→...→0. Free pairs are enumerated by b=0,...,5, then a=b+2,...,7.
Bit k is one exactly when b→a; bit zero means a→b. All adjacent arcs point
down the fixed path. This free-tile order matches the older electric
metagraph scripts. Complement is XOR with 2^21-1. Reflection sends free
pair (b,a) to (7-a,7-b). A mask is blue exactly when fixed by reflection.
The closed endpoint condition is 0→7, equivalently bit 5 is one. Every
complement line has exactly one closed and one open endpoint, including
when the two endpoint nodes coincide after merging.

For a closed source, a=score(7) is the start score and b=score(0) the end
score. The open destination has endpoint scores (8-a,6-b). Reflection
sends source endpoint scores to (7-b,7-a), and its class to the converse
class. These statements are checked on every presentation.

## Why the enumeration is complete

1. Decode all input representatives as literal tournaments. Independently
   minimize every representative over equal-score block permutations,
   recover its automorphism order, and verify distinct canonical codes.
   There are 403,193 permutations in total. The sum of orbit sizes
   8!/|Aut| is exactly 2^28. Since the representatives are pairwise
   nonisomorphic valid tournaments, this certifies the complete ordinary
   class universe without relying on a precomputed class-count claim.
2. Enumerate **every directed Hamiltonian path** of every representative.
   Relabel each path to the descending fixed path and compute its free
   mask. There are 2,143,304 paths altogether. Reject any mask assigned to
   two classes. Every one of the 2^21 masks is reached, and each mask is
   reached exactly |Aut| times from its class. This is also the expected
   free automorphism action on ordered spanning paths: an automorphism
   fixing a full ordered path fixes every vertex.
3. Recover converse classes by separate canonicalization of reversed
   tournaments. Replay every mask, its reflected mask, and its complement.
   Aggregate only the closed endpoint of each line. Preserve color,
   actual source/target classes, and source endpoint scores.
4. Check ordinary Hamiltonian fibres, blue fibres, cyclic-path divisibility,
   all endpoint capacities, all reflection identities, all black twins,
   and all degree/current identities. Only then form graph inventories.

The imported class table is
`04-computation/experiments/tiling_line_metagraph_census_20261020_classes_n8.tsv`,
read from repository snapshot `2130e7633805521d4e24355c37596bd2cd79f044`.
Its exact copied bytes are in `classes_input_n8.tsv`, SHA256
`7d13cd79f9273d903bb8f7e5b34f413cd36a2a98cd4b05e9994fe2baa79c29ab`.
The import's H, Aut, SC and blue-fibre fields are independently checked;
none supplies an unverified target class for the actual transport.

## Conserved electrical and involution data

At each merged node, with t total presentations and s blue presentations:

- literal full degree is t;
- blue degree is s;
- paired-black degree is (t-s)/2;
- fundamental degree is (t+s)/2.

Every degree here counts a loop twice. Every black directed multiplicity,
including every loop multiplicity, is even. For every presentation, if c3
counts cyclic triples, the exact current equation is

    c3(target)-c3(source) = start_score-end_score-1.

Summing all actual ports gives exactly the negative-Laplacian convention
sum_v multiplicity(u,v)[c3(v)-c3(u)], and the total current vanishes.
Thus the complete graph agrees with both the pointwise complement law
and the inherited electrical conservation law over the integers. No
large determinant, approximate resistance, or all-order Smith claim is
needed for this certificate.

## Reusable compressed ledger

All files use deterministic gzip (zero timestamp) and UTF-8 TSV, except the
binary map. A compressed archive is not required: each file can be read
independently using the standard library.

- `class_of_mask.u16le.gz`: exactly 2^21 unsigned 16-bit little-endian class
  IDs in ascending mask order, no header; uncompressed length 4,194,304.
  Its uncompressed SHA256 is
  `636f0ae8f69455d4d39143a7c612260cfbb2a9736c91ca67fe2ac3ea65449be1`.
  This is the complete labelled incidence ledger: together with the stated
  mask decoder it recovers every actual edge and every path position.
- `classes.tsv.gz`: class/converse/node labels, original row ID, H, Aut,
  ordinary and blue fibres, closed fibre, Hamiltonian cycles, c3, scores.
  Hamiltonian cycles are directed cyclic orderings, modulo rotation.
- `nodes.tsv.gz`: one row per merged node, its two ordinary class labels,
  SC flag, total/blue/closed presentations, c3, and exact c3 current.
- `edges.tsv.gz`: `source,target,blue,black`. Each row records directed
  closed→open literal multiplicity. Black is not already divided by two.
  For an undirected nonloop support, add the two possible directions;
  for a loop, use its one diagonal entry once as line multiplicity.
- `ports.tsv.gz`: `source_class,target_class,start_score,end_score,blue,multiplicity`.
  Every row is closed→open; target endpoint scores are determined as above.
- `endpoint_capacities.tsv.gz`:
  `class_id,start_score,end_score,blue,closed,open`. These are capacities of
  actual endpoint types, not inferred feasible capacities.
- `components.tsv.gz`: `graph,component_id,node_id` for literal full,
  blue, literal black, paired-black, and fundamental graphs.
- `receipt.json`: all numerical inventories, check counts, and compressed
  and uncompressed SHA256 hashes.

The compressed ledgers total about 8 MB. The uncompressed binary map is
only a convenience for independent auditors and need not be published.

## Reproduction

    python -B build_full_wiring.py
    python -O -B build_full_wiring.py

Both default commands reconstruct everything, compare every compressed
artifact byte-for-byte and compare the complete receipt, without replacing
frozen artifacts. Use `--write` only to intentionally regenerate them.
There are no removable assertions, external graph-isomorphism libraries,
or numerical dependencies. The full initial run took about one minute in
this environment; elapsed time is deliberately not part of the receipt.
The independent audit uses a different, incidence-based classification
route on every fixed-path mask rather than this Hamiltonian-path route.
