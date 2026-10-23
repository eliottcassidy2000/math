# Full order-eight tournament wiring with incidence-based reconstruction

23 October 2026. **FINITE-EXACT:** every fixed-path tiling at order eight,
its ordinary and converse-merged owners, every complement line, and every
closed/open endpoint-color port. **PROVED:** the interpretation of this
exhaustive certificate and the inherited capacity-peeling implication.
**OPEN:** analogous inexpensive determination at arbitrary order.

## Scope and inheritance

This packet upgrades the class-separation result in
`../metagraph_bucket_obstruction_20261023/README.md` to a **full native wiring
replay at n=8**. Its graph is the converse-merged complement-line multigraph,
with loops and literal line multiplicities. It is not a single-edge-flip
graph. The paired-black graph divides each literal black multiplicity by2.

The old S46 census already established the n=8 numerical totals and the
ordinary-class inventory:

* THM-4731, `01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md`;
* `04-computation/experiments/tiling_line_metagraph_census_20261020_classes_n8.tsv`;
* `04-computation/experiments/tiling_line_metagraph_census_20261020_n8.out`.

Those counts are inherited, not newly discovered here. The additions are
an independently authenticated complete native ownership/port/edge bank,
the full incidence-refinement transportation tests, and a literal
counterfeit showing that the weaker one-round constraints remain
insufficient. The prior class census alone was not this full wiring replay.

## 1. Source identity is proved before aggregation

The fixed Hamiltonian path is0->1->...->7. Its21 non-path pairs (i,j),
i+2<=j, are ordered lexicographically. Bit1 means i->j. A mask m represents
exactly that tournament and that path gauge. Its complement-line partner
is m xor(2^21-1). Reflection sends (i,j) to(7-j,7-i). A tiling is blue iff
reflection fixes its mask; otherwise it is black. It is closed iff the
(0,7) bit is0.

`producer.cpp` starts from the6880 canonical ordinary representatives,
not from a guessed owner table. It independently verifies score-canonical
minimality, automorphism orders, distinct canonical words, and the orbit
mass sum2^28. It enumerates every literal directed Hamiltonian path of each
representative, relabels it to the base path, and assigns the resulting
mask to that ordinary class. Every one of the2^21 masks is reached, and
every mask occurs exactly |Aut(T)| times in its class's path enumeration.
No mask is assigned two different owners.

This uses2,143,304 actual paths, as opposed to performing a fresh canonical
search for each of the2,097,152 tilings. Subset/end-vertex DP independently
checks every class's Hamiltonian count. The producer performs10,626,230
explicit checks. The resulting owner map is the lossless sidecar: later
quotients never infer source identity from a scalar statistic.

The incoming Session37 polynomial checksum is also recovered exactly:
`sum_t H(T_t) = sum_class H(T)^2/|Aut(T)| = 815136768`.

The Python consumer independently generates those paths by appending
vertices and adding only the newly created chord bits. It checks every
leaf against the owner map, the exact automorphism multiplicity of every
mask, all reflection targets, all class fibres, and all blue counts.
It then reconstructs every native port and line directly from the masks.

## 2. Complete wiring results

| Object | Exact count |
|---|---:|
| Ordinary tournament classes |6,880|
| Self-converse classes |176|
| Converse-merged nodes |3,528|
| Fixed-path tilings |2,097,152|
| Closed source ports / complement lines |1,048,576|
| Distinct ordinary source-target-endpoint-color groups |956,339|
| Blue lines |2,048|
| Literal black lines |1,046,528|
| Ordinary unordered class pairs joined |911,027|
| Merged unordered node pairs joined |438,216|
| Ordinary blue / black self-lines |4 /202|
| Merged blue / black self-loops |4 /618|

Every entry in both sparse edge tables is checked independently. Every
merged black multiplicity is even. All inherited ordinary multiplicity
histograms and self-line totals agree exactly. A loop contributes twice
to a node's degree/fibre incidence, and once to the edge count.

For a closed port with start/end outdegrees(a,b), complementing the
non-path arcs gives an open target port with outdegrees(8-a,6-b).
The consumer verifies this on every literal port, separately from the
source/target class assignment and the blue/black test.

## 3. Smaller signatures tested on every actual port

The concurrent `../metagraph_minimal_signatures_20261023/` packet found that
only one initial vertex coordinate is needed in the tested menu:

    v4(v) = number of induced strong four-vertex subtournaments containing v.

Start with r0(v)=(v4(v),). Repeatedly set

    r(k+1)(v) = (rk(v), multiset_{v->w} rk(w)),
    Sk(T) = multiset_v rk(v).

These are invariant profiles, not adjacency encodings asserted complete
at all sizes. Our full replay transports Sk of the **actual complement
tournament**, then imposes the actual source row and target endpoint-color
capacities. The source row includes the ordinary source class, its start
and end scores, its color, and the transported target signature. A column
includes the ordinary target class and its actual open endpoint-color type.

| Rounds k | Target class signatures | Source rows | Allowed entries | Unresolved closed-port mass after capacity peeling |
|---|---:|---:|---:|---:|
|0|2,793|915,746|1,914,275|724,751|
|1|6,878|956,171|956,541|884|
|2|6,880|956,339|956,339|0|
|1 plus guarded wedge bit|6,880|956,339|956,339|0|

The sum of forced and unresolved mass is exactly1,048,576 in every row.
With two rounds, all956,339 aggregated source/target port counts are
forced and coincide with the native table; all ordinary and merged line
multiplicities consequently follow. No Hamiltonian count or blue-count
coordinate is used in this target signature.

A smaller adaptive repair from the concurrent signature packet measures
one typed out-wedge count only on the two ambiguous first-round profiles.
For fixed first-round record types R,B,C it is

    W(R;B,C) = sum_{v:record(v)=R} d_B^+(v) d_C^+(v).

The guarded measurement is0 or1 and separates each residual pair. The
consumer imports and pins that packet's public invariant API, then applies
it to every target in the full native port bank. **One round plus this
guarded one-bit measurement forces exactly the same complete wiring as
two rounds.** This retains within-cell covariance that aggregate edge
counts between record cells lose. The guard is essential; no global
one-bit reconstruction assertion is made.

The one-round residual has two target pairs exchanged by converse:
canonical words(27467784,27467793) and(6427781,6427778). Each pair leaves442
units of closed-port mass,185 source rows, and8 target columns. Thus the
combined residual is884,370,16. The two members of an individual pair are
not converses: their Hamiltonian counts are523 and527.

Getting stuck in capacity peeling alone would not prove ambiguity. The
next section supplies an actual alternative feasible table.

## 4. A same-margin arrangement that is impossible for the native sources

In dense merged-node IDs from `merged_nodes.tsv`, change only these
literal black multiplicities:

    (1238,2246):-2, (1238,2315):+2,
    (2242,2246):+2, (2242,2315):-2.

The receipt expands this into eight ordinary-class port edits with
literal native mask witnesses. It is a2x2 switch together with its converse
reflection. Every resulting count is nonnegative. It preserves:

* every one-round source signature/endpoints/color bucket total;
* every target open endpoints/color bucket total;
* every node's tiling, blue, and black counts;
* the class-level reflection relation and paired-black even multiplicities.

Nevertheless it changes the merged graph. Two-round target incidence
fixes the native table uniquely, so this alternative cannot satisfy that
refined source/target compatibility. This directly realizes the owner's
"same buckets, different arrangement" test: totals and parity allow it,
while native incidence correlations exclude it.

The counterfeit is a finite obstruction to this particular one-round
observer. It is not a counterexample to the tournament metagraph itself,
and it does not imply that every conceivable small observer must fail.

## 5. Portable full-data formats

All IDs are stable IDs from the inherited ordinary-class TSV unless
explicitly marked merged. All compressed content is deterministic gzip;
all emitted TSV bytes use LF and binary integers use little endian.

* `owner_map.u16.gz`: exactly2^21 unsigned16-bit values. Entry m is the
  ordinary owner of literal mask m.
* `ports.bin.gz`: sorted records `(key:uint64,count:uint32)`. Decode key
  by reversing `((((source*6880+target)*8+a)*8+b)*2+blue)`. Both a,b are
  the CLOSED SOURCE endpoint scores. Every underlying literal mask can
  be recovered from owner_map, complement, and these filters.
* `ordinary_edges.tsv.gz`: columns `source,target,blue,black`; endpoints
  are sorted ordinary IDs. Loops and every nonzero pair are retained.
* `merged_nodes.tsv`: dense merged ID, its one/two ordinary member IDs,
  their canonical words, total tiling fibre, and blue fibre. A
  self-converse node repeats the same ordinary ID in its two fields.
* `merged_edges.tsv.gz`: the same sparse columns using dense merged IDs
  0..3527. Black counts are LITERAL counts. Divide by2 for paired-black
  conductances; do not divide blue counts.
* `artifact_hashes.json`: compressed and uncompressed hashes/sizes.
* `receipt.json`: exact source provenance, numerical controls, all
  refinement counts, and the native counterfeit witnesses.

These tables permit electrical and parity computations without copying a
3528-by3528 dense matrix or guessing how loops were counted.

## 6. Reproduction and status boundary

Compile and generate in this packet directory (paths shown from repo root):

    g++ -std=c++17 -O2 05-knowledge/results/metagraph_order8_20261023/producer.cpp -o producer
    ./producer 04-computation/experiments/tiling_line_metagraph_census_20261020_classes_n8.tsv 05-knowledge/results/metagraph_order8_20261023
    python -B 05-knowledge/results/metagraph_order8_20261023/compress.py

The executable path may be any temporary location. The producer emits
redundant raw copies for inspection; only their compressed banks need to
be kept. The consumer requires NumPy and reads the compressed banks:

    python -B 05-knowledge/results/metagraph_order8_20261023/check.py
    python -O -B 05-knowledge/results/metagraph_order8_20261023/check.py

Default runs compare the frozen receipt; `--write` explicitly regenerates
it. No assertions disappear under optimization. The source class hash is
normalized to LF so a checkout's line-ending policy is not mathematical
input. A class source, a literal path, its complement target, and every
integer multiplicity are retained. All conclusions are bounded to order
eight; no all-order signature termination or universal connectivity theorem
is inferred. No theorem ID was reserved here.
