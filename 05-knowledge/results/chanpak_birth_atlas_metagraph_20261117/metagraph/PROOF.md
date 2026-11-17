# Exact first-birth classes, a deletion/flip relay, and a missing unit link

Let b(h) be the least order of a tournament having h directed Hamiltonian paths. A first-birth class is any isomorphism class T with order b(H(T)); there can be several such classes for the same h. The table `witnesses.tsv` records one representative for each of the 1,520 values with birth at most 9, together with its least order. It makes no statement about the birth of values not appearing by order 9.

## Graph and absence certificate

This is the **ordinary edge-edit/deletion graph**, not the path-framed blue/black tiling metagraph. Its level-n vertices are unlabeled tournament isomorphism classes. Horizontal edges mean a single reversed arc in some labeling; omit loops and count each unordered distinct class pair only once. Vertical edges mean that an order-(n−1) class is obtained by deleting a vertex; omit duplicate deletion choices. A birth is *anchored* if it has a vertical neighbor that is itself a first-birth class at level n−1. A relay is an ordinary horizontal edge from a parentless birth class to an anchored birth class, followed by deletion to that earlier birth. It need not preserve or monotonically approach any absolute target. Converse classes are kept separate unless isomorphic.

The exhaustive algorithm certifies least order by induction. Start with the one-vertex class. From each complete order-(n−1) list take all 2^(n−1) outgoing-arc patterns of a new vertex. Every order-n class occurs because deleting any one vertex leaves an order-(n−1) class. Canonical normalization first sorts by vertex outdegree (an isomorphism invariant), then minimizes the packed adjacency code over *all* within-degree permutations; two graphs have the same canonical code iff isomorphic. Thus filtering duplicate codes loses no possible count. A direct subset DP counts paths on each representative. A value new at level n was absent from every previously enumerated complete level. This is the absence certificate; the independent verifier additionally checks that every table row is a valid tournament with the advertised count.

| order | all classes | new values | first-birth classes | no birth parent | simple horizontal arcs | birth↔birth horizontal arcs | birth↔old horizontal arcs |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| 3 | 2 | 1 | 1 | 1 | 1 | 0 | 1 |
| 4 | 4 | 1 | 1 | 0 | 5 | 0 | 3 |
| 5 | 12 | 4 | 6 | 0 | 30 | 9 | 12 |
| 6 | 56 | 12 | 33 | 0 | 290 | 153 | 71 |
| 7 | 456 | 58 | 328 | 0 | 4086 | 2754 | 616 |
| 8 | 6880 | 243 | 5457 | 2 | 91161 | 68454 | 10112 |
| 9 | 191536 | 1200 | 163676 | 0 | 3380751 | 2777815 | 249823 |

The order-3 parentless case is the seed: there are no births at order 2. All observed birth classes of order ≥3 are strong. At orders 3–9 the entire birth-induced ordinary horizontal graph is connected (the first two have one vertex each). Nevertheless a count-increasing flip strategy gets trapped below the top: at order 6 the birth rows (32,17,3,7,13,30) have H=37, with birth-neighbor count set {17,27,31,33,37}, whereas the top count is 45.

## Both order-8 orphans have a one-flip relay

The two H=63 birth representatives are

    (8,9,3,148,7,31,63,119),
    (8,1,3,198,15,31,55,119).

Both are strong; their deletion counts form the same set {1,11,25,33,37}, all born before order 7. Each has exactly two horizontally adjacent anchored birth classes, one of count 163 and one of count 211.

| orphan representative | reversed arc (zero-based) | intermediate H | deleted vertex | final H, born at 7 |
|---|---|---:|---:|---:|
| first | {0,7} | 163 | 1 | 49 |
| first | {0,6} | 211 | 1 | 59 |
| second | {0,7} | 163 | 6 | 49 |
| second | {1,7} | 211 | 6 | 59 |

In the first two rows the unflipped deletion at vertex 1 had H=25. In the last two rows the unflipped deletion at vertex 6 also had H=25. The order-7 anchored edit can therefore be seen before or after deletion, since the removed vertex is not an endpoint of the flipped arc. In general deleting an endpoint of a flipped arc erases that flip, whereas deleting any other vertex keeps exactly the same one-arc edit. Thus the successful relays are literal commuting squares at the edit/deletion level. This is a repair of *genealogy*, not a count-63-preserving construction and not a universal theorem for future levels.

## Adjacent birth values need not have physically adjacent birth classes

At order 5 all three adjacent newly attained value pairs have at least one horizontal unit edge; at orders 6,7,8,9 the counts are respectively 7/8, 31/48, 186/229 and 1104/1167. The first failure is 25 versus 27, both born at order 6. Each value has two classes. The H=25 rows are

    (2,24,3,37,13,23),      (4,9,10,49,7,23),

both with sorted scores (1,2,2,3,3,4). The H=27 rows are

    (32,9,3,5,15,30),       (32,1,19,7,11,30),

with sorted scores (1,2,2,2,4,4) and (1,1,3,3,3,4), respectively. Reversing an arc changes exactly its two endpoint degrees by ±1. Therefore from a H=25 class the only candidates that could give the first H=27 score list reverse the edge between its degree-three vertices; for the second H=27 list one must reverse between its degree-two vertices. These two candidate flips give 19 and 23 in the two orders for either representative. Thus adjacent scalar births do not supply a one-arc +2 source. The exact toy atlas separates these two phenomena without conflating their metagraphs.

## Pinned larger-order assets and their scope

- The [order-11 first-creator audit](https://github.com/eliottcassidy2000/math/blob/eea6f39cc271ab196ca99b422ad86fb320d8b903/05-knowledge/results/first_creator_generation_20261108/audit/aud_o11_cls.out) checks 506,912,096 qualifying **parent/attachment-pattern pairs** and zero no-first-parent candidates; that number is not a count of first-birth isomorphism classes. Its [lineage table through order 10](https://github.com/eliottcassidy2000/math/blob/eea6f39cc271ab196ca99b422ad86fb320d8b903/05-knowledge/results/first_creator_generation_20261108/lineage_le10.txt) covers 7,661 values.
- The audited [order-12 Genesis packet](https://github.com/eliottcassidy2000/math/blob/6f44021adc666340fbe3587c9a8cb1df91bc634f/05-knowledge/results/genesis_20261109/README.md) gives 12 certified non-creator escaper classes in an audited extension family and 193 found overall for a different generation grammar. These are not an order-12 birth inventory or counterexamples to birth-parent existence among order-12 births.

Replay with `g++ -O3 atlas.cpp -o atlas`, then `./atlas 9 --witnesses > witnesses.tsv` for the complete table. Use `./atlas 9` for graph statistics and the raw first-birth-class records through order 8. `python verify.py` checks all 1,520 witnesses and the relay/absent-unit controls independently of the generation program.
