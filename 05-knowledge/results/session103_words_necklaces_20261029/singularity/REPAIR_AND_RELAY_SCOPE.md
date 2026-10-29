# A singular native observer repaired by a unique-cycle vertex extension

**FINITE-EXACT counterexample to unrestricted heredity.** Start with the singular six-vertex child certified in `MARKED_HOSTILE_KERNEL_CERTIFICATE.md`. Add vertex 6 with

    6→{0,1,5},       {2,3,4}→6.

The seven-vertex tournament still has exactly one Hamiltonian cycle,

    0→1→2→3→4→6→5→0,

but its 42×42 native two-run observer has determinant −1. Thus the parent is singular with rank 26 and Smith form 1^26 0^4, while the child is unimodular with Smith form 1^42.

Consequently:

> Singularity of the native M2 observer is **not preserved under arbitrary one-vertex extensions that retain a unique Hamiltonian cycle**.

This is **not** a counterexample to the now-proved be/ee criterion. The repaired child adds a vertex of indegree three, while every actual letter adds a vertex of indegree one. The [last-bad-pair bridge proof](../bridge/BRIDGE_SINK_STAR_LEMMA.md) retains a kernel under exactly the legal suffix operations, not arbitrary unique-cycle extensions.

## Bounded transitive-relay probe

The original six-vertex tournament has two directed triangles, L={0,1,5} and H={2,3,4}. Every interblock edge points H→L except the single bridge 1→2. Insert a transitively ordered relay block R={6,…,5+k}, directed in increasing label order, with all edges H→R→L, retaining the bridge 1→2 and the original two triangles.

Independent subset-DP and direct permutation/path-cut observers give:

| k relay vertices | order | Hamilton cycles | det M2 |
|---:|---:|---:|---:|
| 0 | 6 | 1 | 0 |
| 1 | 7 | 1 | −1 |
| 2 | 8 | 1 | −1 |
| 3 | 9 | 1 | +1 |

This bounded probe suggests the first relay removes the observed kernel, but **no all-order determinant or relay-stability theorem is claimed**. The three positive-relay rows are finite exact data only. The one-vertex repair itself is a complete finite counterexample to unrestricted singularity heredity.

## Reproduction

    python3 -O session103_words_necklaces/singularity/check_repair.py

The checker verifies the repair with three paths: subset-DP observer, direct permutation/path-cut observer, and the proved universal insertion formula. The relay probe compares the first two paths at every listed order. Exact input adjacency, input SHA-256, unique cycle, and bounded results are saved in `repair_receipt.json`. Checks use explicit exceptions. No publication was performed by this lane.
