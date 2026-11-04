# Complete self-converse order eleven: finite reconstruction result

**FINITE-EXACT.** H,S,E,C separates all labelled self-converse tournaments
of order11. The complete279968 isomorphism classes have distinct sidecar
classes, and each sidecar has exactly the underlying tournament's
automorphisms. Thus the full two-port response separates this entire family.
Universal two-port reconstruction remains **OPEN**. This does not test
nonselfconverse order11 tournaments or pairs crossing between the families;
the previous global lower bound11 for a counterexample is unchanged.

## Complete universe, inherited codec and independently retained symmetries

The source is Session89's proved
[odd self-converse lift](../../session89_consecutive_tournament_pairs_20261028/encoding/ODD_SELF_CONVERSE_LIFT.md),
also used by the consecutive metagraph enumeration. Input is a tournament T
on h vertices, an undirected graph R on the same names, and h root bits d.
Create p_i,q_i and a root, orient p_i->q_i, use T on the p vertices and
its converse on the q vertices, orient p_i->q_j iff R_ij=1 for i!=j,
and root->p_i iff q_i->root iff d_i=1. Swapping each p_i,q_i and fixing
the root is an anti-involution, checked on every ordered vertex pair.

Every odd self-converse tournament has such a marking: its automorphism
group has odd order, so its index-two signed automorphism group has
involutions. An anti-involution fixes exactly one vertex. Orienting each
paired arc p_i->q_i picks its tail intrinsically, leaving only the S_h
action on the pairs. All anti-involutions of a tournament are conjugate
by automorphisms. One elementary verification: if r,s are two such maps,
g=sr is an automorphism of odd order m; rgr=g^(-1), and
g^((m+1)/2) r g^(-(m+1)/2)=s. Thus forgetting the marking gives a bijection
of input and unmarked output isomorphism classes.

Use the complete inherited tournament bank through h=5; for each T enumerate
every R,d and keep exactly the lexicographically least member under Aut(T).
All core automorphisms are independently enumerated through5. The source
stabilizer is the centralizer of the chosen anti-involution, not the full
output automorphism group. The census computes the latter separately.
`bank_audit.py` additionally expands every source orbit by literal vertex
permutations and checks disjoint, complete coverage of all labelled
tournaments through5, with normal/-O matching receipts.

| Output order | Classes | Labelled fixed-frame mass |
|---:|---:|---:|
|3|2|2|
|5|8|16|
|7|88|512|
|9|2752|65536|
|11|279968|33554432|

The mass is sum h!/source-stabilizer =2^(h*h), checked literally. At11 the
centralizer histogram is1:279448,3:512,5:8. Full automorphism orders reach81;
the Paley11 control has order55. Seven classes have exact S=E, so ordered
cover data and symmetry tests cannot be replaced by endpoint sorting alone.

## Exact response sidecar and canonicalizer

Subset DP counts every old Hamilton path with its named start/end. For
every nonempty proper support X,

    C_uv = sum_X E_u(T[X]) S_v(T[V\X]).

Absent-arc entries are included. This is the full ordered two-path-cover
matrix, not the number of Hamilton paths using an existing edge. At11,
32-bit unsigned counts are safe: each C_uv<=10!, endpoints<=10!, and
H<=11!. Arithmetic does not use a modular hash for equality.

The canonicalizer starts with (S_u,E_u). It repeatedly refines by the old
colour and the sorted multiset of (colour(v),C_uv,C_vu). If a cell remains
nonsingleton, choose the first such cell and recurse on **every** possible
individualized vertex, then refine again. Terminal leaves retain all22
endpoint coordinates and all121 ordered cover entries, including the zero
diagonal. Exact143-entry vectors are the map keys.

This search tree is equivariant under every named-sidecar isomorphism.
Branches are never heuristically cut off. All minimum-key leaves are
retained. Given one such leaf, sidecar automorphisms act freely and
transitively on the others: equality of their complete labelled matrices
gives the converse bijection. Therefore their count is exactly Aut(sidecar).
The minimum and maximum underlying tournament codes across those leaves
are equal iff every sidecar automorphism also preserves the original arcs.
Different source classes with the same canonical key would be a collision;
no such key collision or extra sidecar symmetry occurred.

At11 there are290128 terminal leaves,295010 total search nodes and at most81
leaves for one source. The output stores full automorphism and centralizer
histograms, their quotient (the number of anti-involutions), and all failure
witnesses (empty). An independent review and bank audit are in
[INDEPENDENT_REVIEW.md](INDEPENDENT_REVIEW.md).

## The full three-path layer is retained and independently checked

The family census needs only a sufficient observable sidecar. It does not
identify the complete response with that sidecar. In [response](../response/),
literal permutations of all11 old names and every admissible cut are compared
entry by entry with squarefree subset multiplication for the transitive11
and Paley11 controls. Both retain S,E, the joint Hamilton endpoints M, the
full C, both three-index two-path endpoint tensors U,V, and all11^4=14641
entries of the ordered three-path tensor D, including singleton middle paths.

For each control every one of the4489 masks of weight at most2 in each port
is evaluated both by the tensor formula and direct Hamilton DP on the actual
13-vertex tournament. Nine further high-weight queries verify interpolation.
The exact total is8996 actual query checks across the two controls.
An independent partial-bijection search
finds all sidecar automorphisms without colour refinement, and64 relabellings
per control preserve the entire canonical object. Normal/-O receipts agree.

## Reproduction and boundaries

Run from this directory:

```text
python3 run.py
python3 run.py --max-half 4
python3 -O run.py --max-half 4
python3 ../response/run.py
python3 -O ../response/run.py
```

`full11.json` is one full compiled-O3 run of the complete family. Bounded
through9 normal/-O controls are separate. Source hashes, the inherited bank
hash and all complete-family cardinalities are retained. Repeating with-O
changes Python assertion semantics only; validation uses explicit exceptions.
The formulas are universal, but these finite results are not a universal
reconstruction proof. If a later sidecar collision appears, the full
three-path response still has to be compared before declaring a two-port
collision.
