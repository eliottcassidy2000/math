# Signed support and the full black core are different graphs

**Status: PROVED IN NOTES after independent root review; see [BLACK_REVIEW.md](BLACK_REVIEW.md).** The contraction and attachment criteria below are elementary all-size conditional graph statements. The listed native tournament consequences are FINITE-EXACT at orders3–7. No all-order connectivity theorem or new theorem ID is claimed.

## Objects and the information loss

Use the native ordinary matrix C, merged dart matrix A, blue part A_B, and signed chiral matrix D from [signed converse transport](../pairings/SIGNED_LIFT.md). Set K=A−A_B. Its diagonal counts twice the number of black loops; parallel black lines retain their multiplicities.

The **full black support** has every merged class as a vertex and an edge u–v for distinct vertices exactly when K(u,v)>0. The **black core** is its induced graph on vertices with positive black dart degree b(u)=sum_v K(u,v). This definition retains a hypothetical vertex with only black loops; it does not mean “delete every isolated vertex.” In the finite census below, however, every isolated vertex has degree0.

The **signed support Sigma** has only the chiral fibres as vertices, including its isolates, and edges where D(u,v)≠0 for distinct u,v. Its components control independent chiral flips of the ordinary matrix. They do not automatically describe components of the full black core.

For a chiral pair of fibres the ordinary block is [[a,b],[b,a]], with a,b≥0. Blue is absent on this sector, so

    K(u,v)=A(u,v)=2(a+b),       D(u,v)=a−b,
    K(u,v) ≥ 2|D(u,v)|.                                  (1)

Thus every signed-support edge is an actual black edge. Equality in the bound holds iff one channel is zero. Conversely a=b>0 erases the signed edge while preserving a positive black edge. Moreover Sigma discards all self-converse vertices and their black connections. These are two distinct losses. Every chiral merged node belongs to the black core: it has no blue darts and has positive total dart degree 2H/Aut.

## Exact contraction criterion

Let S_1,...,S_c be the connected components of Sigma. In the full black core, contract each S_i to one vertex and retain each positive-black-degree self-converse node as an individual vertex. Retain **every** positive K edge between blocks, including edges whose signed difference is zero. A weighted version of the quotient has matrix

    K_bar(B,C)=sum_(u∈B,v∈C) K(u,v),                    (2)

with the same doubled-loop convention. Call its support graph Gamma.

**Conditional theorem.** The full black core is connected iff Gamma is connected. More strongly, contraction preserves the number of connected components.

**Proof.** Every S_i is connected inside the black core by (1); singleton self-converse blocks are also connected. Projecting a black path gives a quotient walk. Conversely a quotient path can be lifted by choosing the black edge witnessing each step and joining successive edge endpoints inside their connected blocks. These two directions identify the component partitions. Equation(2) also preserves each block's total black dart degree. ∎

This uses connected fibres, not merely a partition into chiral and self-converse types. Contracting a disconnected collection in one step would invalidate the reverse implication.

**Attachment corollary.** If Sigma is connected and nonempty, the full black core is connected exactly when every connected component of the black graph induced on positive-black-degree self-converse nodes has at least one edge to the chiral sector.

Indeed the chiral sector contracts to one vertex. An induced self-converse component can reach it only through such an edge, and one edge suffices for all vertices of that component. The condition is vacuous if there are no positive-black-degree self-converse nodes. Connected Sigma alone therefore leaves an explicit attachment obligation; no tournament violating that obligation is asserted here. An abstract black graph consisting of one signed chiral edge and a disjoint positive self-converse-sector edge illustrates the missing hypothesis, without a tournament-realizability claim.

The empty graph is declared **not connected** in this note and checker. The component-count identity remains valid for the empty core, where both counts are0. The attachment corollary explicitly assumes nonempty Sigma.

## Exact finite boundary

The complete frozen matrices give:

| n | full black-support component sizes | positive black-core size | signed component sizes | contraction vertices | positive SC induced components |
|---:|---|---:|---|---:|---:|
|3|1,1|0|empty|0|0|
|4|1,2|2|1|2|1|
|5|1,1,1,7|7|1,1|7|1|
|6|1,1,32|32|22|11|5|
|7|1,1,1,1,268|268|184|85|6|

Every isolated full-support node in these rows has black degree0. The positive black core is connected at4,5,6,7; at3 it is empty. The exact zero-black-degree codes are retained in the receipts.

Order5 decisively refutes the converse “black-core connected implies signed support connected.” Its two chiral fibres, unrestricted codes2 and4, are isolated in Sigma, yet the actual black core is connected. In particular it contains the path

    2 — 8 — 4,

with both edge multiplicities2; the middle class is self-converse. Thus the full black core can connect separate signed components through retained self-converse vertices.

At6 and7 Sigma is connected, but the positive self-converse induced graph has respectively5 and6 components. Every one has a black edge to the chiral sector. Contracting these self-converse components as well therefore produces stars with5 and6 leaves. These stars are finite connectivity certificates, not assertions about later orders. The more detailed quotient in(2), which retains individual self-converse nodes, has respectively11 and85 vertices.

## Verification and scope

[check_black.py](check_black.py) reads only the frozen complete [data.json](data.json). It reconstructs K and D from the ordinary/merged matrices, verifies(1), computes both component systems, verifies that every contraction fibre is black-connected, and checks equality of component counts and quotient dart masses. It checks all self-converse attachment components, the order-five path, the finite component table, and a separately labelled abstract missing-attachment hostile.

The normal and optimized receipts are byte-identical and pass **69,050 always-active checks**. They contain the full quotient matrices and block memberships, all self-converse attachment certificates, source hashes, and the explicit empty-graph convention.

```bash
python3 05-knowledge/results/incidence_closure_20261113/closure/check_black.py
python3 -O 05-knowledge/results/incidence_closure_20261113/closure/check_black.py
```

- Checker SHA-256: `51c4c3245d3bf1db5565d49d64fce0f4b75eea1f7a89428479be05ea130cb2a6`.
- [Normal receipt](black.normal.json) and [optimized receipt](black.optimized.json), each SHA-256: `d737bdea9dca8585d6e62db70f6082214368b713f24de5d0ca4540431e23f819`.
- Frozen matrix input SHA-256: `c0484a99712edff88e8158b5f6000097f32f88874150d059da5aefa20f3de8f7`.

The signed-support connectivity proved in a finite row is a statement about the ordinary symmetry kernel. Full black-core connectivity is a different statement with additional self-converse and zero-current channels. The contraction identifies exactly what still has to be checked to pass from one to the other. Both all-order connectivity questions remain open here.
