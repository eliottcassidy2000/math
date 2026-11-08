# Sharp factor boundary for the entire ordered-triangle deletion class

Status: PROVED by the path-cover decomposition below; independent audit requested. This is a boundary on a specified physical construction class, not on all tournaments with a given deletion count.

Let A be a tournament with a marked root a, and suppose A−a is an ordered sum of t cyclic triangle blocks and any number of singleton blocks. The adjacency between the root and every block is arbitrary. Then

    3^max(t−2,0) divides H(A).                         (1)

For every t≥2, the exponent t−2 is best possible in this class. In particular every member with t≥3 has scalar count divisible by3. Thus unrestricted triangle placements cannot bypass this prime-factor boundary by allowing the interactions excluded from MULTI.md.

Proof of divisibility: delete a from each Hamilton path to obtain an ordered path pair(P,R), where P is empty or ends at an inneighbor of a, and R is empty or begins at an outneighbor of a. Fix the block i containing the first vertex of R, when R is nonempty, and the block j containing the last vertex of P, when P is nonempty. These cases are disjoint and exhaustive.

At most two blocks, i and j, carry root-neighborhood endpoint constraints. Every other block contributes one of the following unrestricted local counts. If only P or only R can meet the block, it contributes H(block). If both can meet it, it contributes Q(block), the ordered two-bin count allowing empty bins. A nonexistent configuration contributes0. To see that these are all cases, before the first R block every vertex must belong to P, after the last P block every vertex must belong to R, and strictly between i<j both paths may receive arbitrary block subsets. If j<i, the paths instead occupy a prefix and suffix, and any gap between them makes the configuration impossible. Coincident endpoint blocks and empty-path cases merely reduce the number of constrained blocks.

Local path pieces concatenate uniquely in block order. A triangle's unrestricted counts are H=3 and Q=12, both divisible by3; singleton counts impose no requirement. Each nonzero case therefore has at least t−2 unrestricted triangle factors when t≥2. Summing proves(1).

Sharpness: use the literal two-triangle graph in MULTI.md's adjacent-interaction hostile: deletion consists of two ordered cyclic triangles, the root beats one vertex in each, and H=41. Prepend t−2 cyclic triangles to this entire rooted graph, directing every edge from the new prefix toward the old graph, including the root. Hamilton paths concatenate uniquely across the ordered sum, so the new count is3^(t−2)·41 and has exactly3-adic valuation t−2. Its root-deletion graph is still an ordered sum of t cyclic triangles. This realizes equality in(1) at every t≥2.

For the more restricted separated family in MULTI.md, formula(1) there gives the stronger common factor3^(t−1). If3 does not divide t, that valuation is exactly t−1, since the remaining factor is3H_r+sum(4^j−2)≡−t(mod3). These facts quantify why exponentially many same-sidecar choices alone do not establish all-odd scalar coverage.
