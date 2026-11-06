# Complete graph structure at the sharp Z-defect floor

**PROVED, independently audited.** This focused packet characterizes every actual positive-start source with D=s, where D=B(ZU)-2B on nontrivial rooted tournaments.

Every such source is strong. Take a nonempty tournament T, add v,a, and orient v->a, T->v, and a->T. There is exactly one additional allowed form: reverse one arc to x->a, where x has both start_x(T)=end_x(T)=0. Equivalently, x lies in an interior strong component of T. These forms are necessary and sufficient, and s=H(T).

The underlying endpoint-cover theorem is exact: L=H if and only if the root has indegree zero, or it has one inneighbor x whose start count after deleting the root is zero. Its proof includes a general two-terminal path-cover lemma.

DENSITY.md counts all physical sources on this equality face. At every fixed start t, only O_t(log²(X+2)) distinct B values up to X occur. The flipped sources reduce to finitely many expressions

    B+t = 2^l (alpha 2^k + beta 3^k - gamma).

This is not merely a finite-seed orbit statement. It excludes the exact equality face as a complete bounded-start strategy. Sources with D>s are not excluded. The separate general bounded-defect extension is not included or promoted here. F(8,137) and the eventual threshold 23 remain OPEN.

Read PROOF.md, DENSITY.md, and audit/REVIEW.md. The original proofs and both status-only deltas are preserved. Immutable published provenance for the inherited equality theorem and counting source appears in dependencies/PUBLISHED_PROVENANCE.json; the earlier singleton correction lineage is retained.

Reproduce the producer controls:

    g++ -O3 -std=c++17 verify.cpp -o verify
    ./verify > replay.json
    cmp replay.json NORMAL.json

The producer checks all 202,013 roots and 502,170 terminal pairs through order 6, all 3,670 sharp equalities in that finite universe, 33,867 baseline cyclic sources, and 10,694 flipped sources. The independent audit also checks literal lost/gained path and cover sets, marked-singleton cases, empty retained cores, and 72 larger graph templates through order 15. These are controls of the all-order proofs.

Every frozen independent-audit manifest entry and its original manifest are preserved. This focused release changes no navigation or theorem-index files.
