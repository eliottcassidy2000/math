# Exact equality at the physical Z-defect floor

**PROVED, independently audited.** For a nontrivial positive-start rooted tournament, put B=h+L-s and let D=B(ZU)-2B=R+b-e-B. Then D>=s. Equality at start s=t is physically realizable if and only if t is the Hamilton-path count of some nonempty tournament. Every attainable t has a strong equality witness.

The necessity proof gives an actual nonempty deletion residual with Hamilton count t. The sufficiency construction takes v->a->T->v and has B=3C_2(T)-t, s=D=t. The packet also gives an exact cover-based test for an individual equality graph.

This does not freely prescribe B and does not establish all-target staggered occupancy. F(8,137) and the eventual threshold 23 remain OPEN. Later complete graph normal-form and sparsity results are separate and are not part of this release.

The singleton boundary is explicit: its formula-defined D is zero, but its physical difference B(ZU)-2B is one. The physical identification is asserted only for nontrivial inputs. The reviewed pre-correction proof, exact boundary correction, approved corrected proof, independent audit and status-only promotion are all preserved.

Read PROOF.md and audit/REVIEW.md. Immutable published dependency provenance is in dependencies/PROVENANCE.json.

Reproduce the producer controls:

    g++ -O3 -std=c++17 verify.cpp -o verify
    ./verify > replay.json
    cmp replay.json NORMAL.json

The normal and debug controls check 202,013 marked inputs through order 6, all 3,670 equalities in that finite universe, 33,867 physical cyclic constructions, and 3,806 physical Z children. The independent audit recounts literal paths and separately checks the exact cover bijections. These controls supplement the all-order combinatorial proof.

The audit directory retains every frozen audit manifest entry and its original manifest. No navigation or theorem-index files are changed by this focused release.
