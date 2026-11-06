# Session 208: exact bridge factorization

Status: elementary all-order theorem, independently audited. No literature-priority claim.

PROOF.md proves exact separation of spanning-tree stretch and fundamental-cut congestion over the components remaining after deleting bridges, and rational coupling with additive rather than multiplicative support. Its metagraph application concerns the ordinary isomorphism-class line multigraph, not labelled degree-three vertices.

The cited attached preprint is Mohammad Hajiaghayi, *Logarithmic Embeddings into Probabilistic Spanning Trees for Distances and Capacities*, arXiv:2611.05189v1, https://arxiv.org/abs/2611.05189. Its main theorem is not independently certified by this release. The exact bridge identities do not depend on it; logarithmic corollaries explicitly use the cited statement.

Run `python3 check.py` for exact finite controls with unequal integer lengths/capacities and parallel edges. The mathematical proof, not these controls, supplies the all-order result.

Publication scope: this directory only, under the session208 result directory. No unrelated checkout changes are included.
