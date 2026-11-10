---
id: THM-5122
title: "Shared-core exchange intervals, negative-inertia capacity and unit-step interaction selection"
status: "PROVED scoped criteria with independently audited finite physical certificate; global coverage OPEN."
session: Session246, 2026-11-10
source: 05-knowledge/results/session246_unifying_principles_20261110/switches/PROOF.md
audit: 05-knowledge/results/session246_unifying_principles_20261110/audit/AUDIT.md
---

# Cancel interactions inside one physical deletion core

For a nonempty tournament C with h=H(C), append a sink and then a new root
beating exactly a cut x of C. The inherited protected-ear formula is

    F(x)=h+2 a^T x−(1/2)x^T W x,

where a records spanning two-path covers by second start and W is the
symmetric endpoint-pair cover matrix with zero diagonal. This identity is
inherited THM-4114 (`THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature.md`), not a new counting formula.

For disjoint pairs (p_i,q_i), a background cut with p_i selected and q_i
unselected, and d_i=e_(q_i)−e_(p_i), every binary vector b gives an actual cut.
Its count is exactly

    F(x0+sum b_i d_i)=F(x0)+sum g_i b_i+sum_(i<j)c_ij b_i b_j,
    g_i=F(x0+d_i)−F(x0), c_ij=−d_i^T W d_j.

**Affine interval criterion.** The cube is affine iff all off-diagonal
pairings d_i^T W d_j vanish. In that case put
B=F(x0)+sum_(g_i<0)g_i and w_i=|g_i|/2. Discard zero weights and sort.
Its exact image fills all odds in [B,B+2 sum w_i] iff

    w_1=1, w_j<=1+sum_(i<j)w_i.

The empty positive-weight list gives a singleton. Descending greedy subset
selection recovers a physical cut for each target in this computed interval.
The baseline B is not freely specified.

**Capacity.** For i!=j the inherited odd-cycle expansion gives
W_ij>=2H(C−{i,j})>0, with H(empty)=1. Thus d_i^T W d_i=−2W_(p_i q_i)<0.
An orthogonal k-exchange family is negative definite, not isotropic, and

    k<=n_−(W), k<=floor(|C|/2).

**A literal shared-core witness.** The six-vertex deletion core has row masks
[10,28,49,36,41,3] and h=43. From cut {0,2,5}, exchanges 2→4 and 5→3
are interaction-orthogonal and each adds two. The four strong eight-vertex
closures have counts 345,347,347,349. Deletion and root-end counts all remain
43, but their response U is respectively 411,423,407,419. Therefore no factor
43 taxes these local increments, while the unequal response still prevents
unjustified continuation-state identification.

**Unit-step interacting extension.** Write g_i=2u_i,c_ij=2v_ij. The normalized
quadratic f(b)=sum u_i b_i+sum_(i<j)v_ij b_i b_j changes by at most one on
every cube edge iff, for every i,

    −1<=u_i+sum_(j:v_ij<0)v_ij,
    u_i+sum_(j:v_ij>0)v_ij<=1.

Its integer image then fills [min f,max f]. These inequalities imply
sum_j |v_ij|<=2, so the interaction graph consists of paths, cycles, isolated
vertices, and isolated weight-two edges. The field inequalities remain
necessary; degree two alone is insufficient. Two-state dynamic programming
computes extrema in O(k) arithmetic operations; a path between extremizers
selects any intermediate integer in at most k flips.

The formal path/cycle optimization is not a theorem that arbitrary such
coefficients embed in tournaments. Growing physical embeddings with controlled
absolute baselines remain the missing global step. Prior THM-4960/5000
already provide arbitrarily long intervals by other moving-core constructions;
no new global range, scalable overlap theorem or Chan–Pak closure is claimed.

[Full proof, all-cut check and independent permutation certificate](../../05-knowledge/results/session246_unifying_principles_20261110/switches/PROOF.md).
