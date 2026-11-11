# A literal electrical network attached to a tournament

**PROVED elementary consequences of the weighted-path construction.**
The finite controls below test the new realization constraints, not inherited
tournament results. No literature-priority claim. Chan–Pak remains OPEN.

## 1. From a path family to an actual Laplacian

For n>=2, let P range over all H directed Hamiltonian paths of T. Forget each
path's direction and give every edge unit conductance. The user's shadow has

    W = sum_P adjacency(P),    L_W = sum_P L_P.

This is an explicit map to an electrical network on the SAME vertices. It is
not a correspondence between unrelated graphs having equal scalar counts.
Converse preserves W and every network observable below. All weights are
ordinary positive conductances; zero-weight edges are omitted.

This creates a literal bridge beyond the analogy in
`01-canon/theorems/THM-4986-the-7-21-omission-proof-has-an-exact-spanning-tree-twin-block-floor-crevasse-lemma-and-the-first-creator-contrast.md`.
That note expressly claims no map between particular graphs and tournaments.
Here the user's W is that map, but its spanning-tree count is NOT H. We do not
transfer the simple-graph spanning-tree spectrum to tournament path counts.

The path-Laplacian mechanism in THM-4609 is a second historical inspiration,
not an imported Collatz assertion. THM-588's Laplacian is instead on the
metagraph of tournament classes; its connectivity4 does not describe L_W.

## 2. Spectral floor with a rigid equality case

**Theorem G1.** The second eigenvalue of L_W obeys

    lambda_2(L_W) >= 2H(1-cos(pi/n)),

with equality iff H=1, equivalently T is transitive.

Proof. For every real x perpendicular to the all-one vector,

    x' L_W x = sum_P x' L_P x
             >= H · 2(1-cos(pi/n)) · ||x||².

The path eigenvalue follows directly by substituting the cosine vector in the
three-term difference equation. If equality holds, a common minimizing vector
is a Fiedler vector for every path. Along a path that vector is a nonzero
multiple of cos((2j-1)pi/(2n)), j=1,...,n, whose coordinates are strictly
ordered. Its coordinate ordering specifies the path up to reversal. Thus all
Hamiltonian paths have the same undirected edge set. A tournament cannot
contain both directed orientations of that set, so H=1. Conversely a unique
path gives equality. The standard SCC factorization and Hamilton-cycle
rotations show H=1 exactly for transitive tournaments. QED.

This uses the path family before taking its average; entrywise signs or a
generic symmetric matrix would not justify the same argument.

## 3. A strict spanning-tree floor

Let tau(W) be the weighted spanning-tree count, the sum over spanning trees
of the product of their W-edge weights.

**Theorem G2.**

    tau(W) >= H^(n-1),

with equality iff H=1. In particular, the normalized tree excess
log(tau(W))-(n-1)log(H) vanishes exactly for transitive tournaments.
This is a structural diversity statistic, not a claimed Shannon entropy.

Proof. Delete the same row and column from every path Laplacian, obtaining
positive-definite matrices R_P. The matrix-tree identity gives det R_P=1.
The function log det is strictly concave on positive-definite matrices: on a
line X+tD its second derivative is

    -tr((X+tD)^(-1) D (X+tD)^(-1) D),

which is strictly negative for nonzero symmetric D. Hence

    log det((1/H)sum_P R_P) >= (1/H)sum_P log det R_P = 0.

Rescaling and using the same matrix-tree identity proves the inequality.
Equality forces all R_P equal. The missing row and column are recovered from
zero row sums, so all path edge sets agree. As above, H=1. QED.

For C3, H=3 and tau(W)=12>9. The four-vertex strong tournament with H=5 has
tau(W)=232>125. These numbers show why the literal map does not preserve the
original counting invariant.

## 4. Electrical response bounds

**Theorem G3.** For distinct u,v, write d_P(u,v) for their distance along P.
The effective resistance in the W-network satisfies

    R_W(u,v) <= (1/H²) sum_P d_P(u,v) <= (n-1)/H.

Proof. On every unit path send one unit of current from u to v. Average those
H flows. The result is a unit flow in the W-network. At each edge,
Cauchy–Schwarz bounds its energy by the sum of the H contributing squared
path currents divided by H². Summing gives sum_P d_P/H². Effective resistance
is the least energy of a unit flow, proving the claim. QED.

This gives a directly composable network sidecar. It does not recover the
oriented rooted responses consumed by the Fano switch or by an apex insertion.
Equal shadows have exactly equal electrical responses while those physical
tournament operations can have different scalar outcomes.

## 5. A concrete inverse-realization barrier

For every path shadow, total edge mass is (n-1)H, each vertex degree lies in
[H,2H], and every nontrivial cut has capacity at least H. The precise tight-cut
reconstruction in `../shadows/PROOF.md` strengthens these constraints.
Yet the formal matrix 3 times the adjacency of the four-vertex path satisfies
all these linear constraints, the spectral lower bound, and the tree lower
bound with H=3. It cannot be a tournament shadow: equality in either theorem
forces H=1. Its positive support also admits only one possible directed
Hamiltonian path in any orientation compatible with a tournament.

Thus an inverse-shadow strategy must preserve the distinct compatible path
family, not merely a point of the Hamiltonian-path polytope or its integral
first moment. This is a small explicit hostile for turning necessary weighted
constraints into a Chan–Pak construction.

## Reproduction

Run `python3 05-knowledge/results/session265_converse_dyadic_20261111/geometry/replay.py`.
The program enumerates actual Hamiltonian paths and uses exact Bareiss
determinants. It checks all1098 labelled tournaments at orders2..5 and128
seeded cases each at orders6,7, for1354 cases total. The tree inequality,
equality classification, total edge mass and degree bounds pass. These finite
checks do not supply the all-order proofs above. No floating eigenvalues are
used as proof evidence.
