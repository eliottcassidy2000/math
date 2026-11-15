# Exact mask energy and a Bézout +2 with the continuation matched

This is selective recontact on the literal same graph: change which existing vertices a fixed apex beats. No deleted path weight, negative copy of a tournament, or modular replacement of absolute counts is used. The construction ends with the full right interface held fixed.

## 1. A quadratic normal form inside a physical continuation fiber

Fix a nonempty tournament R and write h=H(R). For each vertex i let S_i be the sum, over nonempty paths starting at i, of H of their complements; let d_i be its odd-minus-even version. For i≠j set
```
K_ij = 2 sum over paths i→...→j of even vertex order of H(R−vertices(path)),
K_ii=0,    W_ij=K_ij+K_ji.
```
If the apex beats exactly mask X, the [odd-packing identity](https://arxiv.org/abs/2307.05569) gives
```
N(X) = h + sum_(i∈X,j∉X) K_ij
     = h + S(X) − d(X) − sum_(i<j in X) W_ij,
S(X)=sum_(i∈X) S_i,    d(X)=sum_(i∈X) d_i.
```
Indeed an odd cycle using the apex consumes one even-order spine path from i to j; every disjoint residual cycle packing contributes through H of its complement. No packing can use two apex cycles. Also sum_j K_ij=S_i−d_i. Consequently **inside each exact continuation fiber the count is an induced-edge energy**, with positive even edge weights, not an arbitrary lookup table. The right state is the literal (h,S,d) of the [cut transfer](../chanpak_selective_target_20261114/CUT_STATE.md).

Take two templates X,Y with the same (S,d). If Z is a freely chosen common set disjoint from both, the resulting gap is exactly
```
N(X∪Z)−N(Y∪Z) = N(X)−N(Y) − sum_(j∈Z) sum_i ([i∈X]−[i∈Y]) W_ij.
```
Thus for a fixed toggle the gap problem is precisely a signed subset-sum on unused vertices, not an inference from the distribution of marker counts. Different Z choices can have different continuation states; within each pair they match. For several disjoint toggles write each signed difference as a column vector δ_i; their independent additive contributions require δ_i^T W δ_j=0 for every pair i≠j. This is necessary and sufficient, by comparing two-coordinate differences of the quadratic polynomial. It isolates the first physical return obstruction.

## 2. A real crossing interaction and a three-count same-spine window

Use the spine A with outgoing rows (82,0,35,87,6,27,54). Its h is 29. These masks all have (S,d)=(162,0):
```
mask 48: {4,5}, N=175
mask 36: {2,5}, N=177
mask 80: {4,6}, N=169
mask 68: {2,6}, N=179
```
Replacing mark 4 by 2 adds 2 in the first column but adds 10 in the second. Equivalently the mixed second difference is 8. All four masks retain the same continuation, yet the two toggles cannot be treated as independent. Three of the output counts, 175,177,179, are adjacent odds on the **same** spine, no alignment width change needed. Every common transitive tail preserves that three-count gap pattern. If a strong output is needed, first append a common painted 0: the marked apex reaches the first SCC of A, and the end guard returns to it; the initialized counts are 337,339,341 on the three selected arms.

## 3. A real Bézout bridge: two larger corrections cancel to exactly +2

To combine with A, use the regular order-five spine B with outgoing rows (6,24,18,5,9). Masks 3 and 5 give N=43 and 41. They have identical h=15,Q=120, marked starting totals (S,d)=(42,6) and *unmarked* ending totals (T,e)=(63,9). These are exactly the interface needed for inserting a polarized block into an ordinal context. They can be checked independently with path DP or from regular rotational symmetry and equal mask size.

For independent painted choices on A→B with one apex over the entire spine, the physical transfer gives
```
N_AB = h_B N_A + (N_B−h_B)h_A + T_B S_A + e_B d_A .
```
Since both relevant interfaces hold, its two local gap weights are h_B and h_A. More generally a block palette of signed local half-gaps D_i in an ordinal product with total spine count P contributes sums of (P/h_i)D_i when the changed blocks preserve their incoming and outgoing response, or only the outgoing response for the first block. An absolute unit step asks for a literal bounded Bézout solution in those palettes. The gcd is only a necessary condition.

Here use A masks 48 and 68 (local change +4) and B masks 3 and 5 (local change −2):
```
15·4 − 29·2 = 2.
```
The two complete same-spine selective realizations are
```
(A mask 48, B mask 3): N=13643
(A mask 68, B mask 5): N=13645
common (h,S,d)=(435,20658,174).
```
The combined spine is fixed; only its apex contacts change. Both outputs remain exact adjacent counts after **any** common transitive continuation, even though the construction uses two different cyclic reservoirs and opposite changes. A common final 0 yields strong seeds at 34475 and 34477, with (h,S,d)=(435,41316,0). Thus their future exact counts are
```
34475 + 41316 t + 435 p(t) + {0,2},    t≥0,
```
with the usual [canonical suffix address](../chanpak_collatz_exact_addresses_20261114/BINARY_SUFFIX_DECODER.md). It is a unit bridge made of non-unit contributions, not an indefinite additive orbit at one initial target.

## 4. Exact limitation of the method

For independent later ordinal switches, replacing the same local reservoir repeatedly still multiplies an early gap by the ordinary count of the right side. A common transitive tail preserves the gap; a genuinely cyclic common ordinal block multiplies it. Neither this Bézout witness nor zero mixed interactions on a selected switch cube guarantees that arbitrary absolute N can be chosen, nor that the pair can widen indefinitely. The next supply problem is to find cutoff-sized matched palettes whose scaled compositions cover a required integer interval rather than one Bézout identity.

Replay: `python replay_recontact_energy.py` compares the quadratic formula to all masks of both literal small spines by subset DP, checks the mixed second difference, and independently counts both large glues and their common final guards.
