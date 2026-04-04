# THM-299: Recursive Preservation of Multilinear Coefficients

**Status:** PROVED (opus-2026-04-04-S3)
**Proved by:** opus-2026-04-04-S3
**Verified computationally:** n=4→5, 5→6, 6→7 at ALL degrees (0 through max)

## Statement

Let H_n(t_1, ..., t_m) be the multilinear polynomial giving the Hamiltonian path count as a function of the tile variables, where:
- The base path is n → n-1 → ... → 1
- Tiles are non-consecutive arcs (x,y) with x ≥ y+2, labeled in explorer order
- t_k = 0 means forward arc (x→y), t_k = 1 means backward arc (y→x)
- m = C(n-1, 2) is the number of tiles

**Theorem.** For all n ≥ 3 and all tile subsets S involving only "old" tiles (tiles (x,y) with x ≤ n), the multilinear coefficient c_S is EXACTLY preserved when passing from n to n+1:

$$c_S^{(n+1)} = c_S^{(n)}$$

Equivalently: the multilinear polynomial at n is a sub-polynomial of the one at n+1.

## Proof

**Step 1: Source lemma.** At tournament size n+1, the "new" tiles are {(n+1, y) : 1 ≤ y ≤ n-1}. When ALL new tiles are set to 0 (forward), vertex n+1 has:
- Base-path arc: n+1 → n (always present)
- Tile arcs: n+1 → y for all y = 1, ..., n-1 (forward = high→low)

Therefore vertex n+1 beats every other vertex. It is a **source** (no incoming arcs).

**Step 2: HP count of tournaments with a source.** In any directed Hamiltonian path, every vertex except the first must receive an arc from its predecessor. A source (which receives no arcs at all) can therefore only appear as the **first** vertex of any HP.

Given vertex n+1 is the first vertex, the remaining HP visits vertices {1, ..., n} in some order. Since n+1 beats all of them, any HP of T_n can be extended by prepending n+1. Conversely, removing vertex n+1 from any HP of T_{n+1} gives an HP of T_n.

Therefore: **H(T_{n+1}) = H(T_n)** whenever vertex n+1 is a source in T_{n+1}.

**Step 3: Mobius inversion.** The multilinear coefficient for old tile subset S ⊆ {old tiles} is:

$$c_S^{(n+1)} = \sum_{T \subseteq S} (-1)^{|S \setminus T|} \cdot H_{n+1}(\mathbf{1}_T)$$

where **1**_T denotes the tiling with t_k = 1 for k ∈ T and t_k = 0 otherwise.

Since S only involves old tiles, every tiling **1**_T in this sum has all new tiles set to 0. By Step 2, H_{n+1}(**1**_T) = H_n(**1**_T) for each such T. Therefore:

$$c_S^{(n+1)} = \sum_{T \subseteq S} (-1)^{|S \setminus T|} \cdot H_n(\mathbf{1}_T) = c_S^{(n)}$$

**QED.**

## Consequences

1. **The multilinear polynomial grows monotonically**: H_{n+1} extends H_n by adding terms involving new tiles only. No existing coefficient ever changes.

2. **Well-defined infinite limit**: There exists a well-defined formal power series H_∞ = lim_{n→∞} H_n, where each coefficient stabilizes at a finite n.

3. **Linear coefficients are immediate**: c_k = 2^{skip(k)-1} is proved for all n by checking any single n ≥ skip(k)+1 (where the tile first appears).

4. **The quadratic coefficient c_{ij} depends only on the tile pair (i,j), not on n**, as long as n ≥ max(x_1, x_2) (both tiles exist). This means all quadratic formulas (same-end, cross-end, disjoint) are universal.

5. **Layer decomposition**: H_n = H_{n-1} + G_n where G_n involves only terms containing at least one new tile (tile with x = n). The layer increment G_n captures the contribution of the new vertex.

## Computational Verification

| Transition | Old coefficients | Preserved | Changed |
|-----------|-----------------|-----------|---------|
| n=4 → n=5 | 6 | 6 | 0 |
| n=5 → n=6 | 35 | 35 | 0 |
| n=6 → n=7 | 200 | 200 | 0 |

All degrees verified (d=0 through d=max at each n).
