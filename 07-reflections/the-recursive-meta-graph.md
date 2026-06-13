# The Recursive Process in the Merged Meta-Graph

**Session**: kind-pasteur-2026-03-22-S20cf

## The Data

| n | V_merged | E_merged | Blue | Black | Collapsed | SC | NS-pairs |
|---|----------|----------|------|-------|-----------|-----|----------|
| 3 | 2 | 1 | 1 | 0 | 0 | 2 | 0 |
| 4 | 3 | 3 | 1 | 2 | 0 | 2 | 1 |
| 5 | 10 | 21 | 13 | 8 | 0 | 8 | 2 |
| 6 | 34 | 143 | 98 | 45 | 5 | 12 | 22 |

## The Three Growth Mechanisms

As n increases by 1, the merged graph grows through three mechanisms:

### Mechanism 1: New NS Pairs (the bulk)

Most new iso classes at n+1 are NS. They come in complement pairs, each adding ONE merged node. These new nodes connect primarily via BLUE edges to existing NS-pairs (since most of the graph is NS).

At n=6: 22 NS-pairs out of 34 merged nodes (65%). At n=7: this fraction will be even higher.

The NS-pair bulk is where the blue skeleton lives. It grows roughly as (A000568(n) - SC_count(n))/2, which is dominated by A000568(n)/2 since SC_count grows much slower.

### Mechanism 2: New SC Classes (the spine)

SC classes grow much more slowly than NS pairs. The SC count sequence: 2, 2, 8, 12. These form the "spine" of the merged graph -- the axis around which NS-pairs cluster.

SC classes connect to other SC classes via blue edges (forming the blue spine) and to NS-pairs via black edges. At large n, each SC class has MANY black connections (to the NS bulk) and FEW blue connections (to other SC classes).

### Mechanism 3: Collapsed Edges (the new phenomenon at n=6)

At n=6, 5 edges collapse into self-loops when complement pairs merge. These are edges between a class C and its complement C^op. Before n=6: no class is adjacent to its complement. At n=6: some are.

A class C is adjacent to C^op iff there exists an arc flip that sends T to a tournament isomorphic to T^complement. This is a **NEAR-ANTI-AUTOMORPHISM**: a permutation sigma such that sigma(flip(T, e)) = T^complement.

Collapsed edges measure the "proximity" of a tournament to its complement -- how close T is to T^op in the flip graph. At small n (3,4,5): T and T^op are far apart. At n=6: some are adjacent.

## The Asymptotic Picture

As n -> infinity:
1. **Almost all merged nodes are NS-pairs** (SC fraction -> 0)
2. **Almost all edges are blue** (NS-NS connections dominate)
3. **The merged graph approaches a single-color (blue) graph** with a thin layer of black edges connecting the rare SC spine to the NS bulk
4. **Collapsed edges grow** but remain a small fraction of total edges

The merged graph at large n looks like:
```
[SC spine] ---(black)--- [NS-pair bulk]
    |                        |
  (blue,                   (blue,
   sparse)                  dense)
```

The blue SC-SC subgraph is sparse (few SC classes, few connections between them).
The blue NS-NS subgraph is dense (many NS-pairs, many connections).
The black SC-NS connections form a bipartite layer.

## The Recursive Formula Attempt

The MERGED vertex count at n: V_merged(n) = (A000568(n) + SC(n)) / 2.

Since A000568(n) ~ 2^{C(n,2)} / n!, and SC(n) grows much slower:
V_merged(n) ~ A000568(n) / 2 ~ 2^{C(n,2)-1} / n!

The MERGED edge count: E_merged(n) = E_total(n) - collapsed(n) - twin(n).
From the near-formula: E_total ~ V*m*(1-f)/2.
So: E_merged ~ V*m*(1-f)/2 - collapsed - twin.

The collapsed and twin counts are the NEW pieces. If we can predict them:
- collapsed(n) = # edges between complement pairs (0, 0, 0, 5, ?)
- twin(n) = # edges appearing in both C and C^op neighborhoods

Then: E_merged = E_total - collapsed - twin.

## What the Merged Graph MEANS

The merged graph G_n/Z_2 is the **irreducible tournament classification**. Each node represents a genuine structural type -- not distinguishable by any symmetry from its complement.

In the merged graph:
- Blue edges = "same-side" connections (orientation-preserving flips)
- Black edges = "cross-side" connections (orientation-reversing flips)

This is the structure of a **real projective space** -- the quotient of the sphere (G_n on S^2 at n=5) by the antipodal map (complement). RP^2 from S^2 by antipodal identification.

At n=5: G_5 lives on S^2 (genus 0). G_5/Z_2 would live on RP^2 (the real projective plane). RP^2 is non-orientable (has only one side, like a Mobius strip). The "one-sidedness" comes from the complement identification: there's no global way to distinguish a tournament from its complement.

At n=6: G_6 may live on a higher-genus surface. G_6/Z_2 would be the quotient of that surface by the complement involution.

## The n -> n-2 Descent Revisited

In the merged graph:
- G_3/Z_2 = K_2 (2 nodes, 1 edge)
- G_4/Z_2 = K_3 (3 nodes, 3 edges = complete)
- G_5/Z_2 = 10 nodes, 21 edges (density 0.47)
- G_6/Z_2 = 34 nodes, 143 edges (density 0.25)

The density is DECREASING: 1.0, 1.0, 0.47, 0.25. At large n it will approach 0 (sparse graph).

The descent G_n/Z_2 -> G_{n-2}/Z_2 would map:
- G_4/Z_2 = K_3 -> G_2/Z_2 = K_1 (trivial)
- G_5/Z_2 (10 nodes) -> G_3/Z_2 = K_2 (2 nodes)
- G_6/Z_2 (34 nodes) -> G_4/Z_2 = K_3 (3 nodes)

The descent reduces 10 -> 2 and 34 -> 3. The ratio is roughly V_merged(n) / V_merged(n-2) ~ A000568(n) / A000568(n-2). At n=5->3: 12/2 = 6. At n=6->4: 56/4 = 14. Growing.

The PoS classes in G_n/Z_2 should map surjectively onto G_{n-2}/Z_2.

## The Deep Principle

The merged meta-graph G_n/Z_2 is the **moduli space of tournament types** -- the space of all genuinely distinct pairwise comparison structures on n items, up to relabeling and complementation.

This moduli space has:
- A **spine** of self-complementary types (rare, slow-growing)
- A **bulk** of non-SC pairs (dominant, fast-growing)
- **Blue connections** (structure-preserving deformations)
- **Black connections** (structure-reversing deformations)
- **Collapsed points** (where a type touches its complement)

The evolution with n is:
1. The moduli space EXPANDS (more types)
2. The spine THINS (SC fraction -> 0)
3. The connections become BLUER (NS-NS dominate)
4. The density DECREASES (graph becomes sparser)
5. The H-gradient is strong but NOT a perfect DAG — level edges and H-decreasing edges appear at larger n (see MISTAKE-035)

This is the tournament analog of the moduli space of algebraic curves: as the genus increases, the moduli space grows, becomes more complex, but retains its fundamental structure (DAG from low-H to high-H, spine of symmetric types, bulk of generic types).
