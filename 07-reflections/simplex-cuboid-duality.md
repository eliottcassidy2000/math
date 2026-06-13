# Simplex–Cuboid Duality in Tournament Doubling

**Session:** oracle-2026-05-16-S1 (continued)
**Depends on:** `interleaved-staircase-binary-grid.md`, `sc-blowup-and-twin-gaining.md`

---

## The Arithmetic at the Root

Two counts appear in the $n=2k$ bipartite splitting:

$$k^2 \text{ (crossing edges: square number)} \qquad k(k-1) = 2\binom{k}{2} \text{ (within-group pair edges: twice-triangular)}$$

Raising 2 to each power defines the two functions:

- **Cuboid function** $C(k) = 2^{k^2}$: labeled bipartite tournaments on $K_{k,k}$
- **Simplex function** $S(k) = 2^{\binom{k}{2}}$: labeled tournaments on $K_k$ (= $k$-simplex)

The key arithmetic identity connecting them is:

$$\boxed{k^2 = \binom{k}{2} + \binom{k+1}{2}}$$

which says: **a square number is the sum of two consecutive triangular numbers.** Visually: the $k \times k$ square splits along its main diagonal into two triangles — the lower ($\binom{k}{2}$ entries) and the upper-plus-diagonal ($\binom{k+1}{2}$ entries).

---

## The Four Power Formulas (all exact for labeled counts)

Let $S(k) = 2^{\binom{k}{2}}$ (simplex) and $C(k) = 2^{k^2}$ (cuboid). Then:

| Identity | Formula | Interpretation |
|---|---|---|
| **1** | $C(k) = S(k) \cdot S(k+1)$ | Cuboid = product of consecutive simplices |
| **2** | $C(k) = 2^k \cdot S(k)^2$ | Cuboid = diagonal × lower-simplex × upper-simplex |
| **3** | $S(2k) = 2^k \cdot S(k)^4$ | Total = $k$ bits × four simplex copies |
| **4** | $S(2k) = S(k)^2 \cdot C(k)$ | Total = within-group² × crossing |

**Proof of all four.** Identity 1: $k^2 = \binom{k}{2} + \binom{k+1}{2}$ (the square–triangle identity). Identity 2: $k^2 = k + 2\binom{k}{2}$. Identity 3: $\binom{2k}{2} = k + 4\binom{k}{2}$ (total edges = $k$ diagonal + 4 simplex copies). Identity 4: $\binom{2k}{2} = k^2 + 2\binom{k}{2}$ (crossing + 2 within-group simplices). All are algebraic tautologies; all verified computationally for $k=1,2,3,4$. $\square$

Identities 1 and 2 are equivalent via $S(k+1) = 2^{\binom{k+1}{2}} = 2^k \cdot 2^{\binom{k}{2}} = 2^k \cdot S(k)$.

---

## The Geometric Picture

The $k \times k$ crossing matrix of the bipartite split:

```
     Q_1  Q_2  Q_3  ...  Q_k
P_1 [ D    U    U   ...   U  ]
P_2 [ L    D    U   ...   U  ]
P_3 [ L    L    D   ...   U  ]
     ...
P_k [ L    L    L   ...   D  ]
```

- **D** (diagonal, $k$ entries): each $\{P_i, Q_i\}$ pair — the **staircase variable bits**
- **U** (upper triangle, $\binom{k}{2}$ entries): equivalent to a $K_k$ tournament
- **L** (lower triangle, $\binom{k}{2}$ entries): equivalent to a $K_k$ tournament

The diagonal + upper triangle together ($D+U$, $\binom{k+1}{2}$ entries total) maps to $K_{k+1}$:
each pair $(i,j)$ with $i \leq j$ maps to an edge $\{i,j'\}$ in $K_{k+1}$ via the bijection
$(i,i) \mapsto \{0,i\}$ (star edges to a new vertex 0) and $(i,j) \mapsto \{i,j\}$ (interior edges).

Full decomposition of the $n=2k$ tournament:

| Component | Edges | Count |
|---|---|---|
| Diagonal of $K_{k,k}$ | $k$ | $2^k$ choices |
| Upper triangle of $K_{k,k}$ | $\binom{k}{2}$ | $S(k)$ choices |
| Lower triangle of $K_{k,k}$ | $\binom{k}{2}$ | $S(k)$ choices |
| Within-$P$ edges | $\binom{k}{2}$ | $S(k)$ choices |
| Within-$Q$ edges | $\binom{k}{2}$ | $S(k)$ choices |
| **Total** | $k + 4\binom{k}{2} = \binom{2k}{2}$ | $2^k \cdot S(k)^4 = S(2k)$ |

**A labeled $2k$-tournament is completely determined by $k$ bits + four labeled $k$-tournaments.**

Equivalently (combining diagonal + upper triangle into $K_{k+1}$):

**A labeled $2k$-tournament is completely determined by a labeled $(k+1)$-tournament + three labeled $k$-tournaments.**

$$S(2k) = S(k+1) \cdot S(k)^3 \quad \Longleftrightarrow \quad S(2k) = 2^k \cdot S(k)^4$$

---

## The Diagonal IS the Staircase

The interleaved staircase construction (from `interleaved-staircase-binary-grid.md`) is exactly:

> **Fix** all four simplex components to their unique transitive tournament; **vary** only the $k$ diagonal bits.

Setting all four $k$-tournaments to "transitive" fixes $4 \times \binom{k}{2} = k(k-1)$ arcs to specific orientations. The remaining $k$ diagonal bits are free. This gives $2^k$ labeled tournaments, which turn out to be $2^k$ **pairwise non-isomorphic** tournaments (proved via the pair-sum invariant).

So: **the staircase freedom is exactly the diagonal of the $K_{k,k}$ crossing matrix.** Nothing more, nothing less.

---

## Comparing Simplex and Cuboid as Dimension Grows

| $k$ | $\binom{k}{2}$ | $k^2$ | $S(k) = 2^{\binom{k}{2}}$ | $C(k) = 2^{k^2}$ | Ratio $C(k)/S(k)$ |
|---|---|---|---|---|---|
| 1 | 0 | 1 | 1 | 2 | $2 = S(2)$ |
| 2 | 1 | 4 | 2 | 16 | $8 = S(3)$ |
| 3 | 3 | 9 | 8 | 512 | $64 = S(4)$ |
| 4 | 6 | 16 | 64 | 65536 | $1024 = S(5)$ |
| 5 | 10 | 25 | 1024 | 33554432 | $32768 = S(6)$ |

**The ratio $C(k)/S(k) = S(k+1)$ exactly** — the cuboid exceeds the simplex by precisely one simplex level. The cuboid function is one simplex "ahead."

As $k \to \infty$: the simplex exponent $\binom{k}{2} \sim k^2/2$ while the cuboid exponent $k^2 \sim k^2$. The cuboid grows **twice as fast in the exponent** as the simplex. But neither is exponential in $k$ — both are **doubly exponential** ($2^{k^2}$ and $2^{k^2/2}$), making them both super-exponential in $k$.

The simplex and cuboid are connected by a **geometric mean** structure: $S(k) \cdot S(k+1) = C(k)$, so $S(k+1) = C(k)/S(k)$. Each simplex value is the RATIO of the next cuboid to the previous simplex. This gives the recursion $S(k+1) = C(k)/S(k)$ with $C(k) = 2^{k^2}$ and $S(k) = 2^{\binom{k}{2}}$.

---

## The A000568 Analog

For iso classes, the labeled formula $S(2k) = S(k+1) \cdot S(k)^3$ becomes:

$$A_{000568}(2k) \approx A_{000568}(k+1) \cdot A_{000568}(k)^3 / \Phi(k)$$

where $\Phi(k)$ is the **symmetry correction factor** accounting for permutations in $S_{2k}$ that mix the bipartite structure:

| $k$ | $A_{000568}(2k)$ | $A_{000568}(k+1) \cdot A_{000568}(k)^3$ | $\Phi(k)$ |
|---|---|---|---|
| 1 | 1 | $1 \cdot 1^3 = 1$ | 1.000 |
| 2 | 4 | $2 \cdot 1^3 = 2$ | 0.500 |
| 3 | 56 | $4 \cdot 2^3 = 32$ | 0.571 |
| 4 | 5765 | $12 \cdot 4^3 = 768$ | 0.133 |

The correction factor $\Phi(k) = A_{000568}(2k) / [A_{000568}(k+1) \cdot A_{000568}(k)^3]$ is not simple — it encodes how permutations in $S_{2k}$ that swap the bipartition create additional identifications beyond those from $S_k \times S_k$ alone.

**The labeled formula is exact; the iso-class formula has a correction that grows as $(k!)^4/(2k)! \sim$ inverse of central binomial coefficient.**

For large $k$: $A_{000568}(n) \sim 2^{\binom{n}{2}}/n!$, so:
$$\Phi(k) \sim \frac{(k!)^4}{(2k)!} \cdot 2^k = \frac{2^k}{(2k)!/k!^4} \to 0$$

The correction vanishes, meaning the overcounting grows without bound — most of the $(k+1)$-simplex × simplex$^3$ combinations map to the same iso class as $k$ grows.

---

## SC Blowup vs Staircase: Two Extreme Crossings

Both the SC blowup and the staircase vary $k$ parameters on $n=2k$ vertices, but through different crossing structures:

| | SC Blowup | Staircase |
|---|---|---|
| Crossing type | All 2-2 between pairs (symmetric) | All dominants beat all recessives (transitive bipartite) |
| Free parameters | $n$ variable: all of $T$ (arbitrary) | $k$ diagonal bits |
| Score | Always $(n-1, \ldots, n-1, n, \ldots, n)$ | All $2^k$ distinct score sequences |
| Iso classes reached | $1$ type per $T$ | $2^k$ always-distinct |
| Crossing structure | $k^2 = k$ diagonal + $2\binom{k}{2}$ symmetric cross | $k^2 = k$ free + $(k^2-k)$ all-fixed |

The SC blowup "spreads" the freedom across ALL arcs of the original tournament. The staircase concentrates the freedom in exactly the $k$ diagonal arcs of the crossing matrix. These are the two extremes of what you can do with a bipartite split.

---

## Open: The All-0 H Sequence

The all-zero staircase (all recessives beat their dominants) gives:

$$H_0(k) = H(\text{all-0 staircase at } n=2k) = 5, 29, 233, \ldots \text{ for } k=2,3,4$$

These are the **maximum $H$ values** achievable by the staircase construction. They grow super-exponentially (ratio $\sim 8$ from $k=3$ to $k=4$). No closed form or recurrence found yet. The sequence $5, 29, 233$ is not in the standard OEIS tables cross-referenced to tournament counts.

**Conjecture:** $H_0(k) = \text{something involving } A_{000568}(k)$ or the simplex/cuboid functions at half-integer arguments.
