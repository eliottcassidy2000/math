# The Frustration Propagation Formula

*opus-2026-04-04-S19*

## The Formula

**Theorem.** The total frustration of the multilinear polynomial H(t) — defined as the sum of absolute values of all negative quadratic coefficients — satisfies:

$$F(n) = 2^n - 4(n-1)$$

for all $n \geq 3$. The recursion is:

$$F(n) = F(n-1) + 2(2^{n-2} - 2)$$

with $F(3) = 0$.

## What This Means

The frustration grows EXPONENTIALLY ($2^n$ dominates) with a LINEAR correction ($-4n+4$). At $n = 7$: $F = 128 - 24 = 104$. At $n = 10$: $F = 1024 - 36 = 988$.

## The Layer Decomposition

When adding vertex $n$ (going from $n-1$ to $n$), the new frustration has two equal parts:

$$\Delta F(n) = F_{\text{nn}}(n) + F_{\text{no}}(n) = 2(2^{n-2} - 2)$$

where:
- $F_{\text{nn}} = 2^{n-2} - 2$: frustration between the NEW tiles (all share upper vertex $n$)
- $F_{\text{no}} = 2^{n-2} - 2$: frustration between new tiles and OLD tiles (shared lower vertices)

**The exact equality $F_{\text{nn}} = F_{\text{no}}$ is the key symmetry.** It says: the amount of frustration the new vertex creates WITHIN its own tiles equals the frustration it creates BETWEEN its tiles and the existing structure.

## Proof of $F_{\text{nn}} = 2^{n-2} - 2$

The new tiles at step $n$ are $(n, 1), (n, 2), \ldots, (n, n-2)$, with skips $n-1, n-2, \ldots, 2$.

All pairs share upper vertex $n$, so all are same-end with $c_{ij} = -2^{\max(1, |s_i - s_j| - 1)}$.

$$F_{\text{nn}} = \sum_{1 \leq y_1 < y_2 \leq n-2} 2^{\max(1, |s_1 - s_2| - 1)}$$

where $s_i = n - y_i$. Since $|s_1 - s_2| = |y_2 - y_1|$:

$$F_{\text{nn}} = \sum_{d=1}^{n-3} (n-2-d) \cdot 2^{\max(1, d-1)}$$

For $d = 1$: coefficient $= (n-3)$, weight $= 2^0 = 1$. But $\max(1, 0) = 1$, so weight $= 2$.

Actually: for $d = 1$: $\max(1, 0) = 1$, weight = $2^1 = 2$.
For $d \geq 2$: weight = $2^{d-1}$.

$$F_{\text{nn}} = 2(n-3) + \sum_{d=2}^{n-3} (n-2-d) \cdot 2^{d-1}$$

This telescopes to $2^{n-2} - 2$ (verified computationally for $n = 4, \ldots, 8$).

## The Frustration Threshold by Tile Type

Each tile has a frustration-to-energy ratio:

$$\rho_k = \frac{\text{frustration load of tile } k}{c_k} = \frac{\sum_{j: c_{kj} < 0} |c_{kj}|}{2^{\text{skip}(k)-1}}$$

| Tile type | Skip | $\rho$ at $n=7$ | Pattern |
|-----------|------|-----------------|---------|
| Apex $(n,1)$ | $n-1$ | 1.00 | Always exactly 1 |
| Near-apex $(n,2)$ or $(n-1,1)$ | $n-2$ | 1.12 | Slightly above 1 |
| Middle tiles | $\lfloor n/2 \rfloor$ | 2-3 | Moderate |
| Corner tiles $(3,1)$ or $(n,n-2)$ | 2 | $2^{n-3}$ | Exponential! |

**The apex is at the frustration threshold ($\rho = 1$).** Below threshold, the tile's gradient is always positive (flipping it always increases H from transitive). At threshold, the gradient can flip with enough backward neighbors. Above threshold, the gradient is negative for most configurations.

**Corner tiles are exponentially frustrated** ($\rho = 2^{n-3}$). They are the MOST constrained tiles — their tiny energy ($c_k = 2$) is overwhelmed by frustration from many same-end neighbors.

## Connection to the Three Functors

In the three-functor picture:
- Frustration lives in **Functor 2** (tournament → conflict graph). It measures how many same-end cycle-competition edges exist.
- The recursion $F(n) = F(n-1) + \Delta F$ is a direct consequence of **THM-299** (recursive preservation): the old frustration is inherited exactly, and the new frustration comes from the new vertex's tiles.
- The exponential growth $2^n$ in $F(n)$ reflects the exponential growth of the tile count $m = C(n-1,2) \sim n^2/2$, while the linear correction $-4n$ reflects the boundary effect of the staircase.

## Prediction

At $n = 8$: $F(8) = 256 - 28 = 228$. At $n = 10$: $F(10) = 1024 - 36 = 988$. The frustration-to-tiles ratio $F/m$ grows as $F/m \approx 2^n / (n^2/2) = 2^{n+1}/n^2$, which increases exponentially — the system becomes MORE frustrated per tile as $n$ grows.
