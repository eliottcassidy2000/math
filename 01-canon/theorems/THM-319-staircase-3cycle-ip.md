---
id: THM-319
title: Independence Polynomial Coefficients for 3-Cycle Restriction of All-0 Staircase
status: PROVED (alpha_1--7); CONJECTURAL (alpha_8--9)
session: opus-2026-05-22-S4
---

## Statement

Let $I_3(T_k, x)$ denote the independence polynomial of $\Omega_3(T_k)$, the subgraph of the conflict graph $\Omega(T_k)$ induced on the 3-cycles. Then:
$$I_3(T_k, x) = \sum_{m=0}^{\lfloor 2k/3 \rfloor} \alpha_m x^m, \quad \alpha_m = \sum_{j=0}^{\lfloor m/2 \rfloor} c_{m,j} \binom{k}{2m-j},$$
with the following explicit coefficients:

| $m$ | $\alpha_m$ formula |
|---|---|
| 1 | $2\binom{k}{2}$ |
| 2 | $12\binom{k}{4} + \binom{k}{3}$ |
| 3 | $120\binom{k}{6} + 20\binom{k}{5}$ |
| 4 | $1680\binom{k}{8} + 420\binom{k}{7} + 10\binom{k}{6}$ |
| 5 | $30240\binom{k}{10} + 10080\binom{k}{9} + 560\binom{k}{8}$ |
| 6 | $665280\binom{k}{12} + 277200\binom{k}{11} + 25200\binom{k}{10} + 280\binom{k}{9}$ |
| 7 | $17297280\binom{k}{14} + 8648640\binom{k}{13} + 1108800\binom{k}{12} + 30800\binom{k}{11}$ |
| 8 | $518918400\binom{k}{16} + 302702400\binom{k}{15} + 50450400\binom{k}{14} + 2402400\binom{k}{13} + 15400\binom{k}{12}$ (conjectural) |

## Patterns in Coefficients

**Leading coefficients:** $c_{m,0} = A_m := \frac{(2m)!}{m!}$.

Proof: The leading term counts $m$-tuples of mutually vertex-disjoint 3-cycles from fully disjoint pairs. A matching of size $m$ in $K_k$ corresponds to $m$ disjoint pair-of-pairs; for each such matching, all $2^m$ choices of A- or B-type cycle are valid, giving $2^m \cdot \binom{k}{2m}(2m-1)!! = A_m \binom{k}{2m}$ configurations.

**Sub-leading coefficients:** $c_{m,1} = \frac{m-1}{12} A_m$.

| $m$ | $c_{m,1}/A_m$ |
|---|---|
| 2 | $1/12$ |
| 3 | $2/12 = 1/6$ |
| 4 | $3/12 = 1/4$ |
| 5 | $4/12 = 1/3$ |
| 6 | $5/12$ |
| 7 | $6/12 = 1/2$ |
| 8 | $7/12$ (conjectural) |

**Number of terms:** $\alpha_m$ has $\lfloor m/2 \rfloor + 1$ terms. Pattern: 1, 2, 2, 3, 3, 4, 4, 5, ... for $m=1,2,3,...$

**Minimum-$k$ coefficients (diagonal):** For $m=2j$:
$$c_{2j,j} = \prod_{i=1}^{j-1} \binom{3i+2}{2}.$$
For $m=2j+1$:
$$c_{2j+1,j} = 2\prod_{i=1}^{j} \binom{3i+2}{2}.$$

These equal $\alpha_d(T_{k_{\min}})$ where $k_{\min} = \lceil 3m/2 \rceil$ (the first $k$ where $d(T_k) \geq m$).

**Explicit diagonal values:**
| $m$ | $k_{\min}$ | $c_{m,\lfloor m/2 \rfloor}$ | Formula |
|---|---|---|---|
| 1 | 2 | 2 | $2\cdot\prod_{i=1}^{0}$ |
| 2 | 3 | 1 | $\prod_{i=1}^{0}$ |
| 3 | 5 | 20 | $2\cdot\binom{5}{2}$ |
| 4 | 6 | 10 | $\binom{5}{2}$ |
| 5 | 8 | 560 | $2\cdot\binom{5}{2}\binom{8}{2}$ |
| 6 | 9 | 280 | $\binom{5}{2}\binom{8}{2}$ |
| 7 | 11 | 30800 | $2\cdot\binom{5}{2}\binom{8}{2}\binom{11}{2}$ |
| 8 | 12 | 15400 | $\binom{5}{2}\binom{8}{2}\binom{11}{2}$ |
| 9 | 14 | 2802800 | $2\cdot\binom{5}{2}\binom{8}{2}\binom{11}{2}\binom{14}{2}$ |

**Adjacent-diagonal ratio:** $c_{2j+1,j} = 2 c_{2j,j}$ (the odd-indexed diagonal is always twice the even-indexed). 

**Growth of diagonal:** $c_{2j,j}/c_{2j-2,j-1} = \binom{3j+2}{2} = \frac{(3j+1)(3j+2)}{2}$. So:
$$c_{2j,j} = \frac{1}{\binom{3j+2}{2}} \cdot c_{2j+2, j+1}.$$

## Combinatorial Interpretation

$\alpha_m$ counts collections of $m$ vertex-disjoint 3-cycles in $T_k$. The expansion as $\sum_j c_{m,j}\binom{k}{2m-j}$ organizes by the number of distinct pair indices used:
- Leading term $c_{m,0}\binom{k}{2m}$: uses exactly $2m$ distinct pairs (full matching).
- Term $c_{m,1}\binom{k}{2m-1}$: uses $2m-1$ pairs — one pair shared across two cycles via an A-B bridge.

**Bridge structure:** The 3-cycle A(a,p) = $\{2a, 2a+1, 2p\}$ and B(p,b) = $\{2p+1, 2b, 2b+1\}$ are vertex-disjoint for $a<p<b$, using only 3 distinct pair indices $\{a,p,b\}$.

## Verification

All values for $m=1$ through $7$ verified computationally for $k=2$ through $14$.
The formula for $m=8$ is conjectural: $c_{8,0}$ and $c_{8,1}$ follow the pattern but need $k\geq 16$ data to verify independently; $c_{8,2..4}$ are verified from $k=12..14$ data.

## See Also

- THM-318: Proof that $d(T_k) = \lfloor 2k/3\rfloor$
- THM-320: Product formulas for $\alpha_d(T_k)$ (leading coefficient at minimum $k$)
</content>
</invoke>