---
id: THM-317
title: Endpoint Theorem for All-0 Staircase (HPs ending at last dominant = H(k-1))
status: PROVED
session: opus-2026-05-22-S3
---

## Statement

**Theorem:** For the all-0 interleaved staircase $T_k$ on $n=2k$ vertices:
$$\text{ep\_end}(2k-2, T_k) = H(k-1).$$

The number of Hamiltonian paths of $T_k$ ending at vertex $2k-2$ (the **last dominant**) equals the total HP count of $T_{k-1}$.

## Proof

In $T_k$: vertex $2k-1$ (last recessive) beats ONLY vertex $2k-2$ (within pair). So vertex $2k-1$ has out-degree 1.

Any HP of $T_k$ ending at $2k-2$: since $2k-1$ must appear somewhere in the path and its only outgoing arc is $2k-1 \to 2k-2$, vertex $2k-1$ must be the **second-to-last** vertex (immediately preceding $2k-2$ at the end). (If $2k-1$ appeared elsewhere, its successor would need to be $2k-2$, but $2k-2$ is the last vertex — so $2k-1$ must be at position $2k-1$.)

Therefore: HPs ending at $2k-2$ have the form
$$v_1, v_2, \ldots, v_{2k-2}, \;2k-1,\; 2k-2$$
where $(v_1,\ldots,v_{2k-2})$ is a Hamiltonian path of $T_{k-1}$ (on $\{0,1,\ldots,2k-3\}$).

This path can be extended by appending $2k-1$ (valid since ALL vertices of $T_{k-1}$ beat $2k-1$ in $T_k$: $2k-1$ has rank $2k-1$ = worst) and then $2k-2$ (valid since $2k-1 \to 2k-2$ within pair).

Bijection: HP of $T_{k-1}$ $\leftrightarrow$ HP of $T_k$ ending at $2k-2$. Count = $H(k-1)$. □

## Combined with Anti-Palindrome (THM-316)

Since $\text{ep\_start}(1, T_k) = \text{ep\_end}(2k-2, T_k)$ by anti-palindrome:
$$\text{ep\_start}(1, T_k) = H(k-1).$$

The number of HPs starting at vertex $1$ (first recessive) = $H(k-1)$.

Similarly, $\text{ep\_start}(2k-2) = \text{ep\_end}(1)$ and $\text{ep\_start}(2k-1) = \text{ep\_end}(0)$. Combining with the data: $\text{ep\_end}(0) = \text{ep\_end}(1) = a_{k-1}$ (from the endpoint table).

## Verified Data

| $k$ | ep\_end$(2k-2)$ | $H(k-1)$ | Match? |
|-----|------|------|------|
| 2 | $\text{ep\_end}(2) = 1$ | $H(1)=1$ | ✓ |
| 3 | $\text{ep\_end}(4) = 5$ | $H(2)=5$ | ✓ |
| 4 | $\text{ep\_end}(6) = 29$ | $H(3)=29$ | ✓ |
| 5 | $\text{ep\_end}(8) = 233$ | $H(4)=233$ | ✓ |
| 6 | $\text{ep\_end}(10) = 2489$ | $H(5)=2489$ | ✓ |
| 7 | $\text{ep\_end}(12) = 33773$ | $H(6)=33773$ | ✓ |

## See Also

- THM-316: Anti-palindrome theorem (full endpoint symmetry)
- All-0 staircase: H values in `05-knowledge/results/staircase_extended.out`
