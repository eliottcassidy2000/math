---
id: THM-318
title: Independence Polynomial Degree for All-0 Staircase
status: PROVED
session: opus-2026-05-22-S4
---

## Statement

**Theorem:** For the all-0 interleaved staircase $T_k$ on $n=2k$ vertices:
$$d(T_k) := \deg I(\Omega(T_k), x) = \left\lfloor \frac{2k}{3} \right\rfloor.$$

## Proof

**Upper bound:** Every independent set $S$ in $\Omega(T_k)$ is a collection of vertex-disjoint odd cycles in $T_k$. Each odd cycle uses $\geq 3$ vertices, and $T_k$ has $2k$ vertices. So $|S| \leq \lfloor 2k/3 \rfloor$.

**Lower bound:** Construct $\lfloor 2k/3\rfloor$ mutually vertex-disjoint 3-cycles. Group pairs into consecutive triples $(3j, 3j+1, 3j+2)$ for $j=0, 1, \ldots$:
- For each triple, extract 2 vertex-disjoint 3-cycles:
  - Cycle A from pairs $(3j, 3j+1)$: $(2\cdot 3j+1) \to (2\cdot 3j) \to (2\cdot(3j+1)) \to (2\cdot 3j+1)$
  - Cycle B from pairs $(3j+1, 3j+2)$: $(2\cdot(3j+2)+1) \to (2\cdot(3j+2)) \to (2\cdot(3j+1)+1) \to (2\cdot(3j+2)+1)$
- Triple groups provide $2\lfloor k/3\rfloor$ cycles from $2\lfloor k/3\rfloor$ triples.
- Residual pair(s): if $k \equiv 2 \pmod{3}$, one extra pair provides a length-1 partial triple giving 1 additional cycle via A(last-2, last-1) or similar.

Detailed count:
- $k = 3j$: $j$ triples, $2j = 2k/3$ cycles.
- $k = 3j+1$: $j$ triples + 1 leftover pair. Left pair contributes 0 extra cycles. Total $2j = \lfloor 2k/3\rfloor$.
- $k = 3j+2$: $j$ triples + 2 leftover pairs. The 2 leftover pairs give 1 more cycle. Total $2j+1 = \lfloor 2k/3\rfloor$.

In all cases, $\lfloor 2k/3\rfloor$ mutually vertex-disjoint 3-cycles exist, so $\alpha(\Omega(T_k)) \geq \lfloor 2k/3\rfloor$. □

## Structure of 3-Cycles

**Key fact (proved):** ALL directed 3-cycles of $T_k$ arise from exactly 2 pairs. For pairs $a < b$ (0-indexed), exactly two 3-cycles exist:
- **A(a,b):** $(2a+1) \to (2a) \to (2b) \to (2a+1)$, vertices $\{2a, 2a+1, 2b\}$
- **B(a,b):** $(2b+1) \to (2b) \to (2a+1) \to (2b+1)$, vertices $\{2a+1, 2b, 2b+1\}$

No 3-cycle exists using vertices from 3 or more distinct pairs.

Total 3-cycles: $2\binom{k}{2} = k(k-1)$.

## Verified Data

| $k$ | $d = \lfloor 2k/3 \rfloor$ | Verified |
|-----|------|------|
| 2 | 1 | ✓ |
| 3 | 2 | ✓ |
| 4 | 2 | ✓ |
| 5 | 3 | ✓ |
| 6 | 4 | ✓ |
| 7 | 4 | ✓ |
| 8 | 5 | ✓ |
| 9 | 6 | ✓ |
| 10 | 6 | ✓ |
| 11 | 7 | ✓ |
| 12 | 8 | ✓ |
| 13 | 8 | ✓ |
| 14 | 9 | ✓ |

## See Also

- THM-319: Explicit alpha_m formulas for the 3-cycle IP of $T_k$
- THM-316: Anti-palindrome theorem (endpoint symmetry)
