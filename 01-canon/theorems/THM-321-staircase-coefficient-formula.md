---
id: THM-321
title: Complete Closed-Form Formula for 3-Cycle IP Coefficients of All-0 Staircase
status: PROVED
session: opus-2026-05-22-S4
---

## Statement

**Theorem:** The coefficient $c_{m,j}$ in $\alpha_m = \sum_{j=0}^{\lfloor m/2\rfloor} c_{m,j}\binom{k}{2m-j}$ satisfies:
$$c_{m,j} = \frac{A_{m-2j}}{j!}\prod_{i=0}^{j-1}\binom{2m-j-3i}{3},$$
where $A_n = (2n)!/n!$ and empty products equal 1.

In particular:
$$c_{m,0} = A_m = \frac{(2m)!}{m!}, \qquad \frac{c_{m,1}}{A_m} = \frac{m-1}{12}, \qquad \frac{c_{m,2}}{A_m} = \frac{(m-1)(m-2)(m-3)}{144(2m-1)}.$$

## Proof

**Combinatorial interpretation:** $c_{m,j}\binom{k}{2m-j}$ counts collections of $m$ vertex-disjoint 3-cycles in $T_k$ that together use exactly $2m-j$ distinct pair indices.

**Cycle structure** (proved in THM-318): Every 3-cycle of $T_k$ involves exactly 2 pairs. For $a<b$:
- $A(a,b)$: vertices $\{2a, 2a+1, 2b\}$
- $B(a,b)$: vertices $\{2a+1, 2b, 2b+1\}$

**Conflict rules** (vertex-sharing analysis):
- $A(a,b)$ and $A(c,d)$ conflict $\iff \{a,b\}\cap\{c,d\}\neq\emptyset$.
- $B(a,b)$ and $B(c,d)$ conflict $\iff \{a,b\}\cap\{c,d\}\neq\emptyset$.
- $A(a,b)$ and $B(c,d)$ conflict $\iff a=c$ or $a=d$ or $b=d$. (Notably: $b=c$ does NOT cause conflict.)
- $A(a,b)$ and $B(a,b)$ always conflict (same pair).

**Key observation — Bridge:** For $a < p < b$, cycles $A(a,p)$ and $B(p,b)$ are vertex-disjoint. They share pair index $p$ (as second pair of $A$, first pair of $B$) but use distinct vertices: $A(a,p) = \{2a, 2a+1, 2p\}$, $B(p,b) = \{2p+1, 2b, 2b+1\}$.

This is the ONLY way two cycles can share a pair index without conflicting. No other A-B or same-type sharing is vertex-disjoint.

**Counting configurations using exactly $2m-j$ pairs:**

For each $j$: the $j$ "savings" (pairs used twice) come from exactly $j$ independent bridges. Each bridge is an $\{A(a_i, p_i), B(p_i, b_i)\}$ pair sharing one pair index $p_i$.

For a fixed set $S$ of $2m-j$ pair indices:
1. Choose $j$ disjoint bridge triples $(a_1,p_1,b_1),\ldots,(a_j,p_j,b_j)$ from $S$ (no shared pair among the $j$ triples). Number of ways (unordered) = $\frac{1}{j!}\prod_{i=0}^{j-1}\binom{|S|-3i}{3}$.
2. The remaining $|S|-3j = 2m-j-3j = 2m-4j$ pair indices form a matching of size $m-2j$ (one cycle per matched pair-of-pairs). Number of matchings = $(2m-4j-1)!!$ [double factorial].
3. For each matched pair-of-pairs $\{c,d\}$: choose $A(c,d)$ or $B(c,d)$ (2 choices). Total: $2^{m-2j}$.

Total configurations for fixed $S$ with $|S|=2m-j$:
$$\frac{1}{j!}\prod_{i=0}^{j-1}\binom{2m-j-3i}{3} \cdot (2m-4j-1)!! \cdot 2^{m-2j} = \frac{A_{m-2j}}{j!}\prod_{i=0}^{j-1}\binom{2m-j-3i}{3}$$
using $A_{m-2j} = 2^{m-2j}\cdot(2m-4j-1)!!$.

The total count over all $\binom{k}{2m-j}$ choices of $S$ gives $\alpha_m$'s contribution $c_{m,j}\binom{k}{2m-j}$, hence:
$$c_{m,j} = \frac{A_{m-2j}}{j!}\prod_{i=0}^{j-1}\binom{2m-j-3i}{3}. \quad\square$$

## Corollaries

**Leading coefficient:** $c_{m,0} = A_m = (2m)!/m!$ (empty product, no bridges).

**Sub-leading:** $c_{m,1} = \binom{2m-1}{3}A_{m-2}$, and $c_{m,1}/A_m = (m-1)/12$.

*Proof:* $c_{m,1}/A_m = \binom{2m-1}{3}A_{m-2}/A_m = \frac{(2m-1)(2m-2)(2m-3)}{6}\cdot\frac{m(m-1)}{2m(2m-1)(2m-2)(2m-3)} = \frac{m-1}{12}$.

**Third coefficient:** $c_{m,2} = \frac{1}{2}\binom{2m-2}{3}\binom{2m-5}{3}A_{m-4}$, and $c_{m,2}/A_m = \frac{(m-1)(m-2)(m-3)}{144(2m-1)}$.

**Diagonal (minimum $k$):** For $j=\lfloor m/2\rfloor$:
$$c_{2j,j} = \prod_{n=1}^{j}\binom{3n-1}{2}, \qquad c_{2j+1,j} = 2\prod_{n=1}^{j}\binom{3n+2}{2}.$$

(These recover THM-320.)

## Verification Table

| $m$ | $j$ | Formula | Computed |
|---|---|---|---|
| 2 | 0 | $A_2=12$ | 12 ✓ |
| 2 | 1 | $\binom{3}{3}A_0/1! = 1$ | 1 ✓ |
| 3 | 0 | $A_3=120$ | 120 ✓ |
| 3 | 1 | $\binom{5}{3}A_1/1! = 10\cdot2=20$ | 20 ✓ |
| 4 | 0 | $A_4=1680$ | 1680 ✓ |
| 4 | 1 | $\binom{7}{3}A_2/1!=35\cdot12=420$ | 420 ✓ |
| 4 | 2 | $\binom{6}{3}\binom{3}{3}A_0/2!=20\cdot1/2=10$ | 10 ✓ |
| 5 | 2 | $\binom{8}{3}\binom{5}{3}A_1/2!=56\cdot10\cdot2/2=560$ | 560 ✓ |
| 6 | 3 | $\binom{9}{3}\binom{6}{3}\binom{3}{3}A_0/6=84\cdot20\cdot1/6=280$ | 280 ✓ |
| 7 | 3 | $\binom{11}{3}\binom{8}{3}\binom{5}{3}A_1/6=165\cdot56\cdot10\cdot2/6=30800$ | 30800 ✓ |
| 8 | 4 | $\binom{12}{3}\binom{9}{3}\binom{6}{3}\binom{3}{3}A_0/24=220\cdot84\cdot20\cdot1/24=15400$ | 15400 ✓ |

All entries for $k=2$ through $14$ verified computationally.

## See Also

- THM-318: $d(T_k)=\lfloor 2k/3\rfloor$
- THM-319: Explicit coefficient table
- THM-320: Product formula for $\alpha_d(T_k)$
