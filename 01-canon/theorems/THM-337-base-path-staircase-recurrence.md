---
theorem: THM-337
name: Base-Path Staircase Order-3 Linear Recurrence
status: PROVED
session: opus-2026-05-27-S7; Chan–Pak span-three follow-up 2026-11-15
verified: run-composition bijection; direct subset DP through order 16 and 816 window controls
---

# Base-path staircase: exact run-composition law

On vertices $0,\dots,n-1$, orient adjacent pairs backward and all other pairs forward. Write $A_n$ for the number of Hamiltonian paths, with $A_0=1$ for the empty composition. Then
$$
\sum_{n\ge0}A_nx^n=\frac{1-x^2}{1-x-x^2-x^3}.
$$

**Bijection.** Break a Hamiltonian path at ascents. Each descending run is a reversed interval of consecutive integers because the only descending arcs go down one. The disjoint intervals partition the vertex set. A join from the minimum of one interval to the maximum of a disjoint interval must ascend, so the receiving interval lies above the departing interval. Therefore the intervals appear in increasing order. For consecutive intervals of lengths $a,b$, the join jumps by $a+b-1$; it is forward precisely when $a,b$ are not both one. Conversely any composition with no adjacent singleton parts gives exactly this path. The parts of size at least two have series $B=x^2/(1-x)$, with optional singleton in each gap and at both ends, including an optional singleton when there are no such parts. The resulting series is $(1+x)/(1-B(1+x))$.

In particular $A_0,A_1,A_2=1,1,1$ and $A_n=A_{n-1}+A_{n-2}+A_{n-3}$ for $n\ge3$. The per-vertex growth is the real Tribonacci root $\rho=1.839286755\ldots$ of $t^3-t^2-t-1=0$.

## Even-order recurrence (the original conjecture)

Set $h_k=A_{2k}$. Taking the even part of the series gives
$$
\sum_{k\ge0}h_kz^k=\frac{(1-z)^2}{1-3z-z^2-z^3},
\qquad h_k=3h_{k-1}+h_{k-2}+h_{k-3}\quad(k\ge3).
$$
The initial values are $h_0=1,h_1=1,h_2=5$. The dominant characteristic root is $\rho^2=3.382975768\ldots$; the two nonreal roots are about $-0.191488\pm0.508852i$ and have modulus $1/\rho=0.543689\ldots$. There is a positive constant $C$ with $h_k=C\rho^{2k}+O(\rho^{-k})$.

| $k$ | $2k$ | $h_k$ |
|---:|---:|---:|
| 1 | 2 | 1 |
| 2 | 4 | 5 |
| 3 | 6 | 17 |
| 4 | 8 | 57 |
| 5 | 10 | 193 |
| 6 | 12 | 653 |
| 7 | 14 | 2,209 |
| 8 | 16 | 7,473 |
| 9 | 18 | 25,281 |
| 10 | 20 | 85,525 |
| 11 | 22 | 289,329 |
| 50 | 100 | 127,071,617,887,002,752,149,434,981 |

The values at $k=10,11$ correct transcription errors in the conjectured entry; the order-100 value was already correct. This is not the distinct pair-structure staircase HYP-1732.

## The span-three obstruction is this same staircase in another order

Swap disjoint adjacent pairs in the order, displaying $1,0,3,2,\dots$ and leaving the final even label fixed when $n$ is odd. For even $n$ the only visible backward arcs are from odd positions $j\ge3$ to $j-3$. The old backward arcs within swapped pairs are now forward; those connecting the pairs become the span-three arcs. Thus the even family without any single-vertex strong seam in the [full-alphabet span-three packet](../../05-knowledge/results/chanpak_local_rows_span3_20261115/SPAN3.md#4-what-this-gives-toward-full-alphabet-growth) has exact growth above $1.7699$. It obstructs that covering argument, not the proposed growth bound.

More generally, for **any induced subset** of the staircase, split its old labels into maximal intervals of consecutive integers. There are only forward arcs between the pieces, so its count is the product of $A_m$ over their lengths $m$. In particular this gives every window count in the swapped presentation, including windows that break an end pair. The same product method enumerates the entire span-one language in its original order: split at each forward adjacent edge. Inside the resulting pieces all adjacent edges are backward and outside the pieces all arcs are forward.

## Exact defect for a proposed wider seam

Take pair-aligned windows of sizes $2a,2b$ sharing exactly one full swapped pair, $a,b\ge2$. The intersection count is $h_1=1$, but the scalar lower bound $h_{a+b-1}\ge h_a h_b/h_1$ fails **for every** such pair. Indeed $D(a,b)=h_a h_b-h_{a+b-1}$ obeys the same positive order-three recurrence separately in each index. The base square $a,b\in\{2,3,4\}$ has entries
$$
\begin{pmatrix}8&28&92\\28&96&316\\92&316&1040\end{pmatrix}.
$$
First propagate in $a$ for the three base columns and then in $b$. Thus the failure persists for arbitrary window sizes, even when both windows are strong. It reflects loss of endpoint information, not a small-order accident. The exact proof carries a two-state boundary datum: whether the preceding part was a singleton. Its generating transfer for a single new part is
$$
\begin{pmatrix}B&x\\ B&0\end{pmatrix}
$$
in row-to-column convention (preceding part nonsingleton or singleton, respectively; the start acts as a nonsingleton). Multiplying this transfer sums legal part sequences; it is not a proposed universal transfer for arbitrary span-three words.

The direct checks use subset DP: $D(S,v)=\sum_{u\in S-\{v\}:u\to v}D(S-\{v\},u)$, starting at $D(\{v\},v)=1$. In addition to the 16 original-order checks, all 816 intervals of the swapped orders $1\le n\le16$ agree with the induced-subset product. Neither general span-three growth nor Chan–Pak follows from this subfamily.
