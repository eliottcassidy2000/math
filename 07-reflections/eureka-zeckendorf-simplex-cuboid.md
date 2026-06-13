# Eureka, Zeckendorf, and the Simplex–Cuboid Duality

**Session:** oracle-2026-05-16-S1 (continued)
**Depends on:** `simplex-cuboid-duality.md`, `zeckendorf-non-consecutivity-pairing.md`,
`interleaved-staircase-binary-grid.md`
**Computations:** oracle-2026-05-16

---

## The Three Theories, One Object

The same object — the tournament $H$ value and its conflict graph — sits at the
intersection of three classical theorems:

| Theorem | Basis | Fugacity | Output |
|---|---|---|---|
| **Gauss Eureka** ($n = \Delta+\Delta+\Delta$) | Triangular numbers | — | Edge count factorization |
| **Zeckendorf** ($n = \sum F_k$, non-consecutive) | Fibonacci numbers | $x=1$ | $I(\text{path}, 1) = F_{k+2}$ |
| **Tournament OCF** ($H = I(\Omega, 2)$) | Independence poly | $x=2$ | $I(\text{path}, 2) = J_{k+2}$ |

The fugacity $x$ is the bridge: at $x=1$ we get Fibonacci/Lucas numbers; at $x=2$ we get
Jacobsthal numbers and tournament H values. Gauss Eureka governs the EXPONENTS (edge counts),
while Zeckendorf/OCF govern the POLYNOMIAL EVALUATIONS.

---

## Section 1: Eureka and the Simplex–Cuboid Duality

**Gauss's Eureka theorem:** Every positive integer is the sum of at most 3 triangular
numbers $T_k = \binom{k+1}{2}$.

In our framework, $T_{k-1} = \binom{k}{2}$ is the edge count of the $k$-vertex simplex
$K_k$, and $S(k) = 2^{T_{k-1}} = 2^{\binom{k}{2}}$ is the labeled tournament count.

**Consequence:** Every labeled tournament count $2^m$ is a product of at most 3 simplex-
function values:
$$m = T_a + T_b + T_c \implies 2^m = S(a+1) \cdot S(b+1) \cdot S(c+1)$$

**The square-triangle Eureka** ($k^2 = T_{k-1} + T_k$) is the specific case used in the
cuboid identity: the cuboid crossing $K_{k,k}$ has $k^2$ edges, which is the sum of the
upper ($T_k = \binom{k+1}{2}$ edges, $\equiv K_{k+1}$) and lower ($T_{k-1} = \binom{k}{2}$
edges, $\equiv K_k$) triangular halves:
$$C(k) = 2^{k^2} = S(k) \cdot S(k+1) \quad \Longleftrightarrow \quad k^2 = T_{k-1} + T_k$$

**The total edge count** $\binom{2k}{2} = T_{2k-1}$ is itself ONE triangular number — the
trivially minimal Eureka decomposition. The bipartite splitting gives an ALTERNATIVE
decomposition $T_{2k-1} = k + 4T_{k-1}$ (Eureka applied to $k$ itself gives 3 more
triangulars), yielding the 4-simplex formula:
$$S(2k) = 2^k \cdot S(k)^4 \quad \Longleftrightarrow \quad T_{2k-1} = k + 4T_{k-1}$$

The $k$ "diagonal bits" of the staircase construction are the difference between the
one-triangle Eureka ($T_{2k-1}$) and the four-triangle bipartite decomposition
($4T_{k-1}$): they account for the remainder $k = T_{2k-1} - 4T_{k-1}$, which itself
requires at most 3 more triangular numbers by Eureka. So the FULL chain is:

$$T_{2k-1} = \underbrace{4 T_{k-1}}_{\text{4 simplex copies}} + \underbrace{T_a+T_b+T_c}_{\text{diagonal bits, } \leq 3}$$

for some $a,b,c \geq 0$ with $T_a+T_b+T_c = k$. Total: at most 7 triangular summands.

**The polygonal connection.** The $k$-gonal numbers satisfy: square ($4$-gonal) = sum of
consecutive triangular numbers. This is the identity $k^2 = T_{k-1}+T_k$ we've been using —
**the tournament's simplex-cuboid duality is the $k$-gonal identity for $k=4$ (squares).**

---

## Section 2: Zeckendorf and the OCF Fugacity Shift

From `zeckendorf-non-consecutivity-pairing.md` (proved):

**The Zeckendorf basis IS the independence basis** at fugacity $x=1$. The non-consecutive
condition (Zeckendorf's theorem) is exactly the independence condition in the path graph
$P_\infty$ (Fibonacci graph), evaluated at unit weight per element.

**The fugacity shift** $x=1 \to x=2$:

| Graph $G$ | $I(G,1)$ | $I(G,2)$ | Connection |
|---|---|---|---|
| Path $P_k$ | $F_{k+2}$ (Fibonacci) | $J_{k+2}$ (Jacobsthal) | Basis at $x=1$ → Tournament H at $x=2$ |
| Cycle $C_m$ | $L_m$ (Lucas) | $2^m + (-1)^m$ (Mersenne/Fermat-type) | Proved in prior session |
| Complete $K_k$ | $1 + kx\big|_{x=1} = k+1$ | $1 + 2k$ | Linear in $k$ |

**Jacobsthal numbers** $J_n = (2^n - (-1)^n)/3$: the tournament-fugacity analog of
Fibonacci numbers for path conflict graphs. Sequence: $0,1,1,3,5,11,21,43,85,171,\ldots$

The connection: $J_{k+2} = I(P_k, 2)$ is the H value of any tournament whose conflict
graph is a path on $k$ vertices.

**The forbidden values** through this lens:
- $H = 7$ is forbidden ↔ $I(C_3, 2) = 2^3 - 1 = 7$: the cycle $C_3 = K_3$ appears in
  the conflict graph, but the tournament constraints prevent this from occurring.
- $H = 21$ is forbidden ↔ $21 = 3 \times 7$: compounded $C_3$ obstruction.

---

## Section 3: The All-Zero Staircase H Sequence

The all-zero staircase at $n=2k$ (all recessives beat their pair-dominants) gives the
maximum H within the staircase family:

| $k$ | $n$ | $H_0(k)$ | Fibonacci? | Lucas? | Zeckendorf |
|---|---|---|---|---|---|
| 1 | 2 | 1 | $F_1=F_2$ | $L_1$ | trivial |
| 2 | 4 | 5 | $F_5$ | — | $F_5$ (1 term) |
| 3 | 6 | 29 | — | $L_7 = F_6+F_8$ | $F_6+F_8$ (2 terms) |
| 4 | 8 | 233 | $F_{13}$ | — | $F_{13}$ (1 term) |
| 5 | 10 | 2489 | — | — | $F_3+F_7+F_9+F_{13}+F_{15}+F_{17}$ (6 terms) |
| 6 | 12 | 33773 | — | — | $F_4+F_{11}+F_{13}+F_{15}+F_{19}+F_{23}$ (6 terms) |

**Pattern observations:**

1. $k=2,4$: $H_0 = F_5, F_{13}$ — single Fibonacci numbers at indices $5 = F_5$ and
   $13 = F_7$. Fibonacci indices are themselves Fibonacci numbers.

2. $k=3$: $H_0 = 29 = L_7$ (Lucas number). And $29 = F_6 + F_8 = I(P_2, 2) \cdot I(P_4,2)/(?)$...
   Actually $29 = L_7 = F_6 + F_8$ (from the Lucas-Fibonacci identity $L_m = F_{m-1}+F_{m+1}$,
   here $L_7 = F_6+F_8$). The Zeckendorf representation is exactly the Lucas non-consecutive
   Fibonacci pair formula.

3. $k=5,6$: 6 Zeckendorf terms each. The representations share index 13 (and 15 for $k=6$).

4. **The Jacobsthal connection**: $H_0(2) = 5 = J_4$ (Jacobsthal). So the first non-trivial
   all-zero staircase value is a Jacobsthal number, consistent with it being the H value of
   a tournament whose conflict graph has Jacobsthal structure.

**No simple recurrence found** for this sequence: the candidates $(2k-1)a_{k-1}+(2k-2)a_{k-2}$
(= the Cayley-Dickson-flavored recurrence) match for $k=3,4$ but fail at $k=5$.

---

## Section 4: The Unified Picture

The **triangular-number basis** (Eureka) governs EXPONENTS:
$$\text{edge counts} = \binom{k}{2} = T_{k-1}$$
$$S(2k) = 2^{T_{2k-1}}, \quad C(k) = 2^{T_{k-1}+T_k}$$

The **Fibonacci/Lucas basis** (Zeckendorf) governs H VALUES at $x=1$:
$$I(P_k, 1) = F_{k+2}, \quad I(C_m, 1) = L_m$$

The **Jacobsthal basis** (OCF at $x=2$) governs tournament H values:
$$I(P_k, 2) = J_{k+2}, \quad I(C_m, 2) = 2^m + (-1)^m$$

The three layers are connected by the fugacity shift and the simplex-cuboid structure:

```
Triangular numbers (Eureka)
    ↓  [exponent → labeled count]
Simplex/Cuboid functions S(k), C(k)
    ↓  [Burnside quotient]
A000568(n) iso class counts
    ↓  [OCF: I(Ω,x) at x=1,2]
Fibonacci/Lucas (x=1) and Jacobsthal/Mersenne (x=2)
    ↑  [Zeckendorf: independence condition = non-consecutiveness]
```

**The square-triangle identity** $k^2 = T_{k-1}+T_k$ appears at the TOP of this chain
(simplex-cuboid duality), while the Zeckendorf/OCF connection appears at the BOTTOM
(H-value structure). The `eureka formula for k-gonal numbers connects the two:

- **$k=3$ (triangular)**: $T_n = \binom{n+1}{2}$ = simplex edge count. Eureka applies here.
- **$k=4$ (square)**: $k^2 = T_{k-1}+T_k$ = cuboid edge count. Square-triangle identity.
- **$k=5$ (pentagonal)**: $P_n = n(3n-1)/2$. The pentagonal number theorem (Euler) controls
  partition function recurrences — another independence-polynomial structure.

**The forbidden value $H=7$** lives at the BOUNDARY of these layers: it's simultaneously
$I(C_3,2) = 7$ (OCF layer), the forbidden tournament value (iso-class layer), and related to
the triangular number $T_3 = 6 \approx 7$ (edge-count layer). The three layers "conspire"
to make 7 forbidden.

---

## Open Directions

1. **Eureka decomposition of $H$ values**: Given $H = I(\Omega, 2)$, can $\Omega$ always be
   decomposed into at most 3 "simplex components" (cliques, paths, cycles) so that
   $H = \prod$ simplex evaluations? (Analogous to Eureka for the independence polynomial.)

2. **The $k=5,6$ Zeckendorf pattern**: Why do $H_0(5)$ and $H_0(6)$ both have 6 Zeckendorf
   terms? Is there a pattern by Hamming weight of the diagonal bits?

3. **Polygonal H values**: Is there a class of tournaments whose H equals the $k$-th
   $m$-gonal number for each $m$? Triangular $H = T_n$ happens for transitive extensions;
   square $H = n^2$ would correspond to some cuboid-structured conflict graph.

4. **The fugacity $x=\phi$**: At $x = \phi$ (golden ratio), $I(P_k, \phi)$ satisfies a
   clean recursion (since $\phi^2 = \phi+1$). This intermediate fugacity bridges Zeckendorf
   ($x=1$) and OCF ($x=2$) — does it correspond to any natural tournament invariant?
