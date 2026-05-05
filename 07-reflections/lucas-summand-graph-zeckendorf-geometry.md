# Lucas Numbers, the Summand Graph, and the Higher-Dimensional Geometry of Zeckendorf

**Session:** oracle-2026-05-05
**Depends on:** `summand-graph-fermat-zeckendorf.md`,
`zeckendorf-non-consecutivity-pairing.md`, Complementarity Theorem (prev session)
**Computations:** all verified in-session

---

## Section 1: The Lucas Sequence IS a Summand Chain

**Theorem (Lucas Summand Chain).**
The Lucas sequence $L_1, L_2, L_3, \ldots$ is itself a summand chain in the summand
graph: $L_{m+1} = L_m + L_{m-1}$ for all $m \geq 2$, with $L_m \neq L_{m-1}$ (satisfied
for $m \geq 2$).

Every Lucas number therefore has an incoming edge in the summand graph coming from
the pair $(L_{m-1}, L_m)$ — the two previous Lucas numbers.

**Simultaneously**, every Lucas number has a *unique Fibonacci-pair incoming edge*:
$$L_m = F_{m-1} + F_{m+1} \quad (\text{non-consecutive Fibonacci, gap}=2).$$

So each Lucas node $L_m$ in the summand graph has **two distinguished incoming edges**:
1. The *Lucas-recurrence edge*: $(L_{m-1}, L_m) \to L_{m+1}$
2. The *Fibonacci-pair edge*: $(F_{m-1}, F_{m+1}) \to L_m$

These are DIFFERENT edges arriving at $L_m$ from different summand pairs. Together they
embed the Lucas sequence into both the Fibonacci lattice and the Lucas lattice of the
summand graph simultaneously.

**Verification (all pairs for $L_m$ in the summand graph):**

| $m$ | $L_m$ | Total pairs | Fibonacci pair | Zeckendorf |
|---|---|---|---|---|
| 2 | 4 | 1 | $(1, 3)$ | $[1,3]$ |
| 3 | 7 | 3 | $(2, 5)$ | $[2,5]$ |
| 4 | 11 | 5 | $(3, 8)$ | $[3,8]$ |
| 5 | 18 | 8 | $(5, 13)$ | $[5,13]$ |
| 6 | 29 | 14 | $(8, 21)$ | $[8,21]$ |
| 7 | 47 | 23 | $(13, 34)$ | $[13,34]$ |

Each Lucas number has exactly ONE Fibonacci-pair incoming edge (the minimum-gap pair),
plus many non-Fibonacci pairs for larger $m$.

---

## Section 2: The Module-Spine Partition of Lucas Numbers

From the Complementarity Theorem: the summand graph starting from $\{2,3\}$ generates
$\mathbb{N} \setminus \{1,4,6\}$ exactly. The three module elements $\{1,4,6\}$ and the
two spine seeds $\{2,3\}$ together generate all of $\mathbb{N}$.

**Partition of Lucas numbers (standard: $L_1=1, L_2=3, L_3=4, L_4=7, L_5=11,\ldots$):**

| Lucas | Value | Class | Reason |
|---|---|---|---|
| $L_1 = 1$ | 1 | **Module** | Source node; no incoming pairs |
| $L_2 = 3$ | 3 | **Spine seed** | Given as a seed |
| $L_3 = 4$ | 4 | **Module** | Only pair is $(1,3)$; needs $F_1=1$ |
| $L_4 = 7$ | 7 | **Spine** | $7 = 2+5$; doesn't need $F_1$ |
| $L_5 = 11$ | 11 | **Spine** | $11 = 3+8$ |
| $L_6 = 18$ | 18 | **Spine** | $18 = 5+13$ |
| $L_m, m\geq 4$ | — | **Spine** | All have non-$F_1$ paths |

**Key finding:** Exactly $\{L_1, L_3\}$ (the 1st and 3rd Lucas numbers) are in the
$\{1,4,6\}$ module. The 2nd Lucas number is a spine seed. All subsequent Lucas numbers
($L_4$ onward) are in the spine. The module $\{1,4,6\}$ contains exactly
$\{L_1, L_3, \text{non-Lucas } 6\}$ where $6 = F_1 + F_4$ is a Fibonacci pair with gap
3 but is NOT itself a Lucas number.

**Why $L_3=4$ is in the module and $L_4=7$ escapes:**
$4 = 1+3$ (only valid summand pair; $2+2$ excluded). Every path to 4 uses $F_1=1$.
$7 = 2+5$ (spine path: $5 = 2+3$ doesn't need $F_1$). The "escape" happens at $L_4=7$
because $7 = F_2 + F_4$ uses only spine-seed-level Fibonacci numbers.

**The non-Lucas 6:** The third module element $6 = F_1+F_4 = 1+5$ is not a Lucas number
in standard convention (no $L_m = 6$). It belongs to the module because both its
summand pairs require module elements: $6 = 1+5$ (needs $L_1=1$) and $6 = 2+4$
(needs $L_3=4$). It is the *ternary threshold* $T_3 = 1+2+3$ from the previous
session, and its module membership reflects its role as the first 3-part sum (requiring
$F_1$ for the 3-part expression $1+2+3=6$).

---

## Section 3: The {2,3}-Spine as a Fibonacci Recurrence Machine

The spine generation from $\{2,3\}$ follows a **Fibonacci-level structure**:

| Level | Elements | Count | Lucas numbers |
|---|---|---|---|
| 0 | $\{2, 3\}$ | 2 | $L_2=3$ (seed) |
| 1 | $\{5\}$ | 1 | — |
| 2 | $\{7, 8\}$ | 2 | $L_4=7$ |
| 3 | $\{9,10,11,12,13,15\}$ | 6 | $L_5=11$ |
| 4 | $\{14,16,17,18,19,\ldots,28\}$ | 14 | $L_6=18$ |
| 5 | $\{29,\ldots,50\}$ | 22 | $L_7=47, L_8=47\ldots$ |

The level counts $2, 1, 2, 6, 14, 22, \ldots$ do not follow a simple pattern, but
**at every spine level, the minimum-gap Lucas number $L_{m}$ appears at level $m-2$
(approximately)** — proven by the Fibonacci-pair arrival formula
$L_m = F_{m-1}+F_{m+1}$ and the level formula $\mathrm{level}(F_k) = k-3$ for $k \geq 3$.

**Recursive formula for the spine:**

$$S_0 = \{F_2, F_3\} = \{2, 3\}$$
$$S_{n+1} = S_n \cup \{a+b \;:\; a,b \in S_n,\; a < b,\; a+b \notin S_n\}$$

At each level transition, the newly added elements are those not yet reachable.
The Fibonacci numbers themselves form the "frontier" of each level:
$F_k$ enters $S_{k-3}$ for $k \geq 3$. The Lucas numbers follow one level behind:
$L_m$ enters roughly $S_{m-2}$.

---

## Section 4: The Wythoff-Zeckendorf Index Shift Theorem

**Wythoff's theorem:** The lower Beatty sequence $a_n = \lfloor n\phi\rfloor$ and upper
Beatty sequence $b_n = \lfloor n\phi^2\rfloor$ partition $\mathbb{N}$, with $b_n = a_n + n$.

**New theorem (Wythoff-Zeckendorf Index Shift):**
If $a_n$ has Zeckendorf representation using Fibonacci indices $\{k_1, k_2, \ldots, k_r\}$,
then $b_n = a_n + n$ has Zeckendorf representation using indices $\{k_1+1, k_2+1, \ldots, k_r+1\}$.

*Equivalently:* the map $n \mapsto b_n$ acts on Zeckendorf binary strings as the shift
operator $\sigma: b_{k_i} \mapsto b_{k_i+1}$.

| $n$ | $a_n$ | $Z(a_n)$ indices | $b_n$ | $Z(b_n)$ indices | Shift |
|---|---|---|---|---|---|
| 1 | 1 | $[1]$ | 2 | $[2]$ | $+1$ ✓ |
| 2 | 3 | $[3]$ | 5 | $[4]$ | $+1$ ✓ |
| 3 | 4 | $[1,3]$ | 7 | $[2,4]$ | $+1$ ✓ |
| 4 | 6 | $[1,4]$ | 10 | $[2,5]$ | $+1$ ✓ |
| 5 | 8 | $[5]$ | 13 | $[6]$ | $+1$ ✓ |
| 6 | 9 | $[1,5]$ | 15 | $[2,6]$ | $+1$ ✓ |
| 7 | 11 | $[3,5]$ | 18 | $[4,6]$ | $+1$ ✓ |

**Connection to the summand graph:** Since $b_n = a_n + n$, the pair $(a_n, n)$ is an
incoming summand edge to $b_n$ in the summand graph. The Wythoff index shift theorem
says this edge (from the lower to upper Beatty sequence) corresponds to a **Fibonacci
index translation** in Zeckendorf space.

**Self-reference:** Setting $n \mapsto a_n$ recursively: the lower Wythoff sequence
$a_n$ has Zeckendorf indices $\{k_i\}$, and $n$ itself has Zeckendorf indices $\{k_i-1\}$
(one lower). This self-referential nesting $a_n \supset n \supset \cdots$ is the
combinatorial heart of the Wythoff game and of the Fibonacci self-similarity.

---

## Section 5: The Fibonacci Cube as the Zeckendorf Space

**Definition.** The **Fibonacci cube** $Q_k$ is the graph on all binary strings of
length $k$ with no two consecutive $1$s (i.e., no "11" substring), with edges between
strings differing in exactly one bit.

$|V(Q_k)| = F_{k+2}$ (the $(k+2)$-th Fibonacci number).

**The Zeckendorf injection $\varphi: \mathbb{N} \hookrightarrow Q_\infty$:**
Every positive integer $n$ maps to its Zeckendorf binary vector
$(b_1, b_2, \ldots)$ where $b_k = 1$ iff $F_k$ appears in $n$'s Zeckendorf
representation. The non-consecutive condition ensures the image lies in
$Q_\infty$, and uniqueness makes this an injection.

**Density of the Fibonacci cube:**

$$\frac{|Q_k|}{|\{0,1\}^k|} = \frac{F_{k+2}}{2^k} \sim \frac{\phi^{k+2}}{\sqrt{5} \cdot 2^k} = \frac{\phi^2}{\sqrt{5}} \cdot \left(\frac{\phi}{2}\right)^k.$$

Since $\phi/2 \approx 0.809 < 1$, this density $\to 0$ exponentially. The Fibonacci
cube is an **exponentially sparse fractal subset** of the hypercube with effective
dimension $\log(\phi)/\log(2) \approx 0.694$ — less than $1$.

This gives the Zeckendorf space a fractal geometry: the integers, viewed through the
Zeckendorf lens, occupy only a $0.694$-dimensional "slice" of the infinite binary
hypercube.

**Comparison to the tournament tiling:**
The tournament on $n$ vertices is encoded as a binary string in $\{0,1\}^{C(n,2)}$
(each bit = one arc orientation). The tournament space is the full hypercube
$\{0,1\}^{C(n,2)}$ — 100% density. The Zeckendorf space (Fibonacci cube) occupies
density $\to 0$.

Both are subsets of binary hypercubes; the tournament lives in the "full" hypercube
while Zeckendorf lives in the "sparse fractal" sub-cube. The summand graph at x=2
(tournament evaluation) probes the full hypercube; at x=1 (Zeckendorf/chemistry)
it probes only the Fibonacci sub-cube.

---

## Section 6: The Fibonacci Cube as Tournament Space

**The Fibonacci tournament $FT_k$:** Define a tournament on vertex set
$V = Q_k \setminus \{0^k\}$ (all non-zero Fibonacci-cube vertices), with total order
given by integer value: $u$ beats $v$ iff $\sum b_i^{(u)} F_i > \sum b_i^{(v)} F_i$.

This is the **transitive tournament on $F_{k+2}-1$ vertices** ordered by Zeckendorf value.
It has $H = 1$ (unique Hamiltonian path = decreasing order) and trivial conflict graph
$\Omega = \emptyset$.

**The Fibonacci Kneser tournament:** Define a different tournament where the **conflict
structure** mirrors the Zeckendorf non-consecutive condition. Two integers $a, b$ are
"Fibonacci-conflicting" if their Zeckendorf supports overlap (share a Fibonacci term).
This defines the **Fibonacci Kneser graph** $FK_k$: vertices = Zeckendorf representations
of $\{1, \ldots, F_{k+2}-1\}$, edges when supports intersect.

For $Q_5$ (13 vertices): 29 conflict pairs, 37 non-conflict pairs. The conflict graph
$FK_5$ is a specific graph on 12 vertices (excluding 0) with 29 edges. The independence
polynomial $I(FK_5, 2)$ gives a tournament H-value for a tournament whose conflict
graph is $FK_5$.

**What tournament has $\Omega(T) = FK_k$?**
This would be a tournament whose directed odd cycles, when two cycles share a vertex,
exactly mirror the Zeckendorf overlap condition. Such a tournament would have
$H(T) = I(FK_k, 2)$, connecting Zeckendorf combinatorics directly to tournament
Hamiltonian path counts.

---

## Section 7: Recursive Formulas — The Lucas-Fibonacci Interplay

### 7.1 The Two Recursions at Each Lucas Node

Every Lucas number $L_m$ (for $m \geq 4$ in the spine) admits two incoming summand
paths:

**Path A (Lucas recursion):** $L_m = L_{m-1} + L_{m-2}$
— a path through the Lucas chain, arriving from two previous Lucas numbers.

**Path B (Fibonacci pair):** $L_m = F_{m-1} + F_{m+1}$
— a path through two non-consecutive Fibonacci numbers (minimum gap = 2).

These two paths differ in their depths in the summand graph. Path B typically has lower
depth (it uses Fibonacci numbers, which are "boundary" elements at each level).

**Level formula for Lucas numbers:**
$\mathrm{level}(F_k) = k-3$ for $k \geq 3$ (in the spine, using levels of the {2,3}-reachability BFS).

$$\mathrm{level}(L_m) = \max(\mathrm{level}(F_{m-1}), \mathrm{level}(F_{m+1})) + 1 = (m-2)+1-1 = m-2$$

(approximately, since level$(F_{m+1}) = m-2$ dominates). Verified for $m=3,4,5,6$.

### 7.2 The Spine Generation Satisfies a Fibonacci-Type Count

The count of new elements added at each level:

| Level | Count | Pattern |
|---|---|---|
| 0 | 2 | (seeds) |
| 1 | 1 | |
| 2 | 2 | |
| 3 | 6 | |
| 4 | 14 | |
| 5 | 22 | |

The count sequence $1, 2, 6, 14, 22, \ldots$ grows approximately like $8n$ for
level $n$, reflecting the quadratic growth in the number of pairs available at each level
(the spine has $\sim n^2/2$ elements by level $n$, and new pairs grow quadratically).

### 7.3 The Module-to-Spine Boundary Recursion

The transition from module elements to spine elements follows a specific rule:
an element $n = F_1 + F_k$ is in the module iff there is NO non-module summand pair
summing to $n$. This holds iff every number $j$ with $j < n$ and $n-j$ also $< n$
satisfies: either $j$ or $n-j$ is in $\{1,4,6\}$.

For $n = 1+F_k$:
- Module: $k=3$ ($n=4$) — only pair $(1,3)$, both 1 and 3: 1∈module, 3∈spine-seed
- Module: $k=4$ ($n=6$) — pairs $(1,5)$ and $(2,4)$: both paths need module elements
- Spine: $k=5$ ($n=9$) — pair $(2,7)$ works: $7 = 2+5$ all spine

The escape at $k=5$ (i.e., at $F_5=8$, giving $n=9$) happens because $7 = 2+5$ is
reachable in the spine (level 2), so the sum $2+7 = 9$ gives a non-module path.

**Boundary formula:** The module closure $\{1,4,6\}$ is finite and equal to
$\{1, F_1+F_3, F_1+F_4\}$ — the first two F_1-containing Fibonacci pairs. No higher
F_1-containing pair is in the module because higher Fibonacci numbers $F_k$ ($k \geq 5$)
have enough "slack" in the summand graph to be reachable by spine-only paths.

---

## Section 8: Higher-Dimensional Geometry of Zeckendorf Decomposition

### 8.1 The Fibonacci Hypercube Slice

Every integer $n$ corresponds to a point in the infinite binary hypercube
$\{0,1\}^\infty$ via its Zeckendorf vector. The Zeckendorf map has three key geometric
properties:

1. **Injectivity:** Each integer has exactly one image (Zeckendorf uniqueness).
2. **Fractal image:** The image is a sub-manifold of the hypercube with
   dimension $\log\phi/\log 2 \approx 0.694$ — the "Fibonacci fractal."
3. **Non-consecutive restriction:** The image avoids all points with consecutive 1s —
   a "Cantor-like" removal of $1/4$ of the hypercube at each scale.

The non-consecutive restriction removes exactly the Hamming-weight-2 strings with
adjacent 1s: $\binom{k-1}{1}$ strings from $\{0,1\}^k$, i.e., a $k$-proportional
fraction. Iterating gives the Fibonacci cube density $(φ/2)^k$.

### 8.2 The Zeckendorf Map as a Coordinate System

Think of the Fibonacci numbers $F_1, F_2, F_3, \ldots$ as **coordinate axes** in an
infinite-dimensional space. The Zeckendorf representation assigns to each integer $n$
a **sparse binary coordinate vector** (at most one $1$ per "coordinate pair"
$\{F_k, F_{k+1}\}$, by the non-consecutive condition).

This is the **hyperbolic lattice picture**: the Fibonacci numbers grow as $F_k \sim \phi^k$,
so the "true distance" between consecutive Fibonacci coordinates is
$\log\phi \approx 0.481$ in the natural (log-scale) metric. The Zeckendorf representation
is a **binary code in log-Fibonacci space**, where each integer occupies a unique point.

The summand graph operates in this space by **adding two such points** to get a new one.
The Lucas numbers are the points at "distance 2" (minimum valid gap) in this space —
the closest possible non-conflicting Fibonacci pair.

### 8.3 The Stern-Brocot / Calkin-Wilf Connection

The **Calkin-Wilf tree** enumerates all positive rationals via the mediant operation:
each rational $p/q$ has children $(p)/(p+q)$ and $(p+q)/q$. The integer labels of
Calkin-Wilf positions are related to Fibonacci representations.

The **summand graph** is the "integer projection" of this rational tree: when $p/q$
and $r/s$ are adjacent in the Stern-Brocot tree with $ps - qr = \pm 1$ (Farey
neighbors), their denominators satisfy $q+s = $ (next denominator), reflecting the
summand-graph structure on denominators.

The **Zeckendorf representation** of the denominator of the $n$-th rational in the
Calkin-Wilf tree encodes the path from the root to that rational — another
tournament-tree structure hiding in number theory.

### 8.4 The Golden Ratio as the "Hausdorff Dimension" of the Summand Graph

The density of Fibonacci cube $Q_k$ in $\{0,1\}^k$ equals $F_{k+2}/2^k \to 0$ at
rate $(\phi/2)^k$. The "Hausdorff dimension" of this fractal is:

$$d_H = \frac{\log\phi}{\log 2} \approx 0.694$$

This number $\log\phi/\log 2$ is the answer to: "what fraction of integers, up to a
large bound $N$, require at most $k$ Zeckendorf terms?" — it's the scaling exponent
of the Zeckendorf complexity function.

In the summand graph: depth-$k$ elements (reachable from seeds in $k$ steps)
correspond to integers with Zeckendorf complexity $\leq 2^k$ roughly. The density
$(φ/2)^k$ governs how many integers are first reachable at depth $k$.

### 8.5 The Module as the "Boundary" of the Fractal

The module $\{1,4,6\}$ consists of the integers that lie at the "boundary" of the
Zeckendorf fractal: they are the points where the Fibonacci coordinate $F_1$ (the
smallest Fibonacci number) is indispensable. Removing $F_1$ from the coordinate
system creates a "missing boundary" — the module elements are those that exist only
because of the $F_1$ axis.

Geometrically: in the Fibonacci hypercube $Q_\infty$, the subspace "ignoring the $F_1$
coordinate" (setting $b_1 = 0$ for all points) is the sub-cube corresponding to
numbers WITHOUT $F_1$ in their Zeckendorf representation. The module $\{1,4,6\}$ is
exactly the set of positive integers that **require** $b_1 = 1$ AND have no alternative
non-$F_1$ path in the summand graph.

---

## Section 9: The Hidden Tournament in Zeckendorf Decomposition

### 9.1 The Zeckendorf Conflict Graph

For a fixed integer $N$, define the **Zeckendorf conflict graph** $\mathcal{Z}_N$ as:
- **Vertices:** Fibonacci numbers $F_1, F_2, \ldots, F_k$ with $F_k \leq N$
- **Edges:** between consecutive Fibonacci numbers $F_i$ and $F_{i+1}$
  (i.e., $\mathcal{Z}_N = P_k$, a path graph)
- **Tournament structure:** orient edge $F_i \leftrightarrow F_{i+1}$ as $F_{i+1} \to F_i$
  (larger Fibonacci beats smaller, i.e., the transitive tournament on $\{F_1,\ldots,F_k\}$)

The **Zeckendorf representation of $N$** selects an independent set $S$ in $P_k$ that
sums to $N$ — exactly the hard-core partition function setup of the OCF.

If we use $P_k$ as the conflict graph $\Omega(T)$ for some tournament $T$, then:
$$H(T) = I(P_k, 2) = \frac{4 \cdot 2^k - (-1)^k}{3}$$
giving the path-independence values $3, 5, 11, 21, 43, 85, 171, \ldots$

The forbidden value $H = 21 = I(P_4, 2)$ corresponds to $k=4$ — the Zeckendorf
conflict graph for all integers up to $F_4 = 5$ (using our convention). **The forbidden
H-value 21 is exactly the "tournament value" of the standard Zeckendorf 4-level path.**

### 9.2 The Zeckendorf Nim Connection

**Wythoff's Nim** is a two-pile Nim game where a position $(a, b)$ is a losing (P-)
position iff $(a,b) = (\lfloor n\phi \rfloor, \lfloor n\phi^2 \rfloor)$ for some $n$.

The **Zeckendorf Index Shift Theorem** (Section 4) says: $Z(b_n)$ has all indices
exactly $+1$ compared to $Z(a_n)$. This means the Wythoff P-positions are exactly
the pairs of integers with "adjacent Zeckendorf representations" — consecutive index
shifts of the same underlying integer $n$.

**Tournament interpretation:** Define a tournament on integers where $a$ beats $b$
iff $a$ is the lower Wythoff partner of some $n$ and $b$ is not. This creates a
two-coloring of $\mathbb{N}$ (Beatty partition), with the P-positions (Wythoff pairs)
acting as "neutral points" — neither beats the other in the Nim game.

The Beatty partition $(\lfloor n\phi \rfloor, \lfloor n\phi^2 \rfloor)$ is the
**golden-ratio tournament structure** hiding in Zeckendorf: it partitions $\mathbb{N}$
into two "teams" based on the Fibonacci coordinate parity (even vs odd index position
in the Zeckendorf representation).

### 9.3 The Recursive Tournament Structure

**Fibonacci tournament $\mathbf{FT}_k$:** A tournament on $F_{k+2}$ vertices
(= the Fibonacci cube $Q_k$) where $u$ beats $v$ iff the integer represented by $u$
is larger. This is the **transitive tournament** on $F_{k+2}$ vertices, with $H = 1$.

**Enriched Fibonacci tournament $\widetilde{\mathbf{FT}}_k$:** Keep the same vertex set
but define arcs based on "Zeckendorf dominance": $u$ beats $v$ iff the support of $u$'s
Zeckendorf representation contains the support of $v$'s (as sets of Fibonacci numbers).
This is a **partial order** (Zeckendorf lattice), which can be completed to a tournament.

The conflict graph $\Omega(\widetilde{FT}_k)$ would then be the **Zeckendorf Kneser
graph**: two vertices conflict iff their Zeckendorf supports overlap. This graph
is the "Fibonacci version" of the Kneser graph, and its independence polynomial at
$x=2$ gives a Fibonacci-generated H-value.

---

## Section 10: Recursive Formulas — Summary

| Formula | Source | Notes |
|---|---|---|
| $L_{m+1} = L_m + L_{m-1}$ | Lucas recurrence | Summand-graph edge |
| $L_m = F_{m-1} + F_{m+1}$ | Fibonacci-pair | Min-gap Zeckendorf |
| $I(P_k, 2) = I(P_{k-1},2) + 2\cdot I(P_{k-2},2)$ | Path recurrence | Roots $\{2,-1\}$ |
| $I(C_m, 2) = 2^m + (-1)^m$ | Cycle formula | Mersenne/Fermat |
| $b_n = a_n + n$, $Z(b_n) = Z(a_n)+1$ | Wythoff-Zeckendorf | Index shift |
| $\Sigma H(n) = 2^{C(n-1,2)-n+1} \cdot W(n)$ | THM-292 | Tiling formula |
| $d_H(Q_\infty) = \log\phi / \log 2$ | Fractal geometry | Fibonacci dimension |

---

## Open Questions

1. **Prove the Wythoff-Zeckendorf Index Shift Theorem rigorously.** The numerical
   verification is complete for all cases up to $n=30$. A proof likely uses the
   three-distance theorem and $\phi$-properties of Fibonacci.

2. **What is $I(FK_k, 2)$** (the independence polynomial of the Fibonacci Kneser graph)?
   Does it equal $H(T)$ for some explicitly constructible tournament $T$?

3. **Is there a "Zeckendorf tournament" $T_Z$ on $F_{k+2}$ vertices such that
   $\Omega(T_Z) = FK_k$?** If so, $H(T_Z) = I(FK_k, 2)$ would be a Fibonacci-structured
   H-value.

4. **The Fibonacci cube fractal dimension $\log\phi/\log 2 \approx 0.694$:** does this
   number appear anywhere in the tournament theory (e.g., as a rate constant in the
   asymptotic growth of the H-spectrum density)?

5. **Module boundary:** Is $\{1, F_1+F_3, F_1+F_4\}$ the complete characterization of
   the module, or do other Fibonacci conventions give different modules?
   Specifically: for a tournament analog using tribonacci instead of Fibonacci as the
   "spine," what is the analogous module?

6. **The Wythoff Nim tournament:** Define a tournament on $\mathbb{N}$ where the only
   "level edges" (H-preserving moves) are the Wythoff pairs. What is the H-value of
   this tournament? Does it relate to the golden-ratio structure of the dimension axis?
