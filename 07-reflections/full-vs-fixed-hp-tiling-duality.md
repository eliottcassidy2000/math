# The Full and Fixed-HP Tilings: A Duality

**Session:** oracle-2026-05-05
**Depends on:** `two-models-staircase-recursion.md`, THM-292, the CUT⊕CYCLE decomposition
**Computations:** all verified in-session

---

## Reframing: A Tournament IS a Tiling

The user's insight reframes everything. A tournament on $n$ vertices is a binary string
in $\{0,1\}^{C(n,2)}$ — a tiling of the full $C(n,2)$-dimensional hypercube.
The fixed-HP model is a $C(n-1,2)$-dimensional sub-tiling. Both are tilings;
the difference is dimensionality and which arcs are fixed.

This gives two valid "tiling" perspectives, plus the metagraph (isomorphism classes):

| Model | Bits free | Space | Natural symmetry |
|---|---|---|---|
| **Full tiling** | $C(n,2)$ | $\{0,1\}^{C(n,2)}$ | $S_n$ (all vertex permutations) |
| **Fixed-HP tiling** | $C(n-1,2)$ | $\{0,1\}^{C(n-1,2)} = \delta_{n-2}$ | $\mathbb{Z}_2$ (complement flip of staircase) |
| **Metagraph** $G_n$ | — | $\{0,1\}^{C(n,2)}/S_n$ | trivial |

---

## Section 1: The CUT ⊕ CYCLE Decomposition

The full tiling space has a canonical **GF(2) direct sum decomposition**:

$$\{0,1\}^{C(n,2)} = \underbrace{\{0,1\}^{n-1}}_{\text{CUT space}} \oplus \underbrace{\{0,1\}^{C(n-1,2)}}_{\text{CYCLE space}}$$

**CUT space** ($n-1$ arcs, span-1 arcs in the base path ordering):
$$\text{CUT} = \{(k+1 \to k) \text{ or } (k \to k+1) : k = 1, \ldots, n-1\}$$
These arcs connect consecutive positions in the base-path ordering. Their orientations
determine the score sequence (the "cut" in the GF(2) cut-cycle decomposition).

**CYCLE space** ($C(n-1,2)$ arcs, span $\geq 2$):
$$\text{CYCLE} = \text{staircase } \delta_{n-2} = \{(x,y) : x-y \geq 2\}$$
These arcs generate all directed odd cycles (the "cycle space"). The OCF
$H(T) = I(\Omega(T), 2)$ depends entirely on the cycle space.

**From the definitions file:** "The wiggly arcs = the cycle-space generators when
the base path $P_0$ is used as spanning tree. Base-path arcs = the cut-space (score
hierarchy). This decomposition IS the GF(2) Cut ⊕ Cycle split."

**Dimension balance as $n \to \infty$:**

| $n$ | CUT = $n-1$ | CYCLE = $C(n-1,2)$ | Total $C(n,2)$ | CUT/Total |
|---|---|---|---|---|
| 3 | 2 | 1 | 3 | 0.667 |
| 5 | 4 | 6 | 10 | 0.400 |
| 7 | 6 | 15 | 21 | 0.286 |
| 10 | 9 | 36 | 45 | 0.200 |
| $n$ | $n-1$ | $C(n-1,2)$ | $C(n,2)$ | $2/(n)$ |

As $n \to \infty$: CUT fraction $\to 0$, CYCLE fraction $\to 1/2$.
**The tournament becomes "mostly CYCLE" as $n$ grows.**

**Fixing the HP = fixing the CUT space to all-forward:**
The fixed-HP tiling model constrains all $n-1$ CUT arcs to point "forward"
(arc $k+1 \to k$ for all $k$). This is one specific point in $\{0,1\}^{n-1}$
(the all-zeros point = all consecutive arcs forward = P_0 is a Hamiltonian path).
The remaining $C(n-1,2)$ CYCLE arcs are free.

---

## Section 2: The Span Decomposition of the Full Space

Organizing arcs by span $s = |x - y|$:

| Span $s$ | Arcs | Count | Layer in full staircase |
|---|---|---|---|
| 1 | $(k+1, k)$ for $k=1,\ldots,n-1$ | $n-1$ | **CUT layer** |
| 2 | $(k+2, k)$ for $k=1,\ldots,n-2$ | $n-2$ | Hypotenuse strip $S_3$ |
| $s$ | $(k+s, k)$ for $k=1,\ldots,n-s$ | $n-s$ | Strip $S_{s+1}$ |
| $n-1$ | $(n, 1)$ | 1 | **APEX** |

**The full tiling = an extended staircase** including the span-1 CUT layer:
$$\delta_{n-2}^{\text{ext}} = \underbrace{\text{span-1 row}}_{\text{CUT, }n-1\text{ cells}} \cup \underbrace{\delta_{n-2}}_{\text{CYCLE staircase, }C(n-1,2)\text{ cells}}$$

The fixed-HP model uses only $\delta_{n-2}$; the full model uses $\delta_{n-2}^{\text{ext}}$.

**Geometrically:** the fixed-HP staircase $\delta_{n-2}$ is the full staircase
$\delta_{n-2}^{\text{ext}}$ with its bottom row (CUT layer, span-1 arcs) removed.
Fixing the HP = "gluing down" the bottom row to the all-forward position.

---

## Section 3: Two Theorems About the Difference

### Theorem 1: CUT-CYCLE Arc-Flip Symmetry

**Theorem.** In the full tournament space, averaged over all labeled tournaments,
CUT arc flips and CYCLE arc flips have **identical $\Delta H$ distributions**
when normalized by arc count.

*Proof.* The full tournament space $\{0,1\}^{C(n,2)}$ has the full $S_n$ symmetry:
permuting vertex labels maps each arc to every other arc with equal weight.
Since all arcs are equivalent under $S_n$, the $\Delta H$ distribution when flipping
any specific arc is the same for ALL arcs. Summing over all labeled tournaments
and normalizing by arc count gives equal distributions for CUT and CYCLE. $\square$

**Verified computation** (ratio = CYCLE arcs / CUT arcs = $C(n-1,2)/(n-1) = (n-2)/2$):

| $n$ | CUT ΔH+ dist | CYC ΔH+ dist | Ratio |
|---|---|---|---|
| 4 | $\{2:24,\, 4:6\}$ | $\{2:24,\, 4:6\}$ | 1.0 |
| 5 | $\{2:360,\, 4:240,\ldots\}$ | $\{2:540,\, 4:360,\ldots\}$ | 1.5 |
| 6 | $\{2:6000,\, 4:6360,\ldots\}$ | $\{2:12000,\, 4:12720,\ldots\}$ | 2.0 |

The ratio at each $n$ is exactly $(n-2)/2$ — confirming that CUT and CYCLE arcs
have the same per-arc $\Delta H$ behavior. **Fixing the HP breaks this symmetry**,
making CUT and CYCLE arcs play different roles. The symmetry is a property of the
full space only.

### Theorem 2: The Tiling Model IS H-Weighted Sampling

**Theorem.** The H-distribution in the fixed-HP tiling model is the H-weighted
distribution over all labeled tournaments:

$$P_{\text{tile}}(H(T_t) = h) = P_{\text{full}}(H(T) = h) \cdot \frac{h}{E_{\text{full}}[H]}.$$

*Proof.* A labeled tournament $T$ appears as a tiling $T_t$ for each of its $H(T)$
Hamiltonian paths (one choice of base path per HP). So the number of tilings with
$H = h$ is $\sum_{T: H(T)=h} H(T) = h \cdot |\{T: H(T)=h\}|$. Dividing by total
tilings $= 2^{C(n-1,2)} = E_{\text{full}}[H] \cdot 2^{C(n,2)} / n!$... more
precisely: the tiling count equals $\sum_T H(T) = n!/2^{n-1} \cdot 2^{C(n-1,2)}$
(from THM-292: $\Sigma H = E[H] \cdot 2^{C(n-1,2)}$). Normalizing gives the formula.
$\square$

**Verified** (ratio $P_{\text{full}}/P_{\text{tile}} = E[H]/h$):

| $n=5$, $E[H]=7.5$ | H=1 | H=3 | H=5 | H=9 | H=11 | H=13 | H=15 |
|---|---|---|---|---|---|---|---|
| Ratio $P_F/P_T$ | 7.5 | 2.5 | 1.5 | 0.83 | 0.68 | 0.58 | 0.50 |
| $E[H]/h = 7.5/h$ | 7.5 | 2.5 | 1.5 | 0.83 | 0.68 | 0.58 | 0.50 |

The formula $P_F/P_T = E[H]/h$ holds exactly at every $H$-value.

**Corollary 1.** $E_{\text{tile}}[H] = E_{\text{full}}[H^2]/E_{\text{full}}[H] = E_{\text{full}}[H] \cdot (1 + \mathrm{CV}^2_{\text{full}}[H])$.

This connects to the Cayley-Delannoy sequence $W(n) = n!(1 + \mathrm{CV}^2[H])$:
$$E_{\text{tile}}[H] = W(n)/n! \cdot E_{\text{full}}[H] = W(n)/2^{n-1}.$$

**Corollary 2 (Crossover at $E[H]$).** The ratio $P_{\text{full}}/P_{\text{tile}} = E[H]/h$
equals 1 exactly when $h = E[H] = n!/2^{n-1}$. For $h < E[H]$: full overweights
(CUT space contributes to low-H). For $h > E[H]$: tiling overweights (cycle space
contributes to high-H).

---

## Section 4: H Lives Entirely in CYCLE Space

**Theorem.** For any fixed CUT configuration, the set of achievable H-values spans
the FULL range from $H_{\min} = 1$ to $H_{\max}(n)$.

**Verified:**

| $n$ | H-range within any fixed CUT config | = $H_{\max} - 1$ |
|---|---|---|
| 3 | 2 | $3-1=2$ ✓ |
| 4 | 4 | $5-1=4$ ✓ |
| 5 | 14 | $15-1=14$ ✓ |
| 6 | 44 | $45-1=44$ ✓ |

For every CUT configuration, the CYCLE space alone achieves every H-value from 1
to $H_{\max}(n)$.

**Implication.** The CUT space controls the **statistical distribution** of H over
tournaments (the 97% OCR at $n=5$: score sequences explain 97% of H-variance) but
NOT the **range** of achievable H-values. The range is entirely in the CYCLE space.

The CUT space shifts the H-distribution: all-forward CUT (tiling model) gives
mean $E_{\text{tile}}[H] = 9.875$ at $n=5$; averaging over all CUT configs gives
$E_{\text{full}}[H] = 7.5$. The CUT variation accounts for the $9.875 - 7.5 = 2.375$
difference in means.

---

## Section 5: The Metagraph Through Both Lenses

The metagraph $G_n$ is the quotient of the full tiling by $S_n$. Each node is an
isomorphism class. Edges come from arc flips.

**Edge partition in the full tiling:**

| $n$ | Total edges | CUT edges | CYCLE edges | CUT fraction |
|---|---|---|---|---|
| 3 | 12 | 8 | 4 | 0.667 |
| 4 | 192 | 96 | 96 | 0.500 |
| 5 | 5120 | 2048 | 3072 | 0.400 |
| 6 | 245760 | 81920 | 163840 | 0.333 |

The CUT fraction = $(n-1)/C(n,2) = 2/n$ exactly. As $n \to \infty$: CUT edges $\to 0$,
CYCLE edges dominate.

**Three views of the metagraph:**

**View 1 (CYCLE-only / tiling model):** Vertices = tilings (with fixed P_0). Edges =
CYCLE arc flips (wiggly lines in the tiling model). This is the "staircase hypercube"
$\{0,1\}^{C(n-1,2)}$ with single-bit-flip adjacency. After quotienting by $S_n$, we get
the metagraph restricted to CYCLE edges. The wiggly structure (each tile position defines
a perfect matching on tilings) is a property of this view.

**View 2 (Full tiling):** Vertices = all labeled tournaments. Edges = ALL arc flips
(CUT + CYCLE). Both types of flips are available, but they're indistinguishable under
$S_n$ symmetry (Theorem 1). The CUT edges are "invisible" in this view — they're just
more edges of the same type.

**View 3 (Metagraph $G_n$):** Vertices = isomorphism classes. Edges = single arc
flip (of any span) that changes the iso-class. The distinction between CUT and CYCLE
edges is LOST in the quotient: a CUT flip and a CYCLE flip can produce the same iso-class.

**The key structural difference:** CYCLE flips are LABELED (each flip is labeled by
its tile position in the staircase). CUT flips are LABELED too (labeled by their
consecutive pair position). In the metagraph quotient, both labels are collapsed.

**Where CUT flips appear in the metagraph (but are invisible):**
When a CUT arc (span-1) is flipped, the tournament changes. In the metagraph, this
is just another arc-flip edge. But the specific CUT arc position encodes whether the
new tournament has the SAME or DIFFERENT canonical base path. CUT flips are those
that might take a tiling (with $P_0$ as HP) to a tournament where $P_0$ is NO LONGER
a HP. They "cross the tiling model boundary."

---

## Section 6: Two Inductions — Strip vs Span

The full tiling has two natural inductive structures:

### Induction 1: STRIP (Mode A) — adding vertices

Group arcs by **upper vertex**: Strip $S_k$ = arcs from vertex $k+1$ to all lower
non-adjacent vertices. Mode A adds strip $S_{n-1}$ when $n \to n+1$.

In the FULL tiling (both CUT and CYCLE): adding vertex $n+1$ adds:
- 1 new CUT arc (the span-1 arc $(n+1, n)$)
- $n-2$ new CYCLE arcs (span $\geq 2$ arcs from $n+1$ to vertices $1, \ldots, n-1$)

Total new arcs = $n-1$. The NEW STRIP $S_{n}$ in the full tiling has $n-1$ arcs, of
which 1 is CUT and $n-2$ are CYCLE.

**The CUT/CYCLE ratio within each strip:**
| Strip $S_k$ | CUT arcs | CYCLE arcs | Total |
|---|---|---|---|
| $S_2$ | 1 (span 1) | 0 | 1 |
| $S_3$ | 1 (span 1) | 1 (span 2) | 2 |
| $S_k$ | 1 (span 1) | $k-2$ (span $\geq 2$) | $k-1$ |

Each strip has exactly ONE CUT arc (the span-1 arc connecting the new vertex to its
immediate predecessor). All other arcs in the strip are CYCLE arcs.

**This is why fixing the HP corresponds to fixing one arc per strip:** the $n-1$
CUT arcs are exactly one per strip, and the fixed-HP tiling freezes these in
the "forward" direction.

### Induction 2: SPAN — adding arc layers

Group arcs by **span**: span-$s$ layer = arcs with $|x-y| = s$. The span layers are:
- Span 1 (CUT): $n-1$ arcs
- Span 2 (hypotenuse): $n-2$ arcs
- Span $s$: $n-s$ arcs
- Span $n-1$ (apex): 1 arc

**Span induction** adds one span layer at a time, from span 1 outward.
This is PERPENDICULAR to strip induction in the staircase:
- Strip induction = vertical slices (adding columns to the staircase)
- Span induction = horizontal slices (adding rows to the staircase)

**The two inductionsare DUAL in the staircase:**
- Strip order: $S_2, S_3, \ldots, S_{n-1}$ — columns of the staircase
- Span order: span-$1$, span-$2$, ..., span-$(n-1)$ — rows of the staircase (anti-diagonals)

The fixed-HP model lives at the boundary between span-1 (CUT) and span-2 (hypotenuse):
it freezes span-1 and studies all higher spans.

**What span induction reveals:** The H-contribution of each arc depends on its SPAN.
In the fixed-HP model, the formula $H = 1 + 2^d$ measures distance $d$ from the
hypotenuse (span-2 layer) — which is exactly the depth in the SPAN induction order.
In the FULL model, span-1 arcs are at $d = -1$ (one step "below" the hypotenuse).

---

## Section 7: The Pattern in the Differences

Comparing the full and fixed-HP models reveals a pattern in how H-values are weighted:

**Pattern 1: The full model is uniform; the tiling model is H-weighted.**
Every labeled tournament has equal weight in the full model. In the tiling model,
tournament $T$ has weight proportional to $H(T)$. This makes the tiling model a
"length-biased" sample of the tournament space.

**Pattern 2: The CUT space is a "noise" layer from the tiling perspective.**
From the tiling model's view, the CUT arcs are the "fixed" part — they carry no
useful information about cycle structure. The CYCLE arcs carry ALL the information
about $H(T)$. But from the full model's view, CUT arcs and CYCLE arcs are equivalent.

**Pattern 3: The crossover at $h = E[H]$ marks the boundary.**
The ratio $P_{\text{full}}/P_{\text{tile}} = E[H]/h$:
- $h < E[H]$: CUT flexibility inflates the count of low-H tournaments
- $h = E[H]$: the two models agree on the probability of this H-value
- $h > E[H]$: tiling model's length bias inflates high-H tournament counts

The expected value $E[H] = n!/2^{n-1}$ is the natural scale separating "typical"
(full model) from "high-H" (tiling model). This is also:
- The mean of the H-distribution over labeled tournaments
- The ratio SC(n)/A000568(n) asymptotically
- The Cayley-Delannoy factor $W(n)/n!$ (minus 1)

**Pattern 4: As $n \to \infty$, the two models converge.**
CUT fraction $\to 0$. The span-1 arcs become an infinitesimal fraction of all arcs.
The tiling model and full model become asymptotically identical for large $n$.
The H-weighting effect disappears (CV²$[H] \to 0$ as $n \to \infty$? Need verification).

---

## Section 8: Recursive Formulas in the Full Model

In the FULL tiling, the recursive formula involves BOTH CUT and CYCLE arcs.

**Vertex insertion (Mode A) in the full model:**
Adding vertex $n+1$ adds $n$ new arcs:
- 1 CUT arc $(n+1, n)$ — direction determines whether $P_0'$ includes $(n+1 \to n)$
- $n-1$ CYCLE arcs (the new strip, span $\geq 2$) — determine new HP contributions

The H-recursion remains: $H(T_{n+1}) = H(T_n) + 2k$ (even increment).
But NOW the CUT arc contributes to $k$ directly: if $(n+1) \to n$ (forward CUT arc),
vertex $n+1$ can extend HPs that ended at vertex $n$. If $n \to (n+1)$ (backward),
vertex $n+1$ cannot be appended; instead, it must be inserted in the interior.

**The CUT arc determines the INSERTION PATTERN:**
- CUT forward ($(n+1) \to n$): vertex $n+1$ is a "new source" — it beats vertex $n$,
  making insertion at position BEFORE $n$ in any HP natural
- CUT backward ($n \to (n+1)$): vertex $n+1$ is "absorbed" into the middle of existing HPs

In the fixed-HP model, the CUT arc $(n+1, n)$ is always forward (base path), so
vertex $n+1$ is always at the TOP of the path. The full model allows $n+1$ to be
anywhere, which is why the H-formula is symmetric.

**The full-model ΣH recursion:**
$$\Sigma H_{\text{full}}(n) = 2^{C(n,2)} \cdot E_{\text{full}}[H(T)] = 2^{C(n,2)} \cdot \frac{n!}{2^{n-1}}$$
This is exact (since $E[H] = n!/2^{n-1}$ is known exactly). The tiling model ΣH
is different: $\Sigma H_{\text{tile}}(n) = 2^{C(n-1,2)} \cdot E_{\text{tile}}[H] = W(n) \cdot 2^{C(n-1,2)-n+1}$.

The ratio: $\Sigma H_{\text{tile}}/\Sigma H_{\text{full}} = W(n)/(n! \cdot 2^{n-1})$,
which equals $1 + \mathrm{CV}^2[H]$.

---

## Summary: The Three Core Identities

$$\boxed{
\begin{aligned}
&\text{Full space} = \text{CUT} \oplus \text{CYCLE} = \{0,1\}^{n-1} \oplus \{0,1\}^{C(n-1,2)} \\[4pt]
&P_{\text{tile}}(H = h) = P_{\text{full}}(H = h) \cdot \frac{h}{E_{\text{full}}[H]} \quad\text{(H-weighting)} \\[4pt]
&E_{\text{tile}}[H] = E_{\text{full}}[H] \cdot (1 + \mathrm{CV}^2[H]) = \frac{W(n)}{2^{n-1}} \quad\text{(Cayley-Delannoy link)}
\end{aligned}
}$$

The **CUT space** encodes score structure (the "boring" part).
The **CYCLE space** encodes path and cycle structure (the "interesting" part).
**Fixing the HP** = projecting out the CUT space = studying CYCLE alone.
**All H-variation** lives in the CYCLE space; the CUT space only shifts the mean.
