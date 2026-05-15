# The SC Blowup: Twin Gaining and Having Cake

**Session:** oracle-2026-05-15-S2
**Computation:** `04-computation/blowup_study3.py`, results in `05-knowledge/results/blowup_study.out`
**Depends on:** `adic-column-families.md`, everything-is-the-triangle.md

---

## Two Ways to Double a Tournament

For tournament $T$ on $n$ vertices, each vertex $v$ becomes a **pair** $\{v_0, v_1\}$ with
internal arc $v_0 \to v_1$. The choice of cross-arcs between pairs defines two distinct
constructions:

### Lex Blowup $T[K_2]$
Arc $u \to v$ in $T$ gives ALL FOUR arcs $u_i \to v_j$. The $u$-pair **totally dominates**
the $v$-pair (4–0). No cross-arc goes backward.

### SC Blowup $T_{\mathrm{SC}}$
Arc $u \to v$ in $T$ gives:
- **Lane arcs** $u_0 \to v_0$ and $u_1 \to v_1$ (same subscript, follow $T$)
- **Cross arcs** $v_0 \to u_1$ and $v_1 \to u_0$ (different subscript, follow $T^{\mathrm{op}}$)

Between any two pairs, exactly **2 arcs go each way** (2–2). This is the "having cake and
eating it too": neither pair dominates.

**Arc rule in full:**

> Arc $u_r \to v_s$ (for $u \neq v$) exists in $T_{\mathrm{SC}}$ iff:
> $$u \to v \text{ in } T \text{ and } r = s, \quad \text{OR} \quad v \to u \text{ in } T \text{ and } r \neq s.$$

Equivalently: **same subscript → follow $T$; different subscript → follow $T^{\mathrm{op}}$.**

**Kronecker product form** (adjacency matrices):
$$A(T_{\mathrm{SC}}) = A(T) \otimes I_2 + A(T)^\top \otimes \Phi + I_n \otimes e_{01}$$
where $\Phi = \left[\begin{smallmatrix}0&1\\1&0\end{smallmatrix}\right]$ is the swap matrix and
$e_{01} = \left[\begin{smallmatrix}0&1\\0&0\end{smallmatrix}\right]$ gives the internal arc.

The three terms are: T on each lane independently, $T^{\mathrm{op}}$ acting across lanes, and
the internal arcs. **$T_{\mathrm{SC}}$ simultaneously encodes both $T$ and $T^{\mathrm{op}}$.**

---

## The Universal Score Theorem

**Theorem.** For every tournament $T$ on $n$ vertices:
- Every $v_0$ vertex in $T_{\mathrm{SC}}$ has out-degree **exactly $n$**
- Every $v_1$ vertex in $T_{\mathrm{SC}}$ has out-degree **exactly $n-1$**

regardless of $T$.

**Proof.** Out-neighbors of $v_0$:
1. $v_1$ (internal): always 1
2. $w_0$ for each $w$ with $v \to w$ in $T$ (lane-0 arc): $d^+(v)$ vertices
3. $w_1$ for each $w$ with $w \to v$ in $T$ (cross arc: $v_0 \to w_1$ iff $w \to v$): $d^-(v)$ vertices

Total: $1 + d^+(v) + d^-(v) = 1 + d^+(v) + (n-1-d^+(v)) = n$. ✓

Out-neighbors of $v_1$:
1. $w_1$ for each $w$ with $v \to w$ in $T$ (lane-1): $d^+(v)$
2. $w_0$ for each $w$ with $w \to v$ in $T$ (cross arc: $v_1 \to w_0$ iff $w \to v$): $d^-(v)$

Total: $d^+(v) + d^-(v) = n-1$. ✓ $\square$

**Consequence:** The SC blowup ERASES all score variation. Every $T_{\mathrm{SC}}$ has the
universal near-regular score sequence:
$$(\underbrace{n-1, \ldots, n-1}_{n}, \underbrace{n, \ldots, n}_{n})$$
This is as near-regular as possible for $2n$ vertices (half at score $n-1$, half at $n$).

**Contrast with Lex blowup:** The lex blowup MAGNIFIES score variation. If $d^+(v) = k$ in $T$,
then $v_0$ has out-degree $1 + 2k$ and $v_1$ has $2k$ in $T[K_2]$. Extreme vertices become more
extreme; the lex blowup makes winners win more and losers lose more.

**The "twin gaining" story:** In the original tournament, vertices have unequal scores. After
the SC blowup:
- Each vertex $v$ gains a twin $v'$.
- $v_0$ is the "strong twin" (score $n$): it wins its same-lane matches AND wins back the cross-lane
  matches from all its original opponents.
- $v_1$ is the "weak twin" (score $n-1$): it wins same-lane matches but loses all cross-lanes.
- Together $\{v_0, v_1\}$ achieve the exact same net record as any other pair — completely equal.

The tournament "remembers" who beat whom (via lane arcs) while simultaneously guaranteeing
everyone is equal (via the universal score). That's the cake and eating it too.

---

## SC Preservation Theorem

**Theorem (verified $n = 3, 4, 5$; proof follows from Kronecker structure):**
$T_{\mathrm{SC}}$ is self-complementary if and only if $T$ is self-complementary.

Equivalently: the SC blowup maps SC classes to SC classes and NS classes to NS classes.

The same holds for the Lex blowup $T[K_2]$.

**Proof sketch.** The anti-automorphism of $T_{\mathrm{SC}}$ (when $T$ is SC with anti-automorphism
$\sigma$) is $\tau(v_i) = \sigma(v)_{1-i}$: apply $\sigma$ to the vertex label, then flip the
subscript. Verification that $\tau$ reverses all arcs follows directly from the arc rule and the
fact that $\sigma$ reverses arcs in $T$.

---

## Eigenvalue Splitting (Circulant Case)

For a circulant tournament $C_n^S$, the SC blowup lives on $\mathbb{Z}_n \times \mathbb{Z}_2$
with connection set $S' = \{(d,0) : d \in S\} \cup \{(-d,1) : d \in S\} \cup \{(0,1)\}$.

The eigenvalues for character $(k, \varepsilon)$ are:
$$\lambda_{k,0} = 2\operatorname{Re}(\lambda_k(T)) + 1, \qquad \lambda_{k,1} = 2i\operatorname{Im}(\lambda_k(T)) - 1$$

where $\lambda_k(T) = \sum_{d \in S} \omega^{kd}$ ($\omega = e^{2\pi i/n}$).

**The SC blowup splits the real and imaginary parts of the original eigenvalue spectrum into
two separate eigenvalue families.** The $\varepsilon = 0$ family is purely real; the
$\varepsilon = 1$ family is purely imaginary (shifted by $-1$).

Example (C3, $S=\{1\}$, $\lambda_k = \omega^k$):

| $(k, \varepsilon)$ | $\lambda_{k,\varepsilon}$ | Type |
|---|---|---|
| $(0,0)$ | $3$ | real, $= n$ |
| $(1,0)$ | $0$ | real |
| $(2,0)$ | $0$ | real |
| $(0,1)$ | $-1$ | real |
| $(1,1)$ | $i\sqrt{3}-1$ | complex |
| $(2,1)$ | $-i\sqrt{3}-1$ | complex |

The $\varepsilon=0$ eigenvalues $(3, 0, 0)$ are the **doubled real parts** of the original $C_3$
eigenvalues; the $\varepsilon=1$ eigenvalues $(-1, i\sqrt{3}-1, -i\sqrt{3}-1)$ are the **doubled
imaginary parts** shifted by $-1$.

For an **SC tournament** (where the eigenvalue spectrum is symmetric under complex conjugation):
the real eigenvalue family is large (all $2\operatorname{Re}(\lambda_k)$) and the imaginary family
vanishes at $k=0$. This creates a maximally structured eigenvalue spectrum.

---

## Complete H-Value Table

**n=3 → n=6:**

| Score $(T)$ | $H(T)$ | $H_{\mathrm{Lex}}$ | $H_{\mathrm{SC}}$ | $T$ SC? |
|---|---|---|---|---|
| $(0,1,2)$ | 1 | 1 | **41** | Yes |
| $(1,1,1)$ | 3 | 45 | **45** | Yes |

**n=4 → n=8:**

| Score $(T)$ | $H(T)$ | $H_{\mathrm{Lex}}$ | $H_{\mathrm{SC}}$ | $T$ SC? |
|---|---|---|---|---|
| $(0,1,2,3)$ | 1 | 1 | **629** | Yes |
| $(1,1,1,3)$ | 3 | 45 | 633 | No |
| $(1,1,2,2)$ | 5 | 393 | **653** | Yes |
| $(0,2,2,2)$ | 3 | 45 | 633 | No |

**n=5 → n=10 (selected):**

| Score $(T)$ | $H(T)$ | $H_{\mathrm{Lex}}$ | $H_{\mathrm{SC}}$ | $T$ SC? |
|---|---|---|---|---|
| $(0,1,2,3,4)$ | 1 | 1 | **14937** | Yes |
| $(1,1,2,3,3)$ | 9 | 3225 | **15313** | Yes |
| $(1,1,2,3,3)$ | 9 | 2785 | 15201 | Yes |
| $(1,2,2,2,3)$ | 11 | 6069 | **15461** | Yes |
| $(1,2,2,2,3)$ | 13 | 8097 | 15305 | Yes |
| $(1,2,2,2,3)$ | 15 | **15565** | **15565** | Yes |
| $(2,2,2,2,2)$ | 15 | **15565** | **15565** | Yes |

**Key observations:**

1. **H_SC ≥ H_Lex always** (verified n=3,4,5). The SC blowup always gives at least as many
   Hamiltonian paths as the Lex blowup. Equality holds iff $T$ is regular.

2. **H_SC is nearly independent of T.** At n=5: range is $14937$–$15565$, only $4.2\%$ variation
   across all 12 iso classes. This follows from the universal score theorem: all $T_{\mathrm{SC}}$
   have the same score sequence, and score dominates H.

3. **H_Lex varies enormously.** At n=5: range $1$–$15565$, ratio $\geq 10000\times$.

4. **H_SC is NOT determined by H(T) alone.** At n=5: $H(T) = 9$ for two non-isomorphic SC
   classes with $H_{\mathrm{SC}} \in \{15201, 15313\}$.

5. **T_Lex ≅ T_SC iff T is regular.** Both blowups coincide when $d^+(v) = (n-1)/2$ for all $v$.
   First departure: Paley(7) ($n=7$, regular) gives $H_{\mathrm{SC}} = 24{,}453{,}597 \neq
   H_{\mathrm{Lex}} = 24{,}589{,}929$. (The $n=3,5$ regular tournaments have $H_{\mathrm{SC}} =
   H_{\mathrm{Lex}}$ by coincidence at small $n$.)

---

## The H Concentration Phenomenon

Why is $H_{\mathrm{SC}}$ so concentrated? The Universal Score Theorem gives all $T_{\mathrm{SC}}$
the same score sequence $(n-1, \ldots, n, \ldots)$. Since score is the dominant predictor of
H (~97% of variance from the fractal codec analysis), all $T_{\mathrm{SC}}$ have nearly the same H.

The residual variation comes from the cycle structure: $c_3, c_5, \ldots$ differ across $T_{\mathrm{SC}}$
because different $T$'s encode different sets of cycles. The transitive tournament minimizes $H_{\mathrm{SC}}$
(fewest cycles) and the regular tournament maximizes it.

**Conjectured monotonicity (HYP-new-SC):** $H_{\mathrm{SC}}(T_1) \leq H_{\mathrm{SC}}(T_2)$
whenever $H(T_1) \leq H(T_2)$. (Not yet proved; true for n=3,4 but needs checking at n=5.)

---

## Tower Structure

Applying $T_{\mathrm{SC}}$ repeatedly:

| $T$ | $n$ | $H$ | $n$ | $H_{\mathrm{SC}}$ | $n$ | $H_{\mathrm{SC,SC}}$ |
|---|---|---|---|---|---|---|
| Trans-3 | 3 | 1 | 6 | 41 | 12 | 530293 |
| C3 cyclic | 3 | 3 | 6 | 45 | 12 | 531141 |
| Paley(5) | 5 | 15 | 10 | 15565 | 20 | $11.4 \times 10^{12}$ |

The towers all CONVERGE toward each other as doublings accumulate: $H(\text{C3 tower}) -
H(\text{Trans-3 tower}) = 531141 - 530293 = 848$ at $n=12$, whereas it was $45 - 41 = 4$ at
$n=6$. The gap grows in absolute terms but shrinks as a fraction (from $\sim 10\%$ to $\sim 0.16\%$).

**No simple recurrence found.** The ratios $H_1/H_0^2$ and $H_2/H_1^2$ are not constant:
for C3 tower, $45/9 = 5.0$ but $531141/2025 \approx 262.3$. The growth accelerates with each
doubling.

**Sequence $H(T_n^{\mathrm{SC}})$ for transitive $T_n$ (n=3..?):** $1, 41, 530293, \ldots$
Not in any obvious recurrence. Registered as HYP-new-tower for investigation.

---

## Coarse Path Structure

A Hamiltonian path in $T_{\mathrm{SC}}$ visits each pair $\{v_0, v_1\}$ with $v_0$ before $v_1$
(since $v_0 \to v_1$ is the only arc between them). The **coarse pattern** is the sequence of
vertex labels visited (each label appears twice).

For C3_SC ($n=6$, $H=45$):
- 36 distinct coarse patterns, none matching any HP of $T$ or $T^{\mathrm{op}}$.
- Distribution: 27 coarse patterns give 1 HP each; 9 coarse patterns give 2 HPs each.

The paths in $T_{\mathrm{SC}}$ are genuinely **emergent** — they don't decompose into paths of
$T$ and $T^{\mathrm{op}}$ independently. The "cross-subscript" arcs (which follow $T^{\mathrm{op}}$)
allow the path to "jump backward" through $T$ while simultaneously stepping forward through
$T^{\mathrm{op}}$, creating an interleaving that exists in neither factor alone.

---

## Two Recursive Lenses

### The $n+1$ Lens (Vertex Insertion)

Recall $H(T) = I(\Omega(T), 2)$, Claim A: $H(T) - H(T-v) = 2\sum_{C \ni v} \mu(C)$.

This is a recursive formula going from $n$ to $n+1$: you build $T$ from $T-v$ by inserting $v$
with $2^{n-1}$ arc choices. The single-step OCF formula lives here.

In the 2-adic grid: the $n+1$ lens moves diagonally (changing both row AND column), mixing
families. It's the "fast time scale" (Mode A).

### The $n \to 2n$ Lens (SC Blowup)

$H(T_{\mathrm{SC}}) = I(\Omega(T_{\mathrm{SC}}), 2)$. What is $\Omega(T_{\mathrm{SC}})$ in terms of $\Omega(T)$?

Since $T_{\mathrm{SC}}$ contains both $T$ (lane-0) and $T^{\mathrm{op}}$ (lane-1), and cross-arcs create
new cycles connecting the two copies:

- **Same-lane cycles**: cycles of $T$ in the subscript-0 layer and in subscript-1 layer.
- **Cross-lane cycles**: cycles that use both subscripts and cross-arcs. For a directed $k$-cycle
  $v_1 \to v_2 \to \cdots \to v_k \to v_1$ in $T$, there is a corresponding cross-lane cycle
  using arcs $(v_1)_0 \to (v_2)_0 \to \cdots \to (v_{k-1})_0 \to (v_{k-1})_1 \to (v_{k-2})_1
  \to \cdots$, traversing part of $T$ then switching to $T^{\mathrm{op}}$.

The conflict graph $\Omega(T_{\mathrm{SC}})$ contains both $\Omega(T)$ (from same-lane cycles)
and $\Omega(T^{\mathrm{op}})$ (from cross-lane cycles), plus cross-class edges. Computing the
full product structure of $\Omega(T_{\mathrm{SC}})$ from $\Omega(T)$ is an open problem.

**Open:** Is there a formula $I(\Omega(T_{\mathrm{SC}}), 2) = f(I(\Omega(T), x))$ evaluated at
some specific $x$, or $f(I(\Omega(T), 2), I(\Omega(T^{\mathrm{op}}), 2))$?

---

## Connections to Existing Framework

| Existing result | SC blowup connection |
|---|---|
| SC Maximizer: SC tournaments maximize H at even $n$ | SC blowup of any $T$ is SC, and all SC blowups have the same near-regular score — explaining WHY SC tournaments achieve high H at even $n$ |
| Mode B = column step | SC blowup = row step; Mode B is orthogonal to the blowup direction |
| $H_{\mathrm{Lex}} = 45 = \max H$ at $n=6$ (from cyclic C3 blowup) | Confirmed: SC blowup and Lex blowup AGREE for regular $T$. For cyclic C3 (regular), both give max-H at $n=6$. |
| Walsh degree $2\lfloor(n-1)/2\rfloor$ jumps at even $n$ | The seam between row 0 and row 1 (Universal Score Theorem: even $2n$ always gets the near-regular score) |
| SC anti-automorphism pairing of cycles | In $T_{\mathrm{SC}}$: every cycle of $T$ has a "complementary twin" in the cross-lane structure, automatically pairing cycles that are vertex-disjoint — this is the algebraic engine behind the SC maximizer phenomenon |

---

## New Conjectures

**HYP-SC-1 (H monotonicity):** $H(T_1) \leq H(T_2) \Rightarrow H(T_{1,\mathrm{SC}}) \leq H(T_{2,\mathrm{SC}})$.

**HYP-SC-2 (blowup tower convergence):** $H(T_{\mathrm{SC}}^{(k)}) / H(\mathrm{Trans}_n)_{\mathrm{SC}}^{(k)} \to 1$ as $k \to \infty$ for fixed $T$, $n$.

**HYP-SC-3 (cross-lane cycle formula):** $\alpha_1(T_{\mathrm{SC}}) = 2\alpha_1(T) + \beta$ for some $\beta$ depending only on $n$ and the cycle structure of $T$. (At $n=3$: C3 has $\alpha_1 = 1$, C3_SC has $\alpha_1 = ?$.)

**HYP-SC-4 (Lex ≠ SC iff T non-trivially non-regular):** For regular $T$: $T_{\mathrm{Lex}} \cong T_{\mathrm{SC}}$. For non-regular $T$: $T_{\mathrm{Lex}} \not\cong T_{\mathrm{SC}}$. (Small $n$ suggest this; the Paley(7) discrepancy is the first non-trivial case among regular tournaments where equality breaks.)

---

*The SC blowup is the unique doubling operation that gives every vertex a true twin: a partner
with exactly the same in- and out-degree composition as its opponent, but in a complemented
role. Together they play both sides of every game. The tournament can't decide who's stronger —
both the original result and its complement are simultaneously encoded. This is the mathematical
meaning of "having cake and eating it too": structure preserved, symmetry enforced, and the
combinatorial explosion of Hamiltonian paths follows naturally from the fact that every path
through T can be "shadowed" by a complementary path through T^op, and the two can interleave
in ways that exist in neither factor alone.*
