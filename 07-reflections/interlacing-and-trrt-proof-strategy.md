# Interlacing Property and the TRRT Proof Strategy

**Session:** oracle-2026-05-19-S1
**Source:** `interlacing_investigation.py`
**Status:** Observation (computationally verified, proof open)

---

## The Key Observation

**Computational finding (n=6, stride=16 sampling):**

When deleting a single cycle $C^*$ from the conflict graph $\Omega(T)$ reduces the polynomial degree by exactly 1 (i.e., $\deg(I(\Omega \setminus C^*)) = \deg(I(\Omega)) - 1$), the polynomial $I(\Omega \setminus C^*, x)$ **interlaces** $I(\Omega, x)$.

**Result: 444 interlacing cases tested, 0 failures.** ✓

(There are 4468 "skip" cases where deleting $C^*$ changes the degree by more than 1.)

---

## Interlacing Definition

$Q(x)$ **interlaces** $P(x)$ (both real-rooted, $\deg P = \deg Q + 1$) if their roots alternate:
$$\rho_1 \geq \sigma_1 \geq \rho_2 \geq \sigma_2 \geq \cdots \geq \sigma_{d-1} \geq \rho_d$$
where $\{\rho_i\}$ are the roots of $P$ (all positive, since we use $-x$ convention) and $\{\sigma_i\}$ the roots of $Q$.

**Hermite-Biehler theorem:** $P = A + xB$ is real-rooted with alternating roots iff $A$ and $B$ are both real-rooted and $B$ interlaces $A$.

---

## The Deletion-Contraction Recursion

For any cycle $C^*$ in the conflict graph $\Omega$:
$$I(\Omega, x) = \underbrace{I(\Omega \setminus C^*, x)}_{A(x)} + x \cdot \underbrace{I(\Omega - N[C^*], x)}_{B(x)}$$

where:
- $A(x) = I(\Omega \setminus C^*, x)$ — delete $C^*$, keep all other cycles
- $B(x) = I(\Omega - N[C^*], x)$ — delete $C^*$ and all cycles that conflict with it (keep only cycles vertex-disjoint from $C^*$)

**If Hermite-Biehler applies:** $I(\Omega, x) = A(x) + x \cdot B(x)$ is real-rooted iff $A$ and $B$ are both real-rooted AND $B$ interlaces $A$ (i.e., $A$ interlaces $I(\Omega,x)$).

---

## Proof Strategy for TRRT

**Inductive proof sketch (conditional on interlacing):**

*Base case:* If $\Omega$ has no cycles ($\alpha_1 = 0$), then $I(\Omega, x) = 1$ — trivially real-rooted.

*Inductive step:* Given $T$ with $\alpha_1 = m$ cycles. Choose any cycle $C^*$.
1. $A(x) = I(\Omega \setminus C^*, x)$: conflict graph of all OTHER cycles. By induction (fewer cycles), $A$ is real-rooted.
2. $B(x) = I(\Omega - N[C^*], x)$: conflict graph of cycles DISJOINT from $C^*$. This is also a tournament conflict graph (on the same tournament, restricted to cycles not involving $C^*$'s vertices). By induction, $B$ is real-rooted.
3. **Key claim:** $B$ interlaces $A$.

If claim 3 is proved, then by Hermite-Biehler, $I(\Omega, x) = A(x) + xB(x)$ is real-rooted. $\square$

**What's missing:** The proof of claim 3 ("$B$ interlaces $A$"). Computationally verified but not proved.

---

## Why Interlacing Is Plausible

For **matching polynomials** (Heilmann-Lieb 1972): the matching polynomial satisfies the interlacing property, which implies real-rootedness. The proof uses the fact that deleting a vertex from a graph gives a well-controlled polynomial relationship.

For **claw-free graph independence polynomials** (Chudnovsky-Seymour 2007): the multivariate extension is stable (HPP), which implies the univariate interlacing and hence real-rootedness.

For **tournament conflict graphs**: the specific structure of $\Omega(T)$ — coming from the cyclic structure of a tournament — might force a stronger form of the interlacing than what holds for general graphs.

**The "tournament Ramsey" constraint:** Not all graphs can be conflict graphs of tournaments. The tournament structure prevents certain "bad" configurations (like claw in $\bar\Omega$ when $d=3$). This same constraint might force the interlacing.

---

## The Degree-Drop = 1 vs > 1 Cases

The 4468 "skip" cases (degree drops by more than 1) are the interesting ones. These occur when $C^*$ is the **only cycle in some maximal independent set** — deleting it reduces $d$ by more than 1.

For these cases, standard interlacing doesn't apply in the same way. But the polynomial is still real-rooted (by TRRT, computationally). This suggests a more nuanced interlacing structure involving "generalized interlacing" or "partial interlacing."

**Observation:** Even in the skip cases, $A(x) = I(\Omega \setminus C^*, x)$ is real-rooted. The question is whether $xB(x)$ is also real-rooted (it is, by induction) and whether their "combined" polynomial $A + xB = I(\Omega,x)$ is real-rooted.

---

## The Heilmann-Lieb Analogy in Full

For **matching polynomials**: the recursion $\mu(G,x) = \mu(G-v, x) - x \cdot \mu(G-N[v], x)$ has a SIGN change (minus, not plus). This makes the proof simpler: the interlacing follows from the structure of the recursion.

For **independence polynomials**: $I(G,x) = I(G-v,x) + x \cdot I(G-N[v], x)$ (plus). The positive sign makes interlacing harder to prove in general but doesn't rule it out.

For **claw-free graphs**: the plus-sign recursion still works because the neighborhood $N[v]$ has special structure (no claw = the neighbors form a union of cliques).

For **tournament conflict graphs**: the neighborhood $N[C^*]$ = all cycles sharing a vertex with $C^*$. This set has special structure coming from the tournament: cycles sharing a vertex form "fans" around each shared vertex.

---

## The Claw-Freeness Boundary

From our computation (with the claw-detection code, though it timed out before completing):
- At **n ≤ 8**: $\Omega(T)$ is claw-free (need 9 vertices for a claw). Interlacing follows from Chudnovsky-Seymour.
- At **n = 9**: $\Omega(T)$ can have claws (91% of random n=9 tournaments have claws). But still real-rooted (0 complex roots in 50 samples).

This means interlacing at n ≥ 9 (if it holds) must come from tournament-specific structure, not just claw-freeness.

---

## Connection to the Turán-ULC Proof

The Turán-ULC proof (k=1, unconditional) uses:
$$\text{max clique in } \bar\Omega = d \implies \text{Turán density bound} \implies \text{ULC at k=1}$$

The interlacing approach uses:
$$\text{tournament structure of } \Omega \implies \text{interlacing of } I(\Omega\setminus C^*, x) \text{ with } I(\Omega, x) \implies \text{TRRT}$$

Both exploit the "tournament Ramsey" structure — the fact that conflict graphs of tournaments have special properties. The Turán approach extracts the independence number; the interlacing approach would extract a dynamical property.

**If both are proved:** Turán-ULC + Interlacing-TRRT would give a complete, combinatorial proof of both ULC and real-rootedness, without algebraic geometry.

---

## Open Problems

1. **Prove interlacing when deg-drop = 1:** For any tournament T and any cycle $C^*$ where deleting $C^*$ from $\Omega$ reduces the degree by exactly 1, does $I(\Omega \setminus C^*, x)$ interlace $I(\Omega, x)$?

2. **Handle the deg-drop > 1 cases:** Find the right interlacing structure for cases where the degree drops by more than 1.

3. **Prove TRRT by induction:** Once interlacing is established (for both cases), the inductive proof of TRRT would be complete.

4. **Multivariate stability:** Is $\hat{I}(\Omega(T), \mathbf{x}) = \sum_{S \text{ indep}} \prod_{c \in S} x_c$ a stable polynomial for every tournament $T$? This would give interlacing as a corollary.

5. **Connection to matroid theory:** Is $\Omega(T)$ (or $\Delta\Omega(T)$) a matroid or gammoid for all $T$? Matroid independence polynomials are stable (via Choe-Oxley-Sokal-Wagner 2004), which would imply TRRT.
