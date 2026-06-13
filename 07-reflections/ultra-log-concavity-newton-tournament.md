# Ultra-Log-Concavity via Newton's Inequalities: Tournaments as Matroids

**Session:** oracle-2026-05-17-S1
**Depends on:** `real-rootedness-product-formula-erdos.md`, `root-spectrum-n6-computations.md`
**Status:** Theorem (proved conditional on real-rootedness conjecture)

---

## The Central Theorem

**Theorem (Tournament ULC).** If the independence polynomial $I(\Omega(T),x) = \sum_{k=0}^d \alpha_k x^k$ of the conflict graph of tournament $T$ is real-rooted, then the normalized coefficient sequence is log-concave:
$$\left(\frac{\alpha_k}{\binom{d}{k}}\right)^2 \geq \frac{\alpha_{k-1}}{\binom{d}{k-1}} \cdot \frac{\alpha_{k+1}}{\binom{d}{k+1}}, \quad 1 \leq k \leq d-1$$

where $d = \alpha(\Omega(T))$ is the independence number of the conflict graph.

**Proof.** Write $I(\Omega,x) = \alpha_d \prod_{i=1}^d (x + \rho_i)$ with $\rho_i > 0$. Then $\alpha_k = \alpha_d \cdot e_{d-k}(\rho_1,\ldots,\rho_d)$ where $e_j$ denotes the $j$-th elementary symmetric polynomial. The Newton-Maclaurin inequality for positive reals states:
$$\left(\frac{e_j}{\binom{d}{j}}\right)^2 \geq \frac{e_{j-1}}{\binom{d}{j-1}} \cdot \frac{e_{j+1}}{\binom{d}{j+1}}$$
for all $0 \leq j \leq d$. Substituting $j = d-k$ and multiplying by $\alpha_d^2/\alpha_d^2$ gives the result. $\square$

**Computational verification:** Zero violations in exhaustive n=6 and random samples n=7,8,9.

---

## Comparison with the Matroid Theorem

| | Tournament Conflict Graphs | Matroids |
|---|---|---|
| Object | $\Omega(T)$: independence polynomial of cycles | $M$: characteristic polynomial |
| Coefficients | $\alpha_k$ = # size-$k$ indep. cycle collections | $W_k$ = Whitney numbers |
| Log-concavity | $\alpha_k^2 \geq \alpha_{k-1}\alpha_{k+1}$ | $W_k^2 \geq W_{k-1}W_{k+1}$ |
| Ultra-log-concavity | $(\alpha_k/\binom{d}{k})^2 \geq \cdots$ | $(\tilde{W}_k)^2 \geq \cdots$ |
| Proof via | Real-rootedness → Newton's ineqs | Hodge theory (AHK 2018) |
| Key tool | Polynomial with positive real roots | Kähler manifold geometry |
| Status (ULC) | Proved iff TRRT holds | Proved unconditionally |

**The crucial distinction:** The matroid proof (Adiprasito-Huh-Katz) requires algebraic geometry and is far deeper. The tournament proof reduces to Newton's classical inequality — much simpler — but is conditioned on the real-rootedness conjecture.

**If TRRT is proved**, the tournament ULC becomes a clean corollary. The two results would then inhabit the same landscape with different proof strategies: geometry vs. polynomial analysis.

---

## The Three-Level Hierarchy of Log-Concavity

For the tournament independence polynomial with degree $d$:

### Level 1: Log-concavity (from real-rootedness)
$$\alpha_k^2 \geq \alpha_{k-1} \cdot \alpha_{k+1} \quad \text{for all } 1 \leq k \leq d-1$$

### Level 2: Ultra-log-concavity (Newton-Maclaurin)
$$\frac{\alpha_k^2}{\binom{d}{k}^2} \geq \frac{\alpha_{k-1}}{\binom{d}{k-1}} \cdot \frac{\alpha_{k+1}}{\binom{d}{k+1}} \quad \iff \quad \alpha_k^2 \cdot k \cdot (d-k) \geq \alpha_{k-1} \cdot \alpha_{k+1} \cdot (k+1)(d-k+1)$$

### Level 3: Ultra-ultra-log-concavity (Brändén conjecture?)
Define $\beta_k = \alpha_k/\binom{d}{k}$ (the "normalized" sequence). Is $(\beta_k)$ itself ULC?

The Level 2 already says $(\beta_k)$ is log-concave. Level 3 would say $(\beta_k/\binom{d}{k})$ is log-concave, i.e., $(\alpha_k/\binom{d}{k}^2)$ is log-concave:
$$\left(\frac{\alpha_k}{\binom{d}{k}^2}\right)^2 \geq \frac{\alpha_{k-1}}{\binom{d}{k-1}^2} \cdot \frac{\alpha_{k+1}}{\binom{d}{k+1}^2}$$

This would follow if $(\beta_k)$ itself has a real-rooted generating polynomial. Is $\sum \beta_k x^k$ real-rooted? Unknown.

---

## At n=6: The Degree-2 Case

For degree-2 ($d=2$), ULC at $k=1$ is just:
$$\frac{\alpha_1^2}{4} \geq \alpha_0 \cdot \alpha_2 = \alpha_2 \quad \iff \quad \alpha_1^2 \geq 4\alpha_2$$

This is the **discriminant condition** for $I = 1 + \alpha_1 x + \alpha_2 x^2 \geq 0$ on $\mathbb{R}$. It is equivalent to real-rootedness (discriminant ≥ 0).

**Verification:** All 47 n=6 iso classes (including 9 SC classes with $\alpha_2 \geq 1$) satisfy $\alpha_1^2 \geq 4\alpha_2$:
- H=9, α₁=2, α₂=1: $4 \geq 4$ ✓ (equality = double root)
- H=17, α₁=6, α₂=1: $36 \geq 4$ ✓
- H=29, α₁=10, α₂=2: $100 \geq 8$ ✓
- H=45, α₁=20, α₂=1: $400 \geq 4$ ✓
- H=45, α₁=14, α₂=4: $196 \geq 16$ ✓

The only equality case (double root) is H=9, α₁=2, α₂=1. All others have strict inequality.

---

## Mean Convergence: The Geometric Interpretation

By AM-GM, ULC implies:
$$\frac{e_k}{\binom{d}{k}} \geq \left(\frac{e_d}{\binom{d}{d}}\right)^{k/d} \cdot \left(\frac{e_0}{\binom{d}{0}}\right)^{(d-k)/d} = e_d^{k/d}$$

where $e_d = \prod_i \rho_i$ (geometric mean raised to power $d$). This says the normalized elementary symmetric means converge monotonically from $e_0=1$ (top) to the geometric mean (bottom).

For our case $\alpha_0 = 1$, $\alpha_d = \alpha_d$ (leading coefficient):
$$\frac{\alpha_k}{\binom{d}{k}} \geq \sqrt[d]{\alpha_d} \cdot \left(\sqrt[d]{\frac{1}{\alpha_d}} \cdot \sqrt[d]{\alpha_d}\right)^k = \text{interpolation between 1 and } \alpha_d/\alpha_d^{k/d}$$

The sequence $(\alpha_k/\binom{d}{k})$ is not only log-concave but **convex on a log scale**, meaning it forms a concave sequence in the geometric mean sense.

---

## The Hard-Core Gas Interpretation

In statistical physics, $I(\Omega, x)$ is the grand canonical partition function at fugacity $x$. The roots $-\rho_i$ are the **Lee-Yang zeros**. ULC in this context says:

**The free energy $f(x) = \log I(\Omega, x) / d$ is a concave function of $\log x$** when $x > 0$.

This is because for a polynomial with all real negative roots:
$$f(x) = \frac{1}{d}\sum_i \log(x + \rho_i) = \text{average of concave log-functions}$$

The second derivative $f''(\log x) = -\sum_i x^2/(x+\rho_i)^2/d < 0$. So the free energy is always concave in $\log x$ — no "convex dip" that would signal a phase transition.

**Physical meaning:** The hard-core tournament gas at fugacity $x$ has a concave free energy on the positive real fugacity axis. This is the thermodynamic stability condition — no phase transition for $x > 0$.

---

## Why Tournaments Have ULC But General Graphs Don't

For general graphs $G$, $I(G, x)$ need not be real-rooted (e.g., the Petersen graph has complex roots). Tournament conflict graphs have special structure that (conjecturally) forces real-rootedness:

1. **Perfect graph structure of $\Omega(T)$?** If $\Omega(T)$ is always perfect, then $I(\Omega,x)$ is a "clique polynomial" which... may not force real-rootedness directly.

2. **The interlacing structure:** The deletion-contraction recursion $I(\Omega,x) = I(\Omega \setminus C, x) + x \cdot I(\Omega - N[C], x)$ has a tournament-specific form where the two polynomials interlace. This interlacing forces real-rootedness inductively.

3. **The arc-flip symmetry:** Reversing an arc in $T$ changes $\Omega(T)$ in a controlled way. Real-rootedness might be preserved under arc flips, making it an orbit invariant.

---

## Open Problems

1. **Is $\sum \beta_k x^k$ real-rooted?** (Level 3 ULC). If yes, the entire infinite hierarchy is log-concave.

2. **What is the right normalization?** Newton-Maclaurin uses $\binom{d}{k}$ from the "ambient" dimension $d$. But $d$ depends on the tournament. Is there a canonical normalization that gives a clean interpolation?

3. **Does ULC help prove real-rootedness?** Normally ULC follows from RR, but could it also be used to prove RR inductively? E.g., if $I(\Omega \setminus C, x)$ and $I(\Omega - N[C], x)$ are both ULC, does ULC of their sum plus the interlacing give RR of the sum?

4. **The tournament Newton polytope.** The coefficient vector $(\alpha_0,\ldots,\alpha_d)$ lives in a polytope $P_{n,d}$. Real-rootedness is the condition that this vector lies in the "real-rootedness cone." What is the geometry of $P_{n,d} \cap \text{RRC}$?
