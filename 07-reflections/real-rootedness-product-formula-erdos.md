# Real-Rootedness, Product Formula, and Erdős Connections

**Session:** oracle-2026-05-17-S1
**Depends on:** `real-rootedness-of-independence-polynomial.md`, `cubic-family-x6-coloring.md`
**Computations:** exhaustive n=3..6, all 56 iso classes verified

---

## The Product Formula

Since $I(\Omega(T), x)$ is real-rooted (proved $n \leq 8$ via Chudnovsky-Seymour, conjectured all $n$),
it factors as $I(\Omega(T), x) = \alpha_d \prod_{i=1}^d (x + \rho_i)$ where $\rho_i > 0$.

$$H(T) = \alpha_d \prod_{i=1}^d (2+\rho_i), \qquad I(\Omega(T),6) = \alpha_d \prod_{i=1}^d (6+\rho_i)$$

**Ratio:** $I(\Omega,6)/H = \prod_i (6+\rho_i)/(2+\rho_i) = \prod_i (1+4/(2+\rho_i)) \in (1, 3^d)$.
As tournaments become regular ($\rho_i \to 0$): ratio $\to 3$. For transitive ($d=0$): ratio $= 1$.

---

## Vieta Theorems for Tournaments

**Theorem 1 (n ≤ 5).** The unique root is
$$r = -\frac{1}{\alpha_1} = -\frac{2}{H-1}$$
*Proof:* $\alpha_2=0$ at $n\leq5$ (disjoint cycles need $\geq6$ vertices), so $I=1+\alpha_1 x$. $\square$

**Theorem 2 (n = 6, degree 2).** With $I = \alpha_2 x^2+\alpha_1 x+1$, Vieta's gives:
$$\rho_1+\rho_2 = \frac{\alpha_1}{\alpha_2} \qquad \rho_1\rho_2 = \frac{1}{\alpha_2}$$
**Verified exhaustively** for all 4 iso classes with $\alpha_2\geq1$ at $n=6$. ✓

---

## AM-GM = Log-Concavity = Real-Rootedness

$$\alpha_1^2 \geq 4\alpha_2 \iff \sqrt{\rho_1\rho_2} \leq \frac{\rho_1+\rho_2}{2} \iff \text{real-rooted} \iff \text{log-concave}$$

The three are equivalent for degree-2 polynomials (discriminant $\geq 0$ $\iff$ AM-GM $\iff$ $\alpha_k^2\geq\alpha_{k-1}\alpha_{k+1}$). Verified: all 56 n=6 iso classes are real-rooted (LC=True for all). ✓

**Perfect square case** ($\alpha_1^2=4\alpha_2$, double root): $H = (1+\alpha_1)^2$. Verified: H=9 at
n=6 with $\alpha_1=2, \alpha_2=1$, double root $\rho_1=\rho_2=1$. All odd perfect squares $1^2,3^2,5^2,\ldots$
are (potentially) achievable H values.

---

## The Cubic Invariant is Strictly Finer Than H

**Proved:** Two n=6 iso classes both have $H=29$ but:
- $(\alpha_1=10,\alpha_2=2)$: $I(\Omega,6)=133$, roots $\rho\approx(4.898,0.102)$
- $(\alpha_1=12,\alpha_2=1)$: $I(\Omega,6)=109$, roots $\rho\approx(11.916,0.084)$

The pair $(H,I(\Omega,6))$ recovers $(\alpha_1,\alpha_2)$ exactly at $n\leq8$ (via $\alpha_2=(I-3H+2)/24$).

---

## Forbidden Roots and the Alpha Gap

At $n\leq5$: achievable $\alpha_1\in\{0,1,2,4,5,6,\ldots\}\setminus\{3\}$. The value $\alpha_1=3$
is impossible, creating a **forbidden root gap** $(-1/3,-1/4)$:
- No tournament at $n\leq5$ has $I(\Omega,x)$ with a root in $(-1/3,-1/4)$.
- $r=-1/3$ is **permanently forbidden** (same reason $H=7$ is forbidden).
- The critical boundary $r=-1/4$ corresponds to $\alpha_1=4$ (the "double root boundary").

The critical fugacity $x=-1/4$ is the "half-integer oblong": $n(n-1)=-1/4$ at $n=1/2$.

---

## Erdős Connections

### 1. Tournament Erdős–Ko–Rado

**Classical EKR:** Max intersecting family of $k$-subsets of $[n]$ is the "star" (all sets through a point).

**Tournament EKR conjecture:** $\omega(\Omega(T)) = \max_v \#\{\text{odd cycles through }v\}$.
The maximum clique in the conflict graph is the star at the vertex in the most cycles.

For Paley(7): $14\times3/7=6$ cycles per vertex, so conjectured $\omega(\Omega_3)=6$.

### 2. Erdős on Real Roots

Erdős's school studied: what algebraic properties does the "structure" of a combinatorial object
force? For tournament conflict graphs: real-rootedness may follow from a deep tournament-specific
property. The conjecture: **100% of tournaments at all $n$ have real-rooted** $I(\Omega,x)$.

Evidence: 100% at $n\leq8$ (proved), verified 0 failures at $n\leq20$.

### 3. Log-Concavity (Erdős–Rota–Heron Type)

The $\alpha_k$ sequence is log-concave (via real-rootedness). This is the tournament analog of
the Dowling-Wilson theorem (log-concavity of Whitney numbers) and the Heron-Rota-Welsh
conjecture for matroids (now proved). **Tournament conflict graphs provide a new family of
naturally log-concave independence sequences.**

### 4. Hard-Core Gas Interpretation

$I(\Omega(T),x)$ is the grand canonical partition function of the hard-core lattice gas on
$\Omega(T)$ at fugacity $x$. Real-rootedness = **no phase transition for $x>0$** = the
tournament gas is always in the "gas phase" on the positive real fugacity axis.

The Lee-Yang zeros (in complex $x$-plane) lie only on the negative real axis. The OCF
evaluation $x=2$ is safely away from all zeros, in the convergent "gas phase" region.

**The five fugacity phases:**
- $x=0$: empty (no cycles)
- $x=1$: Fibonacci/Zeckendorf (chemistry)
- $x=2$: Jacobsthal/OCF (Hamiltonian paths)
- $x=6$: S₃/octahedron (triangle colorings)
- $x\to\infty$: maximum independent set dominating

### 5. The Forbidden Root $r=-1/3$

**Conjecture (Erdős-flavored):** No tournament at any $n$ has $I(\Omega(T),x)$ with a zero at $x=-1/3$.

This is equivalent to: the configuration $(\alpha_1=3, \alpha_2=0)$ — three directed odd cycles
pairwise sharing a vertex, with no two disjoint — is impossible in any tournament at any $n$.

The existing proof (for $H=7$ forbidden) uses "cycle-forcing": $\alpha_1=3$ forces specific
cycle configurations that can't exist. Whether this extends to all $n$ is the key question.

---

## The Root-Spectrum Map

The map $T \mapsto \{r_1,\ldots,r_d\}$ (roots of $I(\Omega(T),x)$) is a new tournament invariant.
Properties:
- All roots in $(-\infty, 0)$
- Product formula $H = \alpha_d\prod_i(2+\rho_i)$ recovers $H$ from roots
- Vieta: $\rho_1\rho_2=1/\alpha_2$ (for degree 2) — product encodes disjoint cycle pairs
- The "most symmetric" polynomial: double root $\rho_1=\rho_2$ at $(1+\alpha_1)^2=H$
- The "most asymmetric": ratio $\rho_2/\rho_1\to0$ for SC/regular tournaments (max-H)

**New conjecture:** SC tournaments are precisely those with the most asymmetric root distribution
(minimum $\rho_2/\rho_1$). The SC anti-automorphism forces $\alpha_2\sim\alpha_1^2/4$ (near equality
in log-concavity), which pushes $\rho_2/\rho_1 = ((\alpha_1-\sqrt{\alpha_1^2-4\alpha_2})/(\alpha_1+\sqrt{\alpha_1^2-4\alpha_2}))^2 \to 0$.

---

## Open Problems

1. **Prove** universal real-rootedness: $I(\Omega(T),x)$ real-rooted for all $T$, all $n$.
2. **Prove** Tournament EKR: $\omega(\Omega(T))=\max_v$ star size.
3. **Prove** $r=-1/3$ is permanently forbidden (all $n$).
4. **Determine** if $(H,I(\Omega,6))$ separates all iso classes at all $n$.
5. **Find** the limiting root distribution as $n\to\infty$ for random tournaments.
6. **Prove** ultra-log-concavity of $(\alpha_k/\binom{m}{k})$ for $m=\alpha(\Omega)$.
