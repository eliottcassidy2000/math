---
id: THM-313
title: Algebraic Identity Reducing Lemma B to Root Containment of I(Omega,x)
status: PROVED (identity); Lemma B itself computationally verified
session: opus-2026-05-21-S1
---

## Statement

**Theorem (Key Algebraic Identity for Lemma B):** In the Hermite-Biehler decomposition
$$I(\Omega, x) = A(x) + x \cdot B(x)$$
where $A = I(\Omega \setminus C^*, x)$ and $B = I(\Omega - N[C^*], x)$, and $B$ has degree 1 with unique root $x_B = -1/p$ (where $p = \alpha_1^B$):

$$A\!\left(-\tfrac{1}{p}\right) = I\!\left(-\tfrac{1}{p}\right)$$

**Corollary:** $B$ interlaces $A$ (the Hermite-Biehler interlacing condition) iff $I(\Omega, -1/p) \leq 0$.

## Proof of Identity

Since $B(x) = 1 + px$ (degree 1 with leading coefficient $p$):
$$B\!\left(-\tfrac{1}{p}\right) = 1 + p\cdot\!\left(-\tfrac{1}{p}\right) = 1 - 1 = 0.$$

Substituting $x = -1/p$ in $I = A + xB$:
$$I\!\left(-\tfrac{1}{p}\right) = A\!\left(-\tfrac{1}{p}\right) + \left(-\tfrac{1}{p}\right) \cdot \underbrace{B\!\left(-\tfrac{1}{p}\right)}_{=0} = A\!\left(-\tfrac{1}{p}\right). \quad \square$$

## Geometric Interpretation

For $d = 2$ (the case in THM-311 and THM-312): $I(\Omega, x) = 1 + mx + \alpha_2 x^2$ is real-rooted with roots $-1/\rho_1$ and $-1/\rho_2$ (with $\rho_1 \geq \rho_2 > 0$ by Turán-ULC).

The interlacing condition $A(-1/p) \leq 0$ is equivalent to $I(-1/p) \leq 0$, which holds iff:
$$\rho_2 \leq p \leq \rho_1.$$

**In words:** $B$ interlaces $A$ iff the count $p$ (number of directed odd cycles vertex-disjoint from $C^*$) lies **between the two characteristic scales** $\rho_1, \rho_2$ of the polynomial $I(\Omega(T), x)$.

The condition $I(-1/p) \leq 0$ is purely algebraic:
$$1 - \frac{m}{p} + \frac{\alpha_2}{p^2} \leq 0 \quad \Longleftrightarrow \quad \alpha_2 \leq p(m-p).$$

## Special Cases

**$d = 2$, $\alpha_2 = 1$ (Turán base case):** $I = 1 + mx + x^2$. Roots: $\rho_{1,2} = (m \pm \sqrt{m^2-4})/2$. For $p=1$: $I(-1) = 1-m+1 = 2-m \leq 0$ iff $m \geq 2$. Always true. ✓

**Turán extremal ($\alpha_2 = m^2/4$, $K_{m/2,m/2}$ structure):** $\rho_1 = \rho_2 = 2/m$ (double root). The pair-partner C* has $p = m/2$ disjoint cycles. $I(-2/m) = 1 - 2 + 1 = 0$. Equality — $p$ is exactly the double root. ✓

## Equivalent Combinatorial Form

For $d = 2$: Lemma B (for C* with $p = \alpha_1^B$ disjoint cycles) holds iff:
$$\alpha_2(\Omega(T)) \leq p \cdot (m - p)$$

This is a purely combinatorial inequality about the tournament's cycle structure.

**Computationally verified:** $0$ violations in $109$ tests ($n = 6, 7, 8$, $d = 2$, $\alpha_2 \geq 2$).

## Significance for TRRT

This identity means: to prove Lemma B for the pair-partner construction, it suffices to prove $\alpha_2(T) \leq p(m-p)$, a clean inequality about the tournament's cycle counts.

Combined with THM-311 (Lemma A, which finds C* with $d_A = d_B+1 = 2$), a proof of $\alpha_2 \leq p(m-p)$ would complete the proof of TRRT for all $n \leq 8$ via the elementary Hermite-Biehler route (without Chudnovsky-Seymour).

## See Also

- THM-311: Lemma A for $d=2$ (pair-partner construction)
- THM-312: TRRT for $n \leq 8$
- oracle-2026-05-19-S1: Turán-ULC (proves $\alpha_2 \leq m^2/4$, a weaker bound)
- HYP-1729: TRRT conjecture
