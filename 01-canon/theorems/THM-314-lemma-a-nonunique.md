---
id: THM-314
title: Lemma A (Complete Proof) — Non-Unique Maximum Independent Set Case
status: PROVED
session: opus-2026-05-22-S2
---

## Statement

**Theorem (Lemma A, Non-Unique Case):** Let $T$ be any tournament. Let $\Omega = \Omega(T)$ be its conflict graph on directed odd cycles, and let $d = \alpha(\Omega)$ be the maximum independent set size. Let $S$ be any maximum IS of $\Omega$ (size $d$).

If $S$ is **not** the unique maximum IS (i.e., there exists another max IS $S' \neq S$ of size $d$), then there exists a cycle $C^* \in S$ such that:
$$\deg(I(\Omega \setminus C^*, x)) = d = \deg(I(\Omega - N[C^*], x)) + 1$$

In the Hermite-Biehler notation: $d_A = d$ and $d_B = d - 1$.

## Proof

Since $S \neq S'$ and both have size $d$, the set $S \setminus S'$ is non-empty (otherwise $S \subseteq S'$ and since $|S|=|S'|=d$: $S = S'$, contradiction).

**Pick any $C^* \in S \setminus S'$.**

**Step 1 ($d_A = d$):**
$C^* \notin S'$, so $S'$ is a max IS of size $d$ contained in $\Omega \setminus C^*$. Therefore:
$$d_A = \alpha(\Omega \setminus C^*) \geq |S'| = d.$$
Since $d = \alpha(\Omega)$ is the global maximum: $d_A \leq d$. So $d_A = d$. □

**Step 2 ($d_B = d-1$):**

*Upper bound* (by THM-310, Key Inequality): Any IS $S'' \subseteq \Omega - N[C^*]$ satisfies $\{C^*\} \cup S''$ is an IS of $\Omega$ (since each element of $S''$ is non-adjacent to $C^*$ in $\Omega$). So $|S''| + 1 \leq d$, giving $d_B \leq d - 1$.

*Lower bound*: $S \setminus \{C^*\}$ has size $d - 1$. Since $S$ is an IS, all cycles in $S$ are pairwise vertex-disjoint; in particular, every cycle in $S \setminus \{C^*\}$ is vertex-disjoint from $C^*$. So $S \setminus \{C^*\} \subseteq \Omega - N[C^*]$, giving $d_B \geq d - 1$.

Combining: $d_B = d - 1$. □

**Conclusion:** $d_A = d = (d-1)+1 = d_B + 1$. □

## Significance

This theorem proves **Lemma A** for ALL cases where the maximum IS is non-unique. Combined with:
- The **pair-partner argument** (THM-311) for $d=2$, $\alpha_2 \geq 2$ (where the IS is always non-unique since two distinct max ISes exist).
- **Turán-ULC** for $d=2$, $\alpha_2=1$ (no HB needed).

This leaves only the case $d \geq 3$ with a **unique** maximum IS as open for Lemma A.

## Remark on the Unique IS Case

When $S$ is the **unique** max IS ($d \geq 3$):
- For every $C^* \in S$: $d_A = d - 1$ (removing a cycle from the unique max IS reduces $\alpha$), so $d_A \neq d_B + 1 = (d-1)$.
- For $C^* \notin S$: $d_A = d$, and $d_B \leq d-1$ (Key Inequality). Whether $d_B = d-1$ depends on whether the sub-tournament $T[V \setminus V(C^*)]$ supports enough disjoint cycles. Empirically: 0 failures in 1637 tests at n=7..11.

## Connection to the Full TRRT Proof

Given Lemma A (this theorem) + Lemma B (THM-313, computationally verified) + Hermite-Biehler:

For any tournament $T$ with **non-unique** max IS:
- Choose $C^*$ from this theorem: $d_A = d$, $d_B = d-1$.
- $A = I(\Omega \setminus C^*, x)$ is real-rooted by induction (fewer cycles).
- $B = I(\Omega - N[C^*], x)$ is real-rooted by induction (fewer cycles).
- $B$ interlaces $A$ by Lemma B (THM-313, verified 3672+ times, 0 failures).
- By Hermite-Biehler: $I(\Omega, x) = A + xB$ is real-rooted. □

## See Also

- THM-310: Key Inequality (proves $d_B \leq d-1$)
- THM-311: Lemma A for $d=2$, $\alpha_2 \geq 2$ (pair-partner construction; this theorem subsumes it)
- THM-312: TRRT for $n \leq 8$ via elementary proof
- THM-313: Algebraic identity for Lemma B
- oracle-2026-05-19-S1: Turán-ULC (unconditional, handles $d=2$, $\alpha_2=1$)
