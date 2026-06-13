---
id: THM-315
title: Portal-Disjoint Structure of B-B Pairs (HYP-1732 partial result)
status: PROVED (structural), OPEN (full bound)
session: opus-2026-05-22-S2
---

## Setup

Let $T$ be a tournament with $\alpha(\Omega) = 2$. Let $C^*$ be the pair-partner cycle (from THM-311). Define:
- $A = \{a : a \text{ disjoint from } C^*\}$ ($p$ cycles)
- $B = \{b : b \text{ conflicts with } C^*\}$ ($m-1-p$ cycles)

Partition $B$ by **portal set**: $B_S = \{b \in B : V(b) \cap V(C^*) = S\}$ for each non-empty $S \subseteq V(C^*)$.

## Theorem (Portal-Disjoint Structure)

**Theorem:** Two $B$-cycles $b_1, b_2$ can form a $B$-$B$ disjoint pair only if their portal sets are **disjoint**: $V(b_1) \cap V(C^*) \cap V(b_2) = \emptyset$.

**Proof:** If $b_1$ and $b_2$ are vertex-disjoint ($V(b_1) \cap V(b_2) = \emptyset$), then in particular $(V(b_1) \cap V(C^*)) \cap (V(b_2) \cap V(C^*)) \subseteq V(b_1) \cap V(b_2) = \emptyset$. □

## Corollary: No Same-Portal B-B Pairs

**Corollary:** If $b_1, b_2 \in B_S$ (same portal set $S$): they share ALL vertices of $S$ (i.e., $V(b_1) \cap V(b_2) \supseteq S \neq \emptyset$), so they conflict and CANNOT form a disjoint pair.

**Verified:** 0 same-portal B-B pairs in 100 tests at n=9 (p=1 cases).

## Structure of B-B Pairs

$B$-$B$ disjoint pairs arise ONLY between groups $B_{S_1}$ and $B_{S_2}$ with $S_1 \cap S_2 = \emptyset$. This is a **tripartite/multipartite** structure governed by the $2^L$ partition of $V(C^*)$ into subsets.

For $|V(C^*)| = L = 3$ (a 3-cycle as pair-partner):
- Portal sets: singletons $\{u\}, \{v\}, \{w\}$ and pairs $\{u,v\}, \{u,w\}, \{v,w\}$ and $\{u,v,w\}$.
- Disjoint portal pairs: $(\{u\}, \{v\})$, $(\{u\}, \{w\})$, $(\{v\}, \{w\})$, $(\{u\}, \{v,w\})$, $(\{v\}, \{u,w\})$, $(\{w\}, \{u,v\})$, and more.

## Open: Proof of HYP-1732 ($\alpha_2 \leq p(m-p)$)

The portal-disjoint structure constrains B-B pairs but does not immediately give the bound $e_{AB} + e_{BB} \leq p(m-p)$. The full proof remains open.

**Key inequalities proved:**
1. For each B-B pair $(b_1, b_2)$: $e_{AB}(b_1) + e_{AB}(b_2) \leq p$ (from K₃-free structure of $D$).
2. $\Sigma_b e_{AB}(b) \cdot \deg_{D[B]}(b) \leq p \cdot e_{BB}$ (summing over all B-B pairs).
3. Portal disjointness: $B$-$B$ pairs only between groups with disjoint portal sets.

**Verified:** 0 violations of $\alpha_2 \leq p(m-p)$ in 1637 tests at n=7..11.

## See Also

- THM-311: Lemma A for $d=2$ (uses HYP-1732 implicitly via Lemma B)
- THM-313: Algebraic identity (HYP-1732 $\Leftrightarrow$ B interlaces A $\Leftrightarrow$ $\alpha_2 \leq p(m-p)$)
- HYP-1732: The main conjecture
