---
id: THM-310
title: Key Inequality for Lemma A (TRRT)
status: PROVED
session: opus-2026-05-21-S1
---

## Statement

**Theorem (Key Inequality):** Let G be any graph with independence number d = alpha(G). Let S be any maximum independent set of G (|S| = d). Let C* be any vertex of G NOT in S (C* ∉ S). Then:

$$\text{alpha}(G - N[C^*]) \leq d - 1$$

where G - N[C*] = induced subgraph on vertices NOT in the closed neighborhood of C* (i.e., excluding C* itself and all its neighbors in G).

## Proof

Let S' be any independent set in G - N[C*]. By definition, every vertex of S' is not adjacent to C* in G (since S' avoids N[C*]).

Therefore {C*} ∪ S' is an independent set of G of size |S'| + 1.

Since d = alpha(G) is the maximum, |S'| + 1 ≤ d, giving |S'| ≤ d - 1.

Since S' was arbitrary: alpha(G - N[C*]) ≤ d - 1. □

## Corollary (For Deletion-Contraction)

For any tournament T and C* ∉ S (max IS of Omega(T)):
- **d_A = d**: S ⊆ Omega(T) \ C* (since C* ∉ S), so alpha(Omega \ C*) ≥ d; combined with alpha ≤ d: d_A = d.
- **d_B ≤ d-1**: by the theorem above.

Together: deg(I(Omega\C*)) ≥ deg(I(Omega-N[C*])) + 1.

## Significance

This is the **lower half of Lemma A**: the degree drop in the deletion-contraction I(Omega,x) = I(Omega\C*,x) + x · I(Omega-N[C*],x) is at least 1 when C* ∉ S.

Combined with THM-311 (Lemma A for d=2), this gives a complete proof of TRRT for all n ≤ 8.

## See Also

- THM-311: Lemma A for d=2, α₂≥2 (full existence of good C*)
- THM-312: TRRT for n≤8 via elementary proof
- HYP-1729: TRRT conjecture
