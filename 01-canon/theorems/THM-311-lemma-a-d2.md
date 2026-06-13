---
id: THM-311
title: Lemma A for d=2, α₂≥2 (TRRT Inductive Step)
status: PROVED
session: opus-2026-05-21-S1
---

## Statement

**Theorem (Lemma A for d=2):** Let T be any tournament with alpha(Omega(T)) = 2 and α₂(Omega(T)) ≥ 2 (at least two disjoint pairs of vertex-disjoint directed odd cycles). Let S = {C₁, C₂} be any maximum independent set of Omega(T).

Then there exists a cycle C* ∉ S such that:
$$\text{deg}(I(\Omega \setminus C^*, x)) = \text{deg}(I(\Omega - N[C^*], x)) + 1 = 2$$

In the notation of the Hermite-Biehler strategy: d_A = 2 = d_B + 1 = 1 + 1.

## Proof

Since α₂ ≥ 2, there exist at least two distinct pairs of vertex-disjoint cycles. Let {C₃, C₄} be a maximum IS (disjoint pair) distinct from {C₁, C₂}.

Since {C₃, C₄} ≠ {C₁, C₂} as sets, the set {C₃, C₄} \ {C₁, C₂} is non-empty. Pick C* to be any element of {C₃, C₄} \ {C₁, C₂}.

Let C** be the other element of {C₃, C₄} (the "pair-partner" of C*). Then:

**Step 1 (d_A = 2):** C* ∉ S = {C₁, C₂}. Therefore {C₁, C₂} ⊆ Omega \ C*. So alpha(Omega \ C*) ≥ 2. Since d = 2 is the maximum, d_A = alpha(Omega \ C*) = 2.

**Step 2 (d_B ≥ 1):** The pair-partner C** is disjoint from C* (they form a disjoint pair). Therefore C** is NOT in N[C*] (since N[C*] contains only C* and its neighbors in Omega = cycles CONFLICTING with C*). So C** ∈ Omega - N[C*]. Thus alpha(Omega - N[C*]) ≥ 1.

**Step 3 (d_B ≤ 1):** By THM-310 (Key Inequality): alpha(Omega - N[C*]) ≤ d - 1 = 1.

**Step 4:** d_B = alpha(Omega - N[C*]) = 1 = d_A - 1 = 2 - 1. □

## Proof of Existence of {C₃,C₄}\{C₁,C₂} ≠ ∅

Since {C₃,C₄} ≠ {C₁,C₂} as sets (they are distinct pairs) and each is a 2-element set, |{C₃,C₄} \ {C₁,C₂}| ≥ 1. Concretely: if {C₃,C₄} = {C₁,C₄} (shares one element): then C₄ ≠ C₁ (since C₁ ∈ {C₃,C₄}\{C₁,C₂} would mean... actually C₃=C₁ and C₄ ≠ C₂ since the pairs are distinct). Take C* = C₄ ∉ {C₁,C₂}. Its partner C₃ = C₁ ∈ Omega-N[C₄] (since C₃=C₁ and C₄ are a disjoint pair: C₁ disjoint from C₄). ✓

## Verified Computationally

n=7: 33 tests (d=2, α₂≥2), 33 passes, 0 failures.

## Application (Hermite-Biehler)

With C* given by this theorem:
- A(x) = I(Omega \ C*, x): degree 2, real-rooted by induction.
- B(x) = I(Omega - N[C*], x): degree 1, trivially real-rooted.
- By Lemma B (THM-313, computationally verified): B interlaces A.
- By Hermite-Biehler: I(Omega, x) = A(x) + x·B(x) is real-rooted.

## See Also

- THM-310: Key Inequality (proves d_B ≤ d-1)
- THM-312: TRRT for n≤8 (uses this as the main d=2 step)
- THM-313: Lemma B (interlacing, computationally verified)
- oracle-2026-05-19-S1: Turán-ULC (proved unconditionally, covers d=2, α₂=1)
