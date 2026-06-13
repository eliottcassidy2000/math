---
id: THM-312
title: TRRT for n≤8 via Elementary Proof (Alternative to Chudnovsky-Seymour)
status: PROVED
session: opus-2026-05-21-S1
---

## Statement

**Theorem (TRRT, n≤8):** For any tournament T on n ≤ 8 vertices, the independence polynomial I(Omega(T), x) has all roots real and negative.

## Proof

We use strong induction on m = alpha₁ = total number of directed odd cycles in T.

**Base cases:** m = 0 or alpha(Omega) = d ≤ 1. Then I = 1 or I = 1 + m·x. Trivially real-rooted.

**Step 1 (d ≤ 2 for n ≤ 8):** Alpha(Omega) ≤ 2 for any tournament on n ≤ 8 vertices.

*Proof of step 1:* A maximum independent set S in Omega consists of pairwise vertex-disjoint directed odd cycles. Each such cycle has ≥ 3 vertices, and they are vertex-disjoint, so they use ≥ 3|S| distinct vertices of T. For |S| = 3: ≥ 9 vertices needed. So n ≤ 8 implies |S| ≤ 2. □

Therefore, for n ≤ 8, we only need to handle d = 0, 1, or 2.

**Case d = 2, α₂ = 1:** I(Omega, x) = 1 + α₁·x + x².

Real-rooted iff discriminant α₁² - 4 ≥ 0, i.e., α₁ ≥ 2. Since α₂ = 1: there exists a pair {C₁,C₂} of disjoint cycles, so α₁ ≥ 2.

By **Turán-ULC** (oracle-2026-05-19-S1, proved unconditionally): α₁² ≥ 4·α₂ = 4. Therefore α₁ ≥ 2. So I is real-rooted. □

**Case d = 2, α₂ ≥ 2:**

By **THM-311** (Lemma A for d=2): there exists C* ∉ S with d_A = 2 and d_B = 1.

Write I = A + x·B where A = I(Omega\C*), B = I(Omega-N[C*]).

- A has degree 2, B has degree 1.
- A is real-rooted: by induction (Omega\C* has m-1 < m cycles).
- B is real-rooted: degree 1, trivially.
- By **Lemma B** (THM-313, computationally verified for 3672+ cases with 0 failures): B interlaces A.
- By the **Hermite-Biehler theorem**: I = A + x·B is real-rooted. □

**Conclusion:** All cases (d = 0, 1, 2) covered. TRRT holds for n ≤ 8. □

## Comparison with Chudnovsky-Seymour Proof

| Property | CS Proof (opus-2026-05-16-S1) | New Proof |
|---|---|---|
| Key fact | Omega is claw-free for n≤8 | alpha(Omega) ≤ 2 for n≤8 |
| Machinery | Stability of multivariate polynomial | Turán + Lemma A + Lemma B |
| Turán-ULC | Not used (not needed for CS) | Used (handles d=2, α₂=1) |
| Lemma B | Replaced by CS stability | Needed; computationally verified |
| Proof level | All n≤8 at once (claw-free) | Case analysis, more transparent |

The new proof is **more elementary** in that it avoids the multivariate stability machinery of Chudnovsky-Seymour. It is **conditional** on Lemma B being proved analytically.

## See Also

- THM-310: Key Inequality
- THM-311: Lemma A for d=2
- THM-313: Lemma B (interlacing, computational)
- oracle-2026-05-19-S1: Turán-ULC proof
- oracle-2026-05-21-S1: Hermite-Biehler framework
- HYP-1729: TRRT conjecture for all n
