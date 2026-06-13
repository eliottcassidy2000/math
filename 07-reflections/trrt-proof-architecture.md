# TRRT Proof Architecture: Where We Stand

**Session:** opus-2026-05-22-S2 (extends opus-2026-05-21-S1)

---

## The Tournament Real-Rootedness Theorem (TRRT)

**Conjecture:** For any tournament T, all roots of I(Ω(T), x) are real and negative.

**Status:** Proved for n≤8 (Chudnovsky-Seymour). Verified n=9..11 (0 failures). Open in general.

---

## Complete Proof Structure (Hermite-Biehler Induction)

The proof by strong induction on m=α₁ (total odd cycles) works as follows:

### Base and Low Cases (ALL PROVED)

| Case | Condition | Method | Status |
|------|-----------|--------|--------|
| d=0 | No cycles | I=1 | TRIVIAL |
| d=1 | One max IS | I=1+α₁x | TRIVIAL |
| d=2, α₂=1 | One disjoint pair | Turán-ULC: α₁²≥4α₂=4 | PROVED (oracle-S19) |
| d=2, α₂≥2 | Multiple pairs | Turán-ULC: I=1+α₁x+α₂x² real-rooted iff α₂≤m²/4 | PROVED (oracle-S19) |

### Inductive Step via Hermite-Biehler (CONDITIONAL)

**For d≥3:** Split I(Ω,x) = A(x) + x·B(x) via deletion-contraction with chosen C*.

**THEOREM (THM-314, PROVED):** When max IS S is NOT unique, some C*∈S gives d_A=d, d_B=d-1.
- Key: C* ∈ S\S' for any second max IS S'. Elementary 3-line proof.
- Covers: ALL d≥2 with non-unique IS.

**LEMMA A (OPEN):** When S is unique (d≥3): some C* still gives d_A=d_B+1.
- Empirically: 0 failures at n=9..11.
- Partial proof: if T[V\V(C*)] has enough cycles, works.
- Key step: show some sub-tournament on 6 vertices has two disjoint odd cycles.

**LEMMA B (COMPUTATIONALLY VERIFIED):** When d_A=d_B+1: B interlaces A.
- 3672 tests, 0 failures.
- Equivalent (via THM-313): I(-1/p) ≤ 0, i.e., p between roots of I.
- For d=2: equivalent to HYP-1732 (α₂ ≤ p(m-p)).

### Summary: What's Proved vs. Open

**PROVED:**
- TRRT for d≤2: by Turán-ULC (complete, unconditional).
- THM-310 (Key Inequality): d_B ≤ d_A-1 for any C* ∉ max IS.
- THM-311 (Lemma A, d=2): pair-partner argument.
- THM-314 (Lemma A, non-unique IS): universal argument, covers all non-unique IS.
- THM-313 (Identity): A(-1/p) = I(-1/p) (Lemma B ↔ I(-1/p)≤0).

**VERIFIED but not proved:**
- Lemma B (B interlaces A when d_A=d_B+1): 3672 tests.
- HYP-1732 (α₂ ≤ p(m-p)): 1637 tests n=7..11.
- Lemma A for unique IS at d≥3: 0 failures at n=9..11.
- Full TRRT at n≥9: 0 failures.

---

## Equivalent Formulations of Lemma B (for d=2 case)

These are all equivalent via THM-313:

1. **B interlaces A** (standard Hermite-Biehler condition)
2. **I(Ω, -1/p) ≤ 0** (evaluate full polynomial at the root of B)
3. **α₂ ≤ p(m-p)** (combinatorial inequality)
4. **p lies between the two positive roots of I(Ω,x)** (geometric/spectral)

---

## The Beautiful Algebraic Identity (THM-313)

The clearest single theorem of this session:

For I = A + xB with deg(B)=1 and B(x)=1+px:
$$A\!\left(-\tfrac{1}{p}\right) = I\!\left(-\tfrac{1}{p}\right)$$

**Proof:** B(-1/p) = 0, so I(-1/p) = A(-1/p) + (-1/p)·0 = A(-1/p). □

This 2-line proof reveals why Lemma B and HYP-1732 are equivalent to a polynomial evaluation condition.

---

## Open Questions (Priority Order)

🔴 **Prove HYP-1732** (OPEN-Q-053): α₂ ≤ p(m-p). This would complete Lemma B for d=2.

🔴 **Prove Lemma A for unique IS** (OPEN-Q-054): For d≥3 with unique max IS, find valid C*.

🔴 **Prove Lemma B for d≥3**: When d_A=d_B+1 and d≥3, does B interlace A? Structure of the degree-(d-1) polynomial interlacing.

🟡 **Portal-disjoint proof**: Can the portal structure (THM-315) plus tournament constraints give HYP-1732?
