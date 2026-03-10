# CONJ-002: H(T_p) for Paley Tournaments

**Type:** Conjecture
**Certainty:** 0 — REFUTED for p=11; open whether corrected formula exists
**Status:** REFUTED (p=11); open for p=19 and beyond
**Last reviewed:** kind-pasteur-2026-03-05-S2
**Disputes:** none
**Tags:** #paley-tournament #hamiltonian-paths #conjecture #automorphism-group #refuted

---

## Statement

For the Paley tournament T_p (p ≡ 3 mod 4 prime, connection set = QRs mod p):

$$H(T_p) = |\text{Aut}(T_p)| \cdot 3^{(p-3)/2}$$

where |Aut(T_p)| = p(p−1)/2.

---

## Verified Cases

| p | |Aut(T_p)| | 3^((p-3)/2) | Conjectured H(T_p) | Actual H(T_p) | Status |
|---|-----------|------------|---------------------|----------------|--------|
| 3 | 3 | 1 | 3 | 3 | ✓ |
| 7 | 21 | 9 | 189 | 189 | ✓ |
| 11 | 55 | 81 | 4455 | **95095** | REFUTED |
| 19 | 171 | 3^8=6561 | 1,121,931 | 1,172,695,746,915 | REFUTED |

*Note: p=5 gives H(T₅) = 15 = 5·3 = |Aut(T₅)|·3^1. ✓*

**H(T_11) = 95095 = 5 × 7 × 11 × 13 × 19 = 55 × 1729** (directly enumerated by kind-pasteur-2026-03-05-S2)

---

## Status for p = 11 — REFUTED

The conjecture predicted H(T₁₁) = 55 · 81 = 4455.

**Actual value (kind-pasteur-2026-03-05-S2, direct enumeration):**
- **H(T₁₁) = 95095 = 5 × 7 × 11 × 13 × 19 = 55 × 1729**
- h_QR = h_NQR = 201 (directed Ham cycles in T₁₁\{0,1} and T₁₁\{0,2})
- c₉(T₁₁) = (55/2)(402) = 11055 (>> 258 required by conjecture)

**Why the conjecture failed:**
The c₉ constraint (c₉ ≤ 258) was violated by a factor of ~43.

**Complete cycle count table for T_11** (kind-pasteur-S5, confirmed from inbox/other.txt):

| k | c_k(T_11) | C(11,k) | c_k/C(11,k) | integer? |
|---|-----------|---------|-------------|----------|
| 3 | 55 | 165 | 1/3 | no |
| 4 | 165 | 330 | 1/2 | no |
| 5 | 594 | 462 | 9/7 | no |
| 6 | 1595 | 462 | 145/42 | no |
| 7 | 3960 | 330 | 12 | YES |
| 8 | 7425 | 165 | 45 | YES |
| 9 | 11055 | 55 | 201 | YES |
| 10 | 10681 | 11 | 971 | YES |
| 11 | 5505 | 1 | 5505 | YES |

**Structural observation**: c_k(T_11)/C(11,k) is an integer for k >= 7 = (p+3)/2. Note c_6=1595 is NOT divisible by C(11,6)=462 (1595/462 = 145/42). The correct integrality threshold is k >= (p+3)/2, not (p+1)/2.

**p=19 data**: H(T_19) = 1,172,695,746,915, |Aut|=171, H/|Aut|=6,857,869,865.

**OCF verification for T_11** (other.txt, cross-confirmed with H=95095):
  95095 = 1 + 2*(55+594+3960+11055+5505) + 4*10879 + 8*1155
        = 1 + 42338 + 43516 + 9240 = 95095  EXACT

Where alpha_2 = 10879 (pairs of vertex-disjoint odd cycles):
  (3,3)-pairs: 495, (3,5)-pairs: 3630, (3,7)-pairs: 4840, (5,5)-pairs: 1914
And alpha_3 = 1155 (triples):
  (3,3,3)-triples: 550, (3,3,5)-triples: 605
And alpha_4 = 0 (would need >= 12 vertices).

Note: c_7=3960, NOT 1320 as claimed in more.txt (see MISTAKE-007). The trace-method computation in more.txt had errors.

**Note on factorization:** 1729 = H(T₁₁)/|Aut(T₁₁)| is the Hardy-Ramanujan "taxicab number" (7×13×19 = 12³+1³ = 10³+9³). Whether this is coincidence or structure is unclear.

**Computation files:** 03-artifacts/code/compute_H_T11.py, compute_ham_cycles_T11.py

---

## Discarded Conjecture

The original conjecture H(T_p) = p · 3^((p-1)(p-3)/8) is **FALSE** for p=11:
- Formula gives 11 · 3^10 = 649,539, which is not divisible by 55. ✗

---

## Key Results (RESOLVED)

- h_QR = h_NQR = 201 (computed kind-pasteur-2026-03-05-S2)
- c₉(T₁₁) = 11055
- H(T₁₁) = 95095 (refutes CONJ-002 for p=11)

## Open Questions (revised)

- OPEN-Q-013: Find the correct formula for H(T_p) (not 3^((p-3)/2))

## References

- Source: inbox/processed/2026-03-05/new/PALEY_T11_c9_ANALYSIS.md
- LEM-001: c₉ reduction formula
- LEM-002: sub-tournament structure
