# LEM-002: Structure of the Sub-Tournament T₁₁\{0,1}

**Type:** Lemma
**Certainty:** 5 — PROVED (exact computation from definition)
**Status:** PROVED
**Last reviewed:** kind-pasteur-2026-03-05-S1
**Disputes:** none
**Tags:** #paley-tournament #sub-tournament #score-sequence #3-cycles #structure

---

## The Tournament

T₁₁\{0,1}: the Paley tournament T₁₁ (connection set QR = {1,3,4,5,9} mod 11) with vertices 0 and 1 deleted.

**Vertex set:** {2, 3, 4, 5, 6, 7, 8, 9, 10}

**Arc rule:** i → j iff j−i (mod 11) ∈ {1, 3, 4, 5, 9}

---

## Score Sequence (Out-Degrees)

Each vertex originally had out-degree 5 in T₁₁. Deletion of vertices 0 and 1 removes arcs v→0, v→1, 0→v, 1→v for each remaining vertex v.

- v → 0 iff −v (mod 11) ∈ QR iff v ∈ NQR = {2,6,7,8,10}
- v → 1 iff 1−v (mod 11) ∈ QR iff v ∈ {3,7,8,9} ∩ {2,...,10} = {3,7,8,9}

| Vertex | v→0? | v→1? | Out-deg in T₁₁\{0,1} |
|--------|------|------|----------------------|
| 2 | yes | no | 4 |
| 3 | no | yes | 4 |
| 4 | no | no | 5 |
| 5 | no | no | 5 |
| 6 | yes | no | 4 |
| 7 | yes | yes | 3 |
| 8 | yes | yes | 3 |
| 9 | no | yes | 4 |
| 10 | yes | no | 4 |

**Score sequence:** (4, 4, 5, 5, 4, 3, 3, 4, 4) [vertices 2 through 10]

Sum = 36 = C(9,2). ✓

---

## 3-Cycle Count

By Moon's formula: c₃(T) = C(n,3) − Σ_v C(d⁺(v), 2)

$$c_3(T_{11}\setminus\{0,1\}) = \binom{9}{3} - \sum_v \binom{d^+(v)}{2}$$

$$= 84 - \left[\binom{4}{2}+\binom{4}{2}+\binom{5}{2}+\binom{5}{2}+\binom{4}{2}+\binom{3}{2}+\binom{3}{2}+\binom{4}{2}+\binom{4}{2}\right]$$

$$= 84 - [6 + 6 + 10 + 10 + 6 + 3 + 3 + 6 + 6] = 84 - 56 = 28$$

**T₁₁\{0,1} has 28 directed 3-cycles.**

For comparison: T₉ (9-vertex regular tournament) has C(9,3)/3 = 28 directed 3-cycles on average, so this sub-tournament is near the average 3-cycle density.

---

## Hamiltonian Cycle Count h_QR — RESOLVED

**UPDATE (kind-pasteur-2026-03-05-S2):** Direct computation gives:
- h_QR = h({0,1}) = 201
- h_NQR = h({0,2}) = 201

So h_QR + h_NQR = 402, and c₉(T₁₁) = (55/2)(402) = 11055.
The conjecture CONJ-002 is REFUTED: H(T₁₁) = 95095, not 4455.

---

## Sub-Tournament T₁₁\{0,2} (NQR Case)

By a parallel computation with vertex 2 replaced (2 is NQR), the structure of T₁₁\{0,2} can be derived similarly.

**RESOLVED:** h_NQR = h({0,2}) = 201 (same as h_QR). CONJ-002 is refuted.

## References

- Source: inbox/processed/2026-03-05/new/PALEY_T11_c9_ANALYSIS.md (Sections IV, V)
- Depends on: LEM-001
- Relevant to: CONJ-002
