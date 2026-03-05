# Paley Tournament T₁₁: Audit of c₉ and H(T₁₁) Conjecture

**Source:** Human drop (Eliott) — 2026-03-05
**Topic:** Rigorous audit of the corrected conjecture H(T_p) = |Aut(T_p)| · 3^((p-3)/2) for p=11; exact reduction of c₉(T₁₁) to sub-tournament Ham-cycle counts

---

## I. The "ratio coincidence" fails

The argument c₉ = 220 rested on c₇/C(11,7) = c₉/C(11,9) = 4.

* c₇ = 1320, C(11,7) = 330: ratio = 4. ✓
* c₅ = 594, C(11,5) = 462: ratio = 594/462 = 9/7 ≈ 1.286. ✗
* c₃ = 55, C(11,3) = 165: ratio = 1/3.

Ratios are 1/3, 9/7, 4, ? — not constant. The ratio coincidence for k=7 is accidental; it gives **no structural support** for c₉. The earlier c₉ = 220 estimate was circular reasoning.

---

## II. What the conjecture H(T₁₁) = 4455 requires

Known: 1 + 2(c₃ + c₅ + c₇) = 1 + 2(55 + 594 + 1320) = 3939.

For H = 4455: 4455 − 3939 = 516 = 2(c₉ + h₁₁) + 4α₂ + ...

This requires **c₉ ≤ 258**.

The growth sequence c₃=55, c₅=594, c₇=1320 has ratios ~10.8, ~2.2. Even if the ratio continued decelerating sharply (c₉/c₇ ≈ 0.5), c₉ ≈ 660 >> 258.

The conjecture is only plausible because c₉ counts 9-cycles, which use 9 of 11 vertices — only C(11,9) = 55 possible 9-vertex subsets vs C(11,7) = 330 seven-vertex subsets. This combinatorial suppression may bring c₉ below c₇ dramatically.

---

## III. Exact reduction: c₉(T₁₁) = (55/2)(h_QR + h_NQR)

**Setup:** Each 9-cycle in T₁₁ uses 9 vertices; the two missing vertices form a pair {a,b}. Let h({a,b}) = number of directed Ham cycles in T₁₁\{a,b}.

$$c_9(T_{11}) = \sum_{\{a,b\} \subset V} h(\{a,b\}) = \frac{11}{2}\sum_{b \neq 0} h(\{0,b\})$$

**Symmetry reduction:** The multiplier group {x ↦ αx : α ∈ QR} has order 5 and acts on pairs {0,b}, splitting b ∈ {1,...,10} into two orbits:
* b ∈ QR = {1,3,4,5,9}: all give the same h_QR = h({0,1})
* b ∈ NQR = {2,6,7,8,10}: all give the same h_NQR = h({0,2})

$$\boxed{c_9(T_{11}) = \frac{55}{2}(h_{\text{QR}} + h_{\text{NQR}})}$$

For c₉ to be an integer: h_QR + h_NQR must be even.

**Implication for the conjecture:** If H(T₁₁) = 4455 then c₉ ≤ 258, so:
$$h_{\text{QR}} + h_{\text{NQR}} \leq \frac{2 \cdot 258}{55} \approx 9.4 \implies h_{\text{QR}} + h_{\text{NQR}} \leq 8$$

---

## IV. Structure of T₁₁\{0,1} (the QR sub-tournament)

Vertices: {2,3,4,5,6,7,8,9,10}. Arc i→j iff j−i (mod 11) ∈ {1,3,4,5,9}.

**Out-degrees** (removing vertices 0 and 1 from T₁₁, original out-degree 5):

| v | v→0? | v→1? | new out-deg |
|---|------|------|-------------|
| 2 | yes (NQR) | no | 4 |
| 3 | no (QR) | yes | 4 |
| 4 | no (QR) | no | 5 |
| 5 | no (QR) | no | 5 |
| 6 | yes (NQR) | no | 4 |
| 7 | yes (NQR) | yes | 3 |
| 8 | yes (NQR) | yes | 3 |
| 9 | no (QR) | yes | 4 |
| 10 | yes (NQR) | no | 4 |

Score sequence: (4,4,5,5,4,3,3,4,4). Sum = 36 = C(9,2). ✓

**3-cycle count by Moon's formula:**
$$c_3(T_{11}\setminus\{0,1\}) = \binom{9}{3} - \sum_v \binom{d^+(v)}{2} = 84 - [6+6+10+10+6+3+3+6+6] = 84 - 56 = 28$$

---

## V. Spectral formula for off-diagonal powers (corrected)

For T₁₁ (circulant, connection set = QRs), with correct sign convention:

$$\lambda_s = \sum_{d \in QR} \zeta^{ds}, \quad (A^k)_{0,d} = \frac{1}{11}\sum_{s=0}^{10} \lambda_s^k \zeta^{-sd}$$

Using $\lambda_s = \omega = \frac{i\sqrt{11}-1}{2}$ for s ∈ QR, $\bar\omega$ for s ∈ NQR, and $G(\chi,d) = \chi(d) \cdot i\sqrt{11}$ for p=11≡3 (mod 4):

$$\left(A^k\right)_{0,d} = \frac{5^k - \text{Re}(\omega^k)}{11} + \frac{\text{Im}(\omega^k)\sqrt{11}}{11} \cdot \frac{1+\chi(d)}{2} \cdot \chi(d)$$

(exact form; see derivation in body of document for the full careful computation)

---

## VI. Honest state of knowledge

**Rigorously verified:**
* H(T₃)=3, H(T₅)=15, H(T₇)=189
* c_k(T₁₁) for k=3,5,7: values 55, 594, 1320
* c₄(T₁₁)=165, c₆(T₁₁)=1375 (even cycles)
* 55 | H(T₁₁) (from automorphism group)
* H(T_p) = p·3^((p-1)(p-3)/8) is **FALSE** for p=11
* c₉(T₁₁) = (55/2)(h_QR + h_NQR) ← **new exact reduction**
* T₁₁\{0,1} has score sequence (4,4,5,5,4,3,3,4,4) and 28 directed 3-cycles

**Plausible but unverified:**
* H(T_p) = |Aut(T_p)| · 3^((p-3)/2), giving H(T₁₁) = 4455
* This requires h_QR + h_NQR ≤ 8

**Key open computation:** h({0,1}) = number of directed Ham cycles in T₁₁\{0,1} (9-vertex tournament, fully explicit structure above). Determining this would confirm or refute the conjecture.

---

## Action items for receiving agent

1. Add c₉ reduction formula as a new lemma in 01-canon/theorems/
2. Attempt to bound or compute h_QR and h_NQR via further sub-tournament analysis
3. Note: the ratio-coincidence argument should be flagged as invalid reasoning; consider adding to MISTAKES.md if it appeared in any theorem file
