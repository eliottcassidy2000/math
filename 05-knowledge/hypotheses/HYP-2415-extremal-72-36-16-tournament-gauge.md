# HYP-2415 — The extremal [72,36,16] code as a tournament-gauge question

**Status:** OPEN (claimed claudebox-2026-06-11-S4). See OPEN-Q-061 for the full statement.
**Companions:** THM-481 (eQR ladder), THM-484 (involution modulus 24), THM-480/482.

## Claim / program

The eQR tournament-gauge ladder C(I+S(Paley_q)) is extremal Type II at q = 7, 23, 31, 47
(d = 4, 8, 8, 12 = 4⌊(q+1)/24⌋+4, verified) and first fails at q = 71 (eQR(72) has d = 12,
extremal would be 16). Order 72 ≡ 8 mod 16 ⟹ every skew-Hadamard-of-order-72 gauge code is a
Type II [72,36] code. SUFFICIENT route to the famous open extremal [72,36,16]: find an
order-72 skew-Hadamard whose gauge minimum distance is 16.

## Tests / sub-questions

1. Compute gauge minimum-distance LOWER bounds for the catalogued order-72 skew-Hadamard
   classes (Đoković–Kotsireas) via random coset-leader / Brouwer–Zimmermann partial enumeration
   (full 2^36 infeasible); report the best d found and the spectral profile of its H.
2. Which tournament invariant of H controls the gauge d? (Paley = flat character spectrum gives
   d = 12; does a peaked/structured H lift it?) — ties to the flat-vs-peaked dichotomy (THM-441).
3. Honest caveat: NOT an equivalence — a [72,36,16] code need not be a skew-Hadamard gauge; this
   is one construction route among many (cf. the long history of failed constructions).
4. Conjecture (weak): no Paley-type / highest-symmetry order-72 skew-Hadamard reaches d=16
   (extremal codes have trivial-ish automorphism groups; the very symmetry that makes Paley
   computable caps its distance) — the extremal code, if it exists, is asymmetric.
