# THM-153: Geometric Growth of Omega_m for Paley Tournaments

**Status:** VERIFIED (computational, p=7,11,19,23)
**Session:** opus-2026-03-13-S70
**Depends on:** THM-151 (Omega_3 formula), GLMY definitions

## Statement

For the Paley tournament at prime p ≡ 3 mod 4, with Q_k = (p+1)/4:

$$\Omega_m(\text{Paley}_p) = \binom{p}{2} \cdot (Q_k - 1)^{m-1} \quad \text{for } 1 \leq m \leq 4$$

where Q_k = |Ŝ(k)|² is the Fourier magnitude (constant across all k ≠ 0 for Paley).

Equivalently, with q = (p-3)/4:

$$\Omega_m = \frac{p(p-1)}{2} \cdot q^{m-1}$$

## Special case: p = 7

When p = 7: q = (7-3)/4 = 1, so Ω_m = C(7,2) = 21 for ALL 1 ≤ m ≤ 6.
This is the UNIQUE Paley tournament with **constant Omega profile**.

Moreover, the M matrix (regular m-path counts from a to b) has constant
off-diagonal entries: M[a,b] = 1 for all a ≠ b, for ALL m.

## Verification

| p | Q_k | q | Ω₁ | Ω₂ | Ω₃ | Ω₄ | Ω₅ actual | Ω₅ predicted |
|---|-----|---|-----|------|-------|---------|-----------|-------------|
| 7 | 2 | 1 | 21 | 21 | 21 | 21 | 21 | 21 |
| 11 | 3 | 2 | 55 | 110 | 220 | 440 | 715 | 880 |
| 19 | 5 | 4 | 171 | 684 | 2736 | 10944 | 40698 | 43776 |
| 23 | 6 | 5 | 253 | 1265 | 6325 | 31625 | 149017 | 158125 |

Formula exact for m ≤ 4, deviates at m = 5 for p > 7.

## Proof sketch for m = 2

Ω₂ = C(p,3) - t₃. For Paley at p ≡ 3 mod 4:
- t₃ = p(p²-1)/24 (known formula for doubly regular tournaments)
- Ω₂ = p(p-1)(p-2)/6 - p(p²-1)/24 = p(p-1)(p-3)/8 = C(p,2) · (p-3)/4 ✓

## Proof sketch for m = 3

Using THM-151: Ω₃ = C(p,4) - Σ_v[c₃(N⁻(v)) + c₃(N⁺(v))].
By circulant symmetry, all neighborhoods have the same c₃ value C.
So Ω₃ = C(p,4) - 2pC.

For the formula to hold: C = (p-1)(p-3)(p+1)/192.
At p=7: C = 6·4·8/192 = 1 ✓ (each 3-vertex neighborhood is a 3-cycle).

## Key insight: why p = 7 is unique

Q_k = (p+1)/4 = 2 means each Fourier component has magnitude √2.
The "amplification factor" per path step is Q_k - 1 = 1, so there is NO
amplification — the path count remains constant at C(p,2).

This is related to the doubly regular property: |N⁺(i) ∩ N⁺(j)| = 1 for
all i ≠ j. Combined with Q_k = 2, the eigenspace structure is maximally
uniform.

## Notes

- The geometric growth Ω_m ∝ q^{m-1} suggests each path-extension step
  multiplies the count by q = (p-3)/4 ≈ Q_k - 1.
- The deviation at m ≥ 5 is a finite-size effect: paths can't reuse vertices,
  so the geometric model breaks when paths become long relative to n.
- Connection to THM-145: since all Q_k are equal, all eigenspaces contribute
  identically to each Ω_m.
