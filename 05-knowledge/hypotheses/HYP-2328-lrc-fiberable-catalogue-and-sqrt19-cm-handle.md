# HYP-2328 — The fiberable-vs-prime-hard catalogue of LRC(n), and the √−19 CM handle on prime 19

**Session:** S650
**Status:** CONFIRMED (catalogue + Gauss-sum verified; the √−19 arithmetic core formalized)
**Provenance forward:** math-lean `Math/NumberTheory/HeegnerNineteen.lean` (sorry-free)
**Prompt:** (continuing S649) catalogue which `n` are fiberable vs prime-hard, and chase the √−19 / χ=5 CM
handle on `n = 19`.

---

## Part A — The reducibility catalogue: fiberability is a divisor property

The S640 fiber reduction needs an odd prime divisor `p ∣ n` with `LRC(p)` tractable. Since LRC is proven
for `≤ 7` runners, a clean three-tier classification emerges (`lrc_fiberable_catalogue_and_sqrt19_s650.py`,
`n = 8..32`):

| class | `n` | reduces? |
|---|---|---|
| **`n = 2p`, `p ≤ 7`** (and `3·odd`) | `6, 9, 10, 12, 14, 15, 18, 21, 24, 27, 28, 30…` | **YES → a PROVEN base** (fiber over `p`, S640) |
| **`n = 2p`, `p > 7`** | `22 (=2·11), 26 (=2·13)…` | reduces to an **OPEN base** `LRC(p)` |
| **`n` prime** | `11, 13, 17, **19**, 23, 29, 31…` | **NO fiber — the hard end** |
| **`n = 2ᵏ`** | `8, 16, 32` | fully 2-adic (no odd fiber) |

> **Fiberability is a divisor property.** `n = 14 = 2·7` reduces to the *proven* `LRC(7)`. `n = 19` is
> prime: no divisor, no fiber, `2` a primitive root (`ord₁₉(2) = 18`, single cycle) — verified
> *prime-hard*. The tractable frontier is the composite `2p` family (S640); the primes are the wall.

---

## Part B — The √−19 CM handle on prime 19 (the χ=5 frontier)

For prime `n = 19` the leverage is the **CM field**, not a fiber. The 19-runner witnesses `t = j/19` live
in `ℚ(ζ₁₉)`, and:

1. **`√−19` is the quadratic subfield of `ℚ(ζ₁₉)`** (the Gauss sum). Verified: the quadratic Gauss sum
   `g = Σₐ (a∣19) ζ₁₉ᵃ = i√19`, so **`g² = −19`** (exact). Because `19 ≡ 3 (mod 4)`, the sum is `√−19`
   (not `√19`) — the Heegner field lives inside the field where the 19-runner problem's rational
   witnesses live. `ℚ(√−19)` is the **index-2 (Paley/QR) level** of the cyclotomic witness tower
   (opus S704); `Gal(ℚ(ζ₁₉)/ℚ) = (ℤ/19)* = ℤ/18` is abelian, with subfields `↔ {1,2,3,6,9,18}`.
2. **Heegner / class number 1.** `ℚ(√−19)` has class number 1 — the conjectural **`χ = 5` step** in the
   Moser/Heegner tower:
   ```
     √−3  (χ=3, Eisenstein lattice) → √−11 (χ=4, Moser spindle) → √−19 (χ=5, hex(2) / 19-runner)
   ```
   `19 = 4·5 − 1` (the rotation field for Eisenstein norm `N = 5`, HYP-2277); `19 = 1+6+12 = hex(2)`;
   `2n−1 = 37 = hex(3)`.

**Formalized (math-lean, sorry-free): `Math/NumberTheory/HeegnerNineteen.lean`** — the arithmetic core:
- `neg_one_nonsquare_mod19` (`¬ IsSquare (-1 : ZMod 19)`): `19 ≡ 3 mod 4` ⟹ Gauss sum gives `√−19`,
  Paley-19 exists (S638).
- `two_nonsquare_mod19`, `two_pow_nine_mod19` (`2⁹ = −1`): `2` a non-residue ⟹ primitive root ⟹ `⟨2⟩ ≠ QR`
  ⟹ no doubling fiber (prime-hard).
- `card_qr_mod19` (`= 9`): the Paley-19 connection-set size `(19−1)/2`, the `√−19` index-2 level.

> **HONEST.** This *organizes* the 19-runner witnesses (their exact field is `ℚ(ζ₁₉)`, CM-subfield
> `√−19`, the Paley/Gauss/QR structure its arithmetic shadow). It does **not** prove LRC(19) (open). The
> Heegner/class-number-1 property is the rigidity hallmark (conjecturally the χ=5 forcing), not a proof.

---

## Synthesis: the two seams of the arc, one for each frontier number
- **`14` (composite) → the divisor seam.** Its `7` gives the fiber; the doubling is a small cube-root
  cycle; reducible to a proven base. The 2-adic / divisor half.
- **`19` (prime) → the CM seam.** Its `√−19` is the Gauss sum / quadratic subfield / Heegner χ=5 rung; no
  fiber, `2` a primitive root. The cube-root / CM half.

The whole arc has run between these two seams (the `σ`-involution/2-adic vs the `ω`-cube-root/CM); `14`
and `19` are the LRC frontier cases that sit one on each.

## New threads / handoffs
- **Formalize `g² = −19`** via Mathlib `gaussSum`/`quadraticChar` (`gaussSum_sq` for the quadratic
  character; `χ(−1) = −1` here) — would make the "√−19 ⊂ ℚ(ζ₁₉)" statement machine-checked.
- **The `q*` attack at `q = 19, 37`** (S704): does the `√−19` quadratic level give a CM bound on the
  cyclotomic depth, the way `√−11` rigidity forces the Moser spindle's 4th colour?
- **General `n = 2p` reduction theorem** (formalize the S640 fiber for the proven-base family) — the
  tractable frontier; contrast with the prime wall.
