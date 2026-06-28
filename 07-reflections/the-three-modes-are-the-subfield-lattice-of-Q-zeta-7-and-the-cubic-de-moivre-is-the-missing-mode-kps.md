# The three recursion modes ARE the subfield lattice of ℚ(ζ₇) — and the cubic de Moivre is the missing fourth mode

*kind-pasteur-2026-06-28-S31ar. The owner asked to see the similarities among the three signed recursions
(Möbius `A+B+C−D−E−F+G`, Eisenstein `A+B−C` even, Legendre `A+B+D−C−E−F+G` odd) and extend them, plus any
other relevant ones from a deep dive. The deep dive (codex-S114, HYP-2900/2901, `lrc_three_modes_parity_
composition_kps`) plus my own recent threads (HYP-3216 ladder, HYP-3212 Fejér) collapse into one picture:
**the three modes are Gaussian-period / character recursions indexed by the SUBFIELD LATTICE of `ℚ(ζ₇)`
(the divisors of `6 = φ(7)`), times the 2-adic Eisenstein factor — and `14 = 2·7` is exactly that. The
mode the project did NOT name is the CUBIC de Moivre mode (degree 3), which is precisely HYP-3216.**

## The three modes, decoded (VERIFIED)
The recursions reproduce `h(n)=⌊(n−1)²/4⌋` (PRONIC `=k(k−1)` for even `n`, SQUARE `=k²` for odd), via:
| mode | string | reduce by | parity | = subfield of ℚ(ζ₇) | character | Gaussian periods |
|---|---|---|---|---|---|---|
| **Möbius** | `A+B+C−D−E−F+G` | {1,2,3} | all | `ℚ` (degree 1) | trivial / `μ` | — (IE skeleton) |
| **Eisenstein** | `A+B−C` | {1,2} | even | the **2-adic** `ℤ/2` (`14=2·7`) | complement `T↔Tᵒᵖ` | the pronic / 2-fold |
| **Legendre** | `A+B+D−C−E−F+G` | {1,2,3,4} | odd | `ℚ(√−7)` (degree 2) | quadratic `χ₇` | `η_QR−η_NQR=√−7` ✓ |

`η_QR = ζ+ζ²+ζ⁴`, `η_NQR = ζ³+ζ⁵+ζ⁶`; their difference is the Gauss sum `√(−7) = 2.6458i` (verified). The
Legendre 7-term string is a **3-set Venn** (corners `A,B,D` sizes `N−1,N−2,N−1`; overlaps `C,E,F` sizes
`N−2,N−3,N−3`; triple `G` size `N−4`) — and `7 = 2³−1` terms = a **depth-3 Venn**.

## The MISSING fourth mode: the cubic de Moivre (= HYP-3216)
`6 = φ(7)` has divisors `{1,2,3,6}` → four subfields of `ℚ(ζ₇)`. The project named the degree-1 (Möbius)
and degree-2 (Legendre) modes. **Degrees 3 and 6 are exactly my recent threads:**
| degree | subfield | mode | object |
|---|---|---|---|
| 1 | `ℚ` | Möbius | the IE skeleton |
| 2 | `ℚ(√−7)` | Legendre (quadratic) | the Gauss sum `√−7` |
| **3** | **`ℚ(cos 2π/7)`** | **the CUBIC de Moivre mode = HYP-3216** | **the 3 cubic periods `2cos(2πj/7)` = the de Moivre angles** ✓ |
| 6 | `ℚ(ζ₇)` | the full sextic | the Fejér kernel `F₇=(de Moivre)²` (HYP-3212/3214) |

So **HYP-3216's "cyclotomic moment-order ladder" IS the degree-3 (cubic) mode** the three-mode framework
was missing — the cubic Gaussian periods are the de Moivre angles (verified `{−1.8019,−0.4450,1.2470}`),
and the de Moivre cubic `x³+x²−2x−1` is its minimal polynomial.

## The unifying similarity: depth = cyclotomic degree = (p−1)/2
Every face of the structure carries the SAME number, the **cyclotomic degree `(p−1)/2`** of the apex `p`:
- the Legendre **Venn-set depth** (`(p−1)/2`-set, `2^{(p−1)/2}−1` terms): `p=3→1`-set, `p=5→2`-set (3 terms,
  `A+B−C`!), `p=7→3`-set (7 terms) — VERIFIED;
- the **moment-order ladder depth** (HYP-3216): `(p−1)/2`;
- the **de Moivre / cubic degree**: `(p−1)/2`;
- the number of **Gaussian periods** in the maximal real subfield: `(p−1)/2`.
> **`A+B−C` (the even/Eisenstein string) is also the `p=5` Legendre Venn** — the depth-2 case! The
> recursion strings are the SAME object at different `(p−1)/2`: 1-set, 2-set (`A+B−C`), 3-set
> (`A+B+D−C−E−F+G`). The "even vs odd" of the project is the `(p−1)/2`-set Venn growing with the apex.

## The full picture: `14 = 2·7` = Eisenstein(2-adic) × [the ℚ(ζ₇) subfield tower]
```
LRC(2p)  =  (the 2-adic EISENSTEIN fold, the "2")  ×  (the ℚ(ζ_p) subfield tower, the "p"):
                                                        ┌ degree 1: MÖBIUS (μ, skeleton)
                                                        ├ degree 2: LEGENDRE (√−p, quadratic)
                                                        ├ degree 3: CUBIC de MOIVRE (HYP-3216 ladder)   ← apex node for p=7
                                                        └ degree (p−1)/2..(p−1): up to the full Fejér F_p
```
The proof **descends this subfield tower**: the binding (apex) node sits at the deepest realized degree
`(p−1)/2` (for `p=7`, the cubic — the de Moivre/biquadratic node `k=8`), then the 2-adic Eisenstein fold
halves it (`14=2·7`, my HYP-3216 2-recursion), bottoming at the degree-2 (quadratic/pairwise) covariance
statement. **The two "interlocking recursions" (HYP-3216) are the `2`-factor (Eisenstein) and the `p`-tower
(Möbius→Legendre→cubic→sextic); the three project modes are the bottom of that tower.**

## Other relevant recursions found (the deep dive)
- **mac-mini-S75b's Farey/Stern–Brocot cap-kernel recursion** = the **Diophantine** face: the continued-
  fraction / three-gap structure is the *analytic* side of the same tower (the `L`-functions of the
  characters; the Dedekind sums are the Gauss-sum companions). Möbius=`ζ`, Legendre=`L(χ₇)`, cubic=`L(χ₃)`.
- **The Worpitzky/Eulerian odd-dip leg** (codex HYP-3147) = the **odd-degree (cubic)** part of the dip —
  the cubic-mode residue, consistent with the de Moivre cubic being the degree-3 mode.
- **`h(n)=⌊(n−1)²/4⌋`** pronic↔square parity IS the Eisenstein(even)↔Legendre(odd) split: the `2`-adic vs
  `p`-cyclotomic, the same `14=2·7`.

## EXTENSION (new, testable)
1. **The cubic-mode recursion string.** Just as Legendre is the quadratic-character (`χ₇`) IE and is a
   3-set Venn at `p=7`, the **cubic de Moivre mode is the order-3-character (`χ₃ mod 7`, exists since
   `3|6`) IE** — write its explicit sign string (the `χ₃`-twisted `A..G`) and check it reproduces the
   degree-3 invariant (the de Moivre resolvent / the `k=8` dip). This *names* HYP-3216's mode as a signed
   recursion, completing the project's `{Möbius, Eisenstein, Legendre}` to `{…, cubic, sextic}`.
2. **The L-function ladder.** Each mode = a Dirichlet `L(s,χ)` for a character mod `7` (orders 1,2,3,6).
   The cap/dip rationals (`disc 7²`, S31an) are special `L`-values (`L(1−n,χ)`, generalized Bernoulli) —
   one per mode. The Stark lead (HYP-3215) is the `L'(0)` companion of the cubic mode.
3. **Family:** for `LRC(2p)`, the modes are the divisors of `(p−1)`; the binding node is the degree-
   `(p−1)/2` (cubic-or-higher) mode; `n=14` is first-hard because `p=7` is the first apex whose
   `(p−1)/2 = 3` exceeds the quadratic (the periodic-CF wall, S31ak).

## Net
- **SIMILARITY:** the three modes (Möbius/Eisenstein/Legendre) are the **subfield/character lattice of
  `ℚ(ζ₇)` × the 2-adic** — the same `14=2·7` tower as HYP-3216; the strings are one Venn at growing depth
  `(p−1)/2`; `A+B−C` is both the even mode and the `p=5` Legendre.
- **THE MISSING MODE:** the degree-3 **cubic de Moivre mode = HYP-3216** (the 3 cubic Gaussian periods =
  the de Moivre angles); the sextic = the Fejér kernel. The project had the bottom two; my threads are the
  top two.
- **EXTEND:** write the cubic-character (`χ₃ mod 7`) recursion string; the `L`-value ladder (one per mode);
  the family law (binding mode = degree `(p−1)/2`).

→ HYP-3217 (this), HYP-2900/2901 (three modes), HYP-3216 (= the cubic mode), HYP-3212/3214 (de Moivre/Fejér
= sextic), HYP-3147 (Worpitzky = cubic-odd dip), HYP-3162 (cyclotomic), HYP-3215 (Stark `L`), mac-mini-S75b
(Farey = `L`-functions), OPEN-Q-108.
