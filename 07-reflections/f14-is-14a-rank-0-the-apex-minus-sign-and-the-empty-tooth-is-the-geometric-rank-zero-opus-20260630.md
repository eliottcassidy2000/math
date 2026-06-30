# f₁₄ = the elliptic curve 14a: the apex cusp carries Atkin–Lehner w₇ = −1 (the "minus"/hardness), but 14a has RANK 0 so L(14a,1) > 0 — the obstruction's sign is FAVORABLE; the Dirac-comb's irreducible empty tooth is the geometric form of that non-vanishing, and the proved atom 4cos²(3π/7)>0 (THM-590) is its cyclotomic shadow. The wall is the DESCENT, not the sign

*opus-2026-06-30. Owner: dig into f₁₄/14a's cusp coefficient and whether the Dirac-comb/empty-tooth says
anything about its sign. It does — the sign is favorable, and the empty tooth IS the rank-0.*

## f₁₄ = 14a, computed
The genus-1 weight-2 newform of level 14 is the rank-0 elliptic curve **14a** (`y²+xy+y=x³+4x−6`).
Point-counts give the Hecke eigenvalues:
> `f₁₄ = q − q² − 2q³ + 0·q⁵ + q⁷ + … `  (`a₂=−1, a₃=−2, a₅=0, a₇=+1, a₁₁=0, a₁₃=−4, …`).
**Atkin–Lehner** (for `p∥14`, weight 2, `a_p=−w_p`): `w₂ = −a₂ = +1`, `w₇ = −a₇ = −1`, Fricke
`w₁₄ = w₂w₇ = −1`. **Root number `ε = −w₁₄ = +1` ⇒ even rank ⇒ `L(14a,1) ≠ 0` — and 14a has RANK 0, so
`L(14a,1) > 0`.**

## The apex sign, and why it is favorable anyway
> **The apex cusp `d=7` carries `w₇ = −1` — the "minus" cusp (klein's `d=7` APEX-hard/−).** `f₁₄|W₇ =
> −f₁₄`, so the apex-cusp leading coefficient is the `−1`-image of the `q=∞` coefficient `a₁=1`. The
> minus sign is the *hardness* (the apex is where the obstruction concentrates).
> **BUT the curve has rank 0, so `L(14a,1) > 0`.** The cusp form is *non-degenerate* — its central
> L-value does not vanish. So the genus-1 obstruction, while carrying the apex minus-sign, is a
> NON-ZERO POSITIVE mode, not a sink. Empirically (klein HYP-3593): the floor `inf R' = 114382/332563 =
> 0.3439 > 3/π² = 0.3040` — **the cusp-form contribution is net `+0.040`: it does NOT sink the bulk, it
> ADDS.** And the worst LOCAL gap, the apex doublet `4cos²(3π/7) = 0.198 > 0`, is PROVED positive
> (THM-590). **The sign is favorable on every count.**

## The Dirac-comb / empty-tooth IS the geometric rank-0
What does the comb say about the sign? This:
> **The observer's empty tooth in Ш₁₄ is irreducible and non-degenerate — the gap exists and has positive
> width `1/14` — and THAT is the geometric form of `L(14a,1) > 0` (rank 0).** The cusp form `f₁₄` is the
> obstruction *mode*; its non-vanishing (rank 0) means the mode is "really there," a single irreducible
> hole — exactly the comb's one empty tooth. If 14a had positive rank (`L(1)=0`), the obstruction mode
> would be degenerate — the hole could close — and the floor could vanish. **The empty tooth being a
> genuine, single, measure-zero-but-positive-width gap is the real-space picture of "genus 1, rank 0,
> L(1)>0":** one global mode, non-degenerate, favorable.
The dictionary, completed:
| modular (f₁₄ = 14a) | geometric (the comb) | cyclotomic (the atom) |
|---|---|---|
| genus 1 = one cusp form | the ONE empty tooth | one obstruction atom |
| `w₇ = −1` (apex minus) | the apex/mult-of-7 reflection of the tooth | the doublet's `k(b−a)≡3` (angle nearest π) |
| rank 0, `L(1) > 0` | the gap is non-degenerate (width `1/14 > 0`) | `4cos²(3π/7) > 0` (THM-590, PROVED) |
| Eisenstein bulk (dim 3) | the 13 filled teeth | `3/π² = 1/(2ζ(2))` |

## The honest wall: it is the DESCENT, not the sign
The sign is settled — favorable (atom `0.198>0` proved; floor clears bulk; rank 0 ⇒ `L(1)>0`). **The wall
is the REDUCTION**: that the LRC floor for the FULL covering set equals the `Z₇`-core cyclotomic gap
`g(O)` (mac-mini's descent `ρ_j = g(descended core)`) is CONDITIONAL — THM-590 bounds the core gap, but
the descent to the full set is open. So:
> **LRC(14) is not blocked by a bad sign — the obstruction mode `f₁₄` is non-degenerate and favorable
> (rank 0, `L(1)>0`, atom `>0` proved). It is blocked by the DESCENT** — proving the full floor inherits
> the proved core positivity. The empty tooth is provably open; what is unproved is that no covering set
> can route around it.
(klein flagged "floor = L-value/period" as speculative — and it must be subtle, else LRC(14) would
follow from the known rank-0 of 14a. The L-value is the obstruction's *non-degeneracy* (the sign), not the
*floor value*; the descent is the missing link between them.)

## Status
- **Computed (opus):** `f₁₄=14a`: `a₂=−1,a₇=+1`; `w₂=+1, w₇=−1`; root number `+1`; rank 0 ⇒ `L(14a,1)>0`.
- **The sign is FAVORABLE:** apex `w₇=−1` is the hardness, but rank-0 non-degeneracy + the proved atom
  `4cos²(3π/7)>0` + the empirical floor `0.344>0.304` all say the obstruction does not sink the floor.
- **New picture:** the comb's irreducible empty tooth = the geometric `L(14a,1)>0` (rank 0) = one
  non-degenerate global mode; the three languages (modular/geometric/cyclotomic) now have a full
  dictionary.
- **The wall (honest):** the DESCENT (full floor ← `Z₇`-core gap), not the sign. The atom is proved; the
  reduction is open.

Related: the roots-of-unity-convergence reflection (the comb ≡ atom ≡ cusp form); klein THM-590 (atom),
HYP-3586/3587 (X₀(14), cusp form f₁₄=14a, genus); mac-mini HYP-3594 (odd cycle), the descent machinery;
my Dirac-comb + Z₇-SOS + Eisenstein-Φ₆ reflections; HYP-3547 (apex-7 Heegner); OPEN-Q-108.
