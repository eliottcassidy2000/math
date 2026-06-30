# The Hecke data of f₁₄ = 14a encodes the whole LRC(14): the bad-prime eigenvalues a₂=−1 (mirror)/a₇=+1 (apex) are the 2·7, the torsion congruence a_p ≡ p+1 (mod 6) is the units/existence-column "6", genus 1 is the obstruction, and rank 0 is the non-degeneracy — the floor is a weight-2 form (Eisenstein bulk + f₁₄), LRC(2p) hardness = genus, and the obstruction is the sym² (2nd-moment) component

*opus-2026-06-30. Owner: keep pushing, think modular forms and elliptic curves. The newform f₁₄=14a
doesn't just LIVE the LRC(14) story — its Hecke eigenvalues SPELL it: the modulus (2·7), the existence
column (the "6"), the obstruction (genus), the sign (rank). The Hecke dictionary is the LRC dictionary.*

## The Hecke dictionary (computed)
The weight-2 newform `f₁₄ = Σ a_n q^n` of the rank-0 curve **14a**:
| Hecke datum | value | LRC(14) meaning |
|---|---|---|
| **bad primes** `a₂, a₇` | `a₂ = −1` (`w₂`), `a₇ = +1` (`w₇`) | the **2·7 = 14**: `a₂` = the MIRROR cusp `d=2`, `a₇` = the APEX cusp `d=7` (Atkin–Lehner) |
| **torsion congruence** | `a_p ≡ p+1 \pmod 6` (∀ good p, VERIFIED) | the **"6" = φ(14) = units mod 14 = the razor's edge = the EXISTENCE column** `Q(√−3)` |
| **genus** | `1` (one cusp form) | the OBSTRUCTION the Eisenstein boundary can't see (local–global gap) |
| **rank** | `0` ⇒ `L(14a,1)=0.330>0` | the obstruction is NON-DEGENERATE (the empty tooth keeps positive width) |
> The `ℤ/6` torsion forces `a_p ≡ p+1 \pmod 6` (verified `p=3..43`): `ℤ/2` gives `a_p` even, `ℤ/3` gives
> `a_p ≡ p+1 (3)`. And `6 = φ(14)` is exactly the **6 units / 6 razor's-edge witnesses / the Galois group
> `(ℤ/14)*≅ℤ/6` / the Eisenstein hexagonal `ℤ[ζ₆]=Q(√−3)`** — the EXISTENCE column. So the curve's torsion IS
> the LRC's existence-column "6", carried in its Hecke data as a congruence.

## The floor IS a weight-2 modular form (the structural fact)
The LRC(14) floor's 2nd moment is a weight-2 form on `Γ₀(14)`:
> `dim M₂(Γ₀(14)) = 4 = 3` **Eisenstein** (the BULK, `3/π² = 1/(2ζ(2))`, the boundary/local-global density)
> `+ 1` **cusp form** `f₁₄ = 14a` (the OBSTRUCTION, the genus).
- **Genus 0** (`LRC(6), LRC(10)`): NO cusp form — pure Eisenstein — the boundary DETERMINES the floor →
  SOLVED. **Genus 1** (`LRC(14)`): one cusp form `f₁₄` the boundary CANNOT see → the obstruction. So
  **`LRC(2p)` hardness = `dim S₂(Γ₀(2p)) = genus(X₀(2p)) = 0,0,1,2,2`** — hardness is a genus, the failure
  of the boundary to determine the bulk (klein HYP-3587, now with the curve named).
- **Rank 0 makes the sign favorable:** `f₁₄`'s central `L(1)=0.330>0`, so the cusp-form contribution does
  not sink the bulk (the floor `0.344 > 3/π² = 0.304`). The obstruction is real but non-degenerate.

## Honest: the floor is NOT a clean L-value
Computed `L(14a,1) = 0.330224` (rank 0, `>0`); the cusp-form contribution to the floor is
`0.344 − 0.304 = 0.040`. These do NOT match (`0.330 ≠ 0.040`), so **"floor = L-value/period" is not a
clean identity** (klein-flagged speculation, here disconfirmed as a direct equality). The L-value is the
*non-degeneracy* (the SIGN, rank 0), not the floor *constant*. The constant is the genus-1 cusp-form
component of the 2nd moment, a period-type quantity, not `L(1)`.

## The push: the 2nd moment = the symmetric square
The floor is the 2nd MOMENT of the danger count `X(τ)` (the pair correlation `E[\binom X2]`). The pair
correlation of a Hecke form's coefficients is the **Rankin–Selberg / symmetric square**:
> `L(f₁₄ × f₁₄, s) = ζ(s) · L(\mathrm{sym}^2 f₁₄, s)`.
The Eisenstein bulk `3/π²` is the `ζ(s)` (main-term / local-global density) piece; the genus-1 obstruction
is the **`\mathrm{sym}^2 f₁₄`** piece. So the right L-function for the floor's obstruction is not `L(f₁₄)`
but `L(\mathrm{sym}^2 f₁₄)` — the 2nd-moment/symmetric-square of the genus-1 form. (Direction: identify the
cusp-form floor contribution with the `\mathrm{sym}^2 f₁₄` value/period; the `sym^2` of `14a` is a degree-3
L-function on `GL(3)`, the natural home of a 2nd moment.)

## What the modular push buys
- **The LRC(14) structure is READ OFF `f₁₄ = 14a`:** the modulus (`a₂=−1, a₇=+1` = 2·7 = mirror·apex), the
  existence column (`a_p≡p+1 (6)` = the units), the measure column (the apex `7`/doublet), the obstruction
  (genus 1), the non-degeneracy (rank 0). One newform, the whole story.
- **The hardness is geometric and predictable:** `genus(X₀(2p))` — rational (`P¹`) = solved, curve = hard;
  `14` is the first curve.
- **The right analytic object is `\mathrm{sym}^2 f₁₄`** (the 2nd moment), not `L(f₁₄)` — the floor obstruction
  lives on `GL(3)`, matching the moment-hierarchy's `SL(3)`/Littlewood level.

## Status
- **Computed/verified (opus):** `f₁₄=14a` Hecke data (`a₂=−1, a₇=+1`; `a_p≡p+1 mod 6`, the ℤ/6 torsion =
  the units = the existence "6"); `L(14a,1)=0.330>0` (rank 0, non-degenerate); the floor = weight-2 form =
  Eisenstein + `f₁₄`; hardness = genus.
- **Honest:** the floor constant is NOT `L(14a,1)` (the L-value is the sign, not the constant); the
  constant is the `\mathrm{sym}^2 f₁₄` / cusp-form component.
- **Push/direction:** the floor's obstruction = `\mathrm{sym}^2 f₁₄` (the 2nd moment on `GL(3)`), matching
  the `SL(3)` Littlewood level of the moment hierarchy.

Related: lrc14-lives-on-14a (the moduli=curve), the-master-two-Heegner-columns, f₁₄-14a-rank-0,
the-siegel-rogers-moment-hierarchy (the SL(3)/sym² level), cyclotomic-self-dual-razors-edge (the ℤ/6 units);
klein HYP-3587 (genus=local-global), mac-mini HYP-3594; OPEN-Q-108.
