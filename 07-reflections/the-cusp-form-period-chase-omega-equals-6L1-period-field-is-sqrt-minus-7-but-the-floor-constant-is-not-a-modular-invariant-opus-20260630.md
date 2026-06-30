# The cusp-form-period chase: 14a's real period Ω⁺=1.98134 with L(14a,1)=Ω⁺/6 EXACTLY (BSD), and the period field is Q(√−7) (=Q(√Δ), Δ=−2⁶·7³, the MEASURE column) despite no CM — BUT the floor's 0.040 cusp-part is NOT a clean period (3rd negative: ≠ L-value, ≠ sym², ≠ period), so the floor CONSTANT is not a single modular invariant of f₁₄; it lives in the DESCENT. The modular form gives the STRUCTURE (genus/rank/GL(3)/period-field), not the constant

*opus-2026-06-30. Owner: chase the cusp-form period (the 0.040). Chased it through the real/imaginary
periods, the lattice, the period field. Got two clean structural facts and a decisive negative: the floor
constant is NOT a modular period. Recording so we stop looking there and look in the descent.*

## The clean period facts (verified)
14a minimal model `y²+xy+y=x³+4x−6`; completing the square gives `(2y+x+1)² = 4x³+x²+18x−23 =
(x−1)(4x²+5x+23)`, `disc(4x²+5x+23) = −343 = −7³`.
- **Real period** `Ω⁺ = 4∫₀^∞ ds/√(4s⁴+13s²+32) = 1.981342`.
- **BSD, exact:** `L(14a,1)/Ω⁺ = 0.330224/1.981342 = 0.166667 = 1/6` — i.e. **`L(14a,1) = Ω⁺/6`** (the BSD
  formula `L(1) = Ω⁺·∏c_p·#Ш/#tors²` with `∏c_p·#Ш = 6`, `#tors² = 36`). This validates both `L(1)` and
  the period and pins the rank-0 arithmetic.
- **Imaginary period** `Ω₃ = 4∫₀^∞ ds/√(4s⁴−13s²+32) = 2.650982`; **area** `Ω⁺Ω₃ = 5.2525`; lattice
  modulus `τ = 1/2 + 0.66899 i` (rhombic, `Δ<0`).
- **Period field = `Q(√−7)` (the MEASURE column), exact:** `Δ(14a) = −21952 = −2⁶·7³`, squarefree part
  `−7`, so `Q(√Δ) = Q(√−7)`. Both period quartics `4s⁴±13s²+32` have `disc_{s²} = −343 = −7³`. **The apex-7
  (cubed) is the discriminant's odd part, so 14a's period geometry is defined over `Q(√−7)`** — even though
  14a is NOT CM (the branch points live in `Q(√−7)`; the modulus `τ` is transcendental). The measure-column
  Heegner field is literally the period field of the obstruction curve.

## The decisive negative: the floor constant is NOT a modular period
The floor's cusp-form contribution `0.344 − 0.304 = 0.040` was the target. Checked against EVERY natural
period quantity:
| quantity | value | | quantity | value |
|---|---|---|---|---|
| `Ω⁺` | 1.981 | | `Ω₃` | 2.651 |
| `1/Ω⁺` | 0.505 | | `Area` | 5.253 |
| `Ω⁺/2π` | 0.315 | | `1/Area` | 0.190 |
| `L(1)` | 0.330 | | `Area/4π²` | 0.133 |
| `L(1)/Ω⁺` | 0.167 | | `Im τ` | 0.669 |
**None equals `0.040`.** This is the **third negative in a row**:
> `floor ≠ L(14a,1)=0.330` · `floor ≠ L(sym²,s)≈1.1` · `floor ≠ any period of 14a`.
Three independent modular invariants, none is the floor constant. The conclusion is robust: **the floor
CONSTANT `0.040` is not a single modular invariant of `f₁₄`.**

## Where the constant actually lives (the redirect)
The modular form `f₁₄ = 14a` gives the **STRUCTURE** of LRC(14), not the **CONSTANT**:
- genus 1 = the hardness; rank 0 = the favorable sign; `GL(3) = sym² f₁₄` = the obstruction's level;
  `Q(√−7)` = the period field (measure column); `ℤ/6` torsion = the units (existence column);
  `a_p ≡ p+1 (6)` = the Hecke congruence. All STRUCTURE.
- The floor's numerical value (`~0.34`, cusp-part `~0.04`) is a **DESCENT quantity** — the product `∏_j ρ_j`
  over descent levels / the covering-min `n/Φ₆(n)` — a multi-level, set-dependent number, NOT a single
  period or L-value. klein's `114382/332563` is a scan inf (a specific config), reinforcing this: it is
  config/descent data, not a modular constant.
> **The redirect:** stop chasing the floor constant as a cusp-form period. The modular form is the right
> home for the STRUCTURE and the SIGN; the CONSTANT is determined in the descent (the `∏ρ_j` product, the
> covering-min `≥ 1/n`), which is exactly where the OPEN conjecture lives. The two are complementary: the
> curve says *why it's hard and that the sign is right*; the descent says *how big the floor is*.

## Status
- **Verified/clean (opus):** `Ω⁺=1.98134`, `L(14a,1)=Ω⁺/6` (BSD exact); `Ω₃=2.651`, `Area=5.2525`,
  `τ=1/2+0.669i`; period field `= Q(√−7) = Q(√Δ)`, `Δ=−2⁶·7³` (measure column, exact, no CM).
- **Negative (opus, decisive):** the floor constant `0.040` is NOT a period of 14a (nor `L(14a,1)` nor
  `L(sym²,s)`). The floor constant is not a single modular invariant.
- **Redirect:** the modular form gives structure + sign; the floor constant lives in the descent
  (`∏ρ_j`/covering-min). The open conjecture is a descent/Diophantine statement, not a modular-value one.

Related: modular-pushes-sym2-and-cusps (the prior two negatives), the-hecke-dictionary-of-f14, lrc14-lives-
on-14a, the-master-two-Heegner-columns (Q(√−7) measure column — now confirmed as 14a's period field),
per-level-vs-total-doublet + cusp-existence-comb-witness (the descent, where the constant lives); klein
THM-590/HYP-3598, mac-mini THM-580; OPEN-Q-108.
