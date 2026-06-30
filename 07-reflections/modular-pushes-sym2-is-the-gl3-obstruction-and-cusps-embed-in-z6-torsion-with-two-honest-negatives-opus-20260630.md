# Two modular/elliptic pushes on LRC(14), with two honest negatives: (1) the floor's 2nd-moment obstruction is sym² f₁₄ = a degree-3 GL(3) L-function = the Littlewood/SL(3) level (structural ✓; but L(sym²,s)≈1.1, NOT the floor's 0.040 — not a numerical floor constant); (2) the n=4 metagraph embeds in E_tors(14a)=ℤ/6 (Manin–Drinfeld ✓; the transitive class T is a GENERATOR, order 6) — but the complement R is NOT the curve negation [−1] (it swaps an order-3 and an order-2 cusp), so "SC=2-torsion/pair=3-torsion" is FALSE

*opus-2026-06-30. Owner: chase sym² f₁₄ and pin cusps→torsion. Both pushed; both yielded a real structural
fact AND a caught over-reach. Recording the negatives so they don't recur.*

## Push 1: the floor's obstruction is sym² f₁₄ (GL(3)) — structural ✓, numerical ✗
The LRC(14) floor is the 2nd MOMENT (pair correlation) of the danger count. The pair correlation of a
Hecke form factorizes:
> `L(f₁₄ × f₁₄, s) = ζ(s) · L(\mathrm{sym}^2 f₁₄, s)` — Eisenstein bulk = the `ζ(s)` piece, the genus-1
> obstruction = the `\mathrm{sym}^2 f₁₄` piece.
**`\mathrm{sym}^2 f₁₄` is a degree-3 `GL(3)` L-function** (`b_p = a_p^2 − p`). So the obstruction sits one
`GL`-rank ABOVE the curve (`GL(2)→GL(3)`), exactly the **`SL(3)` / Littlewood level** the Siegel–Rogers /
Han–Lee moment-hierarchy thread already named. *Two independent threads (modular, Diophantine-moment) place
the LRC(14) obstruction on `GL(3)`.* — a genuine cross-thread unification.
> **HONEST NEGATIVE.** Computed `L(\mathrm{sym}^2 f₁₄, s) ≈ 1.17, 1.13, 1.10` at `s = 2.5, 3, 3.5` (Euler
> product). The floor's cusp-form contribution is `0.344 − 0.304 = 0.040`. **`L(sym²,s) ≈ 1.1 ≠ 0.040`** at
> any natural point — the floor constant is NOT a clean `\mathrm{sym}^2` L-value (just as it is not `L(14a,1)
> =0.330`). The `GL(3)` identification is STRUCTURAL/heuristic (the right *home* and *rank* of the
> obstruction), not a numerical formula for the floor. The constant remains a cusp-form period, unidentified.

## Push 2: the n=4 metagraph embeds in E_tors(14a) = ℤ/6 — embedding ✓, R=[−1] ✗
The cuspidal group of `X₀(14)` is `ℤ/6` = the full rational torsion of `14a` (Manin–Drinfeld + Ligozat
eta-quotient computation, Smith form gives the 6). So the **4 cusps = the n=4 tournament classes embed in
`E_tors(14a)=ℤ/6`**. Exact orders (base `S=d{=}14`):
| class (cusp) | order in ℤ/6 | as element (T = generator 1) |
|---|---|---|
| **T** (transitive, `d=1`) | **6 — a GENERATOR** | `1` |
| **+** (dom-out, `d=2`) | 3 | `2` |
| **−** (dom-in, `d=7`) | 2 | `3` |
| **S** (regular, `d=14`) | 1 (base) | `0` |
> **The TRANSITIVE class generates the torsion** (order 6) — the "simplest" tournament is a cyclic
> generator of `E_tors(14a)`. The 4 cusps `{0,1,2,3}` **sum to 0** (a cuspidal relation). `6 = φ(14)` = the
> units = the existence column.
> **HONEST NEGATIVE (corrects my prior script narrative).** I claimed the complement `R` acts as the curve
> negation `[−1]`, giving "SC classes = 2-torsion, complement-pair = 3-torsion." **FALSE.** `R` swaps `+`
> (order 3) and `−` (order 2) — elements of DIFFERENT order — which no group automorphism can do. `[−1]`
> fixes `{0,3}` and swaps `{1,5},{2,4}`; `R` fixes `{T=1, S=0}` and swaps `{+=2, −=3}`. So `R ≠ [−1]`. The
> complement is a curve involution fixing the `T–S` (order-6) axis and flipping the `+/−` pair (an affine
> reflection), NOT the group negation. The clean "complement = negation" picture does not hold.

## What survives (the solid modular/elliptic facts)
- **`X₀(14) = 14a`** (genus 1) is the LRC(14) moduli; the obstruction is the cusp form `f₁₄`.
- **The Hecke data spells the structure:** `a₂=−1, a₇=+1` (the `2·7`, Atkin–Lehner mirror/apex);
  `a_p ≡ p+1 \pmod 6` (the ℤ/6 torsion = the units = existence "6"); genus 1 (obstruction); rank 0,
  `L(14a,1)=0.330>0` (non-degenerate).
- **Hardness = genus:** `LRC(2p)` is rational-moduli/solved (`genus 0`) until `14`, the first genus-1.
- **The obstruction lives on `GL(3) = \mathrm{sym}^2 f₁₄`** (the Littlewood/`SL(3)` level) — structural.
- **The n=4 metagraph ↪ `ℤ/6` torsion**, transitive = generator, cusps sum to 0.

## What did NOT pan out (recorded so it doesn't recur)
1. **Floor ≠ a clean L-value.** Neither `L(14a,1)=0.330` nor `L(\mathrm{sym}^2,s)≈1.1` equals the floor
   `0.344` or its cusp-part `0.040`. The L-values give the SIGN (rank 0 = non-degenerate), not the CONSTANT.
2. **Complement ≠ negation.** `R` is not `[−1]` on the torsion (swaps unequal orders). No "SC=2-torsion"
   identity.

## Status
- **Computed/solid (opus):** cuspidal group `X₀(14)=ℤ/6`, exact cusp orders (T generator/6, + 3, − 2, S
  base), cusps sum to 0; `\mathrm{sym}^2 f₁₄` Euler-product values; the structural `GL(3)`/Littlewood match.
- **Honest negatives (opus):** floor is not `L(sym²,s)` nor `L(14a,1)`; complement `R ≠ [−1]`.
- **Open/direction:** identify the floor's `0.040` cusp-form constant as a *period* of `f₁₄` / `sym² f₁₄`
  (not an L-value); the `GL(3)` home is right, the exact period is the remaining target.

Related: the-hecke-dictionary-of-f14, lrc14-lives-on-14a, the-master-two-Heegner-columns,
the-siegel-rogers-moment-hierarchy (the SL(3) level); klein THM-584/HYP-3587 (cusps=classes, genus),
mac-mini HYP-3594; OPEN-Q-108.
