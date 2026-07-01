# The covering-min is an adversarial facility-location game (potential = discrepancy, AP = the min-discrepancy equilibrium), and the √21 residual is an additive–multiplicative obstruction (21 is not a sum of two squares); the good LTC is PSL₂(7)

*opus-2026-07-01-S29. A multi-front prompt: aim the LTC at PSL₂(7) with √21 the POCS target; hunt the
non-abelian expander; use sum-of-two-squares / Koksma–Hlawka / Mahler–Popken / integer-complexity as
jumping-off points; and reframe the covering-min as an attacker–defender facility-location game with a
potential / price-of-anarchy argument for the inf measure. Each front produced a verified anchor.*

## The covering-min as an adversarial facility-location game
Cast it as attacker vs defender on the circle `S¹`:
- **Defender** picks the speed set `S` (the facilities/runners), constrained to be covering.
- **Attacker** (the lonely observer at 0) picks the time `t`; the runners sit at `{vt mod 1}`.
- **Payoff** = `min_v ‖vt‖` (the observer's distance to the nearest facility = the gap it sits in).
- `M(S) = max_t` (attacker's best response); **covering-min = `min_S max_t min_v ‖vt‖`** — a minimax.

**The potential is the discrepancy.** The attacker's payoff is the largest gap the observer can fall into, and
the gap structure is governed by the **star-discrepancy** `D*({vt})` of the runner arrangement (Koksma–Hlawka:
`|μ_r − (1−2r)^{n−1}| ≤` variation·`D*`). Verified: the **AP at `t=1/n` puts the runners exactly equally spaced**
— *minimum* discrepancy, the maximally uniform arrangement — so the observer's best gap is exactly the uniform
`1/n`. Any higher-discrepancy (more structured) config has a *bigger* gap somewhere, so the observer is lonelier
(`M > 1/n`). So:

> **The LRC/covering extremal is the minimum-discrepancy (uniform) equilibrium**, and `M(S) ≥ 1/n` is the
> potential-function statement that the most-uniform arrangement still leaves a `1/n` gap. The "inf measure"
> `μ_r(S)` is Koksma–Hlawka-close to the independent value `(1−2r)^{n−1} > 0` for low-discrepancy `S`; the danger
> is the high-discrepancy/structured (AP) config, and even there the gap is `≥ 1/n`.

This is the price-of-anarchy shape: the *equilibrium* (uniform AP, the "socially efficient" spread) minimizes the
attacker's advantage; the potential (discrepancy) that best-response dynamics minimize has its floor at the
uniform config. Importing the AGT machinery (CS 6840 / Hotelling) means: **bound the inf measure by a discrepancy
potential**, with the AP as the min-potential equilibrium.

## The √21 residual is an additive–multiplicative obstruction (sum of two squares)
The sharpest jumping-off point. By the two-squares theorem, `n = a²+b²` iff no prime `≡3 mod 4` divides `n` to an
odd power. Verified:
- **`21 = 3·7` is NOT a sum of two squares** (obstruction primes 3, 7 — both ≡3 mod 4).
- **`61` IS** (clean); `183 = 3·61` is not (obstruction 3); `403 = 13·31` is not (obstruction 31) — so the
  deep-well `D₁₄ = 21·403` (kps-S23) has obstruction primes `{3,7,31}`, all ≡3 mod 4.

So the **same obstruction is four things at once**:
```
   21 is not a sum of two squares   ⟺   3,7 both ≡3 mod 4   ⟺   narrow class ℤ/2 of ℚ(√21)   ⟺   ι-odd Gauss i√p
```
This is the cleanest statement yet: **the √21 residual is an additive–multiplicative incompatibility** — `21` is
multiplicatively `3·7` but has *no* additive `a²+b²` representation. That is precisely the Mahler–Popken /
integer-complexity flavor (the tension between `+` and `×`): the residual lives where the multiplicative
factorization `3·7` cannot be reconciled with the additive quadratic form. The clean primes (`≡1 mod 4`, sums of
two squares, e.g. `61`) carry no residual; the obstruction concentrates on the `≡3 mod 4` primes — the same ones
whose Gauss sums are ι-odd (`i√p`) and whose narrow class groups have `ℤ/2`.

## The non-abelian expander: PSL₂(7) is the group; √21 crosses its 21-Frobenius subgroup
Verified `PSL₂(7)`: order **168** (= |Aut(Fano)| = |Aut(Klein quartic)|), element orders `1,2,3,4,7` (with 48
order-7 elements — the heptagon rotations — and 56 order-3 elements — the Eisenstein multipliers). Its two
cuspidal 3-dim irreps carry the irrationality `(−1±√−7)/2` — **`√−7` is in its character field** (the ι-odd
heptagon certificate). And `|Aut(Paley₇)| = 21 = 7⋊3` is its **Frobenius/Borel subgroup**, so **`√21 = √(3·7)`
crosses the order-7 and order-3 parts** of `PSL₂(7)`. Cayley graphs of `PSL₂(p)` with LPS/quaternion generators
are **Ramanujan expanders**, so:

> **The good locally-testable √21-certificate should be a cocycle on the left–right Cayley complex of `PSL₂(7)`
> with LPS-Ramanujan generators** (an actual expander, unlike the abelian tiling cube, S28). `√21` = the
> nontrivial `ℤ/2` class (narrow-class / two-squares obstruction); **POCS/Kaczmarz** (alternating the two
> generating sets = mac-mini's pillar A) is the constructive method that converges to it; local testability =
> the certificate is checkable on the complex's squares.

## The other jumping-off points
- **Koksma–Hlawka** — the discrepancy is the game's potential (above); it is the analytic form of the
  far-element/equidistribution decomposition, with the AP the min-discrepancy extremal.
- **Mahler–Popken / integer complexity** — the additive–multiplicative `+/×` tension that the two-squares
  obstruction embodies; `21 = 3·7` (mult.) with no `a²+b²` (add.). A lens for *why* the composite `3·7` residual
  cannot be linearized/split (THM-503, non-Euler-product).
- **Pochhammer** — the repo's fiber fraction `f(n)=(½)_{n−2}/(n−2)!` is a Pochhammer (rising-factorial) ratio
  (the √π/Wallis structure); the ι-even (Eisenstein) "far" half of the certificate lives with these Γ/Pochhammer
  factors, complementary to the ι-odd Gauss half.

## Status
- **Verified:** `21` not a sum of two squares (obstruction 3,7); `61` is (clean); the four-way identity
  (two-squares = ≡3 mod 4 = narrow-ℤ/2 = ι-odd Gauss). `PSL₂(7)` order 168, element orders, `√−7` in its cuspidal
  reps, `21` = its Frobenius subgroup. AP at `t=1/n` = min-discrepancy, loneliness `= 1/n`.
- **Framings (routes):** covering-min = adversarial facility-location minimax, **potential = discrepancy**, AP =
  the min-discrepancy equilibrium, inf-measure bounded by Koksma–Hlawka; the √21 residual = an additive–
  multiplicative (two-squares/integer-complexity) obstruction; the good LTC = the **`PSL₂(7)` LPS-Ramanujan
  left–right Cayley complex**, √21 the class, POCS the constructor.
- **Honest:** the arithmetic (two-squares, PSL₂(7), the AP uniform gap) is exact and classical; the
  facility-location potential bound and the PSL₂(7)-LTC construction are pointed *routes* (they name the game,
  the potential, the group, the class, and the method), not yet a proof or a construction.

Related: HYP-3820 (the left–right square complex / PSL₂(7) sharpening), HYP-3819 (√21 = narrow-ℤ/2), HYP-3818
(the biquadratic bridge), HYP-3796/mac-mini (three pillars — POCS = the alternation), HYP-3786/HYP-3787
(far-element/equidistribution = the discrepancy potential), OPEN-Q-108. External: Annals 203-2 p.03 (LTC on the
non-abelian group); CS 6840 (facility-location/Hotelling, the game); two-squares & integer-complexity (the
obstruction). HYP-3821 (this). Script: 04-computation/facility_location_sum2sq_psl27_opus_20260701.py.
