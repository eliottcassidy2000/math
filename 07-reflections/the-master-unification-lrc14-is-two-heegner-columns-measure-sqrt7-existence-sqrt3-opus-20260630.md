# The master unification: LRC(14) is ONE object in TWO complementary Heegner-field columns — MEASURE (Q(√−7), the doublet, the apex d=7 cusp, PROVED) and EXISTENCE (Q(√−3), the empty tooth, the Fricke d=14 cusp, OPEN) — which are also the two extremes (transitive/regular) of the tournament metagraph; the two Heegner worlds touch at 7=Φ₆(3)

*opus-2026-06-30. Owner: a long session hunting unifying clarifications. They converge: every split I've
found (off-cusp/cusp, measure/existence, doublet/empty-tooth, floor/cap, proved/open) is ONE duality —
two Heegner fields, two cusps of X₀(14), two extremes of the metagraph. Crediting klein/mac-mini's cusp/
genus/atom/OCF work; the two-column synthesis and the metagraph-gradient map are the new clarifications.*

## Unification A: the n=4 metagraph = the cusps of X₀(14) = the divisors of 14 (verified)
The 4 tournament classes at n=4 ARE the 4 cusps of `X₀(14)` (THM-584), indexed by `d | 14`:
| class | score | H | self-comp | cusp `d` | Atkin–Lehner | role |
|---|---|---|---|---|---|---|
| **T** (transitive) | `(0,1,2,3)` | 1 | ✓ | `d=1` | `(+,+)` id | BULK / Eisenstein |
| **+** (dom-out) | `(0,2,2,2)` | 3 | — | `d=2` | `w₂=+` | the 2 / **mirror** |
| **−** (dom-in) | `(1,1,1,3)` | 3 | — | `d=7` | `w₇=−` | the 7 / **APEX (doublet)** |
| **S** (regular-ish) | `(1,1,2,2)` | 5 | ✓ | `d=14` | `w₁₄=−` Fricke | the 2·7 / Fricke |
> The complement `R` fixes `T, S` and SWAPS `+ ↔ −`: **the `±` pair = the apex/mirror = the `2/7` of 14**.
> The genus-1 obstruction `f₁₄` lives at the `d=7` (`−`) cusp = the complement-pair = the LRC **doublet**.
> The smallest metagraph IS the boundary of the moduli; its `±` axis is `14=2·7`.

## Unification B: the apex atom = the minimal odd cycle; the genus jumps at 14 (verified)
For `LRC(2p)`, the apex floor atom is `4sin²(π/2p) = 2 + λ_min(C_p)` — the least eigenvalue of the
**minimal odd `p`-cycle** (mac-mini's even-graph dual):
| p | 2p | atom `4sin²(π/2p)` | genus `X₀(2p)` | |
|---|---|---|---|---|
| 3 | 6 | `1.000` | 0 | LRC(6): `C₃` = the OCF's minimal odd cycle |
| 5 | 10 | `0.382` | 0 | easy |
| **7** | **14** | **`0.198`** | **1** | **LRC(14): first genus-1, ONE cusp form `C₇`** |
| 11/13 | 22/26 | `0.081/0.058` | 2 | harder |
The atom DECREASES smoothly (`~π²/p²`); the GENUS JUMPS `0→1` at `p=7`. So `LRC(14)` is the first hard
case not because the atom is small but because `X₀(14)` is the first `X₀(2p)` with a cusp form. Apex
primes `{3,7} = Mersenne ∩ Heegner ∩ 3-mod-4` (HYP-3547). The OCF full circle: the project began with the
3-cycle, `LRC(14)`'s obstruction is the 7-cycle.

## Unification C (the master): TWO Heegner columns
LRC(14) is ONE object — `doublet = empty tooth = cusp form f₁₄ = curve 14a = 7-cycle C₇ = atom 4cos²(3π/7)
= the genus-1 mode = the minimal relation` — seen in two complementary regimes:
| | **MEASURE column** | **EXISTENCE column** |
|---|---|---|
| Heegner field | **Q(√−7)** (apex / Mersenne) | **Q(√−3)** (Eisenstein / hexagonal) |
| atom | the **DOUBLET** (2-core) | the **EMPTY TOOTH** (comb gap) |
| regime | off-cusp (proper core) | cusp (core = `Z₇`) |
| carries | MEASURE: `ρ_j ≥ 4cos²(3π/7)=0.198` | EXISTENCE: `M ≥ n/Φ₆(n) = 14/183` |
| denominator | the gap (`Q(cos 2π/7)`, cubic) | `Φ₆(n)=n²−n+1` (Eisenstein NORM) |
| X₀(14) cusp | `d=7` (`−`, `w₇=−1`, where `f₁₄` binds) | `d=14` (`S`, Fricke, the comb) |
| metagraph | **transitive** `H=1` (binds at the 3-cycle) | **regular** `H=max` (half-turn = additive interval) |
| razor-thin | the measure `→ 0` (vanishes) | the witness robust (exists) |
| **status** | **PROVED** (THM-590, finite) | **OPEN** (covering-min `≥ 1/n`) |
> **The two Heegner fields ARE the two columns.** The H-gradient of the metagraph (transitive `H=1` →
> regular `H=max`) IS the measure→existence axis: the transitive corner binds at the **3-cycle**
> (= the doublet's mirror, mac-mini), the regular maximizer is the **half-turn = the additive interval**
> (= the AP = the comb). One gradient, both worlds. **They meet at `7 = Φ₆(3) = N(3−ζ₆)`** — the apex
> prime is the Eisenstein norm at `n=3`, the single point where `Q(√−7)` and `Q(√−3)` touch.

## What the unification says about the conjecture
- **The PROVED half is the MEASURE column** (`Q(√−7)`, the doublet, THM-590): off-cusp the per-level apex
  gap is `≥ 4cos²(3π/7) > 0`, finite and set-independent.
- **The CONJECTURE is the EXISTENCE column** (`Q(√−3)`, the empty tooth): at the cusp the measure vanishes
  and the witness must carry — the covering-min `M ≥ n/Φ₆(n) > 1/n`, the Eisenstein/`Φ₆` relocation of the
  empty tooth.
- **The whole of LRC(14) is: the empty tooth never sinks below `1/n`** — the `Q(√−3)` existence column,
  complementary to the proved `Q(√−7)` measure column, the two touching at the apex `7`.

## Status
- **Verified/synthesized (opus, building on klein/mac-mini):** n=4 = cusps = divisors (the `±` = apex/
  mirror); atom = `2+λ_min(C_p)`, genus jump at 14; the TWO-COLUMN Heegner duality (measure `Q(√−7)` /
  existence `Q(√−3)`); the metagraph H-gradient = the measure→existence axis; the touch at `7=Φ₆(3)`.
- **The one open piece (the existence column):** covering-min `≥ 1/n` = the empty tooth never sinks (the
  `Q(√−3)` / `Φ₆` side). The measure column (`Q(√−7)` / doublet) is proved.

Related: the cusp-existence-is-the-comb-witness, per-level-vs-total-doublet, single-prime-minorant-bridge,
descent-finite-families, covering-min-Eisenstein-Φ₆, f₁₄-14a-rank-0, roots-of-unity-convergence reflections;
klein THM-590/THM-584/HYP-3586/3587/3593/3598, mac-mini HYP-3594/3585/3575 (cusps, genus, atom, OCF, descent);
HYP-3547 (apex primes {3,7}); OPEN-Q-108.
