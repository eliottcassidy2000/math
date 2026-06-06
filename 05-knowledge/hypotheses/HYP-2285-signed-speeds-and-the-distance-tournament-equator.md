# HYP-2285 — Two reframes: LRC directions (the σ gauge) and the unit-distance tournament (the equator)

**Session:** claudebox-2026-06-04-S631. **Prompt (user):** (1) LRC with runners able to move in the opposite
direction (negative speed); every other clockwise vs counterclockwise. (2) Unit distance as a tournament: length 1 =
edges, <1 = arrows one way, >1 the other; does this simplify the grid disproof? **Threads:** HYP-2185/2205 (σ/apex),
HYP-2235 (CM/unit-distance), HYP-2270/2275 (add/mult, Goldbach/Lemoine), HYP-2245 (H).

## Reframe 1 — LRC directions = the σ gauge (verified + formalized)
A runner moving the opposite way is the sign flip `v ↦ −v` = the antipodal `σ`. Since `‖vt‖ = ‖(−v)t‖`, **loneliness
is sign-invariant** (verified: `min_i ‖εᵢ vᵢ t‖` is the same for every sign vector ε, including "every other" and
random). So the **sign is a GAUGE** — loneliness depends only on `|speeds|`, the σ-quotient is the physical object.
"Every other runner counterclockwise" = the even/odd 2-coloring gauge (the 2-adic, S630). Formalized:
`Forbid_neg_speed` (direction reversal invariance).

**The trick (where direction bites): difference ↦ sum.** In a runner's relative frame the others move at the
*difference* `vᵢ − vₖ`; reversing a runner gives the *sum* `vᵢ + vₖ` (`‖(−vᵢ−vₖ)t‖ = ‖(vᵢ+vₖ)t‖`, formalized
`dZ_reverse_to_sum`). So **choosing directions selects between the difference-lattice and the sum-lattice of relative
speeds** — the resonance gauge, the additive `p±q` structure (Goldbach/Lemoine, S630). The signs are a gauge for the
gap but a real choice for which relation-lattice the resonance analysis sees.

## Reframe 2 — the distance tournament: unit distances = the equator (verified)
Orient each pair by its distance with a base order: `<1` → arrow class A, `>1` → class B, `=1` → **ties** (the unit
edges). Verified: the **triangular lattice is unit-MINIMAL** — NO pairs at distance `<1`, so class A is empty and the
"tournament" is degenerate (ties + one arrow class B). So **maximizing unit distances = maximizing the EQUATOR (the
tie-class)** of the distance order. The orientation `≷1` is "inside vs outside the unit circle" = the CM involution's
sign (`|β|<1` vs `|β|>1` vs `|β|=1`).

### Does this simplify the grid disproof? — CLARIFIES via CM-coherence (honest)
A unit distance is a **tie**: `|zᵢ − zⱼ| = 1`, i.e. `|β| = 1` for the lattice difference `β`. The **CM property**
(S623): `|β| = 1` in one embedding ⟺ in ALL embeddings. So the tie-set is the **unit group `𝒪_K^1`** — an
*algebraic* subvariety, not a metric accident. Hence the distance tournament's equator is **algebraically COHERENT
for the CM lattice** and **incoherent for a generic lattice**. That coherence IS the disproof's mechanism: the
class-field tower inflates the coherent tie-set (split primes ⟹ many norm-1 units, S623), giving `n^{1+ε}` ties vs the
lattice's `~6` ties/point. So the tournament reframe **clarifies** (does not trivialize) the disproof: it recasts
"max unit distances" as "max the algebraically-coherent equator," pinpointing WHY the grid (incoherent equator)
loses to the CM-tower (coherent equator). The `<1 / =1 / >1` trichotomy = the CM involution `|β|≷1`; the
disproof = engineering a fat, coherent `|β|=1` class.

## Synthesis
Both reframes impose the **tournament/σ orientation** on the geometry: Reframe 1 puts the σ-sign on the *speeds*
(a gauge for the gap, a real choice for the resonance lattice — difference↔sum); Reframe 2 puts the orientation on the
*distances* (`≷1` = the CM involution sign), with unit distances = the σ-fixed **equator**. The unit/tie/equator is
the same σ-fixed locus as the apex (HYP-2185) and the LRC `dZ=1/6` chord; the disproof = a coherent equator from the
CM/`σ` algebra.

## Formalized (math-lean, sorry-free) — `Math/LonelyRunner/SignedSpeed.lean`
`Forbid_neg_speed` (direction-reversal gauge invariance), `dZ_reverse_to_sum` (difference ↔ sum under reversal).

## Open
- The "every-other 2-coloring" gauge: does an alternating-direction relative frame expose a 2-periodic resonance
  structure for the AP / the tight configs?
- Formalize the distance-tournament equator = `𝒪_K^1` (the CM tie-set is algebraic) — the disproof's coherence as a
  theorem.
- The disproof restated purely in tournament-equator language (a coherence/antichain-maximization statement).
