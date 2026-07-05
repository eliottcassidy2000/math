# The Q50 conjecture is false at bound 50 — the height-independent correction

**kind-pasteur-2026-07-05-S11, HYP-4137.** Refutes, as literally stated, the Q50
conjecture (HYP-4119, mac-mini-S55) and `TemplateDichotomy` (HYP-4127, kps-S10).
Reroutes to the correct **height-independent (pinned-only) witness** object.
Everything below is independently verified (fresh implementations of every
predicate); the explicit counterexample is in
`04-computation/lrc_q50_refutation_kps_S11.py` (deterministic, seed-fixed).

---

## 1. The claim that is false

**Q50 / TemplateDichotomy (bound 50).** *Every primitive, compressed, covering
Fin-13 family (at its argmax peel) is either tight-shaped (`|v_i| = c·j`,
`j ∈ [1,12]`) or has a `2/25`-margin rational witness `t = k/s` with `s ≤ 50`
clearing all 12 base runners.*

This is the open predicate `TemplateDichotomy` on which `lrc14_of_template_and_corner`
rests, and mac-mini's Q50 ("every profile survivor witnesses at `q ≤ 50`",
claimed *scale-independent*, "no census-per-height").

## 2. The counterexample (explicit, verified)

Take the NEEDFREE census shape **F1 = {1,3,4,7,10,11,13,16,17,19,23,36}** (a B≤48
census survivor whose only `2/25`-witnesses are at the FREE moduli 37,43,47,49 —
no pinned-only witness `≤ 50`). Lift it to height `~10^22` preserving its residues
mod `L = lcm(2..25)`, choosing the free residues (mod 27,29,31,32,37,41,43,47,49)
by CRT so that **every free modulus in [26,50] is pinned** (no witness there).
Realize as a single-scale family: base multipliers `1..12` (ratio 12, exactly the
spread of the tight AP), plus an `istar ≡ 0 (mod 14)` argmax covering `q = 14`.

The resulting Fin-13 family passes **all** `TemplateDichotomy` hypotheses:

| hypothesis | status |
|---|---|
| all `v_i ≠ 0` | ✔ |
| `CoveringFamily` (`q | some v_i`, `q = 2..14`) | ✔ |
| compressed (`∀i ∃j≠i: |v_i| ≤ 13|v_j|`) | ✔ |
| primitive (`gcd = 1`) | ✔ |
| `istar` is the argmax | ✔ |
| NOT tight-shaped | ✔ |
| **NO `2/25`-witness at any `s ≤ 50`** | ✔ (min witness `s = 53`) |

So `TemplateDichotomy` is **false**, and `lrc14_of_template_and_corner` reduces
LRC(14) to a false hypothesis — a **dead reduction** (kernel-pure but unusable).

**LRC(14) is NOT threatened.** The family has a real witness at `t = 13/53`
(margin `≥ 2/25`), so it satisfies `TightLooseDichotomy`'s loose branch
(`∃ real tstar`). The genuine surface `lrc14_of_dichotomy_and_corner` is
UNAFFECTED — only the bounded-denominator "template" refinement fails.

## 3. Mechanism — MISTAKE-096 at the profile level

The profile pins residues only mod `q ≤ 25`. Split the witness moduli:

- **PINNED-ONLY** `q`: every maximal prime-power factor `≤ 25` ⟺ `q | L`.
  A witness there depends only on residues mod `L` = **HEIGHT-INDEPENDENT**.
  (In [26,50]: 26,28,30,33,34,35,36,38,39,40,42,44,45,46,48,50.)
- **FREE** `q`: has a prime-power factor `> 25` (27,32,49) or a prime `> 25`
  (29,31,37,41,43,47,53,…). Its witness depends on residues the profile does
  **not** control → **KILLABLE** by choosing the free residues (CRT lift).

Every free modulus in [26,50] is *individually* pinnable with 12 profile-shape
runners (verified); being CRT-independent, they are pinnable *simultaneously*.
A **NEEDFREE** shape — no pinned-only witness `≤ 50` — therefore lifts to a family
with no witness `≤ 50` at all. The census's `HARD = 0` is only evidence for
**bounded** height; it never enumerates the `10^22` free-pinned family. mac-mini's
high-scale test (12,095 lifts to `1.16e14`, first-q `≤ 35`) used **diagonal**
lifts (a uniform shift `v_i ↦ v_i + mL`), which do NOT pin the free residues —
so it systematically missed this family. *(Exactly the recurring trap: an
empirical bound is only as uniform as its generator; the aligned/CRT adversary
was not in the sample. Cf. MISTAKE-095/096/098/101/102.)*

## 4. The pinned-only "repair" is ALSO dead (Q₀ = ∞) — the honest conclusion

At first the fix looked clean: every pinned-only `q | L`, so a pinned-only
witness is height-independent, and over all **511,947** B=48 census survivors the
smallest pinned-only witness is `≤ 69` (max 69 at `{1,3,4,5,7,8,9,11,12,19,23,30}`;
F1's exact bound is 55). That suggested "pinned-only witness at `q ≤ Q₀`, `Q₀ ≈ 69`".

**But the census is bounded height, and high-height shapes break it.** A hill-climb
over residue vectors (magnitudes free at high height) finds shapes whose smallest
pinned-only witness is ≥ 156, achieved by **highly-composite runners** (e.g.
`5040 = 2⁴·3²·5·7`, `13650 = 2·3·5²·7·13`) that are `≡ 0 (mod q)` for many
pinned-only `q` and thereby **block** the witness there. Taken to the limit:

> **A single runner `≡ 0 (mod L)` blocks EVERY pinned-only witness.** Since every
> pinned-only `q | L`, that runner is `≡ 0 (mod q)` — always in the danger band —
> so no dilation clears it. `Q₀ = ∞` for that shape.

**Verified** (`04-computation/lrc_q50_refutation_kps_S11.py` sibling checks): the
family `base = {L, 1,3,4,5,7,8,9,11,13,17,19}`, `istar = 2L` satisfies EVERY
`TemplateDichotomy` hypothesis (covering 2..14, compressed, primitive, not tight)
and has **no pinned-only witness at any `q | L`** (scanned `≤ 3000`, all blocked by
the `L`-runner). It is loose only via a FREE witness (`q = 37` here) — height-dependent.

**Both mechanisms combine → the witness modulus is UNBOUNDED.** An `L`-runner makes
`F5`-pinning and covering *vacuous* (`has0` at every `q ≤ 25`, `L` covers 2..14), so
the other 11 runners are unconstrained by the profile. For ANY target `N`: pinned-only
`q ≤ N` are blocked by `L`; free `q ≤ N` are pinned by choosing the 11 runners' residues
(each free modulus is pinnable with 11 free runners — verified). Hence a valid non-tight
covering compressed primitive family with no witness at any `q ≤ N`. `N` arbitrary ⟹
**no fixed modulus bound (pinned-only, free, or any fixed set) closes the loose branch.**

## 5. The only viable surface: real-valued `TightLooseDichotomy`

The loose branch must use a REAL witness (`∃ tstar ∈ ℝ`, as in `TightLooseDichotomy`
and `lrc14_of_dichotomy_and_corner`) or a HEIGHT-DEPENDENT bound `q ≤ Q(height) ~ c·ln(height)`
(the covering-census / analytic-descent route — exactly the two-sided architecture
MISTAKE-096 already forced). The bounded-denominator TEMPLATE program (Q50 at any fixed
modulus set) is a dead end. LRC(14) is untouched: every counterexample here is LOOSE
(has a real witness); only the fixed-bound refinement fails.

## 6. A clean by-product: the parity-split reduction

For `q = 2p` (p odd), a mod-`2p` `2/25`-witness reduces **exactly** to a
two-color avoidance problem mod `p`: `∃ b ∈ (Z/p)*` with every EVEN base runner
avoiding `E_p = {0, ±2}` and every ODD runner avoiding `O_p = {±1}` (μ=3) or
`{±1, ±3}` (μ=4). `E_p ⊔ O_p` = the width-`μ_{2p}` band mod `p`, split by the
parity of the `2p`-lift. Verified exact (100,000 tests, 0 mismatches). Pure
integer arithmetic — Lean-able like the `rational_point_margin` atom. (The
parity dividend: `0 ∈ E_p` always ⟹ no even runner may be `≡ 0 mod p`.)

## 7. Impact / handoffs

- **HYP-4127 (kps-S10)**: `TemplateDichotomy` / `lrc14_of_template_and_corner`
  rests on a FALSE hypothesis — a dead reduction. Neither the `s ≤ 50` form nor the
  pinned-only `q ≤ Q₀` form is provable (Q₀ = ∞). The Lean is kernel-pure and TRUE as
  an implication but can never be discharged. Superseded by the real-valued surface.
- **HYP-4119 (mac-mini-S55)**: Q50 at bound 50 refuted; the "absolute-height /
  no-census-per-height" claim is wrong. The pole-necessity/periodicity result STANDS
  (it is about the profile FILTERS, which are indeed CRT-ray-periodic). But note the
  companion fact: the WITNESS side inherits no such height-independence, because a
  composite `≡ 0 (mod L)` runner blocks all pinned-only witnesses — the census's
  `HARD = 0` relies on FREE (height-dependent) witnesses it never exhibits.
- **Proof map item 1d (mac-mini-S54)**: "finite template at a fixed bound" (any
  flavor) is DEAD. Reroute to the real-witness `TightLooseDichotomy` (already the
  main surface) or the height-dependent covering-census `q ≤ Q(height) ~ c·ln(height)`
  (the MISTAKE-096 two-sided architecture: bounded-magnitude finite census + analytic
  large-magnitude). The "single absolute-height mechanism" the crux asks for does not
  exist as a fixed modulus bound.
- **The real surface is fine.** `lrc14_of_dichotomy_and_corner` (loose = ∃ real tstar)
  is untouched; every counterexample here has a real witness. The lane to drop is the
  bounded-denominator template refinement, not the dichotomy itself.
- ⚠️ **Self-correction**: my initial S11 broadcast suggested repairing at pinned-only
  bound 69 — that is ALSO wrong (Q₀ = ∞ via composite blocking). The honest verdict is
  the real surface, not any fixed template.
