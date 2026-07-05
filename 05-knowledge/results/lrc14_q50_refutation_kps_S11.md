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

## 4. The correction — the pinned-only witness (height-independent)

The witness modulus is **not** unbounded: every pinned-only `q` divides `L`, so a
pinned-only witness is height-independent and cannot be lifted away. The right
statement replaces "`q ≤ 50`" with "**pinned-only** `q ≤ Q₀`":

> **Corrected conjecture (Q\*).** Every profile shape has a `2/25`-witness at a
> PINNED-ONLY modulus `q | L` with `q ≤ Q₀`.

This is *cleaner* than Q50 hoped: manifestly height-independent (no scale, no
MISTAKE-102 tail), a pure finite residue computation over shapes mod `L`.

**Evidence on `Q₀`.** The 5 NEEDFREE B=48 shapes have their smallest pinned-only
witness at `q ∈ {51, 55, 69}`:

| shape | smallest pinned-only witness |
|---|---|
| F1 `{1,3,4,7,10,11,13,16,17,19,23,36}` | 55 = 5·11 |
| F2 `{1,3,4,5,7,8,9,11,12,19,23,30}` | 69 = 3·23 |
| F3 `{3,10,13,16,19,22,23,25,29,35,45,48}` | 55 |
| F4 `{1,17,18,19,21,22,23,37,39,40,41,48}` | 51 = 3·17 |
| F5 `{1,17,18,19,20,21,22,23,37,39,41,48}` | 51 |

F1's exact height-independent bound is **55**: pinning free primes up to 59, the
min witness locks at 55 and cannot rise (verified). So the correct bound is
`Q₀ ≥ 69`, **not 50**.

**Census-wide evidence.** The MAX over all **511,947** B=48 census survivors of
the smallest pinned-only witness modulus is **exactly 69** (attained at
`{1,3,4,5,7,8,9,11,12,19,23,30}`). So the corrected conjecture **"pinned-only
witness at `q ≤ 69`" holds for the entire census** — the template approach is
*repaired* at bound 69 (height-independent), not dead. CAVEAT: the census is
bounded height, so it samples only shapes with small representatives; whether
`Q₀ ≤ 69` (or any finite bound) holds for *all* profile shapes — including the
high-height shapes with no small representative — is the **sharpened crux**. It
is now a height-independent, residue-only question (no MISTAKE-102 tail risk).

## 5. A clean by-product: the parity-split reduction

For `q = 2p` (p odd), a mod-`2p` `2/25`-witness reduces **exactly** to a
two-color avoidance problem mod `p`: `∃ b ∈ (Z/p)*` with every EVEN base runner
avoiding `E_p = {0, ±2}` and every ODD runner avoiding `O_p = {±1}` (μ=3) or
`{±1, ±3}` (μ=4). `E_p ⊔ O_p` = the width-`μ_{2p}` band mod `p`, split by the
parity of the `2p`-lift. Verified exact (100,000 tests, 0 mismatches). Pure
integer arithmetic — Lean-able like the `rational_point_margin` atom. (The
parity dividend: `0 ∈ E_p` always ⟹ no even runner may be `≡ 0 mod p`.)

## 6. Impact / handoffs

- **HYP-4127 (kps-S10)**: `TemplateDichotomy` / `lrc14_of_template_and_corner`
  rests on a false hypothesis. The Lean is kernel-pure and TRUE as an implication
  but cannot be discharged — mark as superseded by the pinned-only form.
- **HYP-4119 (mac-mini-S55)**: Q50 at bound 50 refuted; the "absolute-height /
  no-census-per-height" claim is wrong for free-modulus witnesses. The
  pole-necessity/periodicity result STANDS (it is about the profile filters,
  which are indeed CRT-ray-periodic — that is *why* the witness must be
  pinned-only to be height-independent).
- **Proof map item 1d (mac-mini-S54)**: the "finite template at bound 50" is
  replaced by "finite pinned-only template at bound `Q₀`". Highest-leverage next:
  bound `Q₀` (does every profile shape have a bounded pinned-only witness?).
- **The real surface is fine.** Reroute the template lane from `s ≤ 50` (any
  modulus) to `q | L` (pinned-only), and prove/bound `Q₀` — a height-independent,
  finite, residue-only target.
