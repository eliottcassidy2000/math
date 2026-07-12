---
source: opus-2026-07-11-S232
status: A PICTURE + two PROVED formulas. The remaining hard part of the clean-ruler route, read in the
  summand-graph / multiplicand-graph frame at the natural modulus C = 2n−1 = 27. The clean ruler at C=27
  has an EXACT criterion (verified 0/12000): SHALLOW = the summand shells are not overloaded (additive),
  LIVE = a multiplicative unit-clock or 3-adic clock finds an empty danger-shell (multiplicative). The
  route certifies everything EXCEPT the multiplicand-maximal (AP-coherent) wall — the same hard core every
  other route stalls on.
tags:
  - lrc14
  - summand-graph
  - multiplicand-graph
  - anti-concentration
  - clean-ruler
  - THM-401
  - the-wall
  - visibility-principle
---

# The hard part is the multiplicand-maximal wall — summand shells vs. unit clocks

**opus-2026-07-11-S232.** Owner: "get an actual picture of the remaining hard part; use that deep
understanding to crack statements; look back at the summand graph and the multiplicand graph." Doing so
resolves the S230/S231 "bounded-modulus anti-concentration" gap into its two proper halves and pins the
hard core to a single named family class — with two exact, proved formulas at the **natural** modulus.

## The frame (THM-401 + HYP-2083/S571): the summand graph at C = 2n−1

The clean-ruler route (kps THM-707) asks each residual for a modulus `q` with `liveCount(q) ≥ 1` and
`maxBand(v,q) ≤ 5`. THM-401 (PROVED) says the **pinch/Farey-companion modulus of the floor `1/14` is
`C = 2n−1 = 27`**: `1/14` and `2/27` are Farey-adjacent, and the nonzero residues mod 27 partition into the
**antipodal summand shells** `P_a = {a, 27−a}` (a = 1..13). The danger arc at `q = 27` is
`{0, 1, 26} = {0} ∪ P₁`. The **summand graph** supplies these additive shells; the **multiplicand graph**
enters as the unit group `(ℤ/27)*`, which acts on them by inverse clocks `k = a⁻¹` (HYP-2083/S571). Since
`27 = 3³`, that action has three strata: **9 unit shells** `{1,2,4,5,7,8,10,11,13}`, **3 gcd-3 shells**
`{3,24},{6,21},{12,15}`, **1 gcd-9 shell** `{9,18}` — the 3-adic part is the "blind spot" the inverse
clocks can't reach.

## The crack: two proved formulas (verified 0 / 12000 each)

Writing `z = #{i : 27 | vᵢ}`, `z₃ = #{i : 3 | vᵢ}`, `z₉ = #{i : 9 | vᵢ}`, and
`maxUSL = max_a #{i : vᵢ ∈ unit-shell P_a}`, the multiplier structure at `q = 27` is exactly:

- `p` a **unit**: `bandCount(v,27,p) = z + (load of unit-shell P_{p⁻¹})` — the units *permute* the 9 unit
  shells (S230's `x ↦ px` bijection, now at a composite modulus but on the unit part);
- `p` with `gcd(p,27)=3`: `bandCount = z₉` (danger `∩ 3ℤ = {0}`, and `27 | 3u·vᵢ ⟺ 9 | vᵢ`);
- `p` with `gcd(p,27)=9`: `bandCount = z₃`.

Reading off the max and the zero set gives the two formulas:

> **(1)  `maxBand(v,27) = max( z + maxUSL , z₃ )`.**
>
> **(2)  `q = 27` is a clean ruler for `v`  ⟺  `maxBand(v,27) ≤ 5`  ∧  live, where
> `live ⟺ [ z = 0 ∧ some unit-shell empty ] ∨ [ z₉ = 0 ] ∨ [ z₃ = 0 ]`.**

Both verified with **0 mismatches over 12000** random 13-families (the live half also 0/8000 separately).
This makes HYP-2935's slogan **"addition creates the shells, multiplication tests their visibility"** into a
*formula*: **SHALLOW** (`maxBand ≤ 5`) is the **summand** condition — no shell overloaded and few multiples
of 3; **LIVE** is the **multiplicand** condition — a unit clock, or a 3-adic clock, exposes an empty
danger-shell. This is the additive/multiplicative split of THM-401(v) at the level of the clean ruler, and
it *generalizes the S230 prime criterion* (band `{0,±1}`, all-nonzero shells) to the natural composite
modulus with its 3-adic strata.

## The wall: why the route stops exactly at the AP

By formula (2), `clean-at-27` fails in exactly two ways:

- **(a) SHALLOW-FAIL** — `z₃ ≥ 6` (or a unit-shell carries ≥ 6): a **mod-3 / mod-27 concentration**. These
  are the dilated/coarse-reducible families (`7·{1..13}` and kin) already peeled by the dilation / detuned
  branches — not genuine residuals.
- **(b) LIVE-FAIL** (shallow but not live) — the family **fills all 9 unit shells** *and* has a multiple of
  3 *and* a multiple of 9. This is precisely **multiplicand-maximal** (S560o: "the AP contains a multiple of
  every `q ≤ 13`, escaping only at `q = 14`"). Verified: the **AP `{1..13}`, GW `12→24`, and the sporadic
  `V* = {1..11,13,24}`** are all live-fail — each hits **9/9** unit shells, has `z₃ = 4`, `z₉ = 1`, so
  `maxBand = 4` (shallow) but **not live**. These are exactly the tight (`M = 1/14`) families.

So the clean-ruler route at the natural modulus certifies **everything except the multiplicand-maximal
wall**, and it fails there *only on the live (multiplicative-visibility) half*, never on shallow. The AP
saturates all three visibility routes at once: no empty unit shell (fills them), `z₉ ≥ 1`, `z₃ ≥ 1`. **That
saturation is the wall** — the exact, structural reason the clean-ruler route cannot self-certify the AP.

## The synthesis: all roads lead to the AP

The live-fail (multiplicand-maximal) hard core is the **same family class** every other route stalls on:

- **Covering-moment route** (HYP-3085 / OPEN-Q-108): the crux is that *consec maximizes the pairwise
  co-emptiness `S₂`* — a concentration-extremality. Same extremizer.
- **Additive-energy route** (THM-656 / THM-660): the variance floor is monotone in the reduced energy `R2`,
  minimized (worst) at the AP. Same extremizer.
- **Inverse-additive characterization** (opus-S181, definitive): *tight ⟺ the relation lattice `Λ` is
  1-dimensional AND coherent (a dilate of `{1..13}`)*. Same class.

The summand/multiplicand frame explains *why* they coincide: **multiplicand-maximal = hits every small
modulus = the AP's defining extremal property = 1-dimensional coherent additive structure.** The AP is the
**joint extremum** — summand-minimal (Freiman `|S +̂ S| = 2·13−3 = 23`, fewest pinch denominators) and
multiplicand-maximal (covers every `q ≤ 13`) — and both extremalities are *the same wall* viewed additively
and multiplicatively (the log bridge: the multiplicand graph is the summand graph in log-coordinates).

## What this cracks, and what remains

**Cracked (proved, exact):** the clean ruler at the natural modulus `C = 2n−1` has an explicit
criterion (formulas 1–2); the additive/multiplicative split is a formula, not a slogan; the failure set is
pinned to two named types, one dispatched (mod-3 concentration) and one the AP-coherent wall; and the wall
is characterized as *multiplicand-maximal*, unifying the clean-ruler, covering-moment, and energy routes'
hard cores into a single family class.

**Remains (the genuine hard core, now sharply named):** the multiplicand-maximal / AP-coherent families are
not self-certified by any bounded-modulus clean ruler — their loneliness is the tight-floor content `M = 1/14`
itself, discharged only by the exact census / LRC≤13 structure (not an anti-concentration lemma). This is the
honest terminus: **the clean-ruler route reduces `hB5` to the multiplicand-maximal wall**, which is exactly
where the problem's difficulty has always lived. The "bounded-modulus anti-concentration lemma" of S230/S231
is now resolved: its *shallow* half is the summand condition (true off the mod-3-concentrated dilates), and
its *live* half is the multiplicand condition, which **provably cannot hold at any modulus for the AP** — so
the AP is not an anti-concentration failure to be bounded away, but the fixed point the whole architecture
must hand to the tight-floor census.

→ THM-401 (the modulus `2n−1` + shells, PROVED), HYP-2083/S571 (the unit-clock strata at `C=27`), HYP-2935
(the visibility slogan — now a formula), THM-361 (product-sum defect: `1` is additive slack, the +/× weld),
S560o (AP = joint extremum), opus-S230/S231 (the bounded-modulus route this resolves), opus-S181 (tight ⟺
1-dim coherent), HYP-3085 / OPEN-Q-108 (the covering-route crux — same wall), THM-656/660 (energy floors —
same wall). Files: `lrc14_summand_shell_visibility_opus_S232.py` (+`.out`).
