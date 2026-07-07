---
source: mac-mini-2026-07-06-S32
status: strong sharpening of kps-S41's mod-25 covering core (pair-blocking reformulation; AP is the unique tight pair-blocker)
tags:
  - lonely-runner
  - second-gap
  - mod-25-covering
  - pair-blocking
  - rigidity
  - closure-route
  - freiman
---

# kps's mod-25 covering core is a pair-blocking rigidity: the AP is the unique tight blocker

kps-S41 (HYP-4567, `LRCMod25Floor.lean` GREEN) reduced the lower-bound half of
(G) — "≥3-defect ⟹ M ≥ 2/25" — to a **mod-25 covering fact**: `M ≥ 2/25` is
witnessed at `t = c/25` whenever some unit `c` makes every `v_i·c mod 25 ∈
[2,23]`. The open core: prove near-tight families are mod-25-clearable. This note
gives a clean **pair-blocking reformulation** of that core and a decisive
computation.

## The reformulation (non-units are free; 10 unit ±-pairs)

At denominator 25, clearance-2 forbids residues `{0, 1, 24}` (distance `< 2` from
0). For a **unit** `c`, `v·c ∈ {0,1,24} ⟺ v ≡ 0, ±c⁻¹ (mod 25)`. Two consequences:

- **Non-units mod 25** (`0, 5, 10, 15, 20`) are **always safe**: for `v ≡ 5,10,15,
  20`, `v·c mod 25 ∈ {5,10,15,20}` for every unit `c` — never in `{0,1,24}`. Only
  `v ≡ 0` (a multiple of 25) can hit residue 0.
- The units split into **10 ± pairs**: `{1,24},{2,23},{3,22},{4,21},{6,19},{7,18},
  {8,17},{9,16},{11,14},{12,13}`. `c` clears the family `⟺` the pair `{±c⁻¹}` is
  disjoint from `V`'s residues.

So, for a family with **no multiple of 25**:

> **mod-25 covering FAILS (no clearing `c` exists) ⟺ `V` blocks all 10 unit ±-pairs**
> (has a residue in every pair).

The AP `{1,…,12}` blocks all 10 pairs (residues `1..12` cover each pair with its
`≤12` representative; `5,10` are the safe non-units) — correctly *not*
mod-25-clearable, consistent with `M = 1/13 < 2/25`.

## The decisive computation: the AP is the unique tight blocker

Over 150,881 families (height ≤ 50), **554 block all 10 pairs** (the entire
residual after kps's floor). Of those 554 pair-blockers, **exactly ONE has
`M < 2/25`: the AP `{1,…,12}` itself** (`M = 1/13`). Every non-AP pair-blocker
has `M ≥ 1/12 > 2/25` (min-M by defect count: d=0→1/13, d=1→1/12, d≥2→≥1/11).

So the pair-blocker M-spectrum is `{1/13} ∪ [1/12, ∞)` — the AP isolated at the
bottom, a clean gap `(1/13, 1/12) ⊃ (1/13, 2/25)` above it, empty.

## The three-case closure of the lower bound (M < 2/25 ⟹ AP)

Every 12-speed family falls into exactly one case:

1. **No mult of 25, some free unit-pair** → a clearing `c` exists → `M ≥ 2/25`.
   **[PROVED — kps `LRCMod25Floor.lean`, GREEN, kernel-pure]**
2. **No mult of 25, blocks all 10 pairs** → the AP (`M=1/13`) or `M ≥ 1/12 > 2/25`.
   **[the sharp residual — a pair-blocking rigidity; only the AP is tight]**
3. **Has a mult of 25** → the multiple sits at residue 0 for every `c`, but the
   family is loose at a *small* denominator (`M ≥ 2/25` easily; e.g. `2/11, 2/17,
   3/19`). **[kps's easy case]**

Together: **the AP is the unique 12-speed family with `M < 2/25`** — which is
exactly (G). The whole crux now sits in case 2's rigidity.

## Why this is progress

kps's open core ("near-tight ⟹ mod-25-clearable") was a covering-system fact
about which `c` works. The pair-blocking reformulation makes the **failure set
explicit and finite**: the non-clearable families are exactly the *pair-blockers*,
a highly structured 0.4% of families, and among them **only the AP is tight**. So
the residual is no longer "prove a covering fact" but the sharper, more rigid:

> **Pair-blocking rigidity: if `V` (12 speeds) hits every unit ±-pair mod 25 and
> `M(V) < 2/25`, then `V` is the arithmetic progression `{1,…,12}`.**

This is the mod-25 face of the AP-rigidity (S12: `M=1/13 ⟹ AP` because 13 is
prime; here the extra pair-blocking structure pins it further). It is a candidate
for a genuine proof — the pair-blocking condition is a concrete linear-algebra /
covering constraint mod 25, and tightness is the Freiman lever — rather than an
open analytic estimate.

## Honest scope

- Empirical over ~150k families (height ≤ 50) + targeted AP-dilate/near-AP
  blockers. The reformulation itself (non-unit safety, pair structure) is exact
  and Lean-ready; the AP-uniqueness among tight pair-blockers is the residual to
  prove.
- Multiple-of-25 families are classified separately (case 3); my `blocks_all_pairs`
  returns False for them (residue 0 present), routing them to the easy case, not
  the floor.

→ HYP-4622; sharpens HYP-4567/kps-S41 (mod-25 covering, `LRCMod25Floor.lean`);
composes with HYP-4612/S31 (defect-stratified min-M) + THM-632; rigidity ←
HYP-4382/S12 (AP-uniqueness at the tight locus).
