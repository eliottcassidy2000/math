---
source: opus-2026-07-11-S251
status: BRIDGE (my S248/S249/S250 arc <-> the fleet's Ostrowski ladder). The low M-spectrum is the Ostrowski/
  Farey tree rooted at [0;13,·]; S248's empty gap (1/14, 3/41) is exactly the Farey gap [0;13,1] -> child
  [0;13,1,2]. The tight locus {AP, V*} (= THM-612's {AP, GW}) are the COMPLETE and ONCE-PUNCTURED {k/14}-
  progressions (AP full, g=1; V* punctured at 12/14, g=2) -- so the open core "tight => {k*alpha}-structure"
  HOLDS on the classified locus, with three-gap free (THM-527). The composite 14 = 2·7 is exactly why Ostrowski
  rung 1 has TWO occupants (the punctured progression is integer-realizable-and-tight only via the doubling
  12->24) where the proved prime k=12 has one. Credits prior fleet Ostrowski work; new = the explicit
  {k*alpha}-support placement + the composite two-occupant mechanism.
tags:
  - lrc14
  - ostrowski-ladder
  - k-alpha-structure
  - three-gap
  - tight-locus
  - farey-tree
  - composite-14
  - bridge
---

# The tight locus is the complete and once-punctured {k/14}-progressions (Ostrowski rung 1)

**opus-2026-07-11-S251.** Owner: work the remaining Ostrowski LRC mathematics. My S248–S250 arc turns out to be
the **rung-1 (AP) end** of the fleet's covering-min Ostrowski ladder — and reading the two together pushes the
open core ("tight ⟹ {kα}") onto the *full* tight locus. This is a bridge; I credit the prior work and state
what is genuinely added.

## Prior work I am building on

- **Ostrowski ladder** (mac-mini S38): the covering-min rungs are `M_k = [0;13,k] = k/(13k+1)`; ends are rung
  `k=1 = 1/14` (AP, tight) and rung `k=14 = 14/183` (deep well, DC covering-min).
- **Farey-tree M-spectrum** (mac-mini S65cont54, klein S266/S267): the low spectrum `1/14 < 3/41 < 2/27 <
  14/183 < 1/13` is a Stern-Brocot mediant tree; `14/183 = n/Φ₆(n)` is the crux value.
- **Three-gap rigidity** (THM-527): once the extremal config is known to be a `{kα}`-progression, Steinhaus
  gives `g ≤ 3` for free. The open core is proving the **`{kα}` structure**, not the gap count.
- **Tight locus `{AP, GW}`** (THM-612 / HYP-2621): the `M = 1/14` families are two. *(My S249's "two mod-14
  patterns" re-derives this; `V* = {1..11,13,24}` is `GW`. Credit to THM-612 for the two-element fact; S249's
  addition was the mod-14 residue characterization and the doubling mechanism, not the count.)*

## New: the two tight families are the complete and punctured {k/14}-progressions

At the tight optimum `t = 1/14`, phases are residues mod 14, and the tight locus is exactly the two `{kα}`-
configs (`α = 1/14`):

- **AP `{1..13}`** = the **complete** progression `{k/14 : k=1..13}` — gaps `{1/14}`, `g=1` (uniform).
- **V\* `{1..11,13,24}`** = the progression **punctured at `12/14`** (phase-set `{1..11,13}/14`, since
  `24 ≡ 10` collides on the existing `10/14` and vacates `12/14`) — gaps `{1/14, 2/14}`, `g=2`.

Both are `{kα}`-supported, both three-gap, both closest-approach `= 1/14` (tight). **So the open core "tight ⟹
{kα}-structure" holds on the entire classified tight locus** — S249's two patterns *are* the two `{kα}`-configs
at `n=14` (complete and once-punctured), and the three-gap rigidity is then automatic. This is the explicit
form the ladder frame wanted: the tight configs are `{kα}`-progressions, verified for *both* occupants, not
just the AP.

## New: the composite 14 = 2·7 is why rung 1 has two occupants

Realizing the **punctured** progression by an **integer** family that also **stays at rung 1** (no better
witness at another `t`) is the delicate part. Computing the actual `M` for single moves `12 → m` (all of which
land back on the progression at `t=1/14`, so all are `{kα}`-supported there):

| move | `m mod 14` | `M` | rung |
|---|---|---|---|
| `12→24` | 10 | **`1/14`** | `[0;13,1]` — **stays, rung 1** |
| `12→36` | 8 | `3/41` | `[0;13,1,2]` — lifts to the Farey child |
| `12→26,38,50,64` | 12/10/8/8 | `1/12` | lifts |

So `{kα}`-support at `t=1/14` is **necessary but not sufficient**; only the **doubling `12→24`** keeps `M=1/14`
(the others lift to higher Ostrowski rungs). And the doubling produces a *collision* (`24 ≡ 10`) only because
`14 = 2·7` is **composite** — the map `k ↦ 2k` has a nontrivial kernel mod 14. At a **prime** `n`, `k ↦ 2k` is
a bijection, no collision, so the punctured progression has **no** integer realization and **rung 1 has one
occupant (the AP alone)**. This is the `{kα}`/Ostrowski explanation of the two-element tight locus, and the
precise sense in which the composite `14` separates `k=13` from the proved prime `k=12` case: *rung 1 gains a
second occupant.*

## Net — the remaining Ostrowski math, one rung pushed

- **S248's empty window** `(1/14, 3/41)` = the Farey/Ostrowski gap between `[0;13,1]` and its child
  `[0;13,1,2]` — empty because no simpler continued fraction lies between a value and its own child.
- **S249's tight locus** `{AP, V*}` = the **complete and once-punctured `{k/14}`-progressions** — so "tight ⟹
  `{kα}`" is confirmed *on the classified locus* (both occupants), three-gap free.
- **The composite `14 = 2·7`** is exactly why Ostrowski rung 1 has two occupants (the doubling `12→24`), vs one
  at the proved prime `k=12`.

The still-open piece is unchanged in kind — prove `tight ⟹ {kα}`-support for **arbitrary** families, not just
verify it on the locus — but it is now anchored to the explicit two-config `n=14` picture, and the "why two"
is pinned to the composite doubling. The natural next step in this frame is the **equioscillation / Chebyshev**
route (the repo's `covering-min-is-a-chebyshev-equioscillation` reflection): the tightest configs push all
phases to the band edge `±1/14`, which *forces* near-uniform `{kα}` spacing — the general "tight ⟹ `{kα}`"
implication that would close the classification half for all families.

→ mac-mini S38 (Ostrowski ladder), mac-mini S65cont54 / klein S266-S267 (Farey tree, `14/183`), THM-527
(three-gap rigidity), THM-612 / HYP-2621 (tight locus `{AP, GW}`), opus-S248 (empty window), opus-S249 (mod-14
classification), opus-S250 (base-rigidity), the chebyshev-equioscillation reflection (the general route).
Files: `lrc14_ostrowski_tight_locus_ka_support_opus_S251.py` (+`.out`).
