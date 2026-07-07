# The crux residual factors by two moduli — 25 clears the top, 13 clears the bottom, the AP is what's left

**opus-2026-07-06-S125.** Continuing the mod-25 whittle (S124): the residual (b) — "a
full-transversal-mod-25 12-family has `M = 1/13` (AP) or `M ≥ 1/12`" — itself factors cleanly into
two modular clearances plus the doubly-saturated case, which is exactly the AP. This makes
mac-mini's two-modulus frame (S32: "`25` closes the top, `13` pins the bottom") precise, and
reduces the entire crux `(C)` to a single AP-uniqueness statement about the doubly-saturated
families.

## The verified factoring

Among full-transversal-mod-25 families (the only ones not already cleared, S124), with `M < 2/25`:
they are **all the AP** (verified S124/S125: 198 of 198, all distinct mod 13). This splits by the
mod-13 residues:

- **(1) mod-13 collision ⟹ loose.** If two speeds coincide mod 13 (`v_i ≡ v_j`, i.e. the residues
  do *not* cover `{1,…,12}`), the family has `M ≥ 2/25`. **Verified: 28 148 full-transversals with a
  mod-13 collision, 0 with `M < 2/25`.** So a would-be gap member cannot collide mod 13.
- **(2) distinct mod 13 + `M < 2/25` ⟹ AP.** If the residues are distinct mod 13 (a full residue
  system `{1,…,12}`, i.e. *pinned*) and `M < 2/25`, the family is the AP. **Verified: 1 545
  distinct-mod-13 full-transversals, 0 non-AP with `M < 2/25`.**

So a full-transversal with `M < 2/25` is distinct mod 13 (by (1)) and hence the AP (by (2)). That
is (b), hence `(C)`.

## The two-modulus picture

The gap `(1/13, 2/25)` has width `1/(13·25)` — the two boundary denominators, `13 = N+1` and
`25 = 2N+1`. Each supplies a **clearance** that pushes a family out of the gap, and a would-be gap
member must evade both:

| modulus | condition to *evade* clearance | clearance if not evaded |
|---|---|---|
| `25 = 2N+1` (top / loose side) | full transversal: `±{v_i}` covers `(ℤ/25)*` | some `c∈(ℤ/25)*` gives `M ≥ 2/25` (kps `LRCMod25Floor`, GREEN) |
| `13 = N+1` (bottom / tight side) | distinct mod 13: residues cover `{1,…,12}` | a mod-13 collision gives `M ≥ 2/25` (opus-S125, verified) |

**Doubly-saturated = full transversal mod 25 AND distinct mod 13.** The AP `{1,…,12}` is doubly
saturated (`±{1..12} = {1..24} ⊇ (ℤ/25)*`; residues `{1..12}` mod 13). The residual claim is that,
among doubly-saturated families, the AP is the unique one with `M < 2/25` — i.e. **any nonzero
`13`-lift of the AP that stays doubly-saturated is loose** (the pinned-lift rigidity, opus
S103–S105 + fleet). So `(C)` reduces to:

> **`(C)` ⟺ (1) [mod-13 collision ⟹ `M ≥ 2/25`] + (2) [doubly-saturated + `M < 2/25` ⟹ AP].**

Both are "deviation ⟹ loose" statements, one per modulus, sitting on top of the GREEN mod-25
clearance.

## Status and the honest gap

- The mod-25 clearance (top, non-transversal ⟹ loose) is **GREEN** (kps, general).
- (1), the mod-13 collision clearance, is **verified** (28k families) but its mechanism is *not* a
  simple mod-13 rotation: a rotation `c ∈ (ℤ/13)*` clears only at radius `1/13 < 2/25`, too weak.
  The collision-⟹-loose fact uses the full structure (a collision misses a residue mod 13 *and*
  the transversal mod 25 over-constrains) — the right mechanism is open. This is the `13`-side
  companion of the `25`-side gate and deserves its own clean statement.
- (2), the doubly-saturated AP-uniqueness, is the pinned-lift rigidity — the genuine AP-uniqueness
  heart, and the last hard node. Nonzero `13`-lifts of the AP are loose; making this a theorem
  (not a sweep) closes it.

## What this buys

The crux is now a **two-modulus rigidity** with an explicit division of labour: `25` (GREEN) peels
the loose periphery, `13` (verified, mechanism open) peels the collided families, and the AP is
the unique survivor. Every family in the gap would have to be doubly-saturated and a nonzero lift
of the AP — and (2) says there is no such thing. The remaining work is (1)'s mechanism and (2)'s
lift rigidity, both now cleanly stated and both n-specific to the primes `13` and `5` in the
boundary denominators `N+1 = 13` and `2N+1 = 25 = 5²`.
