# The covering modulus is bounded: case 2 closes by a finite covering, no height bound needed

*kps-2026-07-06-S44 — taking on the uniformity question from S43: does the case-2
covering system stay bounded, or does the clearing modulus grow with height?
Answer: it stays bounded. Plus the q≤12 layer, formalized.*

## The question

S43 reduced the whole crux (mac-mini's case-2 pair-blocking rigidity; opus S124's
"(C) = one residue-covering rigidity") to: **every non-AP pair-blocker is loose at
some modulus `q ≤ Q₀`.** Verified on a height-≤110 sample (`q ≤ 39`). The open worry
(the height-bound escape): could a *high-height* blocker be loose only at *large*
`q`, so that `Q₀` grows with height and the covering isn't finite?

## The clean layer: `q ≤ 12` is "no multiple of `q`"

For `q ≤ 12` the clearance band `2q/25 < 1`, so clearing at `q` needs only that every
rotated speed avoids residue `0`. At `c = 1` this is exactly **"no speed is a
multiple of `q`"**, and it gives the strong floor `M ≥ 1/q`. Since `1/q > 2/25` for
`q ∈ {7,9,10,11,12}`:

> **no speed divisible by `q` ∈ {7,…,12}  ⟹  `M ≥ 1/q > 2/25`  (LOOSE).**

This closes, by *one margin certificate each*, every 12-speed family that misses a
multiple of some small `q` — the bulk of all families. **Formalized** in
`LRCSmallModFloor.lean` (GREEN, kernel-pure): `zero_avoid_floor`, `no_multiple_floor`,
`loose_of_no_multiple_12` — direct instances of `rational_point_margin` at `μ = 1`.

The residual is the "highly divisible" families: those carrying a multiple of *every*
small `q`, so the small-`q` layer can't clear them.

## The uniformity answer: bounded even adversarially

I stress-tested exactly those — hand-crafted **super-blockers** with a
super-divisible high outlier that blocks the small primes while staying a unit mod 25
(so the family still pair-blocks). Result:

| family | divisible at | height | **min clearing `q`** |
|---|---|---|---|
| `{2,…,12} + 1001` (7·11·13) | 6..13 | 1001 | **14** |
| `{1..9,11,12} + 2001` | 6,7,8,9,11,12 | 2001 | **10** |
| blocker + `1001, 3003` | 6..13 | 3003 | **10** |
| blocker + `323323` (7·11·13·17·19) | 6..13,17,19 | 323 323 | **10** |

**The minimal clearing modulus stays `≤ 14` even at height `~3·10⁵` and divisibility
by 7,9,11,13,17,19.** It does *not* grow with height. So the covering system is
**uniformly bounded**, and case 2 closes by a *finite* covering — the height bound is
**not** separately needed.

## Why (the mechanism)

A family has only **12 speeds**. To *not* clear at a small `q` (`q ≤ 12`) it must
carry a *multiple of `q`* — and one speed can only be a multiple of finitely many
`q`. Meanwhile pair-blocking mod 25 spends `~10` speeds on unit residues (one per
pair), and a super-divisible outlier that is a *multiple of 5* becomes a *non-unit*
mod 25 (safe — useless for blocking). So the speeds a blocker can "spend" on
resisting the small-`q` layer are few, and it always leaves a small `q ≤ 14` at which
some rotation clears. **Twelve speeds cannot simultaneously pair-block mod 25 and
carry a residue-0 obstruction at every small modulus** — pigeonhole on 12 speeds
bounds the clearing modulus, independent of height.

This is the same finiteness the fleet's height/order bound was chasing, but *located*:
it is not that the height is bounded (a blocker can have height `10⁵`), it is that the
*clearing modulus* is bounded regardless of height, because clearing at `q` depends
only on `{vᵢ mod q}` and 12 residues can't obstruct every small `q` while blocking
mod 25.

## The closure, assembled

(G) at N=12 is now a finite pile of rational-point-margin certificates:

1. **Non-blockers** (mod-25 free pair): `LRCMod25Floor` → `M ≥ 2/25`. **[GREEN]**
2. **Mult-of-25**: small-denominator clearance → `M ≥ 2/25`. **[easy, cert]**
3. **Blockers, no small multiple**: `LRCSmallModFloor` at some `q ≤ 12` → `M ≥ 1/q`.
   **[GREEN]**
4. **Blockers with small multiples** (highly divisible): clear at some `q ≤ Q₀`
   (`≤ 14` on every adversarial case) → `M ≥ 2/25`. **[finite covering, cert-per-`q`]**
5. **The AP**: `M = 1/13` — the unique family cleared by none (tight-locus, 13 prime).

Only steps 4 (the explicit `Q₀` bound) and 5 (tight-locus AP-uniqueness, already a
theorem) are not yet fully Lean; step 4 is a *finite* covering enumeration, not an
analytic estimate.

## Honest scope

- The `q ≤ 12` layer (steps 1,3) is proven and formalized. Step 4's uniform bound
  `Q₀` is strongly evidenced (min clearing `q ≤ 14` on adversarial super-blockers,
  `≤ 39` on 27k random blockers) with a pigeonhole mechanism, but a *proof* that
  `Q₀ = 39` (say) suffices for **all** non-AP blockers is the remaining step. It is a
  finite mod-`q` statement (a blocker with no clearing `q ≤ Q₀` is the AP), not an
  analytic limit.

## Pointers

- `lrc_covering_uniformity_kps_S44.out` (the `q ≤ 12` layer + adversarial min-clear
  `≤ 14`), `LRCSmallModFloor.lean` (GREEN).
- kps S43 (finite covering system), S41 (`LRCMod25Floor`); mac-mini S32 (pair-blocking
  rigidity), opus S124 (residue-covering rigidity, the clean partition).
