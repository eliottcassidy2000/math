---
source: opus-2026-07-09-S189
status: FORMALIZED the Freiman-equality CORE of `burden = 11 ⟹ AP` (LRCFreimanBurden.lean, kernel-pure):
  restrictedSum_eq_freimanChain -- at the MINIMAL descent burden |A+^A| = 2|A|-3, the restricted sumset is
  EXACTLY the min/max chain (the low chain {m+y} u the high chain {x+M}). Plus the chain infrastructure
  (freimanChain def, _subset, _card) and the lower bound restrictedSum_card_ge (S188, refactored). The
  full AP conclusion (chain => constant differences) is the row-bijection + reflection argument -- complete
  PAPER PROOF below, remaining as the next Lean step (needs order-indexed elements via orderEmbOfFin).
  Also IDENTIFIED the finish-line obligation map: (i) AP step; (ii) burden=12 => 4 shapes; (iii) the
  moment-floor legs hbonf/hB/hsmall/hsize NOW UNBLOCKED by death-star-S4's de-opaquing of witnessG2/shapeOf
  (were opacity-blocked in opus-S186); (iv) burden=13 2-D-GAP routing to the density floor (opus-S187).
tags:
  - lrc14
  - freiman
  - equality
  - formalization
  - obligation-map
  - finish-line
---

# The Freiman-equality core, formalized; and the finish-line obligation map

**opus-2026-07-09-S189.** Owner: formalize `burden = 11 ⟹ AP`, and identify remaining obligations. I
formalized the CORE (kernel-pure) and give the complete proof of the AP step as the next Lean target, then
map every remaining finish-line obligation.

## Formalized (LRCFreimanBurden.lean, kernel-pure)

- `freimanChain s hne` — the min/max chain `{m+y : y∈s, y≠m} ∪ {x+M : x∈s, x∉{m,M}}`.
- `freimanChain_subset` / `freimanChain_card` — the chain is a subset of `restrictedSum s` of card
  exactly `2|s|−3`.
- `restrictedSum_card_ge` — the Freiman lower bound `|A+̂A| ≥ 2|A|−3` (refactored from S188 as a corollary).
- `burden_ge_eleven` — `|s|=7 ⟹ burden ≥ 11` (THM-675 (ii) floor).
- **`restrictedSum_eq_freimanChain`** — the equality CORE: `|A+̂A| = 2|s|−3 ⟹ restrictedSum s =
  freimanChain s`. (Proof: chain ⊆ sumset, both have card `2|s|−3`, so `eq_of_subset_of_card_le`.)

This pins the sumset to the chain at minimal burden — the structural entry point of the AP characterization.

## The AP step (complete paper proof; the remaining Lean target)

**Claim.** For `a₀ < a₁ < ⋯ < a_{k-1}` with `|A+̂A| = 2k−3`, `A` is an AP.

Write `dᵢ = a_{i+1} − aᵢ > 0`. By `restrictedSum_eq_freimanChain`, `A+̂A` is exactly the chain
`a₀+a₁ < a₀+a₂ < ⋯ < a₀+a_{k-1} < a₁+a_{k-1} < ⋯ < a_{k-2}+a_{k-1}` (positions `0 … 2k−4`).

*Row-1 bijection.* Consider the `k−2` sums `a₁+a₂ < a₁+a₃ < ⋯ < a₁+a_{k-1}`, all in `A+̂A` (valid pairs
`1<j`). Each is `> a₀+a₂` (position 1) and `≤ a₁+a_{k-1}` (position `k−1`). For `j ≤ k−2`, `a₁+aⱼ < a₁+a_{k-1}`
equals a chain element that is NOT a high-chain element (those are all `≥ a₁+a_{k-1}`), so `a₁+aⱼ = a₀+a_{l}`
for some `l`, with `l ≥ j+1` (since `a₀+aⱼ < a₁+aⱼ`). The `k−3` sums `{a₁+aⱼ : 2 ≤ j ≤ k−2}` thus inject
(strictly monotonically) into the `k−3` low-chain elements `{a₀+a_l : 3 ≤ l ≤ k−1}` — an order-preserving
injection between two `(k−3)`-sets, hence the identity shift `l = j+1`. So `a₁+aⱼ = a₀+a_{j+1}`, i.e.
`d₀ = a_{j+1}−aⱼ = dⱼ` for every `2 ≤ j ≤ k−2`. Hence `d₀ = d₂ = d₃ = ⋯ = d_{k-2}`.

*Reflection.* Apply the same to `B = {−a_{k-1-i}}` (an AP iff `A` is; `B+̂B` has the same card; `B`'s
differences are `A`'s reversed). Row-1 for `B` gives `d_{k-2} = d_{k-4} = ⋯` down through `d₁`. Concretely
for `k=7`: row-1 gives `d₀=d₂=d₃=d₄=d₅`; reflected row-1 gives `d₅=d₃=d₂=d₁=d₀`; union ⟹ all six equal ⟹
AP. ∎

**Lean status of the AP step.** The proof needs the elements order-indexed (`s.orderEmbOfFin`, `a : Fin k
→ ℤ` monotone), the "`a₁+aⱼ` is a low-chain element" case split, and the order-preserving-injection ⟹
identity-shift lemma. That is a substantial (but routine) formalization — the natural next file, built
directly on `restrictedSum_eq_freimanChain`. It is NOT yet in Lean; flagged as the remaining piece of the
user's request (I formalized the core, not the full AP).

## The finish-line obligation map (identified)

For the LRC(14) grand-assembly residual (THM-671) via THM-675, the remaining obligations are:

1. **AP step** (above) — `burden = 11 ⟹ AP` in Lean, from the formalized equality core. Routine, sizeable.
2. **`burden = 12 ⟹` the 4 shapes** (opus-S187: `{0..7}` minus one interior point). Same row machinery,
   one excess; finite once the AP step's infrastructure exists.
3. **The moment-floor legs — NOW UNBLOCKED.** death-star-S4 de-opaqued `witnessG2`/`shapeOf` in the
   skeleton (opus-S186's flagged blocker). So `lrc14_from_momentfloor_nodes`'s legs `hbonf`
   (`toReal_bonferroni`, kps-S30), `hB` (Lemma B on the concrete `G_P`), `hsmall` (`k≤7` pigeonhole,
   `GOOD = univ`), `hsize` (concrete cluster length) can now be discharged against the concrete
   `witnessG2 = (slowμ (GOOD ∩ G_P)).toReal`, shrinking that route to `{hMoment, hpartA}`. This is the
   most immediately dischargeable obligation and directly closes the density-floor leg.
4. **`burden = 13` 2-D GAPs** (opus-S187): route through the density floor / looseness (THM-661 /
   opus-S186 moment floor), NOT the near-AP check — the S181 dimension boundary.

The two live routes (grand-assembly residual → THM-675 Freiman; moment-floor node → density floor) now BOTH
have their remaining pieces identified and partially formalized; item 3 is the shortest path.

## Ledger

- FORMALIZED (kernel-pure) the Freiman-equality CORE `restrictedSum_eq_freimanChain` (minimal burden ⟹
  sumset = chain) + chain infra + lower bound (LRCFreimanBurden.lean). The AP step (chain ⟹ constant
  differences) has a complete paper proof (row-bijection + reflection) and is the next Lean target.
- Obligation map: AP step; burden=12 shapes; the moment-floor legs (UNBLOCKED by death-star-S4, shortest
  path to the density-floor leg); burden=13 2-D-GAP routing. -> THM-675, THM-671, opus-S187/S186/S181,
  HYP-5682, death-star-S4 (de-opaquing), kps-S30.
