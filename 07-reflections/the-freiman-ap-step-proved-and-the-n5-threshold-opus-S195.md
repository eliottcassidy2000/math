---
source: opus-2026-07-09-S195
status: PROVED the Freiman AP step in Lean (LRCFreimanAP.ap_of_min_burden, kernel-pure): a StrictMono
  sequence whose first n terms have MINIMAL restricted sumset 2n-3 is an AP -- FOR n >= 5. Found + recorded
  that the S189 blueprint OVERSTATED generality: the step is FALSE for n <= 4 (bi-arithmetic non-AP minima
  like {0,1,3,4}); TRUE iff n >= 5 (census: non-AP minimal sets = 44,38,0,0,0,0 for n=3..8). MISTAKE-133.
  The proof uses a CLEANER route than the blueprinted row-bijection+reflection: the interleaved
  consecutive/skip chain. Feeds kps's Freiman ladder (THM-681 rung) and THM-675 (7-classes, n=7>=5).
tags:
  - lrc14
  - freiman
  - ap-step
  - formalization
  - n5-threshold
  - interleaved-chain
  - mistake-133
---

# The Freiman AP step, proved in Lean; and the n ≥ 5 threshold that the blueprint missed

**opus-2026-07-09-S195.** Owner: commit a full session to the AP step (`burden = 2n−3 ⟹ AP`). I proved
it kernel-pure — but the session's real lesson is that **verifying before formalizing caught a false
generalization**, and the verification then handed me a proof *simpler* than the one I had blueprinted.

## The correction: the step is n-dependent (MISTAKE-133)

The opus-S189 blueprint stated the AP step for all `k ≥ 2`. Before writing a Lean proof I ran the census
(`lrc14_ap_step_verify_opus_S195.out`): primitive `n`-sets with `|A +̂ A| = 2n − 3`, non-AP count =

| n | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|
| non-AP minimal sets | 44 | 38 | **0** | **0** | **0** | **0** |

So the equality forces an AP **iff n ≥ 5**. For `n ≤ 4` the minimum is achieved by NON-AP "bi-arithmetic"
sets — `{0,1,3,4}` (differences `1,2,1`) has `A +̂ A = {1,3,4,5,7}`, card `5 = 2·4−3`, yet is not an AP.

**Why the threshold is exactly 5.** The row/interleaved argument pins `d_{i+2} = d_i` (even-index gaps
equal each other, odd-index gaps equal each other) — a *bi-arithmetic* sequence, not yet an AP. Breaking
the even/odd tie needs one more forced collision, and the shortest one available is `a₀ + a₄` — which
requires the 5th element `a₄`. At `n ≤ 4` there is no `a₄`, the tie is never broken, and bi-arithmetic
sets slip through. This is the same "small sumset ≠ near-AP" seam as the S181 2-D-GAP obstruction, in its
smallest incarnation. Had I formalized the blueprint verbatim I would have proved a false statement (or
fought an unclosable goal). The CLAUDE.md rule "check small cases exhaustively before generalizing" is
what caught it — the census was 10 lines and saved the session.

For the LRC application this costs nothing: THM-675 needs majority-parity **7-classes**, `n = 7 ≥ 5`.

## The cleaner route the verification revealed

The S189 blueprint proved the AP step by a **row-bijection + reflection** (order-preserving injection of
`{a₁ + aⱼ}` onto chain positions, then reindex `aᵢ ↦ −a_{k-1-i}` for the missing difference). While
verifying I noticed the `2n − 3` sums split not just as (min-chain ∪ max-chain) but as **consecutive**
`sᵢ = aᵢ + aᵢ₊₁` and **skip** `tᵢ = aᵢ + aᵢ₊₂` sums, which strictly interleave `s₀ < t₀ < s₁ < t₁ < ⋯`.
This gives a proof with no cardinality bijection and no reflection:

- **`Rset = IC`** (`{sᵢ} ∪ {tᵢ}`): the interleaved sums are `2n − 3` distinct elements of the sumset, so
  at minimal burden they *are* it.
- **Step 1** (`diff_two`, `d_{i+2} = d_i`): `aᵢ + aᵢ₊₃` is a restricted sum lying strictly in
  `(tᵢ, tᵢ₊₁)`, whose only chain element is `sᵢ₊₁`, so `aᵢ + aᵢ₊₃ = aᵢ₊₁ + aᵢ₊₂`.
- **Step 2** (`sum04`, `d₀ = d₁`, needs `n ≥ 5`): `a₀ + a₄` is a restricted sum; the strict order rules
  out every `sₖ` (via Step 1 at `i = 0, 1`) and every `tₖ` except `t₁ = a₁ + a₃`, forcing `a₀ + a₄ =
  a₁ + a₃`, hence (with `d₃ = d₁`) `d₀ = d₁`.
- **Assembly** (`ap_of_min_burden`): `d₀ = d₁` and `d_{i+2} = d_i` ⟹ all differences equal, by strong
  induction. Every consecutive difference is `a₁ − a₀`.

Each step is membership + strict monotonicity + `omega` — the kind of reasoning Lean does cleanly. Steps
1, 2, and the assembly each compiled on the first or second try; the whole file is kernel-pure
(`[propext, Classical.choice, Quot.sound]`, no `sorry`, no `native_decide`). The blueprint's reflection
would have cost a painful `ℕ`-reindexing; the interleaving avoided it entirely. **The obstruction to a
clean formalization was a signpost to a cleaner proof** — the same pattern this project keeps rewarding.

## Where it plugs in

kps-S124 names the remaining THM-681 ingredient as "the Freiman ladder [E3 rigidity the top, mac-mini/opus
rungs]." The AP step is a rung: it is the `n ≥ 5` near-AP characterization THM-675 applies to
majority-parity 7-classes (burden 11 ⟹ the class is a dilated interval ⟹ non-covering, klein-S211). The
one packaging task left is the **Finset bridge**: `ap_of_min_burden` is stated for an indexed
`a : ℕ → ℤ`; wrapping it for a `Finset ℤ` (via the sorted enumeration `orderEmbOfFin` extended
`StrictMono`ly to `ℕ`, with `Rset a n = restrictedSum s`) makes it directly citable. That is routine but
fiddly dependent-type work — flagged as the clean next step, not attempted here to avoid risking the
completed core.

## Ledger

- PROVED (kernel-pure) `LRCFreimanAP.ap_of_min_burden`: `StrictMono a` + `|Rset a n| = 2n − 3` +
  `n ≥ 5` ⟹ `a` is an AP. Via the interleaved-chain route (`Rset_eq_IC`, `diff_two`, `sum04`).
- FOUND + recorded MISTAKE-133: the step is FALSE for `n ≤ 4` (bi-arithmetic minima); the S189 blueprint
  overstated it as general. Census in `lrc14_ap_step_verify_opus_S195.out`.
- Updated `LRCFreimanBurden.lean` docstrings (the AP step is proved, `n ≥ 5`, cleaner route).
- Files: `LRCFreimanAP.lean` (+root wire), `lrc14_ap_step_verify_opus_S195.out`, MISTAKE-133.
  → THM-675 (descent burden), THM-681 (kps Freiman ladder), opus-S189 (blueprint, corrected)/S187/S181,
  HYP-5682, LEM-015, klein-S211.
