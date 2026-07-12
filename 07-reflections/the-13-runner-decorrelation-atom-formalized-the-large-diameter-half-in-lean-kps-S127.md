# The 13-runner decorrelation atom, formalized — the large-diameter half in Lean

*kind-pasteur-2026-07-11-S127 cont.48. Owner: "work the 13-runner decorrelation atom bound." The fleet had
the atom's structure (mac-mini cont.49, THM-636) and a numerical verification (opus-S243), both computational.
This session machine-checks it: the 13-runner decorrelation atom at the direct `1/14` threshold, kernel-pure.*

---

## What was formalized

THM-636 (mac-mini-S38) formalized the decorrelation atom for **12 runners** at the old `2/25` gap. The current
endgame is the **13-runner** direct-LRC(14) form (`1/14` threshold). `LRCDecorrelation13.lean` (kernel-pure
`[propext, Classical.choice, Quot.sound]`, builds green):

- **`reach_decorr13`** — the atom itself: for `Vᵢ = bᵢ + L·Kᵢ` (`|bᵢ| ≤ B`, `0 < L`), a 13-speed family,
  `reach(K) − B/L ≤ reach(V)`. The proof is the reverse triangle inequality at the witness `t = t_K/L`
  (`distZ` is 1-Lipschitz), evaluated over `Fin 13` — the infrastructure (`margin`, `exists_max_margin`,
  `le_margin_iff`) is generic over the index type, so the 12-runner proof transferred verbatim.
- **`escape_loose13_le12`** — `≤ 12` distinct lifts ⟹ LRC(≤13) gives `reach(K) ≥ 1/13`; then any `L > 2366`
  forces `reach(V) > 1/14` (since `1/13 − 13/2366 = 1/14` exactly). The family is LOOSE.
- **`escape_loose13_le6`** — the sharp divisor-complete threshold (mac-mini cont.49): DC even-heaviness
  collapses the lifts to `≤ 6` distinct speeds, so LRC(7) gives `reach(K) ≥ 1/7`, and only `L > 182` is needed
  (`1/7 − 13/182 = 1/14`). This is the version that actually fires on large-diameter DC families.

The lift floor enters as a cited `LRC(≤13)` hypothesis — exactly as `LRC(≤13)` is used everywhere in this
project, and exactly as THM-636 cited `LRC(≤12)`.

## Where it sits in the endgame

mac-mini cont.49 showed the large-diameter half of the DC endgame is `[reach(V) ≥ reach(K) − B/L]` +
`[≤ 6 distinct lifts ⟹ reach(K) ≥ 1/7]` + `[bounded-diameter finite check as the descent base]` — "the two
halves are one descent," where the descent base is my cont.47 bounded-diameter finite check. This session
banks the **atom** of that descent as a machine-checked object: given the scale-separation decomposition and
the lift floor, `reach(V) > 1/14` is now proved in the kernel, not just verified numerically. It is the direct
counterpart, at `1/14`, of the `2/25` atom that has been in Lean since S38.

## Honest scope — what is proved and what is cited

What `escape_loose13_le6` proves is: **if** a 13-family decomposes as `bᵢ + L·kᵢ` with `|bᵢ| ≤ 13`, `L > 182`,
and the lift family has reach `≥ 1/7`, **then** it is loose. Two things remain outside the kernel and are the
genuine open pieces of the large-diameter half:

1. **The decomposition exists with a small base and few lifts.** That a large-diameter DC family admits
   `vᵢ = bᵢ + L·kᵢ` with `|bᵢ| ≤ 13` and `≤ 6` distinct lifts is the "lift-collapse" claim — provable, per
   mac-mini/klein-S263, from DC even-heaviness (multiples of `8, 14, …` force even runners, so the odd/lift
   structure collapses to ~6), but not yet a theorem. This is the real remaining mathematics.
2. **`reach(K) ≥ 1/7`** is the LRC(7) citation on the collapsed lift family — settled, but a citation.

So this is the atom made rigorous, not the whole half. The value is that the *analytic* core — the descent
inequality and its `1/14` arithmetic — is now kernel-checked, so the remaining work is purely the
combinatorial decomposition (piece 1) plus the bounded-diameter base, both finite/elementary.

## The convergence

Four independent objects are the same large-diameter looseness (mac-mini cont.49): THM-720 (pair-sum growth),
THM-636 (this atom), LEM-013 (dissociated margin), klein-S263 (~6 odd runners). opus-S243 verified the
13-runner atom numerically; mac-mini has its structure; this adds the Lean. When the same statement keeps
arriving from four directions and now sits kernel-checked, the object is real — the endgame's large-diameter
half is a decorrelation descent, and its atom no longer depends on anyone's arithmetic being right by hand.

*Files: LRCDecorrelation13.lean (kernel-pure; `reach_decorr13`, `escape_loose13_le12`, `escape_loose13_le6`),
root-wired. HYP-6135. Formalizes the 13-runner form of THM-636 (mac-mini-S38); machine-checks mac-mini cont.49
(structure) + opus-S243 (numerical); the descent base is [[the-clearing-window-is-unbounded-the-finite-check-is-only-bounded-diameter-kps-S127]].*
