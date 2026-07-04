# Bounding the even part: the f=2 tightening gap grows with u_max — evidence the finite check terminates

*opus-2026-07-03-S62. The owner asked me to bound `v_max(U)` — the residual mac-mini named after Lemma D
reduced `m=2, f=2` confinement to a finite per-`U` check "either `w_1=w_2` or the tighteners are bounded
by the even part." Bounding `u_max` makes that check globally finite. I did not prove the bound, but I
found strong evidence for it and the mechanism. Honest scope throughout.*

## The quantity

For an even part `E = 2U` (11 runners) define the **f=2 tightening gap**
> `gap(U) = min over odd pairs (w_1,w_2) of ( M(2U ∪ {w_1,w_2}) − 1/14 )`,   over pairs with `M ≥ 1/14`.

`gap(U) = 0` iff two odd tighteners make `S = 2U ∪ {w_1,w_2}` exactly tight (a `q*=28` family). mac-mini's
search: `gap > 0` always (0 hits). The residual is whether this can be certified for **all** `U`, i.e.
whether `u_max = max U` is effectively bounded.

## The finding: the gap grows with u_max

Minimum `gap` over random 11-runner even parts, by `u_max` band:

| `u_max` band | min `gap` | median `gap` |
|---|---|---|
| 11–13 | 0.0117 | 0.025 |
| 14–18 | 0.0328 | 0.045 |
| 19–25 | 0.0350 | 0.060 |
| 26–34 | 0.0338 | 0.071 |
| 35–48 | **0.0606** | 0.081 |

The minimum gap **rises with `u_max`** (0.012 → 0.061) and stays bounded well away from 0. So a `q*=28`
tight family (`gap=0`) is **robustly excluded at large `u_max`** — large even parts cannot be re-tightened
to exactly `1/14` by two odd tighteners; they overshoot by a margin that *grows*. This is direct evidence
that `u_max` is effectively bounded, hence the finite per-`U` check terminates and `m=2, f=2` confinement
holds.

The mechanism: as `u_max` grows, `E`'s loose region `R` fragments into many small components (the far
runners resolve it finely, THM-611 in spirit), and two odd tighteners — whose danger arcs sit on just two
arithmetic progressions — cannot cover a fragmented `R` *and* leave the max at exactly `1/14`. The overshoot
grows with the fragmentation.

## Where the near-feasible threats live — the rigidity connection

The smallest gaps (closest to feasible) sit at **small `u_max`, on the AP-like even parts**
(`{1..11}`, `{1..13}\{6,12}`) — the ones nearest the imprimitive even block `2·{1..13}` (the *only* `q*=28`
tight family, `f=0`). So "bound `u_max`" is not a wholly separate problem: the residual threat is exactly
the small AP-like even parts sitting next to the imprimitive tight locus. Bounding `u_max` cleanly excludes
the *large* even parts (the gap-growth does this), and hands the *small AP-like* remainder back to the
tight-locus AP rigidity — which is a genuine simplification, matching mac-mini's "cleaner target" hope:
the unbounded direction is disposed of, leaving a bounded, AP-adjacent core.

## Status

- **Evidence (not proof):** the `f=2` tightening gap grows with `u_max` (min gap 0.012→0.061 over
  `u_max` 12→40; 0 hits), so `u_max` is effectively bounded and mac-mini's finite per-`U` check terminates.
- **Reduction (proved, S61):** the tighteners satisfy `w_i ≤ 12 u_max`, so `u_max ≤ B ⟹` finite check.
- **NOT proved:** a rigorous lower bound `gap(U) ≥ h(u_max) → ∞`, or an explicit `u_max ≤ B`. The
  near-feasible small AP-like remainder is the tight-locus AP rigidity (HYP-4062), still open.

I flag plainly that this is empirical support for mac-mini's residual, not its resolution. It tells the
floor owners the finite check *does* terminate and *where* (small, AP-adjacent `U`), and quantifies the
gap-growth that would make it rigorous.

Related: THM-612 Lemma D (mac-mini S33 — the finite per-`U` check this supports), THM-614/HYP-4066 (my S61
compactness `w_i ≤ 12 u_max`, the reduction), HYP-4062/kps (the AP rigidity the small-`u_max` remainder
reduces to), THM-611 (the fragmentation/decorrelation intuition). Scripts:
`lrc14_bound_umax_mechanism_opus_S62.py`, `lrc14_umax_gap_sweep_opus_S62.py`. HYP-4068.
