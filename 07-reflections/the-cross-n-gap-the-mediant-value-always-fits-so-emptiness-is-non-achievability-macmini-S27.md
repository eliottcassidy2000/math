# The cross-N gap: the mediant value always fits the window, so emptiness is non-achievability — and q=38 = 2·19 descends by parity

*mac-mini-2026-07-06-S27 (HYP-4572). Owner: understand the LRC at numbers other than
14 runners and leverage it for proof progress. Mapping the second gap across the
number of speeds N, integrated with opus-S117 (complexity is one parameter; crux =
the mediant 3/38 at N=12) and my own parity seed. Verified:
`lrc_mediant_q38_parity_macmini_S27.out`.*

## The cross-N picture

The second gap at `N` speeds is `(1/(N+1), 2/(2N+1))`; a below-`1/(N+1)` value is a
rung `s/(Ns+k)` of **order** `k` (opus-S116/S117). Known gap members (corrected
order `k = q − N·s`):

| speeds `N` | member | `M` | `s` | order `k` |
|---|---|---|---|---|
| 6 | `{1,5,6,11,16,17}` | `5/33` | 5 | **3** |
| 7 | `{1,3,4,5,7,13,18}` | `3/23` | 3 | **2** |
| 12 | — (the crux) | `3/38` mediant | 3 | 2 (unachievable) |

## The key observation: the value always fits — emptiness is non-achievability

The **mediant** `3/(3N+2)` (order `k=2`, `s=3`, the shallowest interior value) sits
at a **scale-invariant** position in the window — always `≈ 0.65` of the way up:

| `N` | mediant | rise above tight edge | fraction of window |
|---|---|---|---|
| 6 | `3/20` | 0.00714 | 0.65 |
| 7 | `3/23` | 0.00543 | 0.65 |
| 10 | `3/32` | 0.00284 | 0.66 |
| 12 | `3/38` | 0.00202 | 0.66 |

So the **value** `3/(3N+2)` is interior for *every* `N` — the window never squeezes
it out. Therefore **(G)'s emptiness at `N=12` is not that the mediant value fails to
fit; it is that no covering family can *achieve* it.** This is the sharp form of
opus-S117's `(O-depth-monotone)`: not "the window outruns the value" but "the window
outruns the *achievability*" — the achievable order drops (`N=6:k=3 → N=7:k=2 →
N=12:` none) even though the target values persist. The whole crux is a pure
**achievability** question, which is why it is n-specific and beyond any value/width
inequality alone.

## Leveraging the cross-N structure: q = 38 = 2·19 descends by parity

opus-S117 localized the crux to the **mediant `3/38` at `N=12`** — a finite
residue-hole-covering system at `q = 38`. And `38 = 2·19`, so the "14 = 2·7"-style
**parity descent** (my E_p/O_p seed) applies directly. At the witness `a/38` a member
clears the hole `{0,±1,±2} = {0,1,2,36,37} (mod 38)`. By CRT (`mod 2 × mod 19`):

- **even speeds** `2w`: residue `2(wa mod 19)`, must avoid the even part of the hole
  `{0,2,36}` = `2·{0,1,18} (mod 38)`, i.e. `wa` avoids `{0,1,18} = {0,±1} (mod 19)` —
  the halved even config is a **clearance-2 system mod 19** (`M' ≥ 2/19`);
- **odd speeds**: avoid the odd part `{1,37}` = `{1,18} (mod 19)`.

So the `N=12` mediant reduces to: the even coverers (multiples of `2,4,6,8,10,12` —
at least the multiple of `8`) **halve to a clearance-2 configuration mod 19**, while
the odd coverers (which may carry `3,5,7,9,11`) avoid `{±1} mod 19`. This is a
concrete finite feasibility problem mod 19 — opus's `(O-mediant)`, now with the prime
`19` and the parity split exposed. The obstruction, if any, lives in whether a
*covering* speed set can meet the halved clearance-2 constraint mod 19 without a
better witness appearing (raising `M` above `3/38`).

## Why this is the right leverage

Understanding other `N` shows the emptiness is n-specific *achievability*, not a
value inequality — so the proof must be a finite feasibility check at each cell, and
the shallowest cell (mediant `3/(3N+2)`) is the last obstruction. At `N=12` that cell
is `q = 38 = 2·19`, and the `2·19` factorization hands us the parity descent for
free: the problem is a **mod-19 clearance-2 covering feasibility**. The cross-N view
also pins *which* cell to attack (the mediant) and *why the prime 19 appears* (it is
the odd part of `2N+2 = 3N+2` at `N=12`, i.e. `38/2`), which no single-`N` analysis
reveals.

## Net

- **The mediant `3/(3N+2)` always fits the window (≈0.65, scale-invariant)** — so
  (G)'s emptiness is **non-achievability**, the achievable order dropping with `N`
  (`k=3→2→` none), not a value/width inequality. This sharpens opus's
  `(O-depth-monotone)`.
- **The `N=12` crux `3/38` = `q=38=2·19` descends by parity** to a mod-19 clearance-2
  covering feasibility (my E_p/O_p seed) — a concrete finite obstruction to attack
  for `(O-mediant)`.
- Honest: the higher-`N` gap members are *hard to construct* (my bordered-AP search,
  like opus's, misses the complex interior-defect ones); the map is anchored on the
  verified `N=6,7` members and the value/order structure, not a full enumeration.

## Pointers

- `lrc_cross_n_gap_depth_macmini_S27.py` (bordered-AP cross-N search — construction
  caveat), `lrc_mediant_q38_parity_macmini_S27.py/.out` (orders + q=38 parity descent).
- opus HYP-4496 (S117 complexity unification, mediant 3/38, O-depth-monotone/O-mediant),
  HYP-4486 (Kravitz spectrum); mac-mini HYP-4562 (kissing=order), HYP-4552 (kissing),
  the E_p/O_p parity seed (owner, S~6), HYP-4432 (q≤2max).
