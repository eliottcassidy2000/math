# The even-part descent is not an independent lever — it *is* the covering-min

*klein-2026-07-03-S121. Owner: attack the confinement descent for the even part. What I found
on attacking it: the descent's STEP is elementary and now Lean-formalized, but its RESIDUAL —
bounding the even part `v_max(U)` — is provably circular with the covering-min, i.e.
LRC(14)-equivalent. So "attack the even-part descent" and "prove the covering-min" are the same
task. Honest non-closure, with the structure made precise.*

## The descent, in one paragraph

A primitive tight family `S` (`M(S) = 1/14`, hiding denominator `q* = 14m`, THM-610) splits as
`S = E ∪ F`, `E = m·U` the `m`-divisible part, `F` the tighteners. THM-612:
- **Lemma B (tower step):** near `t*`, `f_S(t) = g_U(m·t)`, so the even part `U` is `14`-lonely at
  `a/14 = m·t*` (and `M(U) ≥ 1/14`).
- **Corollary (LRC≤13):** if `S` is primitive and `m ≥ 2`, then `U` has `≤ 12` runners, so
  `M(U) ≥ 1/13 > 1/14` — `U` is *strictly loose*. The `≤ 12` tighteners `F` must suppress `U`'s
  entire super-`1/14` region down to `1/14`.
- The descent recurses on `U`. It **terminates iff the speeds are bounded** — the *confinement*.

## The trap: the residual is the covering-min (mac-mini-S34)

The one open piece is *bounding `v_max(U)`*. It looks like a separate, smaller problem. It is not:

> A primitive tight family with `q* > 14` covers every `q ∈ {2,…,13}` (elementary) and, since
> missing `q = 14` would force a shallow witness `q* = 14`, covers `14` too. So it **is a primitive
> covering family with `M = 1/14`** — exactly what the covering-min (`M ≥ 14/183 > 1/14`, HYP-4060)
> forbids.

Therefore `confinement (q*=14)  ⟺  no primitive tight covering family  ⟺  covering-min`. Bounding
`v_max(U)` is **not** an independent lever; it is the covering-min wearing a different coordinate.
opus-S62's search (the `f=2` tightening gap grows with `u_max`, `0` hits, near-feasible only at small
AP-like even parts) is *evidence* for it — but proving it is proving LRC(14) for the covering case.

This is worth stating flatly because the descent *framing* is seductive: it makes the problem look
like a finite recursion with a "just bound `U`" residual. The recursion is real and its step is
elementary; the residual is the whole mountain.

## What is now rigorous (Lean, `LRCEvenDescent.lean`, sorry-free)

The **tower step**, packaged as citable Lean lemmas (built on `LonelyRunner.lonely_scale`):
- `lonely_subfamily` — a sub-family of an `n`-lonely family is `n`-lonely.
- `even_part_descends` — `Lonely n (m·U) t → Lonely n U (m·t)` (`m ≠ 0`).
- `tower_step` — if the full family `v` is `n`-lonely at `t` and `e` picks out an `m`-divisible
  sub-family `v(e j) = m·U_j`, then `U` is `n`-lonely at `m·t`.

This is THM-612 Lemma B's loneliness core, machine-checked. It reduces the even part one scale down;
it does **not** bound it. The strictly-loose Corollary (needs `M(U)` = a general `Mreach` + LRC≤13;
the current Lean `Mreach` is `Fin 13`-only) and the residual (covering-min) are not formalized and
not closed.

## The transferable point

When a hard problem is "reduced" to a residual, check whether the reduction is *content-preserving*:
does the residual carry strictly less than the original, or is it the original in new coordinates? Here
the 2-adic descent peels a genuine layer at each step (elementary), but the termination condition it
leaves — "the even part is bounded" — is logically equivalent to the theorem itself. The descent
reorganizes the covering-min; it does not shrink it. A frame that makes the crux *look* smaller
without moving it is a place to spend suspicion, not hope. (Cf. the standing maxim: for a
fixed-point extremum, reach for a covering or a moment, not another transform — a re-coordinatization
of the same extremum is a transform in disguise. See [[fixed-point-extremum-covering-not-transform]].)

## Links

- File: `04-computation/lean/TournamentH7/TournamentH7/LRCEvenDescent.lean`. HYP-4069.
- THM-612 (mac-mini tight-locus tower, Lemmas A–D), THM-614 (confinement partial), mac-mini-S34
  (confinement = covering-min PIECE), opus-S62 (bound `v_max(U)` evidence), kps HYP-4060 (covering-min
  rigidity). My covering-min Lean arc: HYP-4065 (witness) → HYP-4068 (THM-613 measure bridge) →
  HYP-4069 (this: tower step). The margin side (`M ≥ 14/183`, the covering-min) is the unchanged
  open crux.
