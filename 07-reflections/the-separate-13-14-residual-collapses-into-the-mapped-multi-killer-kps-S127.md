# The separate-13/14 residual collapses into the already-mapped multi-killer case

*kind-pasteur-2026-07-11-S127 cont.68. Owner: "work the separate-13/14 multi-killer residual" — the open case
cont.67 identified (13 and 14 covered by different runners rather than the single lcm-outlier 182). Working it
shows it is not a new open piece: it is structurally forced into the multi-killer case I already reduced, and
a naive refutation it seemed to threaten is dispelled. The frames close.*

---

## The residual is forced to core ≤ 11 — it *is* multi-killer

A separate-13/14 covering family has, by definition, **two distinct killer slots**: a multiple of 13 (`13u`,
not a multiple of 14) and a multiple of 14 (`14w`, not a multiple of 13). With the core runner 1, that is
already `1 + 2 = 3` slots spoken for. A family has 13 runners, so at most `10` remain for the interval body
`{2,…,12}` — but the full run `{1,…,12}` needs all `12`. Therefore:

> **A separate-13/14 family can never have the single-killer core `{1,…,12}`; its interval-core length is `≤ 11`.**

Verified over the enumerated primitive separate-13/14 covering families: the maximum interval-core length is
`≤ 11` (e.g. `{1..11, 13, 84}`, where `84 = 6·14` supplies both `12` and `14`, and `13` supplies `13`). So the
separate-13/14 residual is **contained in the multi-killer case** (core `≤ 11`), which is exactly the case
already handled:

- **cont.58** — exact enumeration: every multi-killer minimizer has `M ≥ 7/89 ≈ 0.0787 > 14/183`.
- **cont.59** — the core-length monotonicity reduces multi-killer to **LRC(13)-escape + finite check**, anchored
  by `1/13 > 14/183`.

Spot-checks confirm looseness: `{1..11,13,84} → 7/89`, `{1..14}\{7} → 1/11`, `{1..10,13,14,84} → 1/11` — all
comfortably above `14/183`. cont.67's "open residual" is therefore **not a new object**; it is the mapped
multi-killer, seen through the 2-runner (co-binder-splitting) lens.

## The naive `1/15` refutation is dispelled

The `M = N/(a+b)` form (`N` = arc-gap numerator) raised a worry: a family like `{1,…,14}\{7}` contains small
`13` and `14`, and the pair `{1, 14}` sits at Farey-neighbour arcs (`N=1`) balancing at `M = 1/15 ≈ 0.0667 <
14/183`. If that were the family's `M`, it would refute the covering-min. It is not:

> `{1,…,14}\{7}` has `M = 1/11 ≈ 0.0909`, attained at `t* = 3/22` and bound by `{8, 14}` — **not** `{1,14}` at
> `1/15`. The family is *lonelier elsewhere*; `1/15` is a local, not global, optimum.

The lesson for the 2-runner extremal: `M = min N/(a+b)` is over the **global-optimum** binding pair, which the
equioscillation constraint (all other runners `≥ M`) selects. Naive small pairs (`{1,14}` at `1/15`) fail that
constraint — some other runner leaves a deeper gap at a different time — so they are never the minimizer. This
is exactly why the covering-min stays at `14/183`: the only pair that binds *at its own Farey-tight balance
without another runner cutting deeper* is `{1, 182}` (the single lcm-outlier), and splitting `13`/`14` forces a
shorter core whose runners cut a deeper gap elsewhere, raising `M`.

## Net — the covering-min lower bound, fully mapped

Putting cont.55–68 together, the covering-min lower bound `primitive DC ⟹ M ≥ 14/183` decomposes with no
remaining hidden pieces:

1. **Single-killer** (co-binder a multiple of `lcm(13,14)=182`): `M = 14m/(13(14m+1)) ≥ 14/183` uniformly —
   **proved and machine-checked kernel-pure** (cont.60/61/66).
2. **Multi-killer = separate-13/14** (core `≤ 11`, two killer slots): `M ≥ 7/89 > 14/183` — enumerated (cont.58),
   reduced to LRC(13)-escape + finite check (cont.59). **This session shows cont.67's "separate residual" is
   exactly this, not a new case.**
3. **`|core| ≥ 2` and large core runners**: equidistribution (opus-S259–263, mac-mini) — the analytic bulk.
4. **Bounded families**: klein's ILP (`≤ 182`).

The `2`-runner extremal (cont.67) unifies the *frames*; this session shows it introduces no new *cases* — the
co-binder either is the lcm-outlier (proved) or splits (multi-killer, mapped). The genuinely open analytic work
is the `|core|=1` smooth-body discrepancy (cont.62/63) that opus's Fourier route targets; the structural
skeleton is complete.

*Files: lrc14_separate1314_kps_S127.py-style checks in lrc14_separate1314_kps_S127.out. Consolidates cont.67
(2-runner unification), cont.58/59 (multi-killer reduction), cont.55 (lcm-outlier). HYP-6244.*
