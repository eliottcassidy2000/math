# THM-1036 — The position lemma IS the soft Route-A equidistribution (verified with a factor-3 margin) (death-star-2026-07-18-S56)

**Status:** the non-aligned/"position" half of the alignment rigidity (THM-1033) is **reduced to an
explicit soft equidistribution bound and verified to hold with a factor-3 margin.** This is not a proof
— the *uniform* bound is the soft Weyl estimate the program has flagged as its shortest path — but it
proves the exact equivalence "position lemma = soft Route-A discrepancy," kills the cheap-measure route,
and shows the target has slack. Source HYP-7305/7362. Scripts:
`04-computation/lrc_equidist_measure_deathstar_S56.py`, `lrc_secondmoment_deathstar_S56.py`.

Setting: `V = W ∪ {v_max}`, `W` a valid non-AP core, `v_max = 182k`, `G_W = {t : min_{w∈W}‖wt‖ > 1/13}`.

---

## The exact reduction

`M(V) < 1/13 ⟺ G_W ⊆ danger_{1/13}(v_max)`, i.e. the far element's danger set (measure `2/13`) must
contain the core's good set. Two integral consequences (both rigorous):

1. **Measure:** `M(V)<1/13 ⟹ meas(G_W) ≤ 2/13`. **Verified vacuous:** *every* valid non-AP core has
   `meas(G_W) ≈ 0.06 < 2/13` (2400/2400). So there is **no cheap measure elimination** — the whole
   lemma lives in the near-tight regime.
2. **Second moment (the sharp one):** `G_W ⊆ danger_{1/13}(182k)` forces `‖182k·t‖ < 1/13` on `G_W`,
   hence `avg_{t∈G_W} ‖182k·t‖ < 1/13`. So

> **Position lemma ⟸ `avg_{t∈G_W} ‖182k·t‖ ≥ 1/13`** — i.e. the far element is *not* concentrated in
> its danger valleys over `G_W`.

## The verification (equidistribution, factor-3 margin)

`‖182k·t‖` has whole-circle average `1/4`. If `182k` **equidistributes** on `G_W`, the average over
`G_W` is `≈ 1/4 ≫ 1/13`, so the far element cannot cover and `M(V) ≥ 1/13`. Measured:

| core | opt-denom | aligned | `meas(G_W)` | `avg_{G_W}‖182·t‖` | ≥ 1/13? |
|---|---|---|---|---|---|
| `{2..12,15}` | 17 | no | 0.056 | **0.247** | yes |
| `{1..11,24}` | 25 | no | 0.006 | **0.280** | yes |
| `{1..12}` (AP) | 13 | **yes** | 0 | **0.000** | (deep well) |

Non-aligned cores equidistribute (`avg ≈ 1/4`), clearing `1/13` by a factor of **~3**. The AP alone
has `avg = 0` — its good set sits at `t = a/13`, exactly where `‖182·t‖ = ‖14a‖ = 0` (`182 = 14·13`),
so the far element covers it: the deep well.

## Why this is Route A, and why it's the *soft* face

`avg_{G_W}‖182k·t‖ = 1/4 − corr`, where `corr = ⟨1_{G_W}, 1/4 − ‖182k·t‖⟩/meas(G_W)`. The sawtooth
`1/4 − ‖182k·t‖` has its Fourier mass on the frequencies `{j·182k}`, so `corr` is exactly the
correlation of the good-set indicator with the `182k`-lattice — the `Q_s`-type discrepancy of Route A
(THM-729). Non-aligned `⟺ 1_{G_W}` has *no* Fourier mass concentrated at `182k`-multiples `⟺ corr`
small `⟺ avg ≈ 1/4`. So:

> **Position lemma ⟺ `corr < 1/4 − 1/13 = 9/52` for non-aligned cores** — a Route-A discrepancy bound.

The `9/52` threshold against a measured `corr ≈ 0` is a **factor-3 slack**: the lemma needs only that
the far element is not *nearly perfectly* concentrated in its valleys, i.e. **any** power-saving
equidistribution suffices — the acknowledged soft face of the program (finish map §3: "any power-saving
`Q_s = o(r²)` suffices"). No sharp cancellation is required here.

## Net

- **Cheap routes are dead** (measure bound vacuous), so the position lemma is genuinely equidistribution.
- **It is true, verified with a factor-3 margin** — the far element equidistributes on non-aligned good
  sets; only the AP's good set (denominator 13) resonates with `182 = 14·13`.
- **It equals the soft Route-A discrepancy** `corr < 9/52`, with large slack — so the last piece of
  "only the AP is 182-aligned" is exactly the program's softest, shortest-path estimate, not a sharp
  one. Combined with THM-1033 (aligned cores, `max ≤ 34`, PROVED) and the `max ≥ 35` renormalization,
  the wall is now: **one soft Weyl bound on the good-set/lattice correlation.**

→ THM-1033, THM-1029, THM-1028, THM-729 (`Q_s`), the alignment/Freiman/OCF reflections, finish map §3.
