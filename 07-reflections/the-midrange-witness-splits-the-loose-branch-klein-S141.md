# The midrange witness splits the loose branch: the tight ratio is closed, the gap lives in the spread region

*klein-2026-07-05-S141 (HYP-4161). Owner: push the open math AND the formalization even further.
Building on opus-S85's midrange bound and my S140 gap, this session (a) FORMALIZES the general
midrange witness kernel-pure, and (b) uses it to PROVE a chunk of the loose branch and localize the
residual precisely. A genuine advance on both fronts.*

## The formalized brick: the midrange witness (general, kernel-pure)

`LRCMidrangeWitness.lean` (Mathlib-only, kernel-pure, corpus-green):
> **`midrange_margin`**: for any speeds `v : ι → ℤ` all in `[m, Mx]` with `0 < m ≤ Mx`, at
> `t = 1/(m+Mx)` every runner clears every integer by `m/(m+Mx)`:
> `∀ i k, m/(m+Mx) ≤ |v_i · t − k|`.
> **`midrange_margin_compressed`**: if additionally `Mx ≤ 13 m`, the margin is `≥ 1/14`.

The proof reduces the real inequality `m/(m+Mx) ≤ |v_i/(m+Mx) − k|` to the pure integer fact
`m ≤ |v_i − k(m+Mx)|` (`midrange_int_core`, by cases `k ≤ 0` / `k ≥ 1` + `omega`). This is the
general `≥`-half of the covering-min — `M(S) ≥ v_min/(v_min+v_max)` — that opus-S85 (THM-526) uses
and that every gap/rigidity argument needs, now a reusable Lean lemma (it was only a canon-paper
statement + a skeleton mention before).

## The math: the midrange witness SPLITS the loose branch

The loose branch needs: a non-tight 12-runner base `B` has `M(B) ≥ 2/25`. The midrange bound
`M(B) ≥ v_min/(v_min+v_max)` clears `2/25` exactly when
```
    v_min/(v_min+v_max) ≥ 2/25  ⟺  v_max ≤ 11.5·v_min.
```
So the loose branch splits at the **ratio 11.5**:

- **Midrange-tight region (`v_max ≤ 11.5 v_min`): PROVEN loose.** Every such base has `M ≥ 2/25`
  directly from `midrange_margin` — no rigidity, no unbounded witness. (Verified: 0/400 sampled
  bases violate.) This whole region is closed by the formalized brick.
- **Spread region (`v_max > 11.5 v_min`): the residual, and where the dichotomy actually lives.**
  The AP `{1,…,12}` has ratio `12 > 11.5`, so **the tight extremizer sits in the spread region** —
  the midrange-tight region contains *no* tight family. The gap `(1/13, 2/25)` is entirely inside the
  spread region, and its threshold `2/25` is exactly the midrange value at the ratio-`11.5` boundary.

This **localizes the hard part precisely**: the loose branch is done except for spread bases
(`v_max > 11.5 v_min`), where the dichotomy is AP (`1/13`, tight) vs everything-else (`≥ 2/25`,
loose). The spread families are exactly the "one (or few) far runner(s)" shapes — the ladder /
killer-offset / far-peel territory (mac-mini THM-618, my S130 two-killer, S126 11-runner peel): a
spread base has a runner `> 11.5×` the others, which peels toward an 11-runner problem.

## Where this leaves the loose branch

Combining the fleet's pieces, the loose branch is now:
1. **midrange-tight (`v_max ≤ 11.5 v_min`)** — `M ≥ 2/25`, PROVEN + FORMALIZED (this session).
2. **spread (`v_max > 11.5 v_min`)** — the AP-rigidity residual (S140): `M = 1/13` iff AP, else
   `≥ 2/25`. Route: peel the far runner (ratio `> 11.5`) → 11-runner analogue (S126 gap `≥ 1/12`,
   gap to `2/23`) + the killer-offset value (THM-618); the spread runner is a "killer" on the
   11-runner base. opus-S85's exact comb formula `M({a,…,a+r−1}) = a/(2a+r−1)` handles the
   *full-comb* spread bases; the non-comb spread bases (e.g. `{1,…,11,24}` at `2/25`) are the
   ladder rungs.

So the open crux has shrunk from "the whole loose branch" to "spread bases only," with a concrete
peel route and the tight extremizer isolated inside it. The `≥ 2/25` for the midrange-tight half is
now a Lean lemma.

## Honest scope

- Formalized: the general midrange witness + compressed corollary (kernel-pure). Genuinely new to
  the corpus.
- Proven (math): the loose branch for `v_max ≤ 11.5 v_min` (via midrange).
- Open: the spread region (`v_max > 11.5 v_min`) — the AP-rigidity / gap, still needing the
  peel-to-11-runner + killer-offset synthesis (or the LRC(13) extremizer-uniqueness citation).
- Credit: opus-S85 (midrange bound THM-526, comb saturation, AP-is-min); my S140 (the gap +
  bounded violators); mac-mini THM-618 (killer-offset for the spread runner); S126 (11-runner gap).
- Convergence (same window): kps-S12 (HYP-4147) independently reframed the loose branch as the
  LRC(13) second-value gap (floor 1/13 cited, σ₂ = 2/25, gap 1/325 empty, "irreducibly analytic") —
  the same picture as S140 — and fixed the stale LRCTemplateSurface banner I flagged in S139.
  opus-S86 (HYP-4146) integrates the combs + S140 (cluster-desert + resonance-immunity q≤6). This
  session's distinct pieces: the **formalized** midrange witness (the Lean brick they use on paper)
  and the **11.5-ratio split** (midrange-tight region proven loose, localizing the residual to
  spread bases).

## Links

- Lean: `04-computation/lean/TournamentH7/TournamentH7/LRCMidrangeWitness.lean`
  (`midrange_margin`, `midrange_margin_compressed`, `midrange_int_core`; kernel-pure). HYP-4161.
- Builds on: opus-S85 HYP-4142 (midrange/comb/AP-min), klein-S140 HYP-4151 (12-runner gap),
  mac-mini THM-618 (killer-offset), klein-S126 (11-runner even-part gap), THM-526 (arc-width).
  Open: the spread region of the loose branch.
