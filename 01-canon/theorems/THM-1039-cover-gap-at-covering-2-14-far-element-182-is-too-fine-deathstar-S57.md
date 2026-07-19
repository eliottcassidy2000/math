# THM-1039 — Conditional cover-gap lenses at far element 182 (death-star-S57; scope corrected codex-S74)

**Status:** the exact cover-gap criterion and two conditional lenses are proved: stability when
`delta>max(W)/(2366k)`, and soft Weyl when `C<=464k*mu`. The scripts explore bounded compact banks,
but no completed exhaustive `max(W)<=34` certificate is stored. The former
claim that every non-AP core has `coverGap=1/2`, and the associated equality/uniqueness
classification, are withdrawn: the scripts provide examples and finite telemetry, not that uniform
theorem. This post-extraction chart does not extract the core, prove AP rigidity, or close LRC(14).
Scripts: `lrc_covergap_214`, `lrc_covergap_uniform`, `lrc_covering214_test`,
`lrc_covering_witness_V` (S57).

Setting: `V = W ∪ {v_max}`, `W` covers 2..12 and misses 13,14 (AP-core-adjacent), `v_max = 182k` (so `V`
covers 2..14). `G_W = {t : min_{w∈W}‖wt‖ ≥ 1/13}` (measure `μ`, `C` components).
`coverGap(W,v_max) := max_{t∈G_W}‖v_max·t‖`. `δ = M(W) − 1/13`.

---

## 0. The correction (MISTAKE-161)

cont22 claimed the binding far element is the smallest multiple of **13** (26, 39, …) "because covering 14
is not required." **WRONG.** For `M<1/14` (the LRC(14) threshold) the sieve forces `V` to cover **2..14**;
the core misses 13,14 ⟹ `13∣v_max` AND `14∣v_max` ⟹ **`182∣v_max`** (boxeph THM-1017, line 27). The
witness that the covering-14 requirement is essential: `V={1,2,3,5,7,8,9,10,11,12,17,19,104}`,
`M=8/105∈(1/14,1/13)`, primitive, covers 2..13 but **misses 14**, non-AP core — a false alarm at level 1/13.
So the analysis runs at `v_max=182k`, not `13k`. THM-1038's original `182` was right.

## 1. Exact criterion

`M(V) < 1/13 ⟺ coverGap(W,v_max) < 1/13` (the far element covers `G_W`). Proof: if some `t∈G_W` has
`‖v_max t‖ ≥ 1/13` then `min_V(t) ≥ 1/13`, so `M(V) ≥ 1/13`; conversely if `coverGap<1/13` then on `G_W`
the far element fails and off `G_W` some `w` fails, so `M(V)<1/13`. ∎

## 2. Finite telemetry: the far element is often very far from binding

The 13-lattice points `a/13 = 14a/182` are **centers** of `v_max`'s danger arcs, each of half-width
`1/(13·182·k) = 1/(2366k)`. A good component around the maximizer `t_0 = p/D_W` has half-width
`>= delta/max(W)` (min_W decreases at rate at most `max(W)`). In the sampled non-AP cores,
`coverGap(182)` lies between `0.44` and `0.50`; the examples `{1..11,24}`,
`{1..10,22,24}`, and `{2..12,15}` give `0.5`. The AP deep-well chart gives `coverGap=0`.
These computations motivate an alignment-rigidity statement, but they do not prove either
`non-AP => coverGap=1/2` or `coverGap=0 => AP` uniformly. The latter kind of classification
contains the open rigidity content.

## 3. Two rigorous lenses (both at `v_max = 182k`)

- **Stability (geometric):** if `δ > max(W)/(2366k)` then `coverGap ≥ 1/13`. *Proof.* The good component
  at `t_0` contains `[t_0 − δ/max, t_0 + δ/max]`. If `t_0` sits in a `182k`-danger arc (`‖v_max t_0‖<1/13`),
  it is within `1/(2366k)` of the arc center; since the half-width `δ/max > 1/(2366k)`, the component
  reaches the arc edge, where `‖v_max t‖ = 1/13`. If `t_0` is outside all arcs, `‖v_max t_0‖ ≥ 1/13`
  already. Either way `coverGap ≥ 1/13`. ∎ (= THM-1028's "empty window" `v_max > max/(13δ)`, geometrically.)
- **Soft Weyl (uniform in max):** if `C ≤ 464k·μ` then `coverGap ≥ avg_{G_W}‖v_max t‖ ≥ 1/13`. *Proof.*
  `avg = 1/4 − (2·7ζ(3)/8)/π³ · Σ_{m odd}|ĉ_{m·182k}|/m² / μ ≥ 1/4 − 2.104C/(π³·182k·μ)`, using
  `|ĉ_N| ≤ C/(πN)`; `≥ 1/13 ⟺ C ≤ (9/52)π³·182k·μ/2.104 = 464.3k·μ`. ∎ The constant is **independent of
  `max(W)`** — this closes every non-fragmented core at every scale (THM-1037, corrected to `v_max=182k`).

The lenses define complementary *conditional regions*: a row not handled by
either one must satisfy

```text
C>464k*mu  and  delta<=max(W)/(2366k).
```

The earlier `99.84%` figure is bounded-search telemetry from THM-1038, not a
uniform assertion that fragmentation forces a quantitative lower bound on
`delta`.

## 4. What is closed, what remains

- **Bounded compact searches:** `lrc_covering214_test` scans AP-adjacent
  substitutions with entries below `35` and a bounded list of `182k` far
  elements. `lrc_covergap_214` contains a parameterized exact scan intended
  to cover all eligible cores up to `MAXW`, but this repository currently
  stores no nonempty completed output for its default `MAXW=34` run.
  Consequently `max(W)<=34` is not claimed closed here.
  The stored `lrc_covergap_uniform` run used far elements `26,39,...`
  determined from coverage only through modulus `13`; by MISTAKE-161 that
  output is historical wrong-class telemetry and is not evidence for the
  Cover14/`182k` statement.
- **Uniform conditional regions:** soft Weyl `C<=464k*mu` and stability
  `delta>max(W)/(2366k)` are proved. Their intersection complement is the
  very-near-tight fragmented boundary; bounded telemetry does not close it.
  HYP-3901 proposes a renormalization route at larger scales.
- **A necessary loose-core inequality, not an AP reduction.** The stability
  lens (§3) runs backwards: at `v_max=182k`,

  ```text
  M(V)<1/13 => coverGap<1/13 => delta<=max(W)/(2366k).
  ```

  For `k=1` and `max(W)<=34`, this gives
  `M(W)<=1/13+34/2366<0.0913`.

  **REFINED (S57, see reflection `the-stability-step-is-not-near-tight-implies-AP...`):** this is NOT quite
  right — "near-tight ⟹ AP" is FALSE as an exact implication (`{1..11,24}`, M=2/25, is near-tight, non-AP).
  The desired stronger sufficient statement is **cover-gap uniqueness**
  `coverGap(W,182)<1/13 => W=dilated AP` (near-tight is
  necessary, not sufficient). THM-1039's two lenses prove it except a **very-near-tight fragmented** residual
  (`C>464μ` AND `δ≤max/2366`). A sequence in this region with `δ→0` enters the
  HYP-7310 near-tight limit, but the two inequalities do not themselves force
  that convergence uniformly. The bounded banks find only AP rows at exact
  equality, but the uniform twelve-speed equality classification is open (HYP-4382); a finite
  verification is not that theorem. Empirically the only covering-2..14 `M<1/13` families in these
  banks are the deep wells `{1..12,182m}` (core = tight AP).
- **The genuine rigidity suppliers remain:** (i) the n=12 equality/stability
  probes; (ii) nested scales (HYP-3901). Cover-gap uniqueness is a stronger
  sufficient post-extraction target, not a theorem established by these two
  lenses.

**Net:** with the correct far element `182`, the proved lenses rule out broad stable and
non-fragmented regions, and the finite bank shows large empirical slack. The very-near-tight
fragmented residual and the AP/cover-gap uniqueness statement remain open. THM-1149 further shows
that this is a post-extraction view: crown collapse and the n=12 equality classification are separate
suppliers, so the position analysis must not be advertised as a global inverse theorem.

→ THM-1017, THM-1028, THM-1033, THM-1037, THM-1038 [restored to `v_max=182`], MISTAKE-161, HYP-7310, HYP-3901.
