# Is LRC(14) proved? — honest end-to-end assessment (mac-mini-2026-07-08-S59)

Written after the density floor closed (S58) and the finite-`Vmax` glue was reduced (S58/S59) and
the shared `𝒲̂`-decay constant proved (LEM-011, S59). Verdict up front:

> **LRC(14) is NOT yet 100% proved, but it is reduced to a SINGLE, sharply-localized, strongly-
> verified analytic lemma — `j* = O(k)` — whose extremal case (arithmetic progressions) is PROVED
> and which has zero counterexamples over 90k+ adversarial clusters. Everything else in the proof is
> closed, a-priori, finite/decidable, or a cited settled result.** Plus Lean transcription.

## The end-to-end chain and each link's status

Let `S ⊂ ℤ`, `|S|=13`, be a putative counterexample (observer = 14th runner; dilation-normalized).

| # | Link | Status |
|---|------|--------|
| 1 | q-witness sieve: counterexample ⟹ covering-saturated `S = P ∪ {Vmax−e_i}` (THM-369) | **PROVED (Lean-checked)** |
| 2 | Non-covering branch ⟹ LRC(≤13) | **SETTLED** (owner directive; cited) |
| 3 | Reformulation: `∃ good period ⟹ M(S) ≥ 1/14` (THM-527/kps-S4) | **PROVED** modulo the finite-`Vmax` glue (link 6) |
| 4 | `ρ*_{1/7}(P,E) ≥ m_P`, `k≤7` branch (pigeonhole `μ=1` a.e.) | **PROVED** (exact, THM-530) |
| 5a | `k≥8` union bound `ρ* ≥ meas(G_P)+μ−1 ≥ m_P` | **PROVED** (elementary Bonferroni) |
| 5b | density floor `μ_{1/7}(E) ≥ bar_k`, `k=8..13`, uniform (the input to 5a) | **CLOSED** (THM-661 `B_d` + LEM-009 tail); tail constant now **a-priori** (LEM-011) |
| 6 | finite-`Vmax` glue THM-527-A: a good period exists on the `Vmax`-grid | **REDUCED** to `j*=O(k)` (LEM-010); AP case proved, general verified |
| 7 | Lean transcription of the finite/decidable pieces | **not done** (engineering) |

## What closed this week (S58–S59)

- **The density floor (link 5b)** — all six legs `k=8..13` at the uniform (min-over-families)
  level, diameter-free: `k=8,9,10` by the degree-4 moment `B_4` (block-minimizer), `k=11,12,13` by
  the degree-3 `D3` + the LEM-009 longest-AP tail. This removed the old "consecutive minimizes `μ`"
  crux entirely (the covering-moment bound needs no extremal lemma).
- **The shared `𝒲̂`-decay constant (LEM-011)** — the EXACT closed form `𝒲̂(n) = (6/7)^z ∏b(n_i)·Q(N)`
  for the uncovered-measure Fourier coefficients, verified. This makes opus-S157's density-floor
  **tail rate a-priori** (`V_j` had been numerically certified) — link 5b loses its last
  certification — AND gives THM-664's Weyl grid-resonance sum a-priori. One object, both bounds.
- **The finite-`Vmax` glue (link 6)** collapsed from "equidistribution at `Vmax ≤ 91^12`" to two
  elementary lemmas (`j=1` wraparound + Dirichlet, LEM-010) plus the single lemma `j*=O(k)`, whose
  worst case (APs, `j* ≤ ⌈7(k−1)/6⌉`) is proved and whose general form has `max j* = k` over 90k+
  adversarial clusters.

## The one remaining mathematical gap

**`j* = O(k)`: for every primitive `k`-set of co-offsets (`k ≤ 13`) and every `Vmax`, some period
`j ≤ C·k` has `maxgap{frac(e_i j/Vmax)} > 1/7`.** Proven for arithmetic progressions (the extremal
family, via Dirichlet on the single step). Open in general — the difficulty is the usual `k > 7`
one (union bound over arcs fails, `k/7 > 1`), i.e. it needs a genuine simultaneous-clustering
argument, not a per-arc count. But: (a) it is a statement about `≤ 13` integer points, purely
Diophantine; (b) its extremal case is proved; (c) `max j* = k` exactly across 90k+ adversarial
clusters, never absent; (d) if it holds, THM-527-A's finite check collapses to `Vmax ≤ O(k)`, inside
kps-S30's exact `M(S)` sweep — closing the covering case unconditionally.

## Formalization readiness

Because link 6 is not fully closed, a sorry-free Lean proof of the *whole* theorem is premature.
But the **ready pieces** are substantial and `native_decide`/elementary-shaped:
- link 1 (sieve) is already Lean-checked; link 2 is a cited hypothesis;
- links 4, 5a, 5b are finite exact-rational computations (the `ρ*` union bound; the `k≤7`
  pigeonhole; the compact `B_d` Farey checks) — `decide`/`native_decide` targets;
- LEM-010(i) (the `j=1` wraparound lemma) and the AP case of `j*=O(k)` are elementary and
  self-contained — good first Lean targets;
- LEM-011's closed form is an elementary Fourier computation.

**Recommended next step:** either (a) prove `j* = O(k)` in general (closes the mathematics), or
(b) begin the Lean assembly with the covering-case skeleton `[sieve] + [LRC≤13] + [ρ*≥m_P] +
[good-period]` sorry-parameterized on link 6, then discharge the finite nodes. This session
formalizes LEM-010(i) as the first concrete node.

## S60 update (mac-mini) — R2 done, R3 refines the glue

- **R2 CLOSED:** the k=12,13 density-floor tail rate is now a-priori. Via LEM-011, `|D3−D3_∞| ≤
  C/(pd)` with `V_j = 0.25/0.16/0.10` (k=12), `0.24/0.16/0.10` (k=13), `C = 19.1/18.7`, `pd₀ =
  101/65` — the k=11 computation reproduces opus-S157's `C=21.2, pd₀=160` exactly. Tail a-priori for
  all of k=11,12,13 (LEM-009 R2 section). Link 5b (density floor) is fully a-priori.
- **R3 REFINES the finite-`Vmax` glue (not just a writeup):** the good period must lie in `G_P ∩
  Good_E`, and LEM-010's `j=1` (which clears `Good_E`) typically FAILS `G_P` (`‖p/Vmax‖ < 1/14` for
  a small observer speed). The full good period is a **spread-out `j`** in the intersection — which
  DOES exist (0 failures / ~6000 admissible spread-dense `(P,E)`), but the elementary `j=1` proves
  only the cluster half. So the finite-`Vmax` glue is the *intersection* `G_P ∩ Good_E` grid-nonempty
  question (measure `ρ*≥m_P`), same structure, and the `j*=O(k)` lemma should be stated for it.
- **`j*=O(k)` still the one gap, now robustly O(k):** hard adversarial search (40k+ structured
  clusters) gives `max j* = 7/9` at k=11/13, all `≤ ⌈7(k−1)/6⌉`; AP case proved. The general
  simultaneous-clustering proof (and its `G_P`-intersection version) is the remaining mathematics.

Net after S60: every analytic constant is a-priori; the covering case rests on the single Diophantine
`j*=O(k)` lemma (for `G_P ∩ Good_E`), AP case proved, `0` counterexamples; plus Lean.

## S61 update (mac-mini) — the HARD branch of j*=O(k) is now ELEMENTARY (klein-S196)

Major advance. The `j*=O(k)` capstone splits on the longest AP `L`:
- **Near-AP (`L ≥ k−5`) — the HARD branch (j* largest, ≈k):** now **ELEMENTARY PROVED** (LEM-012,
  klein-S196): Dirichlet-cluster the `L`-AP into span `< (L−k+6)/7`, its complement is a gap
  `> (m+1)/7` (`m=k−L ≤ 5`), which the `m` stray points split into `≤ m+1` pieces, the largest
  `> 1/7`. `j* ≤ ⌈7(L−1)/(L−k+6)⌉`. No Weyl, no resonance sum. Recovers mac-mini's exact-AP case.
- **Dissociated (`L ≤ k−6`) — the EASY branch:** verified `j* ≤ 3/5` (k=11/13, mac-mini-S61,
  1400+ clusters). Route (a): dissociated ⟹ few small resonances `n·E≡0` ⟹ small partial-sum
  correction `Corr_N` ⟹ small `r_N` ⟹ small `j*` — supported by the mac-mini-S61 `r_N` decomposition
  (`Corr_N` is the near-resonance partial-sum, so few resonances ⟹ small `Corr_N`). Not yet proved.

The two ranges tile all `L`. So **`j*=O(k)` — hence THM-527-A, hence the covering case — is closed
except the dissociated branch (`L ≤ k−6`, `j* ≤ 5`, verified, "easy")**. This is a big improvement
over the S60 state ("one hard Diophantine lemma"): the hard part is done elementarily; only a
tiny-`j*` easy branch (with a clear 𝒲̂-smallness route) and Lean transcription remain.

**Lean:** three good-period nodes now formalized sorry-free (`LRCGoodPeriodJ1.lean`):
`good_period_j1_wraparound` (LEM-010(i)), `good_gap_of_phases_in_interval` (the arc core, shared by
j=1 / Dirichlet / AP / LEM-012), and `goodPeriod_iff_partialSum_pos` (opus-S165's `S_N>0` reduction).

## S62 CONSOLIDATION (mac-mini) — the complete state, math + formalization

**MATH — the covering-case chain, node by node:**
| Link | Status |
|---|---|
| q-witness sieve (THM-369) | PROVED (Lean-checked) |
| non-covering ⟹ LRC(≤13) | SETTLED (cited) |
| reformulation `good period ⟹ M≥1/14` (THM-527) | PROVED (mod the finite-`Vmax` bridge = good-period existence) |
| density floor `μ≥bar_k`, k=8..13 (THM-661) | CLOSED, uniform, a-priori (LEM-011/R2) |
| `ρ*≥m_P` (THM-530/663) | unconditional given the density floor |
| good-period existence, near-AP `L≥k−5` (LEM-012) | PROVED (elementary: Dirichlet + gap-split) |
| good-period existence, dissociated `L≤k−6` (LEM-013 + arc-count) | VERIFIED (exhaustive `s≤22`, adversarial band, large-spread `c≈0.22–0.37 ≪ ρ*≈0.96`) |

**The single math residual:** a closed-form a-priori `#arcs ≤ c(L)·spread` for dissociated
(`c(L) ≤ 0.37` for `L≤7`, verified with margin; the additive-combinatorics closed form is the
fleet's ongoing work) + the finite check on the small-spread band (`s≤22` exhaustive done). Every
OTHER link is proved / a-priori / cited. Both good-period branches close by ELEMENTARY routes
(Dirichlet gap-split + positive arc-count), bypassing the Mertens-hard partial-sum cancellation
(kps-S92/S93). So LRC(14)'s covering case is closed modulo this one verified (not-yet-closed-form)
additive bound + Lean.

**FORMALIZATION — Lean nodes (all sorry-free, axioms = `[propext, Classical.choice, Quot.sound]`):**
- `lonely_of_Mreach_ge` (skeleton): `Mreach ≥ 1/14 ⟹ lonely`.
- `LRCGoodPeriodJ1.lean` — 4 good-period cores:
  - `good_period_j1_wraparound` (LEM-010(i): `spread<6Vmax/7 ⟹ j=1` good);
  - `good_gap_of_phases_in_interval` (the arc core: phases in a `<6/7` interval ⟹ gap `>1/7`);
  - `goodPeriod_iff_partialSum_pos` (opus-S165's `S_N>0 ⟺` good period);
  - `gap_split_pigeonhole` (LEM-012's core: `m+1` gaps summing to `>(m+1)/7 ⟹` one `>1/7`).
- Skeleton `LRCFourteenSkeleton.lean` derives `LRC14Statement` from Prop obligations (`rhoStar>0`,
  `thm527_partA`), which the analytic content (density floor, reformulation) discharges.

**Remaining formalization:** discharge the analytic obligations (`native_decide` the finite density-
floor `B_d` checks and the `ρ*` union bound; transcribe LEM-012's Dirichlet step; the arc-count
existence) — an engineering task, no new mathematics. The elementary good-period cores are done.
