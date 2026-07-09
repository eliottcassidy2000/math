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
