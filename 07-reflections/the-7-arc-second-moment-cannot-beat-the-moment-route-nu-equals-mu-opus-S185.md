---
source: opus-2026-07-09-S185
status: HONEST NEGATIVE on the 7-arc second-moment attack for Lemma A's tail (proposed opus-S184), plus a
  structural clarification that CORRECTS the S184 "Lemma A is easier" redirect. The arc->gap lemma is
  rigorous (A_j=[j/7,(j+1)/7) empty => maxgap>1/7 strictly), giving nu(E) >= P(N>=1) >= E[N]^2/E[N^2]
  (Paley-Zygmund, N=#empty arcs). But it is TOO LOSSY: the PZ bound crosses nuConsec=0.4425 only at spread
  ~60-100 (P(N>=1) at ~24+), leaving an INFEASIBLE core (~C(59,11)~10^11 clusters). ROOT CAUSE (verified):
  nu(E) = mu(E) = P(W>0) EXACTLY (W=Sum(g_i-1/7)_+ the uncovered measure), so ANY tractable lower bound on
  nu IS a moment-route bound; the 7-arc PZ is a weak degree-2 discretized moment bound, strictly weaker than
  the fleet's degree-4 moment LP (THM-661). The moment LP B_d -> mu as d->inf, with min_E B_d -> min_E mu =
  nuConsec=0.4425; the proved D3 gives min 0.3088. So Lemma A (nu>=0.4425) = the d->inf (tight) moment-LP
  minimum, and its tail IS the moment route's coupled-region difficulty -- the S184 "0.55 slack" is
  ILLUSORY for tractable bounds. No shortcut. Redirect: discharge hfloor via the PROVED D3 (THM-661) if the
  G_P-cap arithmetic allows, else prove Lemma A as a HIGH-degree moment LP -- NOT the 7-arc second moment.
tags:
  - lrc14
  - lemma-a
  - second-moment
  - moment-lp
  - honest-negative
  - density-floor
---

# The 7-arc second moment can't beat the moment route — because ν = μ

**opus-2026-07-09-S185.** Owner: prove the tail bound via the 7-arc second moment. I set it up rigorously,
computed it, and it FAILS — too lossy by a wide margin. The failure has a clean structural reason that also
corrects my own opus-S184 optimism, so it is worth recording precisely.

## The attack, made rigorous

Partition the circle into 7 fixed arcs `A_j = [j/7, (j+1)/7)`. **Arc→gap lemma (rigorous):** if `A_j` is
empty of `{frac(eᵢx)}`, the points bracketing it lie `< j/7` and `≥ (j+1)/7`, so the gap covering `A_j` has
length `> 1/7` strictly ⟹ `maxgap > 1/7`. Hence with `N(x) = #{j : A_j empty}`,
```
   ν(E) = meas{maxgap > 1/7} ≥ meas{N ≥ 1} ≥ E[N]² / E[N²]   (Paley–Zygmund).
```
`E[N] = Σ_j P(A_j empty)` and `E[N²] = Σ_{j,l} P(A_j, A_l empty)` are arc-AVOIDANCE moments (theta-sums over
the resonance lattice), tractable unlike `maxgap`. All correct.

## Why it fails: too lossy

Computed (`lrc14_lemmaA_7arc_secondmoment`), decorrelated `E[N] ≈ 7·(6/7)¹³ = 0.944`, `E[N²] ≈ 1.47`, PZ
`≈ 0.605`. But for real clusters the correlations inflate `E[N²]` and shrink `E[N]`, so:

| cluster | spread | PZ bound | `P(N≥1)` | `ν` |
|---|---|---|---|---|
| consecutive | 12 | 0.246 | 0.324 | 0.4425 |
| minhunt | 24 | 0.292 | 0.435 | 0.881 |
| minhunt | 40 | 0.393 | 0.578 | 0.985 |
| minhunt | 60 | 0.440 | 0.605 | 0.983 |
| minhunt | 100 | **0.495** | 0.642 | 0.973 |

The PZ bound clears `nuConsec = 0.4425` only at **spread ≈ 60–100**, leaving a core of spread `≤ 60`
(`~C(59,11) ≈ 10¹¹` primitive clusters) — computationally infeasible. Two losses stack: fixed arcs miss
gaps that straddle a boundary (`P(N≥1) ≪ ν`), and Paley–Zygmund is loose when the empty arcs cluster
(`PZ ≪ P(N≥1)`).

## The root cause: ν = μ, so there is no shortcut

Verified exactly (`lrc14_lemmaA_nu_equals_mu`): `maxgap > 1/7 ⟺ some gap > 1/7 ⟺ W > 0` (`W = Σ(gᵢ−1/7)_+`),
so **`ν(E) = μ(E) = P(W > 0)` identically** (0.44249 = 0.44249 at consecutive, etc.). Therefore:

> Any tractable lower bound on `ν` IS a moment-route bound on `μ`. The 7-arc PZ is a **weak, degree-2,
> discretized** moment bound — strictly weaker than the fleet's degree-4 moment LP (THM-661). The moment LP
> `B_d(E) → μ(E)` as `d → ∞`, so `min_E B_d → min_E μ = νConsec = 0.4425`. Lemma A (`ν ≥ 0.4425`) is exactly
> the `d → ∞` (tight) moment-LP minimum; the proved degree-3 `D3` gives `min 0.3088`. **Lemma A's tail is
> the moment route's coupled-region difficulty** — the same hard analytic step, reached from the other side.

This **refutes the opus-S184 redirect** ("Lemma A plausibly easier, +0.546 slack"). The slack between the
decorrelated `ν ≈ 0.988` and `νConsec = 0.4425` is real but **illusory for a rigorous bound**: every
tractable lower bound (PZ 0.25–0.6, moment LP 0.31–…) lives far below `0.988` and eats the slack. There is
no free lunch; `ν = μ` closes the door.

## The redirect (what actually to do)

- **Lemma A is unnecessary if the moment route suffices.** `hfloor` needs `witnessG2 = μ(GOOD ∩ G_P) ≥ m_P`.
  Bonferroni (PROVED, kps-S30) gives `≥ ν + measGP − 1`; the proved `D3` already gives `ν = μ ≥ 0.3088 ≫
  m_P = 0.0565`. Whether `0.3088 + cap − 1 ≥ m_P` (D3 suffices) or one needs the tight `0.4425` (Lemma A)
  is a `G_P`-cap arithmetic question — but the density floor `μ ≥ m_P` itself is ALREADY PROVED by THM-661
  independent of Lemma A. The cleanest close routes `hfloor` through THM-661, not Lemma A.
- **If Lemma A is genuinely needed**, its home is a **higher-degree moment LP** (`B_d` for larger `d`, whose
  min → 0.4425), executed on the compact core + the decorrelation tail — i.e., exactly the THM-661 / LEM-005
  machinery pushed to a tighter constant, NOT the 7-arc second moment.

## Ledger

- HONEST NEGATIVE: the 7-arc second moment gives a rigorous but too-lossy bound (`ν ≥ E[N]²/E[N²]`, crossover
  spread ~60, infeasible core). Arc→gap lemma is a clean reusable fact.
- ROOT CAUSE: `ν(E) = μ(E) = P(W>0)` exactly ⟹ every tractable `ν`-lower-bound is a moment bound; the 7-arc
  PZ is a weak degree-2 one, `<` the proved degree-4 THM-661. Lemma A = the tight (`d→∞`) moment-LP min;
  its tail = the coupled-region difficulty. CORRECTS the opus-S184 "Lemma A easier" redirect (slack illusory).
- REDIRECT: close `hfloor` via the PROVED THM-661 (`μ ≥ m_P`, Lemma A unnecessary if the cap arithmetic
  allows) or via a higher-degree moment LP — not the 7-arc second moment.
- Files: `lrc14_lemmaA_7arc_secondmoment_opus_S185` (+out), `lrc14_lemmaA_nu_equals_mu_opus_S185` (+out).
  -> HYP-5723/5722 (Lemma A, opus-S184 corrected), LEM-015, THM-661 (moment LP), LEM-005 (coupled region),
  kps-S30 (Bonferroni).
