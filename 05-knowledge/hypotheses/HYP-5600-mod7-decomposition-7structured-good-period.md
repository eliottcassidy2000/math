---
id: HYP-5600
title: The mod-7 / gcd(7,Vmax) decomposition of the 7-structured dissociated good-period residual — the last covering-case gap reframed in the Z/7-coloring framework
status: PARTIAL (verified, not proved). The two resonance mechanisms (missed-residue, ÷7-collapse) are EXACT; the gcd(7,Vmax) split is verified (7|Vmax: m/7 grid point good 83/85; gcd=1: existence 100% of sample); residual = the gcd(7,Vmax)=1 general argument + the "covers-all-residues, |S7|<k-6" corner. Reframes the residual in the mod-7 arithmetic (mac-mini-S6/HYP-2703), NOT closed
source: mac-mini-2026-07-09-S64
depends_on:
  - LEM-013   # the dissociated good-period existence (maxgap margin); this gives its arithmetic structure
  - THM-530   # the k<=7 pigeonhole (maxgap>=1/7 for <=7 phases) -- the ÷7-collapse mechanism
related:
  - MISTAKE-128  # why moments/arc-count fail for 7-structured (this reframes the fix)
  - THM-536   # Sturmian partial sums (the tool the Z/7 coloring needs)
---

# HYP-5600 — the mod-7 decomposition of the 7-structured good-period residual

## The residual and the reframe

LRC(14)'s covering case is one lemma from done: for **dissociated** clusters (longest-AP `≤ k−6`) a
good period exists. The hard sub-case is **7-structured** sets (many co-offsets `≡ 0 mod 7`), where
the arc-count `c` spikes and moments `B_3..B_6` all stay below it (MISTAKE-128, mac-mini-S63) — no
moment `μ`-floor. This hypothesis reframes it in the project's **Z/7 vertex-coloring** framework
(mac-mini-S6 / HYP-2703: `c(e,x)=⌊7·frac(ex)⌋`, slope-band `measS7 = (1/7)Σ_s bandcover_s`): the
7-structured sets are the ones *arithmetically special mod 7*, and the good period is decided by the
mod-7 residues.

## The two exact resonance mechanisms (at `x = m/7`)

At `x = m/7`, `frac(e_i m/7) = (e_i m mod 7)/7` — all phases on the 7-grid. So:
- **Missed residue:** `{e_i mod 7} ≠ ℤ/7` (some residue `r` absent) ⟹ the slot `[r/7,(r+1)/7)` is
  empty ⟹ `maxgap ≥ 2/7 > 1/7`.
- **÷7-collapse ⟹ k≤7 pigeonhole (THM-530):** `|S₇| := #{e_i ≡ 0 mod 7}`; these all collapse to phase
  `0` at `x = m/7`, so the effective count is `k−|S₇|+1`. If `|S₇| ≥ k−6`, effective `≤ 7` ⟹ THM-530
  gives `maxgap ≥ 1/(k−|S₇|+1) > 1/7`.

Verified on the key sets: MISTAKE-128 `{0,7,14,21,26,29,37,44,51,58,67,75,82}` misses `{3,6}`
(`maxgap=2/7` at `x=1/7`); klein's worst has `|S₇|=10 ≥ 7` (`maxgap ≥ 1/5`). Over 652 sampled
7-structured dissociated `k=13` sets, **628 satisfy `misses-residue OR |S₇|≥k−6`**; the 24 that fail
both cover all residues with few `÷7` elements (not the hard high-`c` sets).

## The `gcd(7, Vmax)` split — the good PERIOD

- **`7 | Vmax`:** `x = m/7 = (m·Vmax/7)/Vmax` **is a grid point**, and the mechanism above makes it
  good (margin `≥ 1/5`). Verified: `j = m·Vmax/7` good for some `m` in **83/85** sampled sets.
- **`gcd(7, Vmax) = 1`:** the grid never hits a `1/7`-resonance — it decorrelates and samples the good
  bulk (`μ ≈ 0.94`). Verified: existence in **79/79** sampled sets.

So the 7-structure — which broke the arc-count and the moments — *hands over the good period* when
`7 | Vmax` (a wide good arc sits on the grid at `m/7`), and is invisible when `gcd(7,Vmax)=1`.

## Sub-result: the first-moment-vanishing path is RULED OUT (mac-mini-S64)

A natural hope: since LEM-011's `𝒲̂(n)` vanishes when 7-commensurate (`7|n_i` or `7|N`), maybe for
7-structured sets the surviving Weyl corrections are small/nonnegative, so
`E_grid[W] = (6/7)^k + Σ corrections ≈ (6/7)^k > 0` — forcing a good period with **no moments**
(`first-moment > 0 ⟹ a positive summand ⟹ a good `j`). **Tested decisively and FALSE**
(`04-computation/lrc14_first_moment_vanishing_macmini_S64.{py,out}`, exact `Fraction` grid averages,
300 sampled 7-structured dissociated `k=13` sets):
- The corrections are **genuinely negative**: 71/134 (`7|Vmax`) and 109/162 (`gcd=1`) sets have
  full-grid average **below** `(6/7)^k = 0.1348`. The vanishing does *not* leave only helpful terms;
  it drags the average down to `≈ 0.035` (`7|Vmax`) / `0.089` (`gcd=1`).
- `min E_grid[W]` over real periods `j=1..V−1` is `> 0` (worst `0.0179` at `7|Vmax`, `0.0674` at
  `gcd=1`) — but that is *equivalent* to "a good period exists," so it is no easier a-priori.
- **The mod-7 picture is confirmed, though:** the worst `7|Vmax` case is `V = 49 = 7²`, where the
  single carrying good period is a **multiple of 7** (`x = m/7`, the resonance). The thin first
  moment is *carried by the resonance itself*. So the route is the **pointwise** resonance argument
  (`x = m/7`), NOT any averaged / moment bound — the same lesson as MISTAKE-128, one level up.

## RESOLUTION of the dissociated 7-structured tight case (mac-mini-S64, non-strict)

The 7-structured "hardness" (arc-count spike, moments dead, `|R|/lead → 0.87` on the resonant grid)
was a **strict-vs-non-strict artifact**. LRC(14) loneliness is `M ≥ 1/14 ⟺ maxgap ≥ 1/7`
(**non-strict** — equality gives `M = 1/14` exactly, which satisfies the conjecture). Scoring the
DISSOCIATED (longest-AP `≤ k−6`) 7-structured sets by the exact loneliness margin
`m = max_j(maxgap·7 − Vmax)` on the resonant grid `7∣Vmax`
(`04-computation/lrc14_nonstrict_knife_edge_macmini_S64.{py,out}`) gives, bucketed by spread:
- `spread < 6·Vmax/7`: min margin `= 14` (strictly lonely at `j=1`);
- `spread = 6·Vmax/7`: min margin `= 0` — the **UNIQUE knife-edge** (`j=1`: phases fill `[0,6/7]`,
  wraparound gap `= 1/7` exactly, `M = 1/14`); the extremal sets all have `spread = 6·Vmax/7`, cover
  all 7 residues, `|S₇| = k−6 = 7` (the `S₇` elements form a length-7 AP of step 7);
- `spread > 6·Vmax/7`: min margin `= 77` (**comfortably** lonely at some other `j`).
- **ZERO counterexamples** (`m < 0`) anywhere.

So the dissociated branch's tight case is the wraparound boundary, closed by the **non-strict `j=1`
wraparound lemma** (`good_period_j1_wraparound_nonstrict`, Lean, sorry-free: `7·spread ≤ 6·Vmax ⟹
gapLen ≥ 1/7`). The `spread > 6·Vmax/7` sub-case is comfortable (margin 77 ≫ 0) — kps-S97's kissing
route / opus-S170 smooth route apply there with ample room; the `|R|/lead=0.87` scare was the
knife-edge (`strict-W=0`) leaking into the strict-`W` average. This does NOT contradict klein-S201's
MISTAKE-129 (the tight *AP* `{0..12}` at `V=13` has `maxgap=1/13 < 1/7` genuinely — density-floor
territory, `V ≤ Q`); my knife-edge is `maxgap=1/7` *exactly* — good-period territory via `j=1`.
The earlier `gcd(7,Vmax)` split is superseded by the sharper `spread` vs `6·Vmax/7` threshold.
See reflection `the-nonstrict-criterion-dissolves-the-7structured-hardness-macmini-S64`.

## Honest status / what remains

Not a closed proof. The exact mechanisms + the `7|Vmax` closure are solid; the residual is (a) a clean
a-priori argument for the `gcd(7,Vmax)=1` branch (currently the general bulk `μ`, 100% sampled), and
(b) the `7|Vmax` "covers-all-residues, `|S₇|<k−6`" corner (a low-`c` non-hard set, general argument).
Value: it **reframes the last covering-case residual into the project's Z/7-coloring / slope-band /
`k≤7`-pigeonhole framework**, splits it by the one arithmetic quantity `gcd(7,Vmax)`, and reduces the
hard `7|Vmax` case to a one-line resonance pigeonhole. The moment route died by ignoring the `7`; this
uses it. Files: `04-computation/lrc14_mod7_decomposition_macmini_S64.{py,out}`. See reflection
`the-7-structured-mu-floor-is-a-Z7-coloring-problem-macmini-S64`.
