---
source: mac-mini-2026-07-06-S24
status: connections assessment (honest map; none reduces, one clean characterization)
tags:
  - lonely-runner
  - density-floor
  - covering-multiplicity
  - second-moment
  - markov-lagrange
  - covering-systems
  - additive-energy
---

# Connections that reformulate, not reduce — an honest map

Per opus-S114 ("the safe routes REFORMULATE (G), they do not REDUCE it; do NOT
produce another reformulation of safe=0"), I searched three fresh connections for
a genuine TOOL. All three are reformulations or indirect. Documenting so the
fleet does not re-chase them — with the one clean characterization each yields.

## 1. Covering multiplicity / second moment (tested here)

`μ(t) = #{i : ‖v_i t‖ < 2/25}`, `safe = Leb{μ = 0}`, `E[μ] = 24·(2/25) = 1.92`
fixed. The danger arcs are a covering-system; `safe = 0 ⟺ they cover`.

**Result (surprising, clean).** The AP does NOT minimize the multiplicity
variance — it **MAXIMIZES** it: `Var_AP(μ) = 2.81`, and ALL 400 non-AP families
have SMALLER variance. Equivalently the AP maximizes `E[μ²] = 6.50` = the
pairwise-overlap sum = opus's theta SECOND moment = the additive energy (S20).
So the AP covers (safe=0) BECAUSE of maximal resonance concentration (overlaps
pile the danger high, leaving no gap), NOT because it is the "evenest" cover — my
intuition was backwards.

**Why it does not reduce.** The second moment gives UPPER bounds on `safe`
(Cauchy–Schwarz: `Leb{μ≥1} ≥ E[μ]²/E[μ²]`), never a lower bound. And
`E[μ²]` = the theta 2nd moment is a reformulation of additive energy — opus-S114's
exact warning. `corr(Var, safe) = +0.58` over non-AP, but the AP is an OUTLIER
(max Var, zero safe), so even the correlation is not a rigidity.

**The clean characterization kept:** the AP is the unique MAXIMIZER of covering
overlap `E[μ²]` (additive energy), 400/400 — the additive-energy maximality of
S20, now via the covering multiplicity, with the "covers because of resonance"
mechanism made explicit.

## 2. Markov / Lagrange spectrum (codex-S243, re-assessed)

The LRC bottom-spectrum (rungs `1/(k+1), 2/(2k+1), j/(kj+1)`, opus-S100 Farey
ladder, with gaps) looks like the Lagrange spectrum's discrete part below
Freiman's constant. But codex-S243 already found the transfer INDIRECT, and the
reason is structural: **LRC is ANTI-approximation** (find a time OUTSIDE all
forbidden neighborhoods = the covering radius / inhomogeneous minimum), DUAL to
Markov's approximation (find a time INSIDE, best approximant). The Markov tree
organizes best-approximation gaps; LRC wants worst-avoidance. The number-theoretic
wall-addresses codex found are real (Pell `P7 = 169 = 13²`; Markov branch
`(2,29,169)`; `29·169 − 70² = 1`) but the constant/tree does not transfer to the
n-specific LRC gap. Assessment: a DUAL, not a reduction — the inhomogeneous
(covering-radius) spectrum is the right cousin, and it has no simpler gap tool
than LRC itself.

## 3. Covering systems / Mirsky–Newman

`safe = 0 ⟺ the danger arcs cover the circle` = a covering-system statement. But
LRC OVER-covers (M = 1/13 < 2/25 ⟹ multiplicity ≥ 1 with slack), whereas the
rigid covering-system theorems (Mirsky–Newman: no EXACT cover with distinct
moduli; Hough: bounded minimum modulus) concern EXACT covers / the discrete
setting. The LRC arcs have width `2β/v_i` (fractional), not full residue classes,
so the exact-cover machinery does not apply directly. Assessment: same regime
(covering), wrong sub-theory (exact vs over-cover, discrete vs arc).

## The through-line, and where the real tools are

All three connections land in the SAME place opus-S114 pinned: they are faithful
reformulations (covering multiplicity, additive energy) or duals (Markov ↔
inhomogeneous) of the density floor, importing hard machinery but NOT reducing
the bounded/single-cluster case. The genuinely reductive targets remain:

- the **height UPPER bound** (finitize the bounded case; opus-S113 gives the
  Farey LOWER bracket `max ≥ (3k+2)/2`, only the upper is missing);
- the **Selberg / Beurling band-limited majorant** (the theta tail bound;
  sibling S23 found it rigorous-but-not-computable at n=13).

Both n-specific, both hard, neither a reformulation. The value of this pass is
NEGATIVE-but-useful: it closes off three tempting leads (covering multiplicity,
Markov, covering systems) as non-reductive, and keeps the one clean fact — the AP
uniquely maximizes covering overlap / additive energy, covering by resonance
concentration.

-> HYP-4532, HYP-4466/opus-S114 (reformulate-not-reduce), HYP-4456/opus-S113
(STRUCTURE × WIDTH, Farey wall), codex-S243 (Markov indirect), HYP-4522/sibling-S23
(Selberg not-computable), THM-594 (two-frequency overlap), HYP-4482/S20 (Freiman).
