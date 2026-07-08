---
source: opus-2026-07-07-S147
status: THM-657 lower bound PROVED (uniform, symbolic); upper-bound mechanism identified
  (4-clique tiling, verified 24 instances); Lean bridge upgraded D<=30 -> D<=75 (kernel-pure,
  green); the relation-lattice cross-side connection stated honestly (two axes, one object).
tags:
  - lonely-runner
  - liu-zhu
  - motzkin-density
  - additive-energy
  - relation-lattice
  - divisor-ladder
  - lean
---

# The Liu–Zhu block construction, and the relation lattice seen two ways

**opus-2026-07-07-S147 (HYP-5277).** Owner asked to work the S146 handoffs (general Liu–Zhu
Conj-2 proof, Haralambis past 25, χ_c(G_GW), the Lean skeleton rewiring) and to wire in the
themes the fleet's density-floor side surfaced tonight — additive energy (klein THM-656),
the degree-4 crossover (mac-mini HYP-5267), recursive combinatorics. Two solid results
landed; the honest cross-side statement is the relation lattice.

## 1. THM-657 — Liu–Zhu Conjecture 2, general lower bound (the main result)

Liu–Zhu (2004) conjectured `μ({x,y,y−x,y+x}) = (k+1)m/N` (`x=2k+1, y=2m+1` coprime,
`N=4(k+1)m+1`) and proved only `x=1`. S146 confirmed it exactly on 27 instances via the
window-graph; this session **proves the lower bound in general** by exhibiting the extremal
set in closed form:

> `B = (k+1)` blocks of length `m`, step `2m`, in `[0,(N−1)/2)`; `A = 2x·B (mod N)` avoids
> `M` with density `(k+1)m/N`. So `μ(M) ≥ (k+1)m/N` for **every** both-odd instance.

The proof is fully symbolic (verified 651 instances, and the residues have closed forms for
all `k,m`): `A − A = 2x(B−B)`, and `B−B = [−(mx−1), mx−1]` minus the odd-`m` block gaps; the
four differences map under `(2x)^{-1}` to `−2(k+1)m` (out of range), `−m` (a gap), `mx+1` and
`mx` (out of range), so none is a difference of `A`. The key identity is `y ≡ −2mx (mod N)`,
i.e. `(2x)^{-1}y ≡ −m` — the whole construction hangs on it. For `k=0` this is the classic
`x=1` slab (Liu–Zhu Thm 4.3); for `k≥1` it is `(k+1)` blocks — a genuinely combinatorial
optimum, which is *why* the slab-only methods stalled for 20 years (S146: `μ>κ`, no rotation
realizes it).

The **upper bound** has a matching mechanism: `{0, x, y, x+y}` is a 4-clique of the distance
graph `G_M` (its six differences are `x,y,x+y,y−x,y,x ∈ ±M`) — the general analog of
Liu–Zhu's `x=1` block `{i,i+1,i+2m+1,i+2m+2}`. If `Z_N` tiles into `(k+1)m` such cliques
plus one leftover (`4(k+1)m = N−1`), Haralambis's single-leftover argument gives
`μ ≤ (k+1)m/N`. Verified: the clique property holds symbolically (all `x+y≤60`), the tiling
by exact cover on all 24 instances `x+y≤32`. So **Conjecture 2 holds wherever the 4-clique
tiles `Z_N`** — proved on the tested range, lower half unconditional. The remaining general
step is a uniform tiling construction (recursive combinatorics: the clique is a
`{0,x}×{0,y}` "parallelogram," and the tiling should factor through the CRT structure of
`N`).

## 2. The Lean bridge, D ≤ 30 → D ≤ 75 (turning the AP₇₆ leaf into reach)

The S145 AP₇₆ certificate proved the `muGood(1/7)` floor to diameter 75, kernel-pure. The
skeleton's concrete-event bridge (`LRCGoodSetBridge`, death-star-S3) transferred the floor
to the skeleton's `witnessG2`-shaped `slowμ(goodSet)` object — but only to `D ≤ 30` (the
AP20/AP30 certificates). This session added the `D ≤ 75` variants
(`slowμ_goodSet_ge_mP_diam75`, `slowμ_goodSet_toReal_ge_mP_diam75`), swapping in
`AP76Certificate.hlarge_floor_diam75_unconditional` — kernel-pure, builds green
(`[propext, Classical.choice, Quot.sound]`). So the concrete `witnessG2`-ready floor now
covers the **full k=13 tail-diameter range** for bounded-diameter pure clusters. (The scout
confirmed the remaining gap to a *skeleton* discharge is de-opaquing `witnessG2`/`shapeOf`
and the arbitrary-cluster→bounded-diameter reduction — open work, not a rewiring; the safe
local upgrade is what I landed.)

## 3. The relation lattice, two axes (the honest cross-side statement)

Tonight the density-floor side named **additive energy** as its spread axis (klein THM-656:
`Var(F) ≤ E(A)·V1`; the AP maximizes `E(A)` and minimizes diameter, so it sits at the
crossing of the two floors — retro-explaining why PZ-descent bottoms at the AP). That is
factoid **D1** (the relation lattice `L(AP)` has maximal kissing number = additive energy)
made load-bearing.

The temptation is to unify it with my Motzkin-side **divisor ladder** (S146: `μ=κ` slab
generic, `μ>κ` combinatorial at the small-gcd structure). **They are different axes, and I
will not force them.** The AP is high-additive-energy AND `μ=κ` (a slab); so high energy does
*not* imply the combinatorial (`μ>κ`) regime. The honest common object is the **relation
lattice** `L(S) = {m : Σ mᵢ sᵢ = 0}` (factoid D2), viewed two ways:

- **Additive energy** = the count of its weight-4 vectors (`a−b−c+d=0`) — klein's *variance*
  axis on the density-floor side; extremized by the AP.
- **The 2-adic / small-gcd type** of `L(S)` — my *slab-vs-combinatorial* axis on the Motzkin
  side; the `μ>κ` sets are those whose difference set carries the composite/even structure
  (both-odd `{x,y,y−x,y+x}`, or `y≡0 mod 4`), i.e. a specific 2-adic relation.

This is exactly monad-S13's "one structure seen four ways" (triple atom `θ²·gcd/q`,
resonance ladder, CRT strata, coprime lens) — now with a fifth face (the Motzkin
slab/combinatorial split) and the honest caveat that these are *different projections* of
`L(S)`, not one threshold. The unifying half-page monad asked for should be written around
`L(S)` and its two invariants (short-vector count vs 2-adic type), not around a single number.

The Liu–Zhu construction itself is a relation-lattice statement: `A = 2x·B` avoids `M`
because `(2x)^{-1}M` misses `B−B`, and `B−B`'s gap structure (odd multiples of `m`) is the
2-adic signature of the `{x,y,y−x,y+x}` lattice. The `y ≡ −2mx` identity is the lattice
relation `2mx + y ≡ 0 (mod N)` — the single short vector that pins the whole construction.

## Ledger

- PROVED: THM-657 lower bound `μ(M) ≥ (k+1)m/N`, uniform, symbolic (651 instances +
  closed-form residues); the `gcd(x,N)=1` lemma; `B−B` = range-minus-odd-`m`-gaps.
- MECHANISM: `{0,x,y,x+y}` 4-clique (symbolic, `x+y≤60`) + `Z_N` tiling (24 instances,
  `x+y≤32`) = the upper bound `μ ≤ (k+1)m/N` where the tiling holds.
- LEAN (green, kernel-pure): `slowμ_goodSet_{ge,toReal_ge}_mP_diam75` — the concrete-event
  floor to D≤75 via the AP₇₆ certificate.
- STATED: the relation-lattice two-axis view (additive energy vs 2-adic type); the honest
  non-unification of klein's energy axis and my divisor ladder.
- RESOLVED CONCURRENTLY (kps-S76, THM-658): **χ_c(G_GW) ≤ 27/2 = 13.5 < 14** — the S145
  open rung question is closed by an explicit quasi-periodic two-speed (27,2)-coloring; the
  CIRCULAR rung is *blind* to GW's tightness (only the integer rung χ=14 is faithful), so
  the four rungs separate `χ_f=13 < χ_c≤13.5 < χ=14=1/M`. My coloring agent was stopped as
  redundant. Beautifully, kps's linearization-locus conjecture `χ_c=1/M ⟺ μ=M` has its
  DEFECT side (`μ>M ⟹ χ_c<1/M`) exactly on the Haralambis μ>κ locus — the same
  slab/combinatorial (`μ>κ`) territory as THM-659 and the divisor ladder. THM-652 (χ=14,
  integer faithful) + THM-658 (χ_c<14, circular blind) + THM-659 (the μ>κ construction) are
  one coherent arc at the tight locus.
- PARTIAL (honest): Haralambis |S|=3 μ=κ confirmed to max=20 (0 counterexamples, 1052 sets)
  but did NOT surpass Liu–Robinson's max≤25 — the window graph is ~λ^max; a sub-exponential
  |S|=3 engine is the named open target.
- NOT done: general Conj-2 upper-bound tiling construction; the skeleton de-opaquing.
- Files: `lrc_liu_zhu_{lower_bound,blocks,general_lb,clique_tiling}_opus_S147.py`,
  `lrc_haralambis_extend_opus_S147.py` (+outs); THM-657; the D≤75 Lean bridge lemmas.
- Builds on: HYP-5217 (S146 Conj-2 instances + divisor ladder), THM-652 (GW, |S|=13 member),
  S145 AP₇₆ Lean; fleet: klein THM-656 (additive energy), monad-S13 (the four-way note),
  mac-mini-S56 (small-gcd break). External: Liu–Zhu 2004, Haralambis 1977, Cantor–Gordon 1973.
