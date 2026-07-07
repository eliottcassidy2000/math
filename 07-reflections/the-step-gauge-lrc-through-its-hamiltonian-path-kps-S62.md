---
source: kind-pasteur-2026-07-07-S62
status: FRAME (owner directive) + PROVED (reversal invariance; the Wiener staircase identity)
  + VERIFIED (record-board palindromy; a=2 grid closure at large st; corner sampled 9.8x bar)
  + OPEN (palindromic-extremizer conjecture; the finite a=2 corner; a>=3 strata).
tags:
  - lonely-runner
  - LRC14
  - step-gauge
  - tiling-model
  - complement-symmetry
  - staircase
---

# The step gauge: LRC(14) viewed from its Hamiltonian path

**kind-pasteur-2026-07-07-S62 (HYP-4877).** Owner: *a tree on 8 events has 7 edges and that
is just the Hamiltonian path connecting each element in an 8 player tournament; the tiling
model leverages viewing a tournament from one of its Hamiltonian paths to reveal the symmetry
in the isomorphism class structure forced by the nature of intersubjective binary relation
itself; apply similar analysis to lonely runners.*

Done. The dictionary is exact, two entries are theorems, and the frame reorganizes the
open residual along a symmetry filtration.

## The dictionary

LRC(14) has **14 runners** — 13 moving plus the origin — and 14 nodes need 13 edges: sort the
speeds `0 = v₀ < v₁ < … < v₁₃` and take the path `0 → v₁ → … → v₁₃`. That path is the **base
Hamiltonian path**; its 13 edge weights are the **consecutive differences**
`δᵢ = vᵢ − vᵢ₋₁ ≥ 1`, and the family is their partial-sum sequence. Everything the density
floor sees factors through this gauge:

| tournament (n players) | LRC(14) (14 runners) |
|---|---|
| base Hamiltonian path, n−1 arcs | sorted-speed path, 13 steps `δ` |
| free tiles = the other arcs | the step word `(δ₂,…,δ₁₃)` (12 letters) |
| relabeling gauge | `δ₁` = translation gauge (μ ignores it); common scaling = dilation gauge |
| **complement `T ↦ T^op`** | **step reversal `= E ↦ max(E)+min(E) − E`** |
| self-complementary classes | **palindromic step words** |
| cut ⊕ cycle (GF(2)) | pair-marginals ⊕ weight-≥3 relations (Fourier) |
| the staircase Δ_{n−2} | the Wiener profile `l(14−l)` (below) |

**Theorem (reversal = complement).** `μ_θ` and `E[maxgap]` are invariant under reversing the
step word, i.e. under `E ↦ max(E)+min(E) − E`. *Proof:* the affine map `e ↦ −e + c` composed
with `x ↦ −x` fixes every `frac(ex)` configuration up to reflection of the circle, which
preserves all gap lengths. (Verified exactly on the zoo: diffs 0.00000 including
non-palindromic words.) This is boxeph-S1's `x → −x` reflection, now recognized as the
**complement involution of the step gauge** — the LRC moduli space for the density floor is
**12-letter positive words modulo scaling and reversal**, precisely the staircase-mod-complement
shape of the tournament project.

**Theorem (the crossing staircase).** Pair `(i,j)` of runners (including the origin-runner)
crosses exactly `|vᵢ − vⱼ|` times per period (the S59 difference-uniformity factoid, counted);
the total number of order-cells charged to the step at sorted position `l` carries the
**triangle weight `l(14−l)`** — `[13, 24, 33, 40, 45, 48, 49, 48, 45, 40, 33, 24, 13]`,
summing to `455 = C(15,3)`, the Wiener index of the path on 14 nodes. Middle steps generate
~3.8× the complexity of end steps — the same center-heavy asymmetry as the tournament tiling
model's wiggly classes ("the center vertex is ~5× more neutral"). *The staircase is the
discrete parabola `l(n−l); everything is the triangle, again.*

## The deep parallel the owner pointed at

A tournament is a *complete* system of pairwise binary relations, and the tiling model shows
its isomorphism structure lives not in the pair data (scores/cut) but in the **cycle space**.
LRC's pairwise data is not merely weak — it is **provably featureless**: every pair distance
`‖(vᵢ−vⱼ)x‖` is *exactly* uniform in `x`, and (S59 deficit frame) every deviation of every gap
functional from the iid baseline is carried by **zero-sum relations of weight ≥ 3**. The
elementary such relations in the step gauge are **step equalities**: `δᵢ = δⱼ` ⟺
`vᵢ − vᵢ₋₁ − vⱼ + vⱼ₋₁ = 0`, a zero-sum weight-4 vector. So:

> **The step-alphabet size `a(E)` = #distinct step values is the LRC analogue of the
> cycle-space filtration.** `a = 1` is the AP — all steps equal, the maximal-relation,
> gauge-degenerate point (the analogue of the transitive tournament, and like it, the unique
> extremal class). Small `a` = high symmetry = computable superset (grid) — the ledger side.
> Large `a` = few relations = near-iid tails — the decorrelation side. The honest frontier is
> the middle, exactly as the tournament metagraph's complexity lives in its middle classes.

## What it yields today

1. **Palindromic-extremizer evidence (conjecture, OPEN).** Every family on the exact record
   board is a palindromic word: the AP `(1¹²)`, death-star's prim-sat `(2⁵1²2⁵)`, monad's
   record `(2⁴1⁴2⁴)`; the non-palindromic GW and opus-stretch were both later *beaten* by
   palindromic families. Conjecture: per diameter, `min μ` and `min E[maxgap]` are attained
   on the reversal-symmetric locus. My quick descents were too noisy to confirm (the
   palindromic-restricted search repeatedly found *lower* values than the unrestricted one —
   search variance, reported honestly, not evidence against); the record board is the real
   signal. If true, searches halve and extremal certificates gain a symmetry to exploit —
   the exact echo of the tournament project's SC-maximizer theorem.
2. **The `a = 2` stratum, mostly closed.** A two-letter word (`i` copies of `s`, `12−i` of
   `t`, any order) lies in the `(i+1)×(13−i)` grid — `N = (i+1)(13−i) ≤ 49` for *every*
   order and *every* diameter. The S61/S62 2-torus table gives `μ₂(split) ≥ 0.118 ≥ 2.1·m_P`
   at the worst split (7×7). Assembly: `i ∈ {0,12}` = dilated AP, exact; diam ≤ 75 = S59
   exact; diam ≥ 76 with `st ≳ 80` = grid ledger with geodesic error `C/(st) <` margin.
   **Residue: a finite corner** (`t` small, `st ≲ 80`, diam ≥ 76) — sampled directly at
   `μ ≥ 0.555 = 9.8× m_P` (the grid bound is lossy there, the truth is comfortable); the
   corner is specified for an exact sweep (finitely many `(i,s,t)` × word orders), not closed.
3. **The `a ≥ 3` strata** are the genuinely sparse lane (balanced 3-letter words need
   `>78`-point grids, where the `4.3/N` law crosses the bar) — consistent with S61's
   coverability stratification; the step gauge explains *why* those families were the
   no-cover ones.

## Pointers

- Files: `lrc_step_gauge_kps_S62.py` (+out with corner addendum).
- Builds on: S59 (difference uniformity, deficit frame, diameter floor), S61 (grid ledger),
  boxeph-S1 (x→−x), THM-637/opus-S135 (roof proved), the tournament tiling canon
  (CLAUDE.md: cut⊕cycle, staircase, SC classes).
- OPEN: the palindromic-extremizer conjecture (needs exact per-D sweeps or a symmetrization
  inequality); the finite `a=2` corner; `μ₃` tables for balanced `a=3`; whether reversal
  symmetry can be exploited in the PZ/CE moment machinery (reversal-symmetrize `U`).
- Does NOT prove LRC(14).
