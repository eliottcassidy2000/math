---
source: kind-pasteur-2026-07-07-S72
status: REDUCTION + strong evidence (NOT a proof). Genuine progress on the open route (R2 /
  density-floor AP-minimality): a per-resonance-q decomposition reduces the σ-even floor to
  finite per-q σ-odd residue statements — the first crossing of the S67 grading. Owner: work
  the open route.
tags:
  - lonely-runner
  - LRC14
  - density-floor
  - AP-minimality
  - three-gap
  - residues
  - reduction
---

# The per-q window reduction of the density floor

**kind-pasteur-2026-07-07-S72 (HYP-5117).** The open route is R2 — equivalently the
density-floor **AP-minimality** `μ_{1/7}(E) ≥ μ_{1/7}(AP)`, the one rigidity every proof-lane
bottoms out at (S70/S71). My S67 grading said this is the σ-even measure core and resists the
σ-odd (covering/parity/residue) tools. This session found the first genuine **crossing**: a
per-resonance-q decomposition that reduces the floor to *finite, per-q, residue-structured*
(σ-odd) statements. It is a reduction with strong evidence, not a proof — reported as such.

## The decomposition

A circular gap `> 1/7` requires the phases `{frac(e_i x)}` to cluster, and clustering happens
near `x = p/q` with **`q ≤ 6`** (a gap of `1/q > 1/7` needs `q ≤ 6`). Near `x = p/q` the config
collapses onto the **residues mod `q`**: `frac(e_i x) ≈ (e_i mod q)/q`. Attribute each
`maxgap > 1/7` event to the resonance `q` (the denominator of the nearest low-order rational):

> `μ_{1/7}(E) = Σ_{q=2}^{6} W_q(E)`, `W_q(E) =` measure of `x` whose responsible big gap sits
> at a denominator-`q` resonance.

This is exactly opus-S134's exact fact — `μ_{1/7}(AP_k)` is the `q ≤ 6` Farey-window measure —
extended to general families as an attribution.

## What is verified

**The AP minimizes every `W_q`.** At `k=13` the AP has `W_q(AP) = 0.065, 0.086, 0.054, 0.078,
0.016` (`q = 2..6`), and **every non-AP family has `W_q(E) ≥ W_q(AP)` at each `q`** (μ-excess
`0.27–0.53`). The only rows that read "False" are the AP's *own affine images* (spread-AP,
all-odd): their `μ` is identical to the AP's, and the per-`q` values merely relabel under
dilation — a bookkeeping artifact, not a violation.

**The mechanism, verified.** Missing a residue mod `q ≤ 6` **widens** the `q`-window
(`Δ = +0.27 … +0.33` at `k=13`, `q = 3..6`): the AP hits **all** `q` residues, so near `p/q`
its phases fall into `r = q` groups leaving a gap `≈ 1/q` (the minimal window); a non-AP that
misses residues has `r < q` groups, a gap `≥ 1/r > 1/q`, hence a wider window.

## The reduction (the crossing)

`W_q(E)` is a function of the **residue multiset `{e_i mod q}`** — a finite, algebraic,
**σ-odd** object. So

> **"AP minimizes each `W_q`" reduces the σ-even density-floor AP-minimality to per-`q` σ-odd
> residue-spread statements.**

That is the first genuine crossing of the S67 grading. S67 said the floor resists σ-odd tools
(covering/parity/residue); this does not contradict it — it *decomposes* the floor into per-`q`
pieces each of which *is* σ-odd. The value `μ` remains residue-*invisible* in total (my S65
barrier — lift families share residues, differ in `μ`), but the **extremal decomposition** is
residue-structured: the AP wins each resonance by hitting all residues.

## Honest scope, and the concrete proof program

This is a **reduction plus verified minimality, not a proof.** What remains:

1. **The per-`q` window lemma:** prove `W_q(E) ≥ W_q(AP)` rigorously for each `q ≤ 6`. This needs
   the *exact* window width as a function of the residue multiset — not just the deficit
   (residue-full non-AP families still have `W_q > W_q(AP)` via the *within-residue-class
   spread*, so the lemma must control both the residue miss-count and the intra-class drift as
   `x` leaves `p/q`).
2. **Dilation bookkeeping:** `W_q` is not individually dilation-invariant (only the sum `μ` is).
   Since `μ` is dilation-invariant, fix the *primitive* representative (gcd of differences 1)
   and prove the per-`q` inequalities there.
3. **The `q = 7` borderline:** the AP's gap near `p/7` is `≈ 1/7` exactly — the edge case needs
   its own (finite) treatment.

Concrete program: prove the per-`q` residue-spread window lemma for `q = 2..6` (a finite
residue-optimization: maximally-spread residues minimize the window) plus the `q = 7` edge
`⟹` AP-minimality `⟹` R2 `⟹` `(A')`. Each per-`q` lemma is a *finite* statement about residue
multisets — the kind the σ-odd machinery (my S65 image census, the sieve, covering) can attack,
unlike the monolithic σ-even floor.

## Why this is the right handle

Every prior route (2-anchor stability, 2-anchor decorrelation, μ net-route) hit the *whole*
σ-even rigidity at once and stalled (S71). The per-`q` decomposition is the first that **splits
the rigidity into pieces matched to the tools that work** — finite residue statements per
resonance. It is the mac-mini-S15 three-gap frame ("`μ` is quantitative three-gap rigidity")
turned from a slogan into a decomposition with an explicit extremal mechanism (residue spread)
and a finite proof target per `q`.

## Ledger

- Verified: `μ = Σ_q W_q`; the AP minimizes each `W_q` (k=8, 13); missing a residue mod `q ≤ 6`
  widens `W_q` (`Δ ≈ 0.3`).
- Reduction: σ-even AP-minimality `⟶` per-`q` σ-odd residue-spread lemmas (`q = 2..6`) + `q = 7`
  edge — the first crossing of the S67 grading.
- NOT proved: the per-`q` window lemma (exact widths + within-class spread), dilation
  bookkeeping, the `q = 7` edge. A concrete finite proof program remains.
- Files: `lrc_per_q_window_kps_S72.py` (+out).
- Builds on: mac-mini-S15 (three-gap frame), opus-S134/THM-637 (roof = `q ≤ 6` windows), S65
  (residue image), S67 (σ-grading, now crossed), S70/S71 (R2 = the target), opus-S131 (AP-min
  verified `k ≤ 10`).
- Does NOT prove LRC(14) or R2. It gives the first σ-odd decomposition of the σ-even floor and a
  finite per-`q` proof program.
