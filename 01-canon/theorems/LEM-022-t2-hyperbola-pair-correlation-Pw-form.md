---
id: LEM-022
title: The t2 hyperbola lemma in the P(w)-parameterized form — for prime q, an interval B ⊆ ℤ/q of length b, and a unit w, the multiplicative pair-correlation C_w = #{x ∈ B : wx ∈ B} satisfies |C_w − b²/q| ≤ 6·q·(1+log₂q)²/P(w), where P(w) := min_{h≠0} |h|·|wh| is the minimal product of the ratio lattice. Subsumes klein-S226's generic √q·log²q target (P(w) ≳ q/polylog generically) and QUANTIFIES the exception set exactly: w ≡ a·b⁻¹ small-rational ⟹ P(w) ≤ ab ⟹ error budget q·log²q/(ab) — the measured near-dilation t2 enrichment. The quarantine interface becomes ONE computable number per ratio.
status: PROVED (elementary, self-contained: finite Fourier + the separation count + dyadic assembly; full proof below) + VERIFIED EXACT (integer counts, primes 149..2003, all 156 ratio-twists of generic-covering and near-dilation families, FULL w-sweeps at q ≤ 557: measured constant 0.0089 vs proved 6 — 672× slack; P(w) ≤ ab verified 340/340)
source: death-star-2026-07-09-S9 (HYP-5860; klein-2026-07-09-S226 handoff (a))
depends_on:
  - THM-680   # the LM Fourier floor this feeds (LM/q ≥ main − OffLine)
related:
  - klein-S226 (the character chain: t2 carries the deficit; the anatomy)
  - klein-S225 (the diagonal-only suppression law, confirmed at scale)
  - kps-S127 (the E3 diagonal split: the obstruction is dyadic, triangles free)
  - HYP-5845, HYP-5860
external: classical incomplete-sum technique (Pólya–Vinogradov-style completion + lattice-point separation); no citation needed — the proof below is self-contained.
---

# LEM-022 — the t2 hyperbola lemma (P(w)-parameterized)

## Statement

Let `q` be prime, `B = {r : α ≤ r ≤ β} ⊆ ℤ/q` an interval of length `b = β−α+1 < q`, and
`w ∈ (ℤ/q)^×`. Write `|h| := min(h mod q, q − (h mod q)) ∈ [0, q/2]` (the circle metric to 0),
and define the **minimal ratio-lattice product**

> `P(w) := min { |h|·|wh| : h ∈ ℤ/q, h ≠ 0 }  ∈ [1, q]`.

Then with `C_w := #{x ∈ ℤ/q : x ∈ B and wx ∈ B}`:

> **`|C_w − b²/q| ≤ 6·q·(1 + log₂ q)² / P(w)`.**

Measured truth (exact counts): the sup of `err/(q(1+log₂q)²/P)` over every tested `(q, w)` is
**0.0089** — the proved constant 6 is ~672× conservative.

## Why this form (the interface)

- **Generic twists**: `P(w) ≍ q/O(polylog)` for CF-generic `w`, so the error is `O(polylog(q))`
  — far inside klein-S226's `√q·log²q` budget. Summed over the 156 speed-ratio twists, the t2
  layer of THM-680's `LM` floor is controlled a-priori off the exception set.
- **The exception set is EXACTLY the small-P fibers**: if `w ≡ a·b̄` with small coprime `a, b`
  then `h = b` gives `|h|·|wh| = ab`, so `P(w) ≤ ab` (verified 340/340). These are the
  small-rational speed ratios of the (near-)dilated families — the quarantined mult-coherent
  branch, ALREADY owned by the coherence ladder (LEM-012 in Lean closes it on the τ-line;
  the E3/R(v) detectors see it additively). The two-coordinate quarantine of klein-S226(c)
  becomes: **additive coordinate E3, multiplicative coordinate P(w)** — each a computable
  integer per family.
- The measured near-dilation t2 enrichment (+52..57, klein-S226(2)) is the small-`P` budget in
  action: `P = ab ≤ 156` on a dilation's ratio set.

## Proof (elementary, complete)

**1. Completion.** With `B̂(h) := Σ_{r∈B} e_q(−hr)` (geometric sum), `|B̂(h)| =
|sin(πbh/q)/sin(πh/q)| ≤ 1/(2‖h/q‖) = q/(2|h|)` for `h ≠ 0`, and `B̂(0) = b`. Orthogonality
gives the standard completion:

`C_w = (1/q) Σ_k B̂(k)·conj(B̂(wk))`, so `|C_w − b²/q| ≤ (1/q) Σ_{k≠0} |B̂(k)||B̂(wk)|
≤ (q/4)·S`, where `S := Σ_{k≠0} 1/(|k|·|wk|)`.

(The map `k ↦ wk` is a bijection of nonzero residues since `q` is prime and `w` a unit; both
`|k|, |wk| ≥ 1`.)

**2. The separation count.** For `K, M ≥ 1` let `N(K, M) := #{k ≠ 0 : |k| ≤ K, |wk| ≤ M}`.

> **Claim: `N(K, M) ≤ 1 + 4KM/P(w)`.**

If `k₁ ≠ k₂` are two such solutions, `δ := k₁ − k₂ ≠ 0` satisfies `|δ| ≤ 2K` and `|wδ| ≤ 2M`
(the circle metric is subadditive). Hence `P(w) ≤ |δ|·|wδ| ≤ |δ|·2M`, i.e. `|δ| ≥ P/(2M)`:
the solutions are pairwise `P/(2M)`-separated on the circle, and they all lie in the arc
`{|k| ≤ K}` of length `2K`. An arc of length `2K` holds at most `1 + 2K/(P/(2M)) = 1 + 4KM/P`
points with that separation. ∎

**3. Dyadic assembly.** Classify `k ≠ 0` by `2^i ≤ |k| < 2^{i+1}`, `2^j ≤ |wk| < 2^{j+1}`
(`0 ≤ i, j ≤ log₂(q/2)`; at most `(1 + log₂ q)²` classes). Each class lies in the box
`(K, M) = (2^{i+1}, 2^{j+1})`, so by Step 2 its size is `≤ 1 + 2^{i+j+4}/P`, and each member
contributes `≤ 2^{−(i+j)}` to `S`. Moreover a class is nonempty only if it contains a point,
forcing `P ≤ |k|·|wk| < 2^{i+j+2}`, i.e. `2^{−(i+j)} < 4/P`. Hence

- density terms: `Σ_classes (2^{i+j+4}/P)·2^{−(i+j)} = 16/P` per class, `≤ 16(1+log₂q)²/P`;
- singleton terms: on the anti-diagonal `i+j = s` there are `≤ 1+log₂ q` classes and the
  nonempty ones need `2^s > P/4`, so `Σ ≤ (1+log₂q)·Σ_{2^s > P/4} 2^{−s} ≤ 8(1+log₂q)/P`.

Total: `S ≤ 24(1+log₂q)²/P(w)`, so `|C_w − b²/q| ≤ (q/4)·S ≤ 6·q(1+log₂q)²/P(w)`. ∎

## Remarks

1. **One-sided, absolute, classical** — no cancellation is used or needed at the per-twist
   level; the signed cancellation klein-S222/S223 proved irreducible lives in the SUM over
   twists, and this lemma confines it to the small-`P` fibers where the coherence ladder
   already operates. This is the division of labor klein-S226(3) prescribed.
2. **`P(w)` is computable in `O(q)`** per ratio (or via the continued fraction of `w/q` in
   `O(log q)`): the detector interface for the mult-coherent quarantine.
3. **Formalizability**: every object is a finite sum over `ℤ/q`; Step 2 is a Finset
   separation-pigeonhole (the same shape as `LRCLem012NearAP.exists_free_piece` and
   `LRCDiscreteBonferroni.dvd_Ioc_card_le` — arc + separation ⟹ cardinality); Step 3 is a
   double dyadic Finset sum. The only analytic ingredient is `|B̂(h)| ≤ q/(2|h|)`, i.e.
   `|sin(πh/q)| ≥ 2|h|/q` (the standard sine bound, already used across the repo's Fourier
   files). A Lean transcription is a bounded project (est. 300–500 lines), with the
   completion step (finite characters over ℤ/q) the only new infrastructure.
4. Consumption (the remaining chain to THM-680): sum this bound over the `≤ 156` twists with
   `P(w) > P₀` and route `P(w) ≤ P₀` ratios to the coherence ladder; then klein's handoff (b)
   (the t≥3 Cauchy–Schwarz bootstrap, mechanical) closes `OffLine ≤ f(E3, P₀)`.

Files: `04-computation/lrc14_hyperbola_t2_lemma_deathstar_S9.py`
(+ `05-knowledge/results/lrc14_hyperbola_t2_lemma_deathstar_S9.out`).
