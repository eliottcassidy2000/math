# The Bridge-2 mindset is one rung of a universal moment ladder: identification climbs with n, the families are orthogonal, and the "reconstruction wall" is the ladder saturating

*kind-pasteur-2026-07-01-S18. The owner asked to think about the Bridge-2 mindset — "the 2nd moment separates what the 1st moment can't" — and hunt for it recursively. Quantified across n, it is not a one-off but a **ladder**: invariants come in orthogonal families (count/combinatorial and spectral/moment), each saturates, adding the next buys exactly one more `n`, and the whole thing keeps climbing. The reconstruction wall (S16) and the LRC Lasserre hierarchy (HYP-3789) are the same recursion seen in two worlds.*

## The ladder, quantified (tournament identification by n)

Resolution = `#distinct invariant-values / #iso classes` (`1.000` = identifies all). An increasing ladder of *intrinsic* invariants:

| rung | n=4 | n=5 | n=6 | n=7 |
|---|---|---|---|---|
| score sequence | **1.000** | 0.750 | 0.393 | 0.120 |
| + H (OCF count) | 1.000 | 0.917 | 0.714 | 0.531 |
| + cycle spectrum (c₃,c₅,c₇) | 1.000 | 0.917 | 0.839 | 0.746 |
| **+ d = det(I+S)** | 1.000 | **1.000** | 0.893 | 0.808 |
| + cpA (adjacency spectrum) | 1.000 | 1.000 | 0.893 | 0.822 |
| + cpS (skew spectrum) | 1.000 | 1.000 | 0.893 | 0.822 |

Read down each column, then across:

- **The identification rung climbs with n.** n=4: *score* alone identifies. n=5: the count family saturates at `0.917` (one OCF-cospectral pair) and the **spectral `d` rescues it to `1.000`** — the Bridge-2 move, live. n=6,7: even the full spectrum saturates (`0.893`, `0.822`) — you must climb past intrinsic invariants (to metagraph structure / WL / canonical form).
- **The spectral family buys exactly one `n`.** Adding `d` pushes identification from n=4 (score) to n=5 (spectrum). It does not reach n=6. The 2nd moment is worth *one level of n*, no more.
- **The families are orthogonal.** At n=7, count-only resolves `0.746`, spectral-only `0.376`, but *combined* `0.822` — each family separates classes the other cannot. This is the determinant-lens `d ⊥ H` quantified as *complementary blind spots*, not redundancy.
- **`d` ≈ the whole spectrum (for discrimination).** cpA and cpS barely improve on `d` (`0.893→0.893` at n=6, `0.808→0.822` at n=7). A single number (the determinant) captures almost all the spectral discriminating power — the coarse-but-powerful coordinate of `the-determinant-lens`.
- **The wall is a descent, not a cliff.** Resolution falls monotonically (`1.000, 1.000, 0.893, 0.822`). "The reconstruction wall at n=7" (S16) is really the *continuous saturation* of the intrinsic ladder — n=5 is the last n a bounded-moment invariant fully identifies.

## The recursion (Bridge-2 is one rung)

The pattern generalizes to a filtration by "correlation order," and it recurs:
1. **Order 1 (count / mean / trace):** `H`, `|L_S|`, degree, score — the orbit-symmetric part. *Blind first* (many objects share it).
2. **Order 2 (spectrum / Gram / variance / determinant):** `d`, `cpA`, the danger-Gram `DDᵀ`, `Var(H)`. Separates order-1-cospectral objects. *This is Bridge-2.*
3. **Order ≥3 (higher moments / WL / canonical form / Lasserre level ≥3):** what you climb to when order-2 saturates.

Each rung works until it saturates, then the next rung buys roughly one more level of the size parameter. The *difficulty* of a problem = **the rung you must climb to**.

## The same recursion in two worlds

- **Tournament identification:** score → H → spectrum → WL. Identifies at rung 0 (n=4), rung 2 (n=5), rung ≥3 (n≥6). (This reflection.)
- **LRC covering-min = the Lasserre / trig-moment hierarchy (HYP-3789):** level 1 = union bound (1st moment, `T₁<0`, blind); level 2 = pair correlations (2nd moment, the S75 signed correction, the near/far split); level ∞ = exact (lazy-cut). Converges only at level `m=|S|`. The covering bound *is* the moment ladder; the difficulty is that no *finite* level is exact for an atomic extremum.
- **mac-mini Bridge 2 (LRC champions):** trace of the danger-Gram `= |S|·2r` (1st moment, cospectral for construction/AP/GW) but distinct Gram spectra (2nd moment, top eig `0.5006/0.4992/0.4803`). Exactly rung-1-cospectral / rung-2-distinct — the LRC twin of my tournament twins.
- **THM-589:** `H` mean (rung 1) vs `CV(H)²` = H-variance (rung 2). The 2nd moment of the OCF is a genuine new invariant.
- **the-determinant-lens floor/ceiling/AVERAGE:** the "average" face `E[det(I+S)] = A000085` *is* a moment — the ensemble 1st moment — sitting beside the extremal (spectral) faces.
- **Verblunsky `|α_k|` (HYP-3795):** the OPUC coefficients are literally the moment ladder of the spectral measure; the AP runner-cloud has the harmonic ladder `1/(n-1-k)`.

## The subtle recursive twist

The ladder does not merely climb — **the low rungs' blind spots are structured**, and their structure is the next rung's content. The order-1-cospectral tournament twins are exactly the ones the order-2 spectrum resolves; the order-2-cospectral classes (n≥7) are exactly the ones WL/metagraph resolves; and WL itself is *iterated* order-2 (neighbor-multiset = local 2nd moment, applied recursively) — which is why WL-from-category is already discrete (S14). **So the ladder is self-similar: order-2 applied recursively (WL) is the whole ladder.** The Bridge-2 move, iterated, is reconstruction.

## What to keep hunting (the mindset)

Whenever a **count / trace / mean** invariant is cospectral, reach for the **Gram / spectrum / variance / pair-correlation** — and expect it to buy about one level before it too saturates, at which point iterate (WL / higher Lasserre). Concrete next hunts: the LRC atom pair-correlation (Ramanujan sums, HYP-3793) as the atoms' 2nd moment; the SC-spine's blue-degree (rung 1) vs its blue-graph spectrum (rung 2); and whether the LRC "true spectral twins" (Gram-cospectral champions, mac-mini's next target) force a rung-3 (3rd-moment) LRC separator — the direct analog of the tournament n=7 wall.

## Honest status

- **Computed:** the resolution table (n≤6 exhaustive, n=7 ~sample 433/456); the identification-rung climb; orthogonality of the count and spectral families; `d ≈ full spectrum`; the continuous descent.
- **Synthesis:** the moment-ladder framing unifies the reconstruction wall, HYP-3789 (LRC Lasserre), THM-589, the determinant-lens, and Verblunsky as one recursion — a lens, not a theorem.
- **Refinement of S16:** the intrinsic-invariant ladder starts saturating at **n=6** (50/56), one earlier than the `(OCF,d)`-with-metagraph-degrees fingerprint (which held to n=6) — the metagraph blue/black degrees are a *third* orthogonal family that buys the last bit at n=6.

— Related: `does-H-close-reconstruction-…`, `the-reconstruction-key-is-OCF-and-determinant-…`, `the-quarter-tiling-…` (the reconstruction wall), HYP-3789 (LRC Lasserre/moment hierarchy), HYP-3813 (mac-mini-S87 Bridge 2), THM-589 (H-variance), `the-determinant-lens-sgn-vs-chi-…` (d⊥H, average face), HYP-3795 (Verblunsky ladder), `lrc-difficulty-by-n-…` (S18). Script: `04-computation/moment_ladder_resolution_by_n_kps.py` (+ .out). Not a HYP reservation.
