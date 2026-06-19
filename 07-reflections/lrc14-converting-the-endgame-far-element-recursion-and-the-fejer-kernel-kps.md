---
source: kind-pasteur-2026-06-19-S13
status: reflection on a verified structural reduction (HYP-2644) + a sanity-checked tool lead
tags: [lonely-runner, lrc14, endgame, far-element, recursion, weyl, fejer-kernel, selberg-beurling, bandlimiting, unbounded-to-bounded]
---

# Converting the LRC(14) endgame: the unbounded direction becomes a 1-D Weyl estimate

**Prompt (user):** work on converting the endgame; don't run into a wall; pull ideas from other
agents; mine older repo ideas; think outside the box.

## Where the wall was

LRC(14)-S3 had been reduced to `meas(S7(E)) ≤ cap_k` for all primitive k-sets (k=8,9,10). That splits
into a **bounded** part (the tight `consec`/near-AP, margin 0.0014 at k=9 — handled exactly by the
finite check) and an **unbounded** part (wide/large-element sets). The unbounded part is where every
analytic attempt died: the lattice-Fourier expansion `meas(S7)=M7(k)+Σ_{n∈Λ°(E)}K(n)` has an
**absolute envelope `Σ|K(n)|` that diverges harmonically** (MISTAKE-078), and the signed sum *is* the
quantity (circular). The Minkowski count was flagged but never executed; the binary GAP/stranger
dichotomy was the wrong abstraction. Wall.

## The conversion: the far-element plateau recursion (HYP-2644)

The key move is to stop trying to bound the infinite lattice sum and instead **decorrelate the single
largest element**. If `w = max E` is large, `frac(w x)` equidistributes (Weyl) and acts as an
*independent uniform fill-in* for the bounded-frequency core `E' = E∖{w}`:
```
   meas(S7)(E)  →  Plat(E') = meas(S7)(E') + (1/7)·P(E' misses exactly one sector),   rate O(1/w).
```
(Verified exact: k=9 core `consec_8`, `Plat = 0.36210`, matching the far-`w` data to 5 digits; the error
`≈ 0.5/w`.) The plateau is bounded by a **`(k−1)`-extremal quantity**
`Q(k−1) = max_{(k−1)-sets}[meas(S7)+(1/7)P(miss 1)]`, and — this is the payoff —
```
   Q(7)=0.197 < cap_8;   Q(8)=0.362 < cap_9;   Q(9)=0.448 < cap_10,   margins 0.13–0.18,
```
each attained at the bounded AP-core `consec_{k−1}`. So the unbounded direction is a **recursion DOWN in
k** to a finite check at `k−1` **with an order-of-magnitude more margin** than the tight `consec_k`. The
tightness (0.0014) lives only in the bounded finite check (done); the unbounded part is *loose*.

This is the right reframe of HYP-2610's "stranger contraction" and HYP-2637's "dissociation peel": made
exact, with the plateau `Q(k−1)`, the recursion, and the margin. **Dissociating decreases coverage** (the
AP is the most uniform orbit, hence covers all sectors most often); widening only helps the proof.

## The remaining piece, and the outside-the-box tool

What's left of the unbounded direction is a single **1-D Weyl decorrelation estimate**
`|meas(S7)(E) − Plat(E∖{w})| ≤ C/w` — one fast frequency, margin 0.13, standard flavor. Its lattice-sum
form (relations through `w`) is where the old MISTAKE-078 divergence lurks — but the **old-ideas mining**
agent found the divergence is largely an *artifact*: the per-coordinate factor `ĉ_T(n)` is a **sine
(Fejér) kernel** `~|sin(π|T|n/7)|/(π|n|)`, which is **3.5× smaller** than the crude `0.6973/|n|` envelope
and **vanishes at every multiple of 7** (the apex prime!). The honest catch (I sanity-checked it): the
sharp envelope *still* diverges logarithmically, so it is not by itself enough. The genuine convergence
comes from one of two sound mechanisms:
- **Bandlimiting:** a degree-`D` trigonometric majorant of the sector indicator forces `K(n)=0` unless
  every `|n_j| ≤ D`, **truncating the infinite lattice sum to finitely many relations** plus a
  Selberg–Beurling error `~1/D`. Infinite-divergent → finite + `O(1/D)`. (THM-537 found the *literal*
  Selberg–Beurling for the whole `S7` "doubly blocked"; the per-coordinate version through the
  inclusion–exclusion product is a different application — the block must be re-examined, not assumed.)
- **Signed cancellation** (codex HYP-2640's "signed coset quotient"): the alternating `T`-sum plus the
  relation-lattice constraint, the route codex is pushing from the reciprocal-tail side.

## The shape of the finish now

The endgame is genuinely *converted*: the part that resisted (the unbounded/wide sets) is now a recursion
to a more-margin finite check, with the only analytic residue a standard 1-D equidistribution estimate
whose tool (bandlimiting the Fejér kernel) is identified. The tight part (margin 0.0014) is the finite
check, already done. The two clocks — the apex prime `7` (which makes `ĉ_T` vanish at `7|n`, the source
of the cancellation) and the additive structure — are exactly the levers. What remains is execution
(construct the bandlimited majorant for relations-through-`w`, or the signed-coset bound), not a missing
idea. The wall has become a doorway with a known key still to be cut.
