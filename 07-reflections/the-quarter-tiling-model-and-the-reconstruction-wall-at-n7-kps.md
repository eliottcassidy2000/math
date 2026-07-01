# The quarter-tiling model (Q_m/⟨σ,φ⟩) and the reconstruction wall at n=7: local fingerprints fail, structure recurses

*kind-pasteur-2026-07-01-S15. Working the open next ideas plus the owner's "quarter-tiling model." Three results: (1) the (OCF, determinant) reconstruction key **fails at n=7** — local invariants degrade sharply; (2) the blue/SC spine is **connected** (upgrading the parity min-degree bound); (3) the **quarter-tiling** = the tiling cube folded by BOTH complements `⟨σ, φ⟩`, a clean Klein-4 quotient whose orbits are the blue lines (self-complementary) and the black quads (complement-paired). The lesson: past n=6, reconstruction must lean on the recursive complement structure, not on richer node fingerprints.*

## The reconstruction wall at n=7

At n=6, `(I(Ω,x), d=det(I+S))` was a complete node-fingerprint (the OCF-cospectral twin `{4,6}` was `d`-visible). **This breaks hard at n=7.** Sampling 1600 random tournaments → 427 of the 456 classes (94%):
- **`I(Ω,x)` alone is cospectral for 90%** of sampled classes (385/427 share their OCF poly with another) — unsurprising, since at n=7 the poly is just `1 + i₁x + i₂x²` (`i₃=0`, needs 9 vertices), two integers of limited resolution.
- **`(I(Ω,x), d)` still leaves 142 collision groups** — the determinant no longer rescues the cospectral classes.
- Adding the skew spectrum `cpS` does **not** close it either (the collided classes are also coskew-spectral).

So `(OCF, determinant)` is a complete fingerprint **only through n=6**. The reconstruction obstruction — the twins — proliferate at n=7 into large cospectral families no bounded classical invariant we tried separates. **Local-invariant reconstruction hits a wall at n=7.** (This is the sharpest form of the metric-suite conclusion: the "reconstruction defect" is not just nonzero but grows past the reach of node fingerprints.)

## What still works: the blue/SC spine is connected

Against that backdrop, a *structural* invariant survives. The odd grid-sym parity gave "blue-other-degree ≥ 1 for every SC node" (no blue-isolated SC class). Computing the blue graph on SC nodes directly: **it is connected at n=4,5,6.** So the min-degree bound upgrades to **connectivity** — the SC classes form a single blue-connected backbone (the "spine" is genuinely one component, not merely edge-covered). This is exactly the kind of `n`-robust structural fact that fingerprints lack.

## The quarter-tiling model — an algebraic double-complement quotient, not a second geometric fold

**Correction / convergence (opus-S20, HYP-3811):** the *geometric* `σ`-fold is a **single terminal step** — the half-tiling's folded region has no second complement-reflection, so there is **no geometric quarter-tiling** (and the folded class-count is A051337, cycles ≡ 2 mod 4). That is correct and stands. What follows is therefore **not** a geometric refold but an **algebraic** construction: quotient `Q_m` by the second, *non-geometric* involution `φ` (complement *tiling* = flip-all, the blue-line pairing). It is a valid `Z₂×Z₂` quotient with real structure — just not a "tiling model" in the staircase sense. I keep the name in quotes for continuity but the object is the **double-complement quotient**.

The half-tiling folds the tiling cube `Q_m` by `σ` (grid reflection = tournament complement `φ(T^op)`), giving the `D=⌊(n-1)²/4⌋`-cell folded region and the SC/blue structure. The **double-complement quotient adds `φ`** (complement *tiling*, flip-all) — the two commuting involutions generate the **Klein 4-group `⟨σ, φ⟩ ≅ Z₂×Z₂`**. Its orbits (verified by enumeration n≤6):

| orbit size | when | count | *is* |
|---|---|---|---|
| **2** | `σt=t` (grid-sym) | `2^{D-1}` | a **BLUE LINE** `{t, φt}` |
| **4** | `t` not grid-sym | `(2^m−2^D)/4` | a **BLACK QUAD** `{t, σt, φt, σφt}` |

There are no other sizes (`φ`, `σφ` are fixed-point-free for `n≥4`). Total `|quarter| = (2^m + 2^D)/4 = 2^{m-2} + 2^{D-2}` (n=4..8: `3, 20, 272, 8320, 525312`).

**The structural meaning is clean:**
- A **black quad** `= {t, φt} ∪ {σt, σφt}` = a black line **together with its tournament-complement line** `σ(·)`. So at the quarter level, black lines come in **complement-pairs**.
- A **blue line** is `σ`-self (`σt=t` ⟹ the quad collapses to size 2): **blue lines are exactly the self-complementary lines.** This re-derives the tripartition at the *line* level: blue = self-complementary lines, black = complement-paired lines.

## The descent tower (a complement Cayley–Dickson)

`Full → Half → Quarter` is a descent, each step factoring out one order-2 structure and halving (with a grid-sym correction):

`2^m ⟶ 2^{m-1}+2^{D-1} ⟶ 2^{m-2}+2^{D-2}` (tilings ⟶ half ⟶ quarter),

folding by `σ` (tournament complement) then `φ` (tiling complement). This mirrors the project's Mode-B / Cayley–Dickson descent (`the-determinant-lens` R→C→H→O): each fold loses a complement symmetry. The blue subgraph = the quarter-tiling's own size-2 stratum = the SC/blue spine, and its connectivity + odd parity are the descent's invariants.

## Synthesis: reconstruct from structure, not fingerprints

The two results point the same way. Node fingerprints `(category, degree, H, d, cpS, …)` **degrade** — complete at n≤6, failing at n=7 — so pushing reconstruction by adding invariants is a losing race against cospectral families. But the **complement structure recurses cleanly**: `σ` folds to the half-tiling (SC/blue), `φ` folds to the quarter (blue lines = self-complementary lines), the blue spine is connected, and the odd grid-sym parity is n-robust. **Reconstruction of the metagraph should be organized by the descent tower** — assemble the quarter-tiling's fundamental units (self-complementary blue lines + complement-paired black quads) under the pairing rules and the recursive parity — rather than by ever-finer node labels. The quarter-tiling is the natural home for that assembly: it is the metagraph with *both* complement redundancies removed.

## Honest status & next

- **Computed:** `(I(Ω,x),d)` non-injective at n=7 (94% coverage, 142 collision groups); blue-spine connected n≤6; quarter-tiling = Klein-4 quotient with size-2 (blue, self-complementary) + size-4 (black, complement-paired) orbits, count `2^{m-2}+2^{D-2}`, verified n≤6.
- **Conjectural:** blue-spine connected for all n; the quarter-tiling descent gives the blue counts a THM-550-style recurrence; a genuinely-complete reconstruction from the quarter-tiling units + parity (the structural program).
- **Reframe:** "does H (or (H,d)) close reconstruction?" → complete only through n=6; at n=7 local fingerprints fail and one must recurse on the complement structure (half → quarter tiling).

— Related: `the-reconstruction-key-is-OCF-and-determinant-…` (S15, the n=6 (OCF,det) key), `does-H-close-reconstruction-…` (S14, metric suite), `buckets-and-pairs-…` (blue=half-tiling, #SC even), THM-549/550 (half-tiling), THM-281 (SC sizes odd), `the-determinant-lens-…` (d⊥H, Cayley–Dickson), HYP-3809/3810 (mac-mini, SC twin-pairing, κ_SC). Scripts: `04-computation/n7_ocf_det_injectivity_kps.py`, `quarter_tiling_model_kps.py` (+ .out). Not a HYP reservation.
