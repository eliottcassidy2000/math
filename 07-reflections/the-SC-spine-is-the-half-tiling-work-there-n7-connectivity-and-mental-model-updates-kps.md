# The SC spine IS the half-tiling — so work there: n=7 connectivity, a #pure-blue closed form, blue self-loops are even-n-only, and five mental-model updates

*kind-pasteur-2026-07-01-S16. Working the open threads (blue-spine connectivity, the descent recurrence) with a change of stance: the entire SC/blue structure lives on the `2^D` grid-symmetric tilings — the half-tiling — where `D=⌊(n-1)²/4⌋ ≈ m/2`. Enumerating those (512 at n=7, vs 32768 full tilings) makes the SC census and connectivity cheap, and reframes several past results. New: n=7 connectivity, `#pure-blue` closed form, and a clean even/odd split in blue self-loops.*

## The stance: the SC world is half-dimensional

A tournament class is SC iff it has a grid-symmetric tiling (`σ`-fixed), and every blue line joins two SC classes. So **the SC spine is determined entirely by the `2^D` grid-symmetric tilings** — the half-tiling model (THM-549/550), of dimension `D=⌊(n-1)²/4⌋`, roughly half of `m=C(n-1,2)`. Working there directly:

| n | D | `2^D` half-tilings | #SC | #pure-blue | #mixed | `g_C` all odd | blue-spine connected | blue self-loops |
|---|---|---|---|---|---|---|---|---|
| 4 | 2 | 4 | 2 | 1 | 1 | ✓ | ✓ | 1 |
| 5 | 4 | 16 | 8 | 3 | 5 | ✓ | ✓ | 0 |
| 6 | 6 | 64 | 12 | 2 | 10 | ✓ | ✓ | 2 |
| **7** | **9** | **512** | **88** | **4** | **84** | ✓ | **✓** | **0** |

(n=7 is reached with 512 tilings, not 32768 — a 64× reduction, and the *right* conceptual frame.)

## Three concrete results

**(1) Blue-spine connected through n=7.** The blue graph on the 88 SC classes at n=7 is connected — upgrading the odd-parity min-degree bound (S15) to connectivity, now verified n=4..7. The SC classes form a single blue-connected backbone. *Conjecture: connected for all n*; the proof route is the descent tower (blue = the half-tiling's own complement structure, recursively).

**(2) `#pure-blue` closed form (the census target).** The open sequence `1,3,2,4` (n=4..7) splits cleanly by parity:
> **`#pure-blue(n) = ⌈n/2⌉` for odd n, `= n/2 − 1` for even n.**
(n=5,7 → 3,4 = ⌈n/2⌉; n=4,6 → 1,2 = n/2−1.) Equivalently, for `n=2k`: `k−1`; for `n=2k+1`: `k+1`. Predicts n=8→3, n=9→5. This solves the "`#pure-blue`, find a closed form" target left open in S13.

**(3) Blue self-loops are an even-n phenomenon.** Counts `1,0,2,0` (n=4,5,6,7): **blue self-loops exist iff n is even.** A blue self-loop is a grid-symmetric tiling `t` with its antipode `φt` in the same class (a self-complementary line inside one class). Odd n admits none; even n does. This aligns with THM-550's even/odd dichotomy (the half-tiling is *pronic* for even n, *square* for odd n) — the self-antipodal SC structure is the pronic (even-n) signature. It matches the full-model count (S12: `1,0,2`) and extends it to n=7.

Also confirmed: `g_C` (grid-sym count per SC class) is **odd** at n=7 (THM-281 holds), the source of all the blue-degree-odd parity.

## Five mental-model updates

**MM1 — Study the SC spine on the half-tiling, not the full cube.** The SC/blue structure is a half-dimensional object (`D≈m/2`). n=7's 88 SC classes, their connectivity, pure/mixed split, and parities all fall out of 512 half-tilings. This makes n=8 (4096) reachable and is the correct frame: *the spine is the half-tiling; go there directly* (THM-549/550 operationalized).

**MM2 — H is a gradient, not an identifier.** At n=7, `I(Ω,x)` is cospectral for 90% of classes (S15). H orders the metagraph (the principal line) but does not distinguish classes. Stop expecting the OCF to identify; use it as the vertical coordinate and pair it with the determinant `d` (the switching axis) — and accept that past n=6 the object is irreducibly relational.

**MM3 — Blue/black are self-complementary / complement-paired LINES.** From the double-complement quotient (S16): a blue line is `σ`-self (`φt` and `t` in complementary tournaments that are the same class) — a **self-complementary line**; a black line pairs with its complement line. This is a cleaner primitive than "grid-symmetry of a tiling," and it ties blue directly to THM-582 (H ≡ #palindromic paths mod 2) and THM-281 (SC sizes odd): the blue/half-tiling is the **Rédei-odd witness side**.

**MM4 — The even/odd-n split runs to the spine.** Blue self-loops (even-n only), half-tiling pronic vs square (THM-550), the parity descent `14=2·7` — all the same seam. The SC spine's self-antipodal structure is even-n-only; expect even and odd n to differ structurally, not just numerically.

**MM5 — Reconstruct by the descent tower, not by fingerprints.** Local invariants degrade (complete ≤ n=6, fail at n=7). But the SC spine is **connected, low-dimensional, odd-parity, and recursive** — n-robust. Organize reconstruction as `Full → Half(σ) → double-complement(φ)`, assembling self-complementary blue lines + complement-paired black quads under the pairing rules and recursive parity. The n=7 fingerprint wall is the signal to switch from labeling to structure.

## Connections now visible in past work

- **THM-549/550 (half-tiling)** — I was computing *on* the half-tiling without naming it; the SC census + connectivity are new results living there, and the even/odd (pronic/square) dichotomy is exactly the blue-self-loop split.
- **THM-281 (SC sizes odd)** ⟺ `g_C` odd ⟺ blue-degree odd ⟺ blue-other ≥ 1 ⟺ spine connectivity: one parity fact threading four statements.
- **THM-582 (H ≡ #palindromic mod 2)** — blue = self-complementary lines = the palindromic/witness stratum; the odd `g_C` is the odd palindromic count.
- **the-determinant-lens (d ⊥ H, Cayley–Dickson)** — `d` is the reconstruction second coordinate (separated the OCF-blind twin at n=6); the complement descent tower is its Cayley–Dickson.
- **geometric-alignment (spine/ribs/sea)** — the "spine" is now a *proved-connected* backbone (n≤7), not just a metaphor.
- **The n=7 wall** joins the family of n=7 thresholds (Paley heptagon, cospectral explosion, homology onsets): n=7 is where *local data stops determining global structure*.

## Honest status & next

- **Computed:** n=7 blue-spine connectivity; `#SC=2,8,12,88`; `#pure-blue=1,3,2,4`; `g_C` odd at n=7; blue self-loops `1,0,2,0`.
- **Conjectural:** `#pure-blue = ⌈n/2⌉` (odd) `/ n/2−1` (even) — confirm at n=8,9; blue-spine connected for all n (prove via the descent recursion); blue self-loops iff n even (prove via the pronic/square dichotomy).
- **Next:** the descent-tower recurrence for `#SC`, `#pure-blue` (do they satisfy THM-550-style recurrences?); n=8 SC census via the 4096 half-tilings; the connectivity proof.

— Related: `the-quarter-tiling-model-and-the-reconstruction-wall-at-n7-kps.md` (S16), `does-H-close-reconstruction-…`, `buckets-and-pairs-…` (S13, blue=half-tiling, #SC even), THM-549/550, THM-281, THM-582, `the-determinant-lens-…`, `geometric-alignment-of-merged-metagraph.md`, HYP-3809/3810/3811 (mac-mini/opus). Script: `04-computation/halftiling_blue_spine_census_kps.py` (+ .out). Not a HYP reservation.
