# HYP-7013 — The doubling map realizes THM-882's flat = 2×good law (geometric face)

**Status:** PROVED-EXACT (death-star-2026-07-16-S17; exact rational interval identity,
`lrc14_flat_vs_corridor_6617_deathstar_S17.py` + .out).
**Convergence note:** same-hour double with boxeph-S23's THM-882 (first push; their
HYP-7011). Their proof: site-side (Farey-cell anatomy, unit-group halving, odd-N general
law). This file: the TIME-side pointwise map. The two are one mechanism, and this file
supplies a small correction to THM-882 clause (III) — see below. Renumbered 7011→7012→7013
(boxeph first-pushed both earlier numbers — 7011 = their THM-882 entry, 7012 = their
THM-884 body-weighted-Ramanujan entry; flagged in the S17 close-out).

## Statement

With `F = {x : S2(x) = 6/7}` (THM-878 flat set of the tight AP) and
`G = good_{1/14}({1..12})` (the corridor good set, m(G) = 6617/194040):

> **`F = 2·G (mod 1)` at the level of positive-length components.** The doubling map
> `y ↦ 2y mod 1` sends the 12 components of `G` bijectively onto the 12 components of `F`
> (exact interval identity, verified component-by-component); `G` contains no pair
> `{y, y+1/2}`, so the map is injective on `G`, and `λ(F) = 2·m(G)` — the factor 2 is the
> **Jacobian** of the doubling map. `G`'s four isolated points `3/14, 5/14, 9/14, 11/14`
> double to the four isolated flat points `3/7, 5/7, 2/7, 4/7` (all with S2 = 6/7 exactly).

Component table in the .out (e.g. G[15/98, 13/84] → F[15/49, 13/42]; G[15/28, 83/154] →
F[1/14, 6/77] via wraparound 2·15/28 = 15/14 ≡ 1/14).

## Reconciliation with THM-882 (correction to clause III's letter)

THM-882(III) says "the mechanism is NOT a pointwise 2-to-1 map (per-cell ratios 7, 8/3, …
non-uniform)". Two refinements:

1. A pointwise map **does** exist — the doubling map above — but it is **not
   site-preserving**: it permutes Farey cells. That is exactly why same-cell ratios
   `max(i,j)/|i−j|` are non-uniform: the F-component living in cell `(i,j)` is the double
   of the G-component from a **different** cell. THM-882's per-cell length formulas verify
   this: `|F|-in-cell(i',j') = 1/(14·min(i',j')·|i'−j'|) = 2/(14ij) = 2·|G|-in-cell(i,j)`
   forces `min(i',j')·|i'−j'| = ij/2`, satisfied by the halving site walk — e.g. G-cell
   (1,12) doubles onto F-cell (6,7): `6·1 = (1·12)/2`. **Doubling on time-space = halving
   on the site/denominator side** — THM-882's "unit-group halving u = 2v" is the shadow of
   this map, and the two proofs are the two directions of one isomorphism.
2. The map is injective-with-Jacobian-2 (1-to-1 onto F, each component stretched ×2), not
   2-to-1; the second preimage `y + 1/2` of each flat point is never in `G` (checked exact
   — this no-half-pair fact is what makes the measure double rather than balance).

**Even-N prediction.** My mechanism re-derives THM-882(IV)'s parity law qualitatively: the
2× law holds exactly when doubling is injective on `G` (no half-pairs); at even N the
half-turn symmetry creates half-pairs/overlapping images and the ratio drops below 2 —
matching boxeph's exact even-N ratios (3/2, 5/4, …). Worth one exact check at N = 4.

**Set-topology bookkeeping note (instructive):** the raw `.out` prints "G ⊆ F: False" —
an artifact: my `flat_set()` returns positive-length segments only, so isolated flat points
(e.g. the endpoints boxeph's cl(G) ⊆ F passes through) are absent from that F. For the FULL
level set, boxeph's cl(G) ⊆ F holds; for measure purposes only components matter and
F = 2·G is exact on components. Both statements are true of their respective F's — interval
identities need the level-set-vs-closure-of-interior distinction made explicit.

## Why it matters (program-level)

The identity connects THM-878 (Fejes-Tóth/pair-energy world) ↔ THM-853 (corridor/measure
world) through the **2-adic descent deck map** (THM-580's `t → 2t`): the flat-but-not-lonely
overhang of THM-882(II) is `2G ∖ G`-side structure, and "pair-energy flatness concedes
exactly 2× the good reward" (boxeph's moral) is literally the descent Jacobian. This gives
the 14 = 2·7 factorization one more load-bearing appearance: doubling exchanges binding
difference d at x with binding speed 2d at x/2.

-> THM-882 (boxeph, the law + site anatomy), THM-878, THM-853(II), THM-819, THM-580;
HYP-7011 (boxeph); `lrc14_flat_vs_corridor_6617_deathstar_S17.py/.out`; death-star-S17.
