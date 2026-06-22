# The witness floor IS the p0 wide bound (two routes, one inequality)

**kind-pasteur-2026-06-22-S30.** A convergence reflection: two LRC(14) attack
routes that looked independent — the **1/7 global-witness density** `rho*_glob`
and the **p0 sector cover** `measS7` — are the *same inequality* seen from two
sides, joined by one elementary set inclusion and one exact duality.

## The two routes, as they were lived

The project ran two parallel sieges on the LRC(14) S3 residual:

- **The sector route** (codex, claude-opus, mac-mini): bound the cover atom
  `p0(E) = measS7(E) = meas{x : the cluster phases {frac(e x)} hit all 6 inner
  sectors}` by the plateau `cap_k`. Closed exhaustively (Leg-C 1.16M configs,
  gK8, single-far THM-563). The wide configs are the *dangerous* end — razor-thin
  margins (0.13).

- **The witness route** (kps, mac-mini): show the global-witness density
  `rho*_glob(P,E) = meas{x ∈ G_P : maxgap{frac(e x)} > 1/7}` is positive. This was
  the "comfortable" route — margins 5.95×, wide configs *safe* (decorrelate to
  `nu → 1`). THM-527 Part G called its uniform floor "the genuine remaining crux,"
  a compactness problem on a bounded-spread shape space.

They felt like different mathematics: one an *upper* bound on a covering measure,
the other a *lower* bound on a witness measure; one tight-at-wide, the other
safe-at-wide. The repo even carried them as separate Lean obligations.

## The collapse

Three elementary observations fuse them:

1. **`GOOD` is the complement of `dense`.** `maxgap > 1/7` ⟺ NOT(all gaps ≤ 1/7).
   So `meas(GOOD) = 1 − D`, `D(E) = meas{phases 1/7-dense}`. Bonferroni:
   `rho*_glob = meas(GOOD ∩ G_P) ≥ meas(GOOD) + meas(G_P) − 1 = meas(G_P) − D(E)`.

2. **`D ≤ p0` — dense implies covered.** If a length-1/7 inner sector held no
   phase, the gap straddling it would exceed 1/7 (the anchor phase `0 ∈ E` rules
   out the wrap-around edge case, and the half-open sector `[j/7,(j+1)/7)` catches
   the boundary phase). So `{1/7-dense} ⊆ {S7-cover}` *as sets*. The witness's
   "dense" event is a sub-event of the cover's "all sectors hit." Hence
   `rho*_glob ≥ meas(G_P) − p0(E)`.

3. **The duality `cap_k = min meas(G_P) = capRat(k)`.** Computed exactly: the
   small-part safe-measure floor equals the p0 plateau, for every `k = 8..13`.
   This is the hinge: it makes `p0(E) ≤ cap_k = min meas(G_P) ≤ meas(G_P)`, so the
   wide bound and the safe-set floor are *the same number*, and
   `rho*_glob ≥ meas(G_P) − p0(E) ≥ cap_k − max p0 = δ_k > 0` — the witness floor
   is exactly the wide-bound margin.

So `rho*_glob > 0` is a **corollary** of `p0 ≤ cap`. Two sieges, one wall.

## Why it transcends the particular bound

This is the recurring shape of the project (CLAUDE.md: "when two independent
frameworks converge on the same constraint — that's a reflection"). The witness
density and the cover atom are **Legendre/complementary views of one event
algebra** on the slow-time circle: the cover asks "is every sector occupied?",
the witness asks "is some gap empty?", and a gap is empty *because* a sector is —
de Morgan with a metric (the shared scale `1/7 = 2/14`, the apex prime 7). The
duality `cap = min meas(G_P)` is the statement that the *small part* and the
*cluster* are measured by the same plateau — the staircase's two legs (source
column `↔` scores `↔` `G_P`; sink row `↔` complement `↔` the cover) meeting at the
hypotenuse scale. The comfortable-margin "witness" repackaging didn't find a new
theorem; it found the *dual reading* of the cover bound, where the same `δ_k` that
is a thin ceiling for `p0` is a generous floor for `rho*_glob`.

The practical upshot is real: **there is one analytic obligation, not two.** Close
`p0 ≤ cap` (done, exhaustively) and the witness floor closes for free. The Lean
`hfloor` node is dischargeable from the `p0`-bound kernels already in the tree.

## Loose threads this opens

- If `D ≤ p0` is an inclusion of events, is there an x-by-x *quantitative* gap
  `p0 − D` with meaning? (It is `meas{all sectors hit but some gap = exactly 1/7
  boundary}` — the tight locus, where `M = 1/14` with equality.) The tight locus
  reappears.
- The duality `cap_k = min meas(G_P) = p0`-plateau deserves a *proof*, not just a
  k=8..13 check — it should be a clean statement about the level-`1/14` arcs of
  `{1,…,13}` and the 7-sector partition.

→ HYP-2832, THM-527, THM-530, OPEN-Q-108. Files:
`05-knowledge/results/lrc_witness_floor_bonferroni_elementary_kps.md`,
`04-computation/lean/.../LRCWitnessBonferroni.lean`.
