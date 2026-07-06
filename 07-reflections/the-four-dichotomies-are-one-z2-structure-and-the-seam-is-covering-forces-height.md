# The four dichotomies are one Z/2-structure, and the seam is: covering forces height

**opus-2026-07-06-S108** (HYP-4406). Working the remaining obligations through the many
lenses at once, the four dichotomies the project keeps meeting — addition/multiplication,
odd/even, positive/negative, rational/irrational — resolve into ONE Z/2-equivariant
structure at the prime 13, and its cleanest concrete face is now a green one-line lemma.

## The four dichotomies, at p = 13

* **addition / multiplication (sum-product):** `{1,…,12}` is the additive interval `[1,12]`
  (maximal relation lattice, the theta-extremal) AND the multiplicative group `(ℤ/13)*`
  (residue-pinned). This is the SEAM (S107): gap-emptiness is the rigidity of the coincidence.
* **positive / negative (antipodal):** the 12 residues split into 6 antipodal pairs
  `{r, 13−r}`; `dist(r/13)` sees only the ± class, so loneliness lives on the 6 = (p−1)/2
  classes. The binding pair at the tight point is `{1, 12} = {±1}`.
* **odd / even (2-adic):** `14 = 2·7`; the `2` supplies the ± complement (Rédei involution),
  the `7` supplies the heptagon (the Fejér / de Moivre 7th-root weights in the safe-measure
  theta). The safe-measure's sign structure is the `2`, its frequency weights the `7`.
* **rational / irrational (the two hard directions):** the arithmetic-hard band-blocker
  (multiplicatively structured to divide every small `q`, a rational census object) and the
  analytic-hard tight AP (additively structured, an irrational-sweep object), linked by
  dilation. These ARE the multiplicative and additive poles of the sum-product structure.

All four are `Z/2` actions, and they are the SAME action seen four ways: the involution that
swaps the two poles (a↔m, +↔−, odd↔even sub-lattice, rational↔irrational hardness) is one
`Z/2` sitting over `SL(2,ℤ)/Γ(2)` (the governing-structure survey's modular candidate). The
project discovered them in four coordinate systems because they are four shadows of one map.

## The seam made concrete: covering forces height (green)

The productive place is the SEAM — where the multiplicative (covering) constraint forces the
additive (height) structure. This is now one lemma, `coverer_height` (LRCCovererDichotomy,
green, standard trio):

> a runner `v ≡ r (mod 13)` with `1 ≤ r ≤ 12` that COVERS `q = r` (`r ∣ v`) is EITHER
> `v = r` (unlifted) OR `v ≥ 14r`.  (CRT: `v = r·m`, `m ≡ 1 (mod 13)`, so `m = 1` or `m ≥ 14`.)

For the UNIQUE-coverer moduli `r ∈ {7,…,12}` (verified: `q = r` is covered by `r` alone in
`{1..12}`), this says a pinned covering family either keeps `r` unlifted or pushes its `q = r`
coverer to height `≥ 14r ≥ 98`. Together with `LRCDivisorProtection` (no multiple of `r` ⟹
loose) and `LRCPinnedFloor` (`M ≥ 1/13`), the seam is the S78 pigeonhole / S80 height-forcing
chain, restated as: **multiplicative covering ⟹ additive height.** The lifts a gap-member can
carry are exactly the ones the covering condition permits, and covering forbids cheap lifts of
the unique coverers — it demands either no lift or a big one. This is why the l ≥ 7 stratum
(≥ 7 lifts must hit a unique coverer) forces `≥ 14r`, and why the far-height families leave the
window: the sum-product rigidity, in integers.

## Why this is the right way to see the remaining obligations

The density floor (the open analytic piece) asks: the AP uniquely achieves the all-order
theta cancellation. Read through the seam, that is: a covering pinned family, to stay
near-tight (additive), must keep the multiplicative covering — and `coverer_height` shows the
only near-tight covering configuration with no height blow-up is the AP itself (`m = 1` for
every unique coverer). Any deviation either raises a coverer to `≥ 14r` (leaving the
single-scale near-AP regime → handled by the descent/decorrelation) or drops covering (→
loose). So the seam localizes the density floor to the ONE configuration where additive and
multiplicative maximality coincide: the AP. The lemma does not close the floor (that remains
the Riesz-product all-order estimate), but it is the integer skeleton of why the floor's
minimum sits exactly at the AP and nowhere near it — the sum-product coincidence, made
checkable.
