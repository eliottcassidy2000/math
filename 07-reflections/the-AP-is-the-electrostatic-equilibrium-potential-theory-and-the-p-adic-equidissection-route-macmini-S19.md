# The AP is the electrostatic equilibrium — potential theory joins the equi-family, and the p-adic equidissection route to the floor

*mac-mini-2026-07-06-S19 (HYP-4472). Owner: keep seeing how equicontinuity relates
to the equi- family; keep investigating, integrating, reframing. This note adds two
unexplored classical layers (0 prior repo mentions each): **potential theory**
(Fekete / equilibrium / Stolarsky) — a verified new face that *unifies* the metric
equi-notions — and **equidissection** (Monsky / p-adic valuation) — a candidate
*proof route* to the equidecomposability floor the fleet just isolated. Verified:
`lrc_fekete_equilibrium_energy_v2…out`. Integrates opus-S111 (combinatorial/metric
split), kps-S27 (dodge/patch), my S18 (equicontinuity = regularity axis).*

## Where the fleet arrived (this same day)

- **opus-S111:** the equi- prefix splits **combinatorial** (equinumerosity,
  equidecomposability, equipartition — structure, universal, n-blind, *necessary*)
  vs **metric** (equicontinuity=Arzelà–Ascoli, equidistribution=Weyl,
  equioscillation=Chebyshev — quantity, n-specific, *decisive*). (G) lives on the
  metric side. `LRCLeaveOneOut.lean` GREEN.
- **kps-S27:** the floor is an **equidecomposability residual** — *dodge ≠ patch*:
  equidistribution/equinumerosity give blind uniform coverage (a far runner "rides
  along", growing the safe set — verified: safe *doesn't degrade* with height), but
  only the AP's roots-of-unity fiber is a **rigid scissors class** that can *patch*
  a concentrated lonely set. Unbounded case = equidistribution (resolved by my
  decorrelation); bounded case = the equidecomposability floor.

Both integrate my S18 (equicontinuity's failure = the floor, dual to
equidistribution). The picture is complete on the metric axis. This note adds the
two layers underneath it.

## New face: the AP is the electrostatic equilibrium (potential theory)

Include the **observer** at `0`. Then the 13 points `{0} ∪ {v_i·(1/13)}` for the AP
are exactly the **13th roots of unity** — and (verified) they are the discrete
**Fekete / equilibrium** configuration:

> **min logarithmic energy** `−Σ_{i<j} log|z_i−z_j|` (AP `−16.672`; all 160
> perturbations strictly higher) **= max sum of chordal distances** (Stolarsky) **=
> min star-discrepancy** (opus-S65) — one config, the AP.

By the **Stolarsky invariance principle**, `L²`-discrepancy + distance-sum is
constant, so these are literally the same extremal. So potential theory *unifies*
the metric equi-family: **equidistribution (uniform equilibrium measure) =
equilibrium (Fekete points) = min energy = min discrepancy = the AP**, and the LRC
tight value `1/n = 1/13` is the **equilibrium spacing** — the observer's gap when
observer + runners sit at electrostatic equilibrium. This is the cleanest statement
yet of *why* the AP is tight: it is the ground state of `n` unit charges on the
circle with one pinned at the observer.

**The floor in this frame — an energy/gap trade-off.** `M(S)` is the observer's
*maximal openable gap*. At equilibrium (AP) it is `1/13`. Opening a larger gap forces
runners away from `0`, i.e. *off* equilibrium, which **strictly raises the energy**.
Covering (`safe=0`) bounds the energy from above (the arcs must tile). The floor is
the assertion that no covering configuration can open the observer's gap into
`(1/13, 2/25)`: the energy cost of a gap of width `g ∈ (1/13, 2/25)` exceeds the
covering budget, except at the exact equilibrium `g = 1/13`. This is the
potential-theoretic form of the equidecomposability rigidity (kps) — and, like it,
it is n-specific (the energy-vs-budget balance depends on `n`).

## New route: the p-adic equidissection obstruction (Monsky)

The floor is a **tiling-impossibility**: only the AP's arc-lattice equi-covers
(patches) the circle at `2/25`. The canonical theorem of exactly this shape is
**Monsky's theorem** — *a square cannot be equidissected into an odd number of
equal-area triangles* — proved not by geometry but by a **`2`-adic valuation**
(extended to `ℝ`) plus **Sperner's lemma**: color points by the valuations of their
coordinates, and a combinatorial fixed point forces a triangle whose area has odd
2-adic valuation, contradicting an equidissection. The moral: **equidissection
(equidecomposition into equal pieces) is obstructed by a p-adic valuation.**

This is a genuinely new candidate route to the LRC floor, and the structure lines
up:

- the floor is an *equal-piece covering* rigidity (kps's patch) — a Monsky-type
  statement;
- the witness lives at `q ∣ v_i±v_j` (my HYP-4432) with covering primes
  `5,7,8,9,11` forced (kps) — a rich **p-adic** profile of the residues;
- the **n-specificity** carries a p-adic scent: `2n−1 = 13` (prime) at `n=7` where
  the gap is nonempty, but `2n−1 = 25 = 5²` at `n=13` where it is empty — the target
  denominator `2/25 = 2/5²` has a **prime-square** (higher `5`-adic valuation)
  structure absent at n=7,8 (`13`, `3·5`).

**The programme this suggests:** assign a `p`-adic valuation (for `p = 5`, or `2`
via the `14=2·7` parity) to the runner residues at the witness; a Sperner/parity
count over the clearance band should show that an *exact* covering (safe = 0) at
`c/q ∈ (1/13, 2/25)` forces a valuation contradiction — that only the AP's fiber
(valuation-balanced roots of unity) survives. This would give a **combinatorial-
topological** proof of the metric floor, exactly as Monsky gives a p-adic proof of a
measure equidissection. It also explains why the structural (n-blind) lenses fail:
the obstruction is a valuation that *sees* `2n−1`'s factorization, hence is
n-specific by construction.

## Almost periodicity (the equicontinuity closure)

One more equi-adjacent fact tying it together: an **equicontinuous** flow is
**almost periodic** (Bohr). The runner maps `t ↦ v_i t` are individually periodic;
their min `f_S(t)` is periodic. The AP's `f_AP` is the *most* symmetric almost-
periodic profile — equioscillating at the `φ(n)` units. `M`'s non-equicontinuity
(kps-S26) is the statement that the *family over directions* is **not** uniformly
almost periodic — the almost-periods grow with height. So "almost periodicity holds
per-config but not uniformly over directions" is the same `L ~ height` obstruction in
Bohr's language.

## Net — the equi-atlas, three floors deep

| layer | the AP is… | role |
|---|---|---|
| **combinatorial** (equinum/equidecomp/equipartition) | max relations / rigid scissors class | structure, universal, *necessary* |
| **metric** (equicont/equidistrib/equiosc) | Chebyshev extremal at φ(n) units; non-equicont resolution | quantity, n-specific, *decisive* |
| **potential-theoretic** (equilibrium/Fekete/Stolarsky) | electrostatic ground state; energy = discrepancy = distance-sum | *unifies* the metric layer; tight value = equilibrium spacing |

- **New verified face:** the AP (with observer) is the 13-point Fekete/equilibrium
  config; potential theory unifies equidistribution = energy = discrepancy = distance-sum.
- **New proof route:** the floor is a Monsky-type equidissection rigidity; a p-adic
  valuation (`p=5` from `2n−1=5²`, or `p=2` from `14=2·7`) + Sperner is the candidate
  obstruction — and its `n`-specificity is built in (it sees `2n−1`'s factorization).
- **Closure:** `M`'s non-equicontinuity = non-uniform almost-periodicity = the floor.

## Pointers

- `lrc_fekete_equilibrium_energy_v2_macmini_S19.py/.out` (AP = min log-energy / max
  chord-sum, observer included).
- opus-S111 (combinatorial/metric split), kps-S27 (dodge/patch), mac-mini S18
  (equicontinuity axis), HYP-4432 (witness denominator), HYP-4452 (leave-one-out /
  n-specific), opus HYP-4074 (discrepancy inversion). New leads: Stolarsky invariance,
  Fekete points, Monsky's theorem / p-adic equidissection, Bohr almost periodicity.
