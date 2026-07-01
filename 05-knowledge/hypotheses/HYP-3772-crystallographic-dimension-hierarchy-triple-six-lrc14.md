---
id: HYP-3772
title: MORE crystallographic connections to LRC14 -- the generalized restriction dim>=phi(n) places LRC14 at phi=6, and the TRIPLE 6. Extending S65's (2,3,n) spine: an n-fold rotation needs an ambient lattice of dimension >= phi(n) (cyclotomic degree; Z[zeta_n] minimal; a planar n-fold quasicrystal is a cut-and-project from R^phi(n)). Hierarchy: phi=2 {3,4,6} = 2D crystals (5 Bravais lattices); phi=4 {5,8,10,12} = 4D quasicrystals (Penrose-5/octagonal-8/decagonal-10/dodecagonal-12); phi=6 {7,9,14,18} = 6D quasicrystals. LRC14 sits at phi=6: phi(7)=phi(14)=6, so the 7-/14-fold apex symmetry first lives in SIX dimensions (Z[zeta_14]). THE TRIPLE 6: 6 = phi(14) (cyclotomic degree) = the quasicrystal embedding dimension = |(Z/14)^*|=6, the LRC14 extremal LONELY SET (units {1,3,5,9,11,13}, HYP-3571/3615) -- so the lonely runners ARE the 6 independent directions of the 14-fold cyclotomic star, and the 3 antipodal unit-pairs (phi(14)/2) are the 'physical' half of the 6D=3+3 cut-and-project splitting. ALSO: 3D crystallography has 7 crystal systems and 14 Bravais lattices (=2*7=LRC14; hexagonal system = the order-6 covering-min side). Coxeter/Weyl: crystallographic <=> branch labels m in {2,3,4,6} (rank-2: A1xA1(2), A2(3, hexagonal), B2=C2(4, square), G2(6)); the apex is I2(7) (dihedral heptagon), the smallest NON-crystallographic dihedral after H2=I2(5) -- LRC covering-min = A2/G2 (crystallographic), apex = I2(7) (non-crys) = the 4cos^2(3pi/7) floor gap. 14-fold soft-matter quasicrystals exist, so the apex is genuinely quasicrystalline; the LRC covering may be a cut-and-project from the 6D cyclotomic lattice (a research direction).
status: SYNTHESIS / unifying map (extends S65 HYP-3771). Identities exact/verified: dim>=phi(n) hierarchy; phi(7)=phi(14)=6=|(Z/14)^*|=#lonely-units (the TRIPLE 6); 3D = 7 systems/14 Bravais; Coxeter labels {2,3,4,6}, apex=I2(7) non-crys. These are classical facts assembled + the triple-6 identification (the units-as-quasicrystal-dimension is the sharp new observation). The cut-and-project view of the LRC covering is a DIRECTION, not established. NOT a new proof step.
source: mac-mini-2026-06-30-S66
related:
  - HYP-3771   # my S65: the (2,3,n) angle-defect spine; 5 tilings/5 solids; apex-prime genus walk
  - HYP-3768   # my S64 + klein-S56: the order-6 hexagonal Dedekind margin; B2 square anomaly-free
  - HYP-3571   # the extremal lonely set = (Z/14)^* units = phi(14)/2 antipodal pairs (the '6'/'3')
  - HYP-3615   # the units {k/n: gcd=1} are the measure-0 lonely touch-points
  - HYP-3586   # X0(14) cusps / apex cusp d=7 / genus obstruction
---

# HYP-3772 -- the crystallographic dimension hierarchy and the triple 6 of LRC14

Extending S65 (HYP-3771, the `(2,3,n)` angle-defect spine): more crystallographic connections, and the sharpest
new one -- **the "6" of LRC14 is three things at once**.

## The generalized crystallographic restriction: `dim >= phi(n)`
An `n`-fold rotation preserving a lattice needs the ambient dimension to be at least `phi(n)` (Euler totient =
degree of the cyclotomic field `Q(zeta_n)`; `Z[zeta_n]` is the minimal such lattice). A **planar** `n`-fold
quasicrystal is a **cut-and-project** from `R^{phi(n)}`. This grades all `n`:

| `phi(n)` | `n` | realization |
|----------|-----|-------------|
| 2 | 3, 4, 6 | **2D crystals** (the 5 Bravais lattices, S65) |
| 4 | 5, 8, 10, 12 | **4D** quasicrystals (Penrose-5, octagonal-8, decagonal-10, dodecagonal-12) |
| 6 | **7, 9, 14, 18** | **6D** quasicrystals (7-, 14-fold -- exist in soft matter) |
| 10 | 11, 22 | 10D | 12 | 13, 21 | 12D |

The `phi(n) <= 2` row is exactly the crystallographic orders `{1,2,3,4,6}`; everything else is quasicrystalline.

## The triple 6 (the sharp new identity)
`LRC14` has apex prime 7, and `phi(7) = phi(14) = 6`. So the 7-/14-fold symmetry of the apex first lives in
**six dimensions** (`Z[zeta_14]`). And

> `|(Z/14)^*| = 6 = phi(14)`, with `(Z/14)^* = {1,3,5,9,11,13}` **= the LRC14 extremal lonely set**
> (the units / touch-points at the cusp `{1,...,13}`, HYP-3571/3615).

So the **`6`** of LRC14 is simultaneously (i) the cyclotomic degree `phi(14)`, (ii) the quasicrystal embedding
dimension, and (iii) the number of lonely units. **The lonely runners are the 6 independent directions of the
14-fold cyclotomic star.** The `phi(14)/2 = 3` antipodal unit-pairs are the "physical" half of the `6D = 3 + 3`
cut-and-project splitting (the same `3+3` as the icosahedral quasicrystal's physical+internal spaces). This
turns the repo's "lonely set = units" (measure-zero, counting-not-measure) into a **crystallographic-dimension**
statement: the lonely set *is* the cyclotomic star of the forbidden 14-fold symmetry.

## More resonances
- **3D crystallography: 7 crystal systems, 14 Bravais lattices** (`= 2*7 = LRC14`). The 7 systems span the
  crystallographic orders `{1,2,3,4,6}`; the **hexagonal** system (order 6) is the covering-min side. (2D: 4
  systems, 5 lattices -- the S65 five.)
- **Coxeter/Weyl:** a Coxeter group is crystallographic (a Weyl group, preserves a lattice) iff all branch
  labels `m in {2,3,4,6}`. Rank-2: `A1xA1 (2)`, `A2 (3, hexagonal, |W|=6)`, `B2=C2 (4, square, |W|=8)`,
  `G2 (6, |W|=12)`. The **apex is `I2(7)`** (the dihedral group of the heptagon) -- the smallest
  **non-crystallographic** dihedral after `H2 = I2(5)` (the pentagon). LRC covering-min = `A2/G2` (order 6,
  crystallographic); the apex = `I2(7)` (non-crystallographic) = the `4cos^2(3pi/7) > 0` floor gap.

## What it connects to / the direction
Every crystallographic guise says the same thing about the apex 7: no 7-fold crystal, `I2(7)` not a Weyl group,
`2cos(2pi/7)` irrational, minimal dimension `phi(7)=6` -- all `=` the LRC apex floor gap. The one genuinely new
*direction*: since 14-fold quasicrystals exist and are cut-and-projects from the 6D `Z[zeta_14]`, the **LRC14
covering may itself be a cut-and-project**: the 13 circle-speeds as a projection of a 6D cyclotomic lattice
covering, with the lonely set = the 6 star-directions. If made precise, the LRC14 floor would become a covering
statement about `Z[zeta_14]` in 6D -- the quasicrystalline face of the apex.

## Honest scope
Extends S65. The `dim >= phi(n)` hierarchy, the triple-6 (`phi(14) = |(Z/14)^*| = 6`), the 7-systems/14-Bravais,
and the `I2(7)`-non-crystallographic apex are exact/classical facts assembled; the **triple-6 identification**
(lonely set = quasicrystal dimension = cyclotomic star) is the sharp new observation. The cut-and-project view
of the LRC covering is a research direction, not established, and this is a synthesis, not a new proof step.
