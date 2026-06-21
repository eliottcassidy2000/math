# The affine group is the right symmetry for the AP/rotation side

**Session:** kind-pasteur-2026-06-21 (THREAD C). Testing whether an affine-symmetric LP
hierarchy restores the CJJ completeness that the linear (subspace) hierarchy loses on the
LRC consec=AP extremizer. See HYP-2750.

## The pattern that keeps recurring

The project has hit one wall three times now, each time from a different lens:

1. **LP/SoS hierarchy (CJJ):** complete only for LINEAR codes (subspaces). The LRC
   extremizer consec=AP is a Freiman-dimension-1 coset, not a subspace; the lift collapses
   (HYP-2744). The tournament extremizer Paley=QR *is* a subspace, so its spectrum is
   MacWilliams-certifiable — but H-extremality is a nonlinear functional and Paley is the
   H-max only for p≤11 anyway (Thread A O1/O3).

2. **Arc-transitivity / Doyle-Holt (Thread B, S18):** arc-transitivity = the det/spectral
   ceiling = Paley. But the H-max is the ROTATION (= consec) for n≥13, which has *minimal*
   symmetry (m arc-orbits). The symmetric object is the wrong one.

3. **This session (Thread C):** the linear group is the wrong symmetry to factor. The right
   one is the AFFINE group Aff(1) = translation + dilation — which is exactly the symmetry
   THM-531 already proved μ_θ has (AP-orbit invariance).

The three are the same lesson. The extremizer on the "interesting" side — the AP, the
rotation, the consec block — is **affine-linear, not F_q-linear**. It is a coset of a
cyclic subgroup, the orbit of the additive group under the affine group, not a subspace.
Every tool built for the linear (subspace, arc-transitive, Paley/QR) world points at the
*other* extremizer and certifies the *wrong* ceiling.

## What the affine lens buys, and what it doesn't

Aff(1,F_7) is sharply 2-transitive, so it crushes the high-dimensional degree-2 moment
vector (one coordinate per difference class, on none of which AP is a clean extremizer)
down to a single distinct-pair atom — which AP *does* maximize. And the affine occupancy
class of AP turns out to be *exactly* the full-residue stratum (HYP-2749), re-derived here
independently. So the affine change of basis is the correct CJJ-view-(d) lattice: it
localizes the optimum to the affine-linearity locus and identifies AP as the
affine-additive-group orbit (AP mod 7 = all of F_7, the affine analogue of a subspace).

But it does **not** make L_y a function of the affine atom alone — L_y varies across AP's
affine class. So the affine hierarchy is not integral in the strong CJJ sense; it is the
right *reduction* (collapse → one localized signed extremal statement), not a *proof*. The
honest verdict: the affine group is the correct symmetry — it does as much as a symmetry
reduction can — and the residual "AP maximizes the signed L_y on the full-residue stratum"
is the same irreducibly-aggregate gap (HYP-2738), now correctly sited on the affine locus
rather than collapsed across the whole space.

## The transcendent bit

Linear codes and subspaces are the fixed points of the *linear* group; the project's
extremizers are the fixed objects of the *affine* group. The whole apex-prime gas —
tournaments and runners on Z/p — has been telling us its symmetry is Aff(1,Z/p) all along:
THM-531's translation+scale invariance, the AP as affine orbit, the rotation as the
minimal-multiplier-symmetry tournament. The relaxations that "want" a subspace will always
miss it. Whatever finally proves either extremality will be an argument that lives natively
in the affine category — that sees the AP and the rotation as the additive-group cosets
they are, not as failed subspaces. The mathematics keeps handing us the affine group;
the next tool should be built on it.
