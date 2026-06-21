---
id: HYP-2750
title: An AFFINE-symmetric (translation+dilation) CJJ-style LP hierarchy makes consec=AP the integral optimizer where the LINEAR (subspace) hierarchy collapses
status: OPEN target / scout in progress (kind-pasteur-2026-06-21 THREAD C)
source: kind-pasteur-2026-06-21-S?? (THREAD C: linear-vs-nonlinear obstruction + affine fix)
depends_on:
  - THM-531   # AP-orbit invariance: mu_theta is translation- AND scale-invariant (the affine group is the right symmetry)
  - THM-534   # L_y moment-LP dual; consec maximizes the SIGNED L_y (the open scalar extremality)
  - HYP-2744  # Boolean-Mobius hierarchy (linear/subspace lift COLLAPSES for the extremality)
  - HYP-2738  # consec-max is irreducibly aggregate (signed L_y, no nonneg monotone cert)
related:
  - HYP-2726  # cover bound IS a Delsarte LP; consec saturates it (LP-tight)
  - HYP-2749  # stratum-localization: consec-max reduces to the full-Z/7-residue stratum
  - THM-126   # Paley H-maximizer (tournament side, the LINEAR/QR optimizer)
external: Coregliano-Jeronimo-Jones complete LP hierarchy (arXiv:2211.01248, 2112.09221)
---

# HYP-2750 — Affine-symmetric CJJ hierarchy for the AP extremizer

## Claim (the THREAD C target)

CJJ completeness (level O(n^2) retrieves the optimum; the SoS/Mobius lift is exact) is
proved for codes whose optimizer is a LINEAR code (closed under linear combination =
F_q-subspace). On the tournament side the extremizer Paley = QR is a genuine linear cyclic
code, so the hierarchy CAN in principle certify it. On the LRC side the extremizer
consec = AP is NON-linear (Freiman dimension 1, a coset/translate, not a subgroup), which
is the structural reason the Boolean-Mobius / subspace lift COLLAPSES (HYP-2744): consec is
not the linear optimizer of the lifted LP.

**Proposed fix.** Replace the linear group (subspace symmetry, S_n permutation) by the
AFFINE group Aff(1, Z/p) = {x -> a x + b : a in (Z/p)^*, b in Z/p} = translation + dilation.
By THM-531 the LRC functional mu_theta (hence L_y, S7, the whole cover object) is invariant
under exactly this affine group: translation-invariance (b) AND scale/dilation-invariance
(a). An AP is precisely an AFFINE orbit of {0,1,...,k-1} (an affine-linear object: closed
under affine combination / a coset of a cyclic subgroup). The conjecture: a CJJ-style
hierarchy whose moment variables are indexed by AFFINE-invariant types (not linear/subspace
types) has consec=AP as its INTEGRAL (linear-in-the-affine-sense) optimizer, restoring
completeness for the AP extremizer.

## Honest status at reservation
STUB. The precise obstruction (which CJJ view; which closure property AP lacks and QR has)
and the level-2 affine-symmetric LP test are in progress this session. Scripts:
04-computation/cjj_linearity_obstruction_kpswf6.py,
04-computation/affine_symmetric_lp_kpswf6.py. Outputs in 05-knowledge/results/.
