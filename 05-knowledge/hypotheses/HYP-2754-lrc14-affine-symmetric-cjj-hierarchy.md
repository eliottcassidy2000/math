---
id: HYP-2754
title: An AFFINE-symmetric (translation+dilation) CJJ-style LP hierarchy makes consec=AP the integral optimizer where the LINEAR (subspace) hierarchy collapses
status: SUPPORTED (partial), exact k=8,9,10 -- the affine fix gives the right REDUCTION (localize to full-residue stratum + AP unique extremizer there) but NOT full CJJ-integrality, NOT a proof (kind-pasteur-2026-06-21 THREAD C)
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

## Status: SUPPORTED (partial), exact verification k=8,9,10. NOT a full proof.

### (1) The precise obstruction (verified)
CJJ (arXiv:2211.01248) build the hierarchy for LINEAR binary codes A_2^Lin(n,d), exact at
level n via pseudoprobabilities; the Loyfer-Linial polytope is INTEGRAL at level n. The
structural property the optimizer needs (across the four views): it must be a LINEAR code =
F_q-subspace = closed under linear combination (view c); equivalently its higher-order
Krawtchouk/interaction moments are DETERMINED by the lower (views a,d: the Mobius/SoS lift
is self-tightening around a subspace, vacuous beyond level 1).
- Paley/QR HAS it: the QR cyclic code is an ideal in F_p[x]/(x^p-1), an F_p-subspace; its
  tournament spectrum is FLAT, |lambda|^2=(p+1)/4 (Gauss sum), MacWilliams-exact. PASS (L1).
- Consec/AP LACKS it: Freiman dimension 1, an additive COSET/translate, NOT a subgroup
  (verified: AP={0..k-1} is not additively closed in Z/14). FAILS (L1). In the LINEAR
  (difference-mod-14) degree-2 basis AP is a per-atom extremizer on NONE of the 7 atoms with
  a single sign (verified: 0/8/148/708/1408/1688/1460 shapes beat AP per difference class),
  so no single-atom certificate exists -> the lift collapses to the signed aggregate L_y
  (HYP-2738/2744). Script cjj_linearity_obstruction_kpswf6.py.

### (2) The affine fix -- WHAT IT BUYS (exact, k=8,9,10)
Aff(1,F_7)={x->ax+b} is SHARPLY 2-TRANSITIVE on F_7 (verified: 42 maps, 42 ordered distinct
pairs, bijective). So the degree-2 AFFINE-type quotient has only TWO atoms: distinct-residue
pair vs collision pair (mod 7). Exact breakpoint-engine results (affine_symmetric_lp +
affine_integrality_check kpswf6):
- AP MAXIMIZES the single affine atom (distinct-residue-pair count) at all of k=8,9,10. YES.
- L_y(AP) is the EXACT global max with ZERO beaters at k=8,9,10 (independent re-verification
  of THM-534 via the exact engine: 2633/7350, 26083/52920, 45253/79380).
- KEY: AP's affine OCCUPANCY class = EXACTLY the full-residue stratum. The affine signature
  class of AP has size 256/432/400 at k=8/9/10, identically equal to the count of
  all-7-residues-occupied shapes. AP is the stratum-MAX of L_y on it (0 beaters). This is an
  INDEPENDENT computational re-derivation of HYP-2749's stratum-localization, via the affine
  group rather than the resonant-arc bound.

### (3) Honest limit -- why this is SUPPORTED, not PROVED
L_y is NOT a function of the affine atom alone: it VARIES across AP's affine signature class
(0.10..0.358 at k=8). So the affine factoring does NOT make consec a trivially-integral
vertex ("only point of its type"). It does the correct REDUCTION: (a) localizes the optimum
to one affine occupancy class = the full-residue stratum (= HYP-2749's CJJ-linearity locus),
and (b) within it AP -- the affine orbit of the additive group F_7 (the affine analogue of a
subspace) -- is the unique L_y-maximizer. The residual "AP maximizes L_y on the full-residue
stratum" is the SAME signed THM-534 extremality, now correctly localized to the affine
locus, not collapsed. Full CJJ-integrality would require L_y constant on the affine optimal
face; it is not. So the affine hierarchy RESTORES the right reduction (collapse -> localized
single extremal statement) but does NOT by itself prove the AP extremality.

Tournament cross-check (sibling Thread A, threadA_cjj_paley_certify_kpswf6): the linearity
split is real at the CERTIFICATE level (Paley's QR spectrum has a MacWilliams/Delsarte dual;
AP's does not) but does NOT close H-extremality, because (O1) H=I(Omega,2) is nonlinear in
the spectrum and Paley is the H-max only for p=3mod4, p<=11 (at p=13 the NON-linear AP-half
wins, THM-128), and (O2) theta_LP bounds alpha(Omega) not H, and does not even separate
Paley at p=7. So linearity SPLITS the certificates but hands neither side a free proof.

## Scripts / outputs
- 04-computation/cjj_linearity_obstruction_kpswf6.py  (+ .out)
- 04-computation/affine_symmetric_lp_kpswf6.py        (+ .out)
- 04-computation/affine_integrality_check_kpswf6.py   (+ .out)
