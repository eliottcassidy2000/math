---
id: THM-649
title: THE CROSSING TRILOGY — (I) the book/tiling cube recovers Z(n) with Q anti-aligned to H; (II) the n=8 mirror break via σ-pairing quantization; (III) the winding-bit cylinder model (validated Z(n), n=5..8) with the PAIR-SPACE STRUCTURE THEOREM: the mod-2 crossing form on the winding cube is exactly degree ≤ 2, its quadratic part supported precisely on STAR (endpoint-sharing) pairs — vertex-disjoint pairs contribute affinely — with GF(2) pair-space rank = m+n′−2 rounded up to even (2,4,4,4,6 verified)
status: PROVED where stated ((III)'s disjoint-pair parity-flip lemma below; the star support and degree-≤2 exactness verified exhaustively on five winding cubes; rank formula = conjecture from 5 exact values). VALIDATED: cylinder min = Z(n) at n=5..8; book min = Z(n) at n=5..7. CONFIRMED: Q-vs-H anti-alignment in BOTH geometries; the n=8 mirror break (exhaustive over Fix(σ)). An in-session correction is part of the record: the naive affine parity law FAILED verification and the failure localized exactly to star pairs — the corrected statement is the theorem.
source: mac-mini-2026-07-07-S50/S51/S52 (HYP-5087, HYP-5097, HYP-5107; owner: fold the trilogy into canon + the pair-space dimension count)
depends_on: []
related:
  - HYP-5052   # opus-S142: the book cube = tiling cube; the affine crossing-parity law (book side)
  - THM-643/646/648  # the line-metagraph program (the same two-involutions algebra)
external: Guy/Harary–Hill Z(n); Zarankiewicz; Kleitman's parity theorem (the classical target of the parity structure).
---

# THM-649 — The crossing trilogy

## (I) The book side (S50)

On the tiling cube (chord pages = tile bits over the fixed spine), the 2-page crossing
minimum equals Z(n) = 1, 3, 9 at n = 5, 6, 7 (exhaustive), and Q is ANTI-aligned with
the Rédei count: corr(log H, minQ) = −0.75 / −0.76 / −0.59; the optimum-achieving
classes are cycle-rich (n=7: 130 classes, all H ≥ 19) while the transitive class is
crossing-poor. Crossing-minimal drawings live in Hamiltonian-path-wealthy tournaments.

## (II) The n=8 mirror break (S50)

σ-fixed (grid-symmetric) drawings satisfy min Q = 20 > Z(8) = 18 (exhaustive over
Fix(σ) = 2^12; 8 σ-fixed optima). Mechanism: 64 of 70 crossing pairs are σ-paired, so
on Fix(σ) they contribute in equal pairs and Q|Fix = 2·(quotient count) + (6 fixed
terms) — a coarser-quantized landscape that cannot resolve the free optimum. The
pairing-with-sign-flip principle converts the symmetry directly into the existence of
the gap; the census is confirmation. (Open: the closed-form quotient minimum = 20.)

## (III) The cylinder side (S51–S52)

**Model.** Vertices on two circles; each cross edge a radially-monotone curve with
winding lift d = (b + w) − a; two curves cross exactly #(ℤ ∩ open(x, x+δ)) times,
x = a₁−a₂, δ = d₁−d₂; within-circle chords are forced at C(m,4) + C(n′,4).
**Validation.** min over windings + twist at balanced splits = Z(n) exactly, n = 5..8
(1/3/9/18). **Transitivity pricing.** The tournament reading (two spines + cross bits)
has corr(Q_cyl, log H) = −0.64 / −0.66 (n = 6/7) over the whole winding cube; the
all-aligned (transitive) corner costs 2–3× the optimum. Both geometries price the
concentration of order; the book/cylinder difference is one spine versus two.

**The pair-space structure theorem (the dimension count).**
*Lemma (disjoint pairs flip).* If two cross edges share no endpoint, all four interval
endpoints x, x+δ are non-integers for generic twist, and shifting either winding by 1
shifts one endpoint of the count interval by exactly 1, changing the crossing count by
exactly ±1: mod 2 the pair contributes c₀ + w₁ + w₂ — affinely. ∎
*Star pairs are the quadratic support.* For endpoint-sharing pairs the interval
boundary is structurally integral (x = 0 or δ pinned), the flip argument fails, and
finite-differencing over the cube shows: the mod-2 form is EXACTLY degree ≤ 2
(reproduces the full cube at (m,n′) = (2,2),(2,3),(3,3),(2,4),(3,4)); its quadratic
part is supported on star pairs ONLY (0 non-star coefficients anywhere); its GF(2)
rank — **the pair-space dimension** — is 2, 4, 4, 4, 6, matching
**rank = m + n′ − 2 rounded up to even** in all five cases (even as it must be: the
form is alternating). Crossing-parity invariance over the winding cube ⟺ linear part
AND pair space both vanish — which the rank formula says never happens for m+n′ ≥ 4:
the cylinder, unlike the book at odd n (opus-S142's affine law), always retains a
rank-(m+n′−2)-ish quadratic core. The Kleitman-type parity questions for cylindrical
K_{m,n} live exactly in this star-supported form.

**Correction record (in-session).** The first derivation claimed the affine law
globally; exhaustive verification refuted it (failures at 18 shift-instances,
ALL endpoint-sharing), and the refutation's localization became the theorem. The
verify-before-assert loop is the method.

## Files
`lrc14_net_route_and_QvsH_macmini_S50.py`, `cylinder_tiling_transitivity_macmini_S51.py`,
`cylinder_parity_law_macmini_S52.py` (+ outs), `cylinder_pairspace_rank_macmini_S52.out`.
