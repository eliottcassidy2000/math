---
id: THM-888
title: CLOSED-FORM COMB CERTIFICATES FOR Q_s — the csc² multiplication identity Σ_{k mod q} csc²(π(x+k/q)) = q²csc²(πqx) telescopes the resonance law's k̂-masses EXACTLY: (A) THE COMB DIAGONAL IS A FUNCTION OF w MOD 7 with full closed form — for coprime w, owner e's comb contributes 2π²[Σ_{r=1}^{6}|ν̂_e(rw)|²csc²(πrw̄/7)/(2P²)·(P/7e)²-collapsed + |ν̂_e(0)|²((P/7e)²−1)/(6P²)], no sums, O(1) per owner (referee: exact); (B) the off-comb tail collapses per distance-slice to EXACT masses csc²(π‖dw̄/e‖)/(2e²) (referee: exact), giving a SCAN-FREE O(Σ_e e) certificate Q_s(w) ≤ (Σ_e√(D_comb+D_tail))² — sound on the full battery; (C) HONEST VERDICT: off-resonance this closed form is NOT sharper than THM-887's crude ℓ-sum (106–235×M vs 39–89×M) — the loss is structural: the envelope² per slice carries an M_e²-floor (Σ_d csc²(π‖dw̄/e‖) is w-independent in total), while the true slice-interior mean square is ~M_e; THE NAMED SHARPENING is the slice-interior Parseval (replace envelope² by the slice's true ℓ²-mass) — where the real factor ~M_e/6 lives
status: (MI)/(MI0), the slice collapse, and the comb closed form all PROVED (partial-fraction expansion of csc²; three lines) and machine-refereed EXACT; the assembled certificate sound (covers on every battery point); the sharpness claim honestly NEGATIVE off-resonance with the loss mechanism identified and the sharpening named
source: boxeph-2026-07-16-S27 (owner directive: take the sharp off-resonance constant)
depends_on: [THM-887 (the profiles being telescoped), THM-886(III) (the k̂ csc² bound — exactly the form (MI) collapses)]
related: [klein cont.6's CRT audit factorization — the same divide-the-period spirit; their per-prime independence is the natural tool for THIS file's cross-term/compound-hit refinement]
script: 04-computation/lrc14_sharp_offresonance_boxeph_S27.py -> 05-knowledge/results/lrc14_sharp_offresonance_boxeph_S27.out
---

# THM-888 — closed-form comb certificates; the honest sharpness verdict

**(MI) The multiplication identity.** Σ_{k=0}^{q−1} csc²(π(x + k/q)) = q²csc²(πqx),
and Σ_{k=1}^{q−1} csc²(πk/q) = (q²−1)/3. *Proof:* csc²(πx) = π⁻²Σ_{n∈Z}(x+n)⁻²;
summing over k regroups the double sum over (1/q)Z as q²·the same expansion at qx. ∎
This is EXACTLY the shape of THM-886(III)'s k̂-bound (1/(2P²))csc²(πn/P) — so every
AP of frequency indices carries an exactly-collapsible k̂-mass.

**(A) The comb diagonal in closed form.** For gcd(w,P) = 1, the weight indices n whose
image nw lies on owner e's comb form eZ_P; the comb point em carries tooth |ν̂_e(m mod 7)|.
Slicing m by residue mod 7 and applying (MI) with q = P/(7e):
> the r-class total k̂-mass = (P/(7e))²·csc²(π·‖rw̄/7‖)/(2P²), r = 1..6; the r = 0 class
> totals ((P/(7e))²−1)/(6P²).
So D_comb^{(e)}(w) is a FUNCTION OF w MOD 7 ONLY, in closed form, O(1) per owner.
(Referee: exact to 1e-9 against direct summation, both owners, two w.) In particular the
comb-diagonal's w-dependence is a permutation of six tooth-weight pairings — the finite
"resonance types" of a cluster are the 6 unit residues mod 7.

**(B) The tail slices.** The off-comb indices at integer distance d from eZ split into
two signed APs whose total k̂-mass is EXACTLY csc²(π‖dw̄/e‖)/(2e²) each (MI, q = P/e;
referee: exact). With THM-887(II)'s envelope, the whole certificate
> Q_s(w) ≤ (Σ_e √(D_comb^{(e)} + D_tail^{(e)}))²,
> D_tail^{(e)} = (π²/e²)·Σ_{d≤e/2} c_d·min(M_e², R_e²csc²(πd/e))·csc²(π‖dw̄/e‖)
is SCAN-FREE: O(Σ_e e) arithmetic, no Z_P pass, no spectral evaluation. Sound on the
full five-cluster × w battery.

**(C) The honest sharpness verdict.** Off-resonance this bound lands at 106–235×M —
WORSE than THM-887's crude ℓ-truncated law (39–89×M). The mechanism: within a distance-d
slice the certificate pays envelope², but Σ_d csc²(π‖dw̄/e‖) is w-INDEPENDENT in total
((MI0)), so the M_e²-capped slices contribute ≥ ~M_e²/12·(2π²/M)-scale REGARDLESS of w —
an envelope floor. The true slice-interior ℓ²-mass is ~M_e (Parseval), a factor ~M_e/6
below the envelope². **The named sharpening: slice-interior Parseval** — replace
envelope² by the slice's true mean square (computable per owner from the run structure,
or boundable by |A_c|-autocorrelations). That is where the genuinely sharp constant
lives; this file contributes the exact collapse machinery it will ride on.

## What survives as value now

1. **(A) is a little theorem with content:** the comb diagonal of Q_s depends on w only
   through w mod 7, with explicit csc²(πr/7)/(98e²)-type weights. Resonance-type
   classification of dilations is FINITE (6 classes) at the comb level.
2. **Scan-free certificates:** O(Σe) a-priori Q_s bounds (engineering-grade: the comb
   bookkeeping (M_e, N^{(e)}, R_e) + closed forms = a certificate library candidate).
3. **The exact slice masses** are the input the CRT compound-hit refinement (klein
   cont.6's factorization spirit) and the slice-Parseval sharpening both need.

## Evidence log
- [x] (MI)/(MI0) refereed; slice collapse exact; comb closed form exact; battery covered
- [ ] slice-interior Parseval (the real sharp constant)
- [ ] CRT compound-hit refinement of the cross terms
