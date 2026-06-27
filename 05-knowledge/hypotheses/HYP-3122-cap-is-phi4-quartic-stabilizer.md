---
id: HYP-3122
title: The covering-bound cap is a φ⁴ partition function — cap_k = C(k+1,2)/91 (quadratic/S2/b) − dip_k (quartic/S4/λ); the miss-count 4th cumulant κ₄ changes SIGN, going negative exactly at the hard binding row k=8 (the λ>0 Lee-Yang/φ⁴ stabilized regime); so coverage extremality reduces to a single sign-guaranteed quartic-cumulant bound
status: VERIFIED (cap split; κ₄ sign-change at k=8; #real=0 throughout) + PROOF-RELEVANT reframe + BOLD (the dip-bound = a uniform κ₄ bound; ear/κ₃ bridge). Not a proof.
source: mac-mini-2026-06-27-S67
reflections:
  - the-cap-is-a-phi4-field-theory-and-the-quartic-marks-the-hard-row
extends:
  - HYP-3103   # the miss-count PGF zeros (S66) — this adds the quartic/cumulant layer
  - HYP-3113   # codex two-map: quartic_cumulant_stabilizer (proposed; here computed)
  - HYP-3111   # codex Ising/Lee-Yang carriers (zeros on |z|=1 for toy ferromagnets; LRC is φ⁴ not Ising)
related:
  - HYP-3085   # gK8/S2 = the quadratic (b) term; the +6 S4 is the quartic
  - HYP-3092   # cap = pair-Pascal mass (the quadratic), dip = the higher-Pascal (the quartic)
  - THM-577    # the cap value + the dip
  - OPEN-Q-108
external: (φ⁴)₂ Euclidean QFT, exp(−λS⁴−bS²)dS; Simon–Griffiths / Lee–Yang / Lieb–Sokal
---

# HYP-3122 — The cap is a φ⁴ partition function; the quartic marks the hard row

## The φ⁴ split (VERIFIED)
```
cap_k = C(k+1,2)/91  −  dip_k          ( = quadratic b·S²  −  quartic λ·S⁴ , the exp(−λS⁴−bS²) measure )
```
- quadratic = the pair-normalized Pascal mass (`S2`), EXACT for k≥10.
- quartic = the dip, nonzero only at the sparse binding rows k=8,9 (the gK8 `+6 S4` term, which appears only
  at k=8; k=9,10 stop at `S3`, k≥11 at `S2`).

## The new signal (VERIFIED, `lrc_phi4_quartic_stabilizer_macmini_S67.py`)
The 4th cumulant `κ₄` of the miss-count `N` for `consec_k`:
```
 k:      8       9      10      11      12      13
 κ₄:  -0.79  +1.61  +3.92  +6.36  +8.19  +9.80
 dip: .0141  .00025    0       0       0       0
```
**κ₄ changes sign, going NEGATIVE exactly at k=8** — the unique row with the largest dip. `κ₄<0 ⟺ sub-Gaussian
⟺ a genuine λ>0 φ⁴ measure` (Simon–Griffiths/Lee–Yang). So **k=8 is the unique φ⁴-stabilized binding row**: the
quartic stabilizer engages precisely where the cap dips below the quadratic pair-Pascal. (`#real roots = 0`
for all k — Lee–Yang zero-confinement holds throughout; the quartic SIGN is the finer hard-row signal.)

## Proof-relevant reframe
Coverage extremality reduces (S63/S64) to **bounding the dip** — the only non-pairwise content. The φ⁴ reframe
makes it **a single 4th-cumulant inequality with a guaranteed sign**: the dip = the quartic `S4`; φ⁴/Lee–Yang
says the quartic is a *stabilizer* (`λ>0`, `κ₄<0` at the binding row), so the correction is bounded and the
right sign. Target: a uniform bound on `κ₄/κ₂²` over the binding family ⟹ the dip is bounded ⟹ the cap closes.
Complementary half (codex HYP-3108/3111): `corr(p0, nearest-zero)=+0.899`, `corr(p0,#real)=−0.48` — high
coverage ⟺ zeros off the real axis ⟺ the φ⁴-stabilized config. **Coverage extremality = φ⁴ stabilization =
Lee–Yang confinement.**

## Bold / creative
- **φ⁴ row := a binding row with κ₄<0** — the φ⁴ row is the hard row (k=8, verified).
- **Ear bridge:** the odd cumulant `κ₃ (3.7→5.6)` ↔ the **odd-ear / odd-cycle / OCF** structure (factor-critical
  ⟺ odd ear decomposition; the odd cycles that make the PGF complex-rooted and define `H=I(Ω,2)`); the even
  quartic `κ₄` ↔ the φ⁴ stabilizer. The cumulant odd/even split mirrors the ear odd/nested split.
- honest: the miss-PGF zeros are NOT on `|z|=1`, so the sector model is φ⁴, not plain Ising — the right
  single-spin measure that yields the zeros is the open modeling question.

## Next
1. Test the dip-bound = κ₄/κ₂² bound across the bounded binding bank (does a uniform quartic bound give the
   dip uniformly?). 2. Fit the effective `(λ_k,b_k)`; verify `λ_k≥0`. 3. Test the odd-ear(winding tournament)
   ↔ κ₃ bridge. 4. Derive the zero-confinement from a φ⁴ sector model (the Lee–Yang half).
