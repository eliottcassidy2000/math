---
id: HYP-3132
title: The single remaining hard node (k=8 bounded-core dip) is a quintic whose gK8/Delsarte dual is the resolvent QUARTIC (t-1)(t-2)(t-4)(t-5); under the reflection s↦6-s it is a BIQUADRATIC u⁴-5u²+4 (a φ⁴ potential), De Moivre-solvable by radicals with perfect-square discriminant 9 (=why cap_8 is rational); so the open dip bound reduces to a solvable degree-2 (biquadratic-resolvent) bound
status: VERIFIED structure (resolvent quartic; biquadratic fold; perfect-square discriminant; rational collapse) + REDUCTION/IMPROVED-ARGUMENT (k=8 dip bound ⟸ degree-2 biquadratic bound via reflection). Not a proof.
source: mac-mini-2026-06-27-S70
extends:
  - HYP-3122   # the φ⁴ cap (S67) -- the φ⁴ potential IS this biquadratic resolvent
  - HYP-3085   # gK8 dual = the resolvent quartic
related:
  - THM-577    # the rational cap = the De Moivre solvability collapse to Q
  - HYP-3131   # far subsumes into bounded core (so k=8 is THE node)
  - HYP-3129   # multi-far SPEC bound proved elementarily (kps, c>=0.642)
  - THM-573    # apex-majority (proved), THM-576 (k>=10 proved)
  - OPEN-Q-108
external: De Moivre solvable quintics; resolvent quartic; (φ⁴)₂; Abel-Ruffini
---

# HYP-3132 — The k=8 hard node is a solvable biquadratic resolvent

## Comprehensive state: ONE node remains
Covering-bound tree (scout-verified): apex-majority THM-573 (PROVED), single-far/Node-3 THM-546/547 (PROVED),
multi-far SPEC bound HYP-3129 (PROVED ELEMENTARILY, `R′≥0.642`, not EH-dependent), far-pushes-out HYP-3131
(far subsumes into bounded core), bounded-core k≥10 THM-577 (PROVED). **The entire open core = the bounded-core
extremality at k=8** (the dip `dip_8=1081/76440` / the φ⁴ `κ₄/κ₂²` bound, HYP-3122).

## The k=8 node is a solvable quintic resolvent (the 120/320 / De Moivre lens)
`k=8 ⟹ |P|=5` — a QUINTIC (minimizer `{1,5,7,8,9}`). The gK8/Delsarte dual (THM-534) is its **resolvent
QUARTIC** `g(t)=(t−1)(t−2)(t−4)(t−5)=t⁴−12t³+49t²−78t+40`. (Dual degree stratifies the rows: quadratic k≥11,
cubic k=9,10, quartic k=8 — the deepest, the hard row.) Under the reflection `s↦6−s` (center `t=3`, the order-2
antipodal of S60), with `u=t−3`:
```
(t−1)(t−2)(t−4)(t−5) = (u²−4)(u²−1) = u⁴ − 5u² + 4      (a BIQUADRATIC = an EVEN/φ⁴ potential).
```
- **The φ⁴ structure of S67 IS the biquadratic resolvent** (φ⁴-evenness = the reflection symmetry =
  De Moivre-solvability). Discriminant `25−16=9=3²` is a **perfect square** ⟹ radicals collapse to `ℚ`
  (`t∈{1,2,4,5}`) — **why `cap_8`,`dip_8` are exact rationals** (THM-577), not surds.
- `1024=2¹⁰` (the owner-example resolvent constant = root product) = the tilings at `n=6` (`Q_{10}`); the
  `2^{1..4}` root tower = the dyadic H4 face.

## The improved argument
Bound the k=8 dip uniformly = a **solvable degree-2** problem, by the reflection fold:
1. the gK8 quartic obligation is `s↦6−s`-symmetric (S60: the pairwise covariance folds via `Z/2`) ⟹ it is a
   **biquadratic `u⁴+b u²`**, degree-2 in `u²` (always solvable);
2. the φ⁴ stabilizer sign (`λ>0`, `κ₄<0`, S67) fixes the direction (bounded, right sign);
3. so "bound `κ₄/κ₂²`" (open) ⟹ "bound the biquadratic resolvent's `u²`-coefficient" — a degree-2 bound, not a
   general quartic.

**Meta:** the bounded-core dual never exceeds degree 4 (k=8 the deepest), staying below the Abel–Ruffini
quintic wall — the same n≤7 tameness as the A000568 sandwich (S69). LRC(14)'s hard core is *solvable* because
its resolvent degree is ≤4.

## Next
1. Make the reflection-fold rigorous: the k=8 dip obligation = the biquadratic `u²`-coefficient bound over the
   binding family (the degree-2 reduction of HYP-3122's `κ₄/κ₂²`).
2. Bound that single coefficient (now solvable) + the φ⁴ sign ⟹ close the k=8 dip ⟹ close coverage extremality
   ⟹ (with the proved branches) close the covering bound.
