---
id: HYP-+2888
title: LRC(14) strict-threshold realizability is a BOUNDARY / rational-witness structure -- the extremal APs cover EXACTLY to measure 1 (max additive energy A=1469), witness on the measure-ZERO boundary; positive-measure realizability is REFUTED at the strict threshold
status: GROUNDED (exact-fraction coverage). Complements kps HYP-2885 (additive-energy extremality) via EXACT coverage; redirects my own forced-overlap idea (REFUTED for extremal sets).
source: mac-mini-2026-06-22-S39 (user: creative realizability arguments; "missing structure slightly different")
related:
  - HYP-2885   # kps: LRC missing realizability = additive-energy extremality (interval is global p0-max)
  - HYP-2876   # my finite rational-witness certificate (every 13-set has witness D<=41)
  - HYP-2878   # strong-atom + covering-system route
  - HYP-2873   # additive energy A(E)=int|Ehat|^4, AP maximizes (mac-mini)
---

# HYP-+2888: LRC strict-threshold realizability is a BOUNDARY / rational-witness structure

## The exact-coverage computation (GROUNDED, lrc14_exact_coverage_macmini_S39.py)
LRC(14) <=> the safe set {t: ||s t|| >= 1/14 for all s in S} is NON-EMPTY <=> the OPEN unsafe arcs
U_s = {||s t||<1/14} (s arcs of width 1/(7s), total meas Sum=13/7~1.857) do not cover the CLOSED circle.
Computed meas(safe) EXACTLY (fractions) for 13-sets:
| set | meas(safe) | meas(union) | A(E) |
| consec {1..13}    | **0** | **1 (exact)** | **1469 (MAX)** |
| AP d=2, AP d=3    | **0** | **1 (exact)** | **1469** |
| Sidon-ish | 0.140 | 0.860 | 405 |
| random    | 0.090 | 0.910 | 685 |
| powers-of-2-ish | 0.278 | 0.722 | 325 |

## What this means (and what it REFUTES)
- The extremal sets are exactly the **APs** (max additive energy A=1469, affine-invariant), and they
  cover [0,1) to measure **EXACTLY 1**: the open arcs TILE, leaving only the measure-ZERO boundary safe.
- **REFUTED (discipline):** my forced-overlap / positive-measure realizability (meas(safe)>0 forced)
  FAILS -- the extremal APs have meas(safe)=0. At the STRICT threshold 1/14, there is NO positive-
  measure floor (the universal Farey floor 3/pi^2 lives at the RELAXED threshold 1/7). The strict-
  threshold witness is a measure-ZERO BOUNDARY point.
- For consec, the witness is **t = 1/14**: ||s/14|| >= 1/14 for all s in {1..13} because no runner is
  == 0 mod 14 (no runner = 14). The apex-7/D=14 structure: the tight witness is the boundary rational.

## The "slightly different" realizability (the user's hint)
Tournament realizability: combinatorial forcing (Omega=K_3 forces a C_5). LRC strict-threshold
realizability: a **BOUNDARY / rational-witness** structure -- the safe time is a bounded-denominator
RATIONAL a/D (D<=41, HYP-2876), found by the residue/covering structure (HYP-2878), NOT by measure.
This **BYPASSES extremality** (a point for any set), unlike kps's additive-energy route (HYP-2885,
which must prove the AP maximizes coverage). 
- **Two complementary routes to finish:** (a) kps -- AP maximizes coverage AND covers exactly 1 (not
  over), boundary safe [extremality]; (b) me -- every 13-set has a safe rational a/D, D bounded
  [rational-witness, extremality-free]. My S39 EXACT-coverage grounds the "covers exactly 1, not over"
  half of (a): the max-additive-energy AP tiles to precisely measure 1 with the boundary surviving.
- **The realizability obstruction to a counterexample:** to OVER-cover (kill even the boundary), a set
  must exceed the AP's coverage; but the AP already achieves the maximum (additive-energy extremality)
  at exactly 1, so no set over-covers => the boundary always survives => LRC(14). The open crux is the
  same extremality (a); my contribution pins (b) the boundary value is exactly 1 and the witness is
  the rational D=14, making the rational-witness route the extremality-free alternative.
-> kps HYP-2885 (the other route), HYP-2876 (bounded D), HYP-2878 (covering system).


## REFINEMENT (S39): the exact-coverage extremal is SCALING-invariant {1..13}*d, NOT additive energy
Verified (lrc14_exact_coverage): meas(safe)=0 (exact coverage) holds for **d*{1,...,13}** (consecutive
multiples) ONLY -- {1..13}, 2*{1..13}, 5*{1..13} all =0. But TRANSLATES/general APs have the SAME max
additive energy A=1469 yet POSITIVE safe measure: {2..14}->0.061, {3..15}->0.098, AP(a=2,d=3)->0.118.
- **Coverage != additive energy:** A is TRANSLATION-invariant (all length-13 APs: A=1469); strict-1/14
  coverage is SCALING-invariant (d*{1..13} only). They disagree -- {2..14} is max-A but not extremal.
  (Consistent with my S29 'L_y not a function of A'; sharpens kps HYP-2885 from a trend to the exact
  scaling-invariant extremal.) The right extremal invariant is scaling-, not translation-, invariant.

## THE FINISHING STRUCTURE (clean reduction)
LRC(14) <= the conjunction of:
- **(2) VERIFIED, RIGOROUS:** every consecutive-multiple set {d,2d,...,13d} has the explicit safe
  witness **t = 1/(14d)**: ||jd/(14d)|| = ||j/14|| = min(j,14-j)/14 >= 1/14 for all j=1..13. So the
  tight sets are explicitly safe (no measure needed).
- **(1) THE OPEN CRUX (strongly evidenced):** every 13-set that is NOT d*{1,...,13} has meas(safe)>0.
  Random+perturbation search: 0/5668 non-(consec-multiple) sets reached meas(union)=1 (max 0.998 for
  near-consec [1..5,7..13,69]). So d*{1..13} is conjecturally the UNIQUE exact-coverer.
- Combine: the unique tight sets are the consecutive-multiples, resolved explicitly by (2); all other
  sets have a positive-measure safe set by (1). => safe set always non-empty => LRC(14).
- **This is the 'slightly different' realizability:** the obstruction is a SCALING-invariant exact-
  tiling extremality (only {1..13}*d tiles [0,1) by its U_s arcs), with the tight case explicitly
  witnessed -- complementing kps's additive-energy route and bypassing the translation-invariant A.
  The crux (1) is now CLEANLY stated: characterize the exact tilers of [0,1) by arc-systems U_s.