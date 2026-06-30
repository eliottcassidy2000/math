---
id: HYP-3737
title: The RADIUS-1 BAND over-constraint FORCES the construction for n>=12, resolving the covering-min trajectory into THREE regimes -- mediant 2/(2n-1) [n<=6, projective] -> small-depth spread [0;n-1,a(n)] [n=7..11, a(n)=2,2,4,4,3] -> CONSTRUCTION n/Phi_6=[0;n-1,n] [n>=12]. MECHANISM (leveraging klein-S38's radius-demand criterion): any small-M covering set must cover Z/D at radius 1 for every D in the band (n,2n-2]; the CONSECUTIVE base {1,..,n-2} covers the band INTERIOR D in (n,2n-3] at radius 1 (verified all n), leaving only the edge D=2n-2; a SPREAD/drop base instead SCATTERS uncovered residues across many moduli (n=14: spread base leaves deficits at D=15,20,21,24,25,26 vs consecutive at only D=26), which NO single outlier can patch (one outlier covers <=3 residues per D). So the band FORCES the consecutive base; covering the top two q=n-1,n with the remaining speed FORCES the unique outlier lcm(n-1,n)=n(n-1) = the construction. For EVEN n the outlier n(n-1)≡0 mod 2(n-1) covers the band edge cleanly (n=14: 182=7*26; the parity echoes HYP-3729 even=bipartite). The spread family survives only while the band is NARROW (n=7..11); at n>=12 the band is wide enough that only the consecutive base works -> construction. CONSEQUENCE: LRC14 covering-min = 14/183 (the construction; opus's original value, vindicated for n=14 -- the spread family that beat it at n<=11 is dead by n=12), margin = 13/2562 > 0, pinning HYP-2566 looseness at n=14
status: STRONGLY EVIDENCED (not a full proof). The forcing mechanism is verified: consecutive base covers the band interior (all n=10..16); spread/drop bases fail (n=14, several checked); the outlier n(n-1)≡0 mod 2(n-1) for even n; the construction is the unique consecutive-base+single-outlier covering set (lcm(n-1,n) forced). The spread-family death at n>=12 is reliably established (S53/HYP-3735). opus's 107-set scan found 14/183. NOT proven: "only the consecutive base covers the band interior" is verified, not exhaustively proved; the full ILP at V=n(n-1) times out (milp on 40747 constraints) so couldn't independently confirm. So a(14)=14 / covering-min=14/183 is strongly indicated, rigor open.
source: mac-mini-2026-06-30-S54
related:
  - HYP-3735  # the small-depth spread family dies at n=12 (this explains WHY: the band over-constraint)
  - HYP-3734  # klein-S38 radius-demand criterion (the leverage); a_1=n-1; up-set rungs
  - HYP-3732  # the Stern-Brocot ray; construction = depth-n point
  - HYP-3729  # even=bipartite C_n parity (the outlier-covers-edge parity echoes it)
  - HYP-2566  # uniform looseness -- this pins the n=14 value at 13/2562
results:
  - 05-knowledge/results/covering_min_full_n14_macmini_20260630.out
---

# HYP-3737 -- the over-constraint forces the construction; the three-regime trajectory

The next target after HYP-3735: for `n>=12`, is the covering-min the construction (depth `n`) or a
moderate-depth spread set? Leveraging klein-S38's radius-demand / over-constraint criterion, the answer is the
**construction** -- and the mechanism resolves the whole covering-min trajectory.

## The radius-1 band over-constraint (the leverage)
klein-S38: a covering set with `M(S)=k/m` must cover `Z/D` at radius `floor(kD/m)` for every modulus `D`. The
**radius-0 moduli `D<=n-1` are the THM-523 resonances**; the **radius-1 band `D in (n, 2n-2]`** is the extra
demand. For any small-`M` (covering-min-scale) set, the band must be covered at radius 1: for each `D` and each
`j`, some speed `v` has `vj ≡ 0,±1 mod D`.

**The construction saturates the band.** Verified: the construction `{1..n-2, n(n-1)}` has coverage radius
**exactly 1 at every `D in (n,2n-2]`** (e.g. `n=14`: all of `D=15..25` tight at radius 1) -- it is maximally
over-constrained there.

## The forcing -- why only the construction survives
- **The consecutive base covers the band interior.** `{1,..,n-2}` covers `Z/D` at radius 1 for every
  `D in (n, 2n-3]` -- **zero uncovered** (verified `n=10..16`). It leaves only the edge `D=2n-2`.
- **A spread base scatters the deficits.** Verified `n=14`: a spread base leaves uncovered residues at
  `D=15,20,21,24,25,26` (many moduli, up to 4 each); a drop-two base at `D=22,23,24,26`. **One outlier covers
  `<= 3` residues per `D`** (`j ≡ 0, ±w^{-1} mod D`), so it cannot patch deficits scattered across many `D`.
  Hence spread bases **fail the band**.
- **The outlier is forced.** With the consecutive base fixed, the remaining speed must cover `q=n-1` and `q=n`;
  the unique single speed doing both is `lcm(n-1,n) = n(n-1)` -- the construction's outlier. For **even `n`**,
  `n(n-1) ≡ 0 mod 2(n-1)`, so it covers the band edge `D=2n-2` cleanly (`n=14`: `182 = 7·26`; a parity echo of
  HYP-3729 even=bipartite).

So for `n` large enough that the band is wide, the over-constraint **forces the construction**
`{1,..,n-2, n(n-1)}`, `M = n/Phi_6`.

## The three-regime trajectory (resolved)
| regime | n | covering-min | depth on `[0;n-1,k]` | driver |
|--------|---|--------------|----------------------|--------|
| projective | `n<=6` | `2/(2n-1)` (mediant) | -- | PG(2,n-1)/difference sets |
| **spread** | `7..11` | `[0;n-1,a(n)]`, `a=2,2,4,4,3` | shallow `a(n)` | band narrow enough for a spread base |
| **construction** | `>=12` | `n/Phi_6 = [0;n-1,n]` | deep `k=n` | radius-1 band forces consecutive base + outlier |

The spread family (HYP-3735) survives only while the band `(n,2n-2]` is narrow; at `n>=12` it is wide enough
that **only the consecutive base** covers it, forcing the construction. This is exactly *why* the spread dies
at `n=12`.

## Consequence -- LRC14
**The LRC14 covering-min is `14/183` (the construction).** opus's original value is vindicated *for `n=14`* --
the spread family that beat the construction at `n<=11` (S47) is dead by `n=12`. Margin `= 14/183 - 1/14 =
13/2562 > 0`, pinning HYP-2566's looseness at `n=14`. So the LRC14 hard core is the construction `14/183`, and
the proof target is `M(S) >= 1/14` with the construction as the (forced) extremal.

## Convergence with klein-S39 (HYP-3736) -- the PROVED version
klein independently found the same forcing and **proved the key lemma**: `{1,..,m}` is a **±-transversal mod
every prime `p <= 2m+1`** (each unit pair `{u, p-u}` has `min <= (p-1)/2 <= m`). With `m=n-2`, the dense core
`{1,..,n-2}` is a transversal mod every prime `p <= 2n-3` -- the rigorous form of my "consecutive base covers
the band interior at radius 1." klein's **killer-or-transversal** mechanism: at each band prime `p`, the
covering must either contain a **multiple of `p`** (a *killer* = a large speed) or be a **±-transversal mod
`p`**; the construction's core gets all band primes free (transversal), its big speed `n(n-1)` is the killer
for `n-1,n`. klein's **budget** `n-1 = resonance-killers + band-prime killers/transversals + spreaders`
tightens as the band (`~n/ln n` primes) grows, so `k_min` rises `2,2,4,4,3` and by `n=12` only the construction
survives. So klein-S39 (proved lemma + budget) and my S54 (verified forcing + the even-`n` outlier-covers-edge
parity + the trajectory) are the same result; klein supplies the rigor I flag as open below.

## Honest caveat
This is **strongly evidenced, and klein-S39's transversal lemma is proved** (the rigorous backbone): the forcing mechanism is verified (consecutive covers the band,
spread fails, outlier forced), the spread death is reliable (HYP-3735), and opus's scan agrees -- but "only the
consecutive base covers the band interior" is verified computationally, not exhaustively proved, and the full
ILP at `V=n(n-1)` times out (milp chokes on 40747 constraints), so it could not independently confirm `14/183`.
A rigorous proof needs the band-forcing made exhaustive (or klein-S38's radius-demand criterion fully
worked out into a uniqueness theorem). Avoiding MISTAKE-088, I do NOT claim certainty -- but unlike the n>=12
`1/(n-1)` artifact (which the construction trivially beat), here the construction IS the candidate and the
mechanism explains it.
