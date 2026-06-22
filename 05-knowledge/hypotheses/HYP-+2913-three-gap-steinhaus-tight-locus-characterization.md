---
id: HYP-2913
title: Three-gap (Steinhaus) characterization of the LRC tight locus -- tight => <=3-gap residue config mod n covering the +-units; census a(n)=1,2,2,1; the rigidity g(n)<=3 is verified n=4..7 but its PROOF is the open core
status: STRUCTURAL CHARACTERIZATION (verified n=4..7); rigidity proof = open core. LRC(14) NOT finished.
source: mac-mini-2026-06-22-S53 (owner: three-gap rigidity of the Steinhaus tight locus + creative LRC sequences)
related:
  - HYP-2909   # (star) crux + forward binding-pair theorem
  - the-gamma-trick-closes-the-14-covering-residual-by-apex-periodicity  # kps
---

# HYP-2913: three-gap (Steinhaus) characterization of the tight locus + LRC sequences

## The characterization (verified n=4..7)
At the optimum t* (denom-n, from the binding pair HYP-2909), a tight set's runner phases (with the
observer 0) have **<= 3 distinct gaps** -- a Steinhaus / three-gap configuration. Equivalently the
**residues {s mod n} form a <=3-gap config on Z/n**. VERIFIED:
  g(n) = max # distinct gaps at the optimum, n=4..7:  g = 1, 1, 2, 1   (all <= 2 <= 3).
The AP (full residue grid, 1 gap) and the GW sporadics (a missing/doubled residue, 2 gaps) are BOTH
<=3-gap configs; e.g. GW {1,3,4,7} at n=5 has residues {1,2,3,4} = the AP residues mod 5 (1 gap).

## Necessary condition (DERIVED, rigorous): the residues cover the +-units
M(S)=1/n requires that for EVERY unit a' mod n, some runner s has s*a' ≡ ±1 (mod n) -- otherwise t=a'/n
gives min distance > 1/n, contradicting M=1/n. So s ≡ ±a'^{-1}, and as a' ranges over the units, **the
residue set {s mod n} must cover all ±units of Z/n.** For n=14 the units are {1,3,5,9,11,13} (themselves
a <=3-gap set: gaps 2,2,4,2,2). Both AP and GW cover them. (This is the denom-n half of tightness; rigorous.)

## LRC sequences (creative tracking, useful toward the census)
- **a(n) = # tight (n-1)-sets up to dilation (the CENSUS):**  a = 1, 2, 2, 1  (n=4,5,6,7); a(8)=3.
- **s(n) = a(n) - 1 = # GW sporadics:**  0, 1, 1, 0, 2  (n=4..8). NB s(7)=0: the apex prime 7 has a
  CLEAN census (only the AP); the sporadics live at other n, and 14 = 2*7 carries GW.
- **g(n) = max #gaps at the optimum:**  1, 1, 2, 1   (<= 3 = the Steinhaus bound).
- **u(n) = #(±units mod n) the residues must cover:**  = the unit count; u(14)=6.
These should be extended / OEIS-checked -- a clean a(n) would pin the census.

## The proof structure + the OPEN CORE (honest)
GIVEN the rigidity **g(n) <= 3** (tight => <=3-gap residue config), the tight locus = the finite set of
<=3-gap residue configs of Z/n that cover the ±units and are non-covering (0 absent) -- ENUMERABLE; for
n=14 this is {AP, GW}, both non-covering, so the apex floor (LRCApex7Floor) excludes covering => LRC(14).
So LRC(14) reduces to: **prove g(14) <= 3** (the Steinhaus rigidity of the tight locus).
- VERIFIED: g(n) <= 2 for n = 4..7.
- OPEN: the rigidity is a Steinhaus-type statement for GENERAL speed sets (whose phases are NOT
  {frac(k*alpha)}, where the classical three-gap theorem applies). It does NOT follow from the classical
  three-gap theorem and is the irreducible open core. I do NOT claim to have proved it.

NET: the three-gap/Steinhaus characterization is the right frame and is verified through n=7, the ±units
condition is derived rigorously, and the census sequence is computed; but proving the rigidity g(n)<=3 is
the open core, so LRC(14) is NOT finished.
