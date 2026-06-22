---
id: HYP-2914
title: The LRC tight "kind" -- a convoluted battery of necessary conditions (covers{2..n-1}, survives n, residues co-finite missing <=1, zsum=MAX binding pairs); the kind = AP-skip-double with divisibility compensation. Parity/char-energy/apex-count REFUTED as common.
status: NECESSARY conditions (rigorous + verified n=5,6,8,14); the convoluted characterization of the tight kind. (Sufficiency = the open core.)
source: mac-mini-2026-06-22-S54 (owner: derive more necessary conditions, what AP & GW & any of their kind share; abstract/convoluted/creative)
related:
  - HYP-2913   # three-gap/Steinhaus characterization
  - HYP-2909   # binding-pair forward theorem
---

# HYP-2914: the tight "kind" -- a battery of necessary conditions

What do AP and GW (and any tight set) share? Testing a battery across n=5,6,8,14:

## RIGOROUS necessary conditions
1. **Almost-covering: kills every b in {2,...,n-1}** (a runner ≡ 0 mod b for each). Else t=a/b gives
   min distance >= 1/b > 1/n, so M > 1/n. (Resonance-killing; rigorous.)
2. **Survives n:** no runner ≡ 0 mod n (non-covering). Else t=a/n is unsafe; M cannot be the 1/n witness.
3. **Residues cover the ±units mod n:** for every unit a', some runner ≡ ±a'^{-1} (else t=a'/n beats 1/n).

## VERIFIED-common conditions (n=5,6,8,14)
4. **The residue SET is co-finite, missing <= 1:** the distinct residues mod n equal {1,...,n-1} (AP)
   or {1,...,n-1}\{k} for a single k (GW). The residue-GAP multiset is {1: n} (AP) or {1: n-2, 2: 1}
   (GW) -- only gap-sizes 1 and 2, with AT MOST ONE size-2 gap. (Far stronger than the bare <=3-gap.)
5. **zsum = MAX binding pairs:** the number of pairs s_i+s_j ≡ 0 (mod n) equals (n-2)/2 (n even) /
   (n-1)/2 (n odd) -- ALL antipodal residue classes are realized. When GW skips residue k, the loss is
   EXACTLY compensated by a DOUBLED residue j (e.g. n=14: skip 12, residues 10 & 24 both ≡10 give the
   pair (4,24) replacing the lost (2,12)). So zsum is conserved at the maximum.

## The convoluted synthesis: the kind = "AP-skip-double with divisibility compensation"
A tight set is the complete residue system {1,...,n-1} (the AP), OR the AP with one residue k SKIPPED
and one residue j DOUBLED, where the doubled runner's LIFT re-kills the b's that the skip lost. n=14:
GW skips residue 12 and doubles residue 10 via the lift **24 = 2*12** -- 24 ≡ 10 (mod 14) (doubles the
residue) AND 12 | 24 (re-kills b=12, which only residue 12 killed in the AP). The "kind" is thus the
orbit of the AP under **divisibility-preserving residue-twists** (skip a residue, double another via a
higher multiple that carries the lost divisibility). The apex 7 = n/2 is the self-antipodal residue
(the only one without a distinct partner), which is why 14 = 2*7 is the critical modulus.

## REFUTED as common (discipline -- tested, do NOT hold)
- **Parity** #odd - #even: varies (n6: AP=(3,2), GW=(4,1); n8: GW1=(4,3), GW2=(5,2)).
- **Apex-runner count** (#{s ≡ n/2}): varies (1 or 2).
- **Character-energy** Σ_k |Σ_s e(sk/n)|^2: varies (n14: AP=13, GW=41).
- **Sum/product of residues mod n:** vary.

## Status (honest)
These are NECESSARY conditions (1-3 rigorous; 4-5 verified n=5,6,8,14), tightly constraining the kind to
"AP-skip-double with divisibility compensation." They are not yet shown SUFFICIENT, and proving that the
list FORCES the tight locus to be exactly {AP, GW} for n=14 (no other skip-double survives all of 1-5)
is the remaining open core (the consec-maximizes / Steinhaus rigidity). The battery is a sharp net; the
catch is the open part.
