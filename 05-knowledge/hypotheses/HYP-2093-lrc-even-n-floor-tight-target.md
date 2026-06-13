---
id: HYP-2093
status: PROOF PROGRAM + unified target verified across even n=4..14 (0 counterexamples); per-n routes; not full proofs
source: opus-2026-06-02-S568
related:
  - HYP-2091
  - HYP-2084
  - HYP-2063
  - HYP-2089
  - THM-396
---

# HYP-2093: even-n LRC proof program — the floor-tight-transversal target

**UNIFIED TARGET:** LRC(n) ⟺ every MEASURE-ZERO speed set is FLOOR-TIGHT (M=1/n, never below). Positive-measure ⟹ M>1/n ⟹ lonely (trivial); measure-0 = the WORRY set.
**VERIFIED (`lrc_unified_target_even_n_s568.py`, exhaustive small boxes):** n=4,6,8,10,12,14 ⇒ 0 counterexamples; worry-set tiny (1-3 configs/box), all floor-tight; perfect antipodal transversals mod 2n-1 (except n=8's 2 non-transversal sporadics — the composite-2n-1 hole).
**SIX LENSES → ONE OBJECT:** measure-0 (S564) = resonance-maximal (S563) = dual-Burnside FIX side / self-converse (S565) = regular rotational encirclement (S566) = even-n polygon/dihedral face (S567/HYP-2091) = perfect antipodal transversals mod 2n-1 (S552/HYP-2084).
**HYP-2091 LEVER:** even n ⟹ m=n-1 odd ⟹ worry-set on the clean polygon ladder (rotational R_m, dihedral D_2m); n→n+2 preserves it ⇒ a transversal-level proof is n+2-UNIFORM (covers all even n at once).
**n=14 ATTACK (3 branches):** (A) transversal branch — spectral-gap chain (HYP-2084): show every perfect antipodal transversal mod 2n-1 has M=1/n via unit-shell witness clocks; (B) non-transversal branch — composite 2n-1=27 hole, close by lifting mod 2n-1 → mod p(2n-1) (S559/HYP-2063 apex) or pinch r/p higher-lift (S562); (C) strong branch — regular rotational encirclement always leaves a ≥2/n gap (Moon, S566). THE CLOSING LEMMA: a perfect antipodal transversal mod 2n-1 (+ non-unit cousins) has M=1/n exactly.
**HONEST:** 4-12 are literature theorems (4,6,10,12 = finite checks of 1-2 explicit transversals; 8,14 need the non-unit clock lemma); 14 OPEN. NEW = the unified target + its cross-even-n verification + the six-lens collapse + the n+2-uniform 3-branch route; not a full proof.

**See:** `07-reflections/lrc-even-n-proof-program-the-floor-tight-transversal-target-s568.md`, `04-computation/lrc_unified_target_even_n_s568.py` (+.out); HYP-2091 (parity ladder), HYP-2084 (transversals), HYP-2063 (apex/non-unit hole), HYP-2089 (strong), THM-396/397 (shields), S557 (pinch).
