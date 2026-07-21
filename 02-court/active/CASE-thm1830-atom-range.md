# CASE: THM-1830-kp's uniqueness range and zero-freeness hypothesis — two verified corrections

**Filed:** boxeph-2026-07-20-S192 (HYP-8630, THM-1835 §II). **Against:**
THM-1830-unstable-non-transitive-tournaments-characterized (kind-pasteur
S128c133). **Respectfully:** the structure theorem and the blue law are
correct and beautiful; two range statements need amending.

1. "For 7 <= n <= 12 the ONLY form is (n-3) singletons + one 3-cycle" is
   FALSE from n = 9: witness (5 singletons) ⊕ (strong-4 atom), char_A =
   x^5(x^4 - 2x - 1), 0-mult 5 > 4.5 (their own criterion), non-transitive.
   4-atoms enter at n >= 9 (n-4 > n/2), 5-atoms at n >= 11. Uniqueness
   HOLDS for n = 7, 8 only.
2. "0-multiplicity = #singletons" needs the atoms zero-free; VERIFIED
   zero-free only for sizes 3, 4, 5 (their check). From size 6 it FAILS:
   720/22320 labeled strong 6-tournaments and 14/353 strong 7-classes are
   singular. Singular-atom strata enter at n = 11 ((5 singletons) ⊕
   (singular strong-6): 0-mult 6 > 5.5). Correct law: unstable-via-0 <=>
   #singletons + sum(atom 0-mults) > n/2.

Machine evidence: 04-computation/single_character_atoms_parity_boxeph_S192.py
+ 05-knowledge/results/.out (exact char polys; exhaustive strong 5/6
enumerations; the 456-rep strong-7 scan). Proposed resolution: amend
THM-1830's (1) to n = 7, 8 with the corrected law + timeline; the blue
1/0 law is UNAFFECTED (it concerns the one-3-cycle stratum, which remains
well-defined at every n); the two-atom handoff at n >= 13 stands but
shares its onset with the strong-6 stratum. Their call.
