# Message: oracle-2026-06-01-S526: proving small LRC with only the new methodology — n=3 PROVED (covering + mod-3 character), n>=4 reduced (HYP-2004)

**From:** oracle-2026-06-01-S?
**To:** all
**Sent:** 2026-06-01 14:57

---

User: try proving cases < n=14 using only the new methodology. RESULT: the covering reformulation (S525) + roots-of-unity characters (S522/S523) give |SAFE| = (1-2/n)^{n-1} + resonance corrections; the main term (1-2/n)^{n-1} is EXACTLY opus-S524's independence value ((6/7)^13=0.1348 at n=14). n=3 PROVED (complete, elementary): 2 runners a<b, gcd=1; m=2 => only 2-term resonances; the mod-3 character is the Legendre symbol chi; closed form |B_a ∩ B_b| = 4/9 + (2/9)chi(a)chi(b)/(ab) (verified vs numeric integration). Then |SAFE| = 1/9 + (2/9)chi(a)chi(b)/(ab) >= 1/9 - (2/9)/(ab) >= 0 (ab>=2), equality iff {a,b}={1,2}; at {1,2}, t=1/3 is lonely on the boundary. So LRC(n=3) holds for ALL speed pairs, {1,2} the unique tight case = the AP / regular-polygon extreme (S522). The proof is nothing but covering + the A_2/Paley-3 quadratic-residue character. WHY n>=4 is open: m>=3 => |SAFE| is an m-fold intersection whose resonance sum has 3-term and higher resonances sum k_i s_i=0 that do NOT factor into a single character = opus's multi-way correlation; the main term stays positive (->e^-2), so a proof needs only a SIGN/lower bound on the higher-resonance correction = the same wall as n=14. VERIFIED (lrc_small_n_covering_proof_s526.py): n=4..7 exact scans, 0 LRC failures, AP always tight; HONEST caveat: AP is NOT the unique tight set at n=5,6 ({1,3,4,7}; {1,3,4,5,9} also tight). Scorecard: n=2 trivial, n=3 PROVED, n>=4 reduced to bounding the higher-resonance character sum. New HYP-2004. HANDOFF: bound the 3-term resonance correction for n=4 (first case past pairwise; look for a mod-4 odd-harmonic character analogue of the n=3 Legendre trick). Files: 04-computation/lrc_small_n_covering_proof_s526.py (+.out); HYP-2004-lrc-covering-harmonic-n3-proof.md; reflection 07-reflections/lrc-small-n-covering-proof-s526.md.

---

*Reply by writing to `agents/oracle/inbox/` or run `python3 agents/processor.py --send --to oracle`*
