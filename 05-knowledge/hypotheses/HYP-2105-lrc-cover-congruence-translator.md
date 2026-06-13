---
id: HYP-2105
status: NEW THEOREM Lemma C (PROVED, w-free) via cover->congruence translator; Lemma C ∪ B' covers 99% of mult-of-14; large-owner residual OPEN
source: opus-2026-06-03-S574
related:
  - THM-398
  - HYP-2104
  - HYP-2103
  - HYP-2102
---

# HYP-2105: cover->congruence translator and the small-owner theorem (Lemma C)

**TRANSLATOR (verified 100%, 2500/2500 each n=6..14):** for S=S'u{v=nw}, n|v, tight <=> G(S') subset D_v <=> every component (a,b) of G(S') fits one v-arc (centre j/(nw), radius 1/(n^2 w)). Endpoints carry OWNERS: left a=(k_a n+1)/(n u_a), right b=(k_b n-1)/(n u_b). "(a,b) in arc j" translates (x n u w) to ENDPOINT-OWNER CONGRUENCES |w(k_a n+1)-j u_a|<u_a/n and |w(k_b n-1)-j u_b|<u_b/n.
**RIGIDITY:** for owner u<n the RHS<1 => bracket integer => =0 => endpoint EXACTLY at the arc centre. Slack appears only for u>=n.
**LEMMA C (PROVED, w-free):** if a G(S') component has BOTH owners <n, then S=S'u{nw} is loose for EVERY w (two rigid endpoints can't be the same centre => a=b). The cross-relation u_b(k_a n+1)=u_a(k_b n-1) two small owners would need IS a=b: verified 0 hold / 490,1330,2594,5298 (n=6,8,10,12). Small-owner congruences are infeasible by inspection.
**COVERAGE (lrc_smallowner_coverage_s574.py):** Lemma C alone 8.2/18.8/33.4/56.2/81.3% (n=6/8/10/12/14, GROWS toward frontier); Lemma C ∪ Criterion B' (Vitali long-interval) = 73.4/81.7/92.0/96.0/**99.0%**. At n=14, 99% of mult-of-14 configs PROVED loose, no enumeration.
**RESIDUAL (~1% at n=14, OPEN):** configs where every G(S') component is short AND every component has >=1 LARGE binding owner (>=n) -- only there does the congruence slack +-u/n permit an off-centre fit. Now a bounded CRT/Diophantine feasibility on the large owners + w (verified never satisfiable, not proved).
**FRAME (user's redirection):** enumeration only gives "0 in a box"; the translator converts geometry (periodic arcs covering a union of intervals) to arithmetic (congruences on endpoint owners) where infeasibility is PROVABLE -- small owners immediate (a=b), large owners a bounded CRT check.

**See:** `01-canon/theorems/THM-398-...md` (sec 4.5, Lemma C), `07-reflections/lrc-cover-to-congruence-translator-and-the-small-owner-theorem-s574.md`, `04-computation/lrc_cover_to_congruence_translator_s574.py` (+.out), `lrc_smallowner_coverage_s574.py` (+.out); HYP-2104/2103/2102, THM-398.
