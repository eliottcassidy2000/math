---
id: HYP-2120
status: SYNTHESIS+CORRECTION -- LRC's slice of the perspectives curiosity is the SOURCE-perspective recursion (exact, THM-381), not the full count (fails at n=6); worry-set = source-less regular perspective
source: opus-2026-06-03-S582
related: [THM-381, THM-385, THM-400, HYP-1977, HYP-1824]
---

# HYP-2120: the perspectives curiosity -- the SOURCE slice is the LRC key

See `07-reflections/lrc-the-perspectives-curiosity-source-slice-not-full-count-s582.md` for the full account. Summary:
- Perspective Conjecture (T075): #vertex-orbits(n) = #iso-classes(n+1). Works n=3 (4=3+1=T(4)), n=4 (12=4+4+2+2=T(5)); FAILS n=5 (full perspectives 48 != T(6)=56, gap 8 = chirality => the 56-bridge/HYP-1824/1825 detour).
- VERIFIED: SOURCE-perspectives (vertex-orbits that are sources) = 1,2,4 = A000568(n-1) EXACTLY (= THM-381). The LRC-relevant recursion is source(n)=T(n-1), clean & gapless; deleting the source is canonical.
- WORRY-SET = source-LESS regular perspective (n=5: the (1,0) regular class, observer never a source = M=1/n boundary). Source-having classes => observer strict source => strictly lonely (easy).
- UNIFIES with THM-400 (observer-coupled eps!=0 = rooted = perspective; observer-blind eps=0 = unrooted = structure) and THM-381/S509 (root at observer; unrooted class doesn't determine loneliness = projection defect).
- MISINTERPRETATION: LRC was analyzed unrooted (round/self-converse iso-classes) and chased the failing full-perspective gap; the key is the source slice + its source-less worry fiber. WHERE/WHY works-and-doesn't: skeleton source(n)=T(n-1) exact (no gap); content = source-less fiber, degenerate at composite seam 14=2*7 (S580).
