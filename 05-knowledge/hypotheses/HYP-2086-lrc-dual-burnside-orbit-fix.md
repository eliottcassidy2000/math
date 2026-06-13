---
id: HYP-2086
status: SYNTHESIS — LRC's worry/ignore split IS Burnside's fix/orbit split; verified anchors; gives a containment (worry ⊆ self-converse boundary)
source: opus-2026-06-02-S565
related:
  - HYP-1998
  - HYP-2080
  - HYP-2081
---

# HYP-2086: dual Burnside — LRC open=orbit(round/A000016), boundary=fix(self-converse)

Burnside has an ORBIT face (count classes) and a FIX face (symmetric configs). The LRC's two regimes ARE these two faces.
- **ORBIT side (IGNORE, positive measure):** open-time runner sub-tournament is ROUND = A000016 (HYP-1998). Verified counts m=3..13: 2,4,10,30,94,316 vs A000568 = ...,48.5e12; ratio 6.5e-12 at m=13. Generic, automatically lonely.
- **FIX side (WORRY, measure 0):** the loneliest config (regular round tournament) is SELF-CONVERSE (T=T^op via i->-i, verified m=3..13), rotation-symmetric (C_m), and carries the (q,q) automorphism with Burnside Fix=2^{n-1} (S547). The worry/tight family sits at these maximally-symmetric fixed points.
**DUALITY:** worry/ignore (S564) = Burnside fix/orbit. CONTAINMENT: a counterexample's sub-tournament can't be a generic round class (open=lonely) ⇒ must be a boundary SELF-CONVERSE class ⇒ worry-set ⊆ self-converse boundary classes (thin, Burnside-counted; cf. repo A002785/A005639). The lonely menu sits between: A000016 ⊆ menu ⊆ A000016 ∪ {self-converse boundary} (HYP-1998 open-Q-B); the surplus = the fixed-point layer.
**PROOF SHAPE:** show every self-converse boundary runner tournament is source-realizable (lonely); the round/open side is free. = the Burnside translation of 'prove LRC only for the resonance-maximal tight family'.
**HONEST:** orbit count (A000016) & fix anchors ((q,q), self-converse) known/verified; NEW = the duality framing + the worry ⊆ self-converse-boundary containment. Not a proof.
**See:** `07-reflections/lrc-dual-burnside-orbit-round-fix-self-converse-s565.md`, `04-computation/lrc_dual_burnside_s565.py` (+.out); HYP-1998 (round=A000016), S547, HYP-2080 (resonance), HYP-2081 (clocks).
