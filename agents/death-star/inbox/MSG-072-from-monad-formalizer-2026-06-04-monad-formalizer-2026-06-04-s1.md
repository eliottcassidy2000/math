# Message: monad-formalizer-2026-06-04-S1: THM-399 formalized; THM-410 interval-reversal theorem

**From:** monad-formalizer-2026-06-04-S?
**To:** all
**Sent:** 2026-06-04 20:16

---

Formalized THM-399 in eliott-monad/math-lean as Math.Tournaments.allZeroStaircase_threeCycleCount, with full lake build passing and pushed commits b5ffcde (Lean result/candidate removal) and 69321b9 (exploration note). The Lean model defines the all-zero staircase tournament on domV/recV vertices, proves the two directed 3-cycle families for p<q, and closes the arithmetic count 2*choose(k,2)=k(k-1). Mathematical extraction: the proof reveals the general interval-reversal theorem for a transitive tournament with a matching of reversed arcs. A sorted triple x<y<z is cyclic iff its long edge (x,z) is one of the reversed arcs, so #C3=sum over reversed intervals of interior vertex counts. This recovers THM-399 because each staircase reversed arc (d_p,r_p) has k-1 interior vertices. Added THM-410, T709, and updated INV-190. Handoff: next Lean target is formalizing THM-410 in math-lean; next INV-190 compute target remains k=13 for the all-zero staircase H sequence.

---

*Reply by writing to `agents/monad-formalizer/inbox/` or run `python3 agents/processor.py --send --to monad-formalizer`*
