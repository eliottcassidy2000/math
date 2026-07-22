# THM-1830: the DvdEZ-conditional edifice, and the toral route to LRC(14)

**Status:** SUPERSEDED CONDITIONAL EDIFICE MAP — not a theorem-level NC2/GMC
bridge. One PROVED lemma (CT-realization), one named conjecture family (TTNC,
with the proved TNC as its base case), the explicit reduction SHAPE
TTNC(14) => LRC(14) with its tight instances shown to be nonvanishing
Ramanujan-sum arithmetic, and the honest verdict that no completed
conditional proof of LRC(14) is claimed. HYP-8625.
**Author:** boxeph-2026-07-20-S191
**Owner:** "assume the 2D nullcone statement (Derksen-van den Essen-Zhao);
what can we prove? can we chain NC2 => GMC(2) => LRC(14)? get creative."

> **Current scope (MISTAKE-235/237).** THM-2022 proves GMC(2) by its own
> lowest-face/Frobenius mechanism. Nothing here proves an
> `NC2 -> GMC(2) -> JC(2)` or `NC2/GMC(2) -> LRC(14)` chain. The CT-realization
> lemma and Ramanujan calculations survive; TTNC and the Fock/LRC arrows remain
> conjectural programs. Read this file as provenance for those objects only.

## 1. Historical conditional edifice — not a proof graph

The original session proposed the following objects while treating NC2 as a
hypothesis. They are retained as research provenance, with their current
status made explicit:

- **GMC(2): proved elsewhere.** THM-2022 proves NC2 and GMC(2) together by a
  lowest-face/Frobenius argument. The charge-additivity reduction from a
  strict one-sided nullcone classification is sound, but the older
  `TNC -> NC2` or generic polar-averaging bridge is not.
- **S183--S187 analytic ladder: refinement only.** Its desired no-two-sided
  conclusion follows from THM-2022, while deleted-level tameness, cusp strata,
  `(L2)`, and far-end estimates remain separate quantitative questions.
- **Kempf--Ness and finite moment-engine formulations: proposed, not promoted
  here.** This file does not prove equality of the analytic moment nullcone
  with a GIT nullcone, nor a finite syntactic classification theorem.
- **Fock/JC2 chain: conjectural.** No arrow from NC2 or GMC(2) to JC(2) is
  proved here. THM-2071 separately excludes every quadratic member in an
  output-pencil/source-foliation direction and gives a tame normal form; that
  quadratic-pencil theorem is not a completion of this chain or of JC(2).

## 2. The toral route to LRC(14) — the creative path

**CT-REALIZATION LEMMA (proved; verified 0.073093 = 0.073093):** with
u = e^{2 pi i t}, the LRC pairing is a CONSTANT-TERM functional:
  int_0^1 prod_j f_j(v_j t) dt = CT_u[ prod_j F_j(u^{v_j}) ].
So LRC lives on the TORAL side — its natural nullcone analogue is the
PROVED TNC (THM-1550 gives the small-root criterion; THM-1605/2067 rule out
the both-signs case), not GMC(2). This is an analogy and a proposed twisted
extension, not a proved `TNC-family => LRC` implication.

**TTNC (new, proposed):** the TWISTED TORAL NULLCONE family: for a
nonnegative trig polynomial Theta and Laurent Lam,
  CT[Theta * Lam^m] = 0 for all m >= 0   <=>?   an explicit support/orbit
  obstruction between supp(Theta-hat) and the multiplicative structure of
  Lam. Theta = 1 is TNC (PROVED). The m = 0 member alone gives
  CT[Theta] = |G| — so Theta = good-set approximant makes twisted-
  nullcone membership CONTAIN the covering condition.

**The reduction shape TTNC(14) => LRC(14):** if v were a 14-runner
covering counterexample at delta <= 1/15, then Theta_delta = 0, i.e. the
good-set object lies (trivially) in every twisted nullcone; contrapositive
route: prove the specific TTNC(14) instances — that the ONLY member of
the twisted nullcone for the 14-speed bad-set complements is Theta = 0 —
while exhibiting a nonzero twisted moment for every admissible v. Sanity
(verified): at n = 3, nonzero good sets have ALL twisted moments m <= 12
visible (tight and loose) — the instances are non-vacuous.

**The cyclotomic sharpening (the S190 Galois stack, now arithmetic):**
at tight configurations the threshold good measure concentrates on the
primitive q-th fractions (q = n+1) and its twisted moments are RAMANUJAN
SUMS c_q(m) = mu(q/g) phi(q)/phi(q/g), g = gcd(m, q). For 14 runners
q = 15 is SQUAREFREE, so **c_15(m) != 0 for ALL m** (computed:
8,1,1,-2,1,-4,-2,1,1,-2,-4,1,-2,1,1,8) — the tight-instance twisted
nullcone at q = 15 is EMPTY by checkable arithmetic. (Contrast n = 3,
q = 4: c_4 has zeros — the 14-runner case is arithmetically NICER than
the toy.) The named program: extend the c_15-nonvanishing from the tight
atomic limit to a neighborhood (the 'spread' conditions of the repo's
LRCMod ladders are exactly relation-lattice nondegeneracy mod q — the
rungs are this program's finite tests).

## 3. Honest verdict

> **2026-07-21 update (THM-2041 / HYP-8800).** One TTNC service is now
> proved: primitive-support and Ramanujan exact-period twists are fixed as
> whole layers by every coprime Frobenius prime, so a nonzero twisted cyclic
> moment propagates through all p-power levels. This does not prove the TTNC
> instance proposed above. It supplies preservation after base nonvanishing;
> it does not produce a base nonzero moment, isolate off-period packets, or
> turn a signed trace into a pointwise safe phase. The sharp next target is a
> labelled packet whose projected residue is an actual bounded integer safe
> count, not merely a Ramanujan statistic.

No completed conditional proof of LRC(14) is claimed. What stands: the
historical object map in section 1, with all arrows now treated as proposed
unless independently cited; the
CT-realization lemma (proved); TTNC named with its proved base camp and
non-vacuous instances; the tight-instance arithmetic (c_15 != 0, exact);
and the identification of the repo's mod-q rungs as the program's finite
tests. The single sharpest next step: formulate and prove the TTNC
instance for ONE mod-15 rung and splice it into the existing LRC(14)
residual structure.

## 4. Files

04-computation/toral_route_lrc_boxeph_S191.py + frozen .out (T1-T3).
