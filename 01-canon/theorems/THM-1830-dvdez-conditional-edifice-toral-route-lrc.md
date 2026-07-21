# THM-1830: the DvdEZ-conditional edifice, and the toral route to LRC(14)

**Status:** DELIVERED — an EDIFICE MAP (what NC2 buys, each arrow graded),
one PROVED lemma (CT-realization), one NEW NAMED CONJECTURE FAMILY (TTNC,
with the proved TNC as its base case), the explicit reduction SHAPE
TTNC(14) => LRC(14) with its tight instances shown to be nonvanishing
Ramanujan-sum arithmetic, and the honest verdict that no completed
conditional proof of LRC(14) is claimed. HYP-8625.
**Author:** boxeph-2026-07-20-S191
**Owner:** "assume the 2D nullcone statement (Derksen-van den Essen-Zhao);
what can we prove? can we chain NC2 => GMC(2) => LRC(14)? get creative."

## 1. The edifice: what NC2 (assumed) makes unconditional

- GMC(2): via the repo's charge/polar reduction (klein THM-1645 + DvdK
  angular closure): NC2's algebraic classification is exactly the missing
  classification input — GMC(2) becomes a THEOREM. All span/degree
  partial results become instances.
- The S183-S187 analytic ladder: its CONCLUSION (no two-sided moment-
  nullcone member) holds outright; the remaining obligations (deleted-
  level tameness, cusp strata, (L2), far-end) become independent-interest
  analysis of the integral representation, no longer load-bearing.
- The Kempf-Ness form (S188/THM-1800): "analytic moment-nullcone = GIT
  nullcone of the hyperbolic torus w.r.t. the Fischer pairing" — becomes
  a theorem; DvdEZ = Kempf-Ness for the Bargmann pairing (quotable form).
- The moment engine (THM-1800): finite tests upgrade from decidability to
  syntactic CLASSIFICATION (nullcone membership <=> visible one-sidedness).
- The Fock chain toward JC2 (S166/HYP-8350, arrows still conjectural
  beyond the first): NC2 supplies arrow one of
  NC2 => GMC(2) => A_1-vacuum rigidity => {P_top,Q_top} obstruction =>
  DC_1 => JC2 — the only JC-adjacent chain that stays in N = 2 (the
  standard de Bondt-van den Essen route inflates dimension into the FALSE
  GMC(N>=3) zone). DIMENSION AUDIT: with ¬JC external, Mathieu(n_cx)
  fails at the counterexample's dimension and VC/IC fail with it; NC2 is
  2D-specific and UNTOUCHED — assuming it is consistent with everything
  known, and the Fock chain is where its truth would bite hardest.

## 2. The toral route to LRC(14) — the creative path

**CT-REALIZATION LEMMA (proved; verified 0.073093 = 0.073093):** with
u = e^{2 pi i t}, the LRC pairing is a CONSTANT-TERM functional:
  int_0^1 prod_j f_j(v_j t) dt = CT_u[ prod_j F_j(u^{v_j}) ].
So LRC lives on the TORAL side — its natural nullcone partner is the
PROVED TNC (THM-1605: CT(Lam^m) = 0 for all m <=> the small-root product
is ct), NOT GMC(2). The owner's chain re-routes: not NC2 => GMC(2) =>
LRC, but TNC-family => LRC, with NC2's role indirect (the edifice above).

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

No completed conditional proof of LRC(14) is claimed. What stands: the
edifice map (section 1, each arrow graded PROVED/CITED/CONJECTURAL); the
CT-realization lemma (proved); TTNC named with its proved base camp and
non-vacuous instances; the tight-instance arithmetic (c_15 != 0, exact);
and the identification of the repo's mod-q rungs as the program's finite
tests. The single sharpest next step: formulate and prove the TTNC
instance for ONE mod-15 rung and splice it into the existing LRC(14)
residual structure.

## 4. Files

04-computation/toral_route_lrc_boxeph_S191.py + frozen .out (T1-T3).
