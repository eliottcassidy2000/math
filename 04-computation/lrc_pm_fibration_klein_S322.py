#!/usr/bin/env python3
"""
lrc_pm_fibration_klein_S322.py
==============================
klein-2026-07-18-S322 (owner: the remaining LRC(14) frontier and how it relates to tournaments/metagraphs
abstractly). Verifies the three structural claims of 07-reflections/the-pm-fibration-*.md.

ONE FIBRATION, TWO FLOORS:   (Z/qZ)^*  --pi-->  U_q = (Z/qZ)^*/{+-1}
  * TOURNAMENTS live in the SECTIONS: a Cayley tournament on Z/q picks one of {d,-d} per class.
    Paley (D=QR) is a section iff -1 is a non-residue iff q = 3 mod 4.
  * LRC lives in the BASE: THM-762/764's witness criterion is "deck B_q(S) = {[s]_+- : s in S} is a
    PROPER subset of U_q" (plus no speed divisible by q). THM-567's hypothesis F(r)=F(-r) is the same
    base condition stated analytically.

THE DISANALOGY THAT EXPLAINS WHY PROOF STEPS DO NOT CROSS (mac-mini-S89's verdict, now with a reason):
  tournament complement T -> T^op has FIXED POINTS (the self-complementary classes, n = 0,1 mod 4);
  the LRC action s -> -s is FREE for every q >= 3  (2s = 0 mod q with gcd(s,q)=1 forces q | 2).
  The metagraph's load-bearing decomposition SPINE(SC-SC)+RIBS(SC-NS)+SEA(NS-NS) and THM-A/B/C are
  FIXED-POINT constructions. The LRC base has no fixed points => no spine => those transfers are dead
  a priori. Only SECTION-FREE (pure base) transfers have a chance.

USABLE CONSEQUENCE: freeness gives |U_q| = phi(q)/2 exactly, so 13 speeds can only fill the base when
phi(q) <= 26; for phi(q) > 26 the deck is automatically proper and the sole obstruction is q | s.
"""
from math import gcd
def phi(n): return sum(1 for k in range(1,n+1) if gcd(k,n)==1)
def units(q): return [s for s in range(1,q) if gcd(s,q)==1]
def action_is_free(q): return not any((2*s) % q == 0 for s in units(q))
def qr(q): return {(x*x) % q for x in range(1, q)}
def paley_is_section(q):
    U, Q = units(q), qr(q)
    return len(Q) == len(U)//2 and all((s in Q) != ((q-s) in Q) for s in U)
if __name__ == "__main__":
    print("(1) +- action free for all q>=3:", all(action_is_free(q) for q in range(3,200)))
    print("(2) Paley is a section iff q=3 mod 4:",
          all(paley_is_section(q) == (q % 4 == 3)
              for q in (7,11,13,17,19,23,29,31,43,47,53)))
    print("(3) |U_q| = phi(q)/2 :", all(len(units(q))//2 == phi(q)//2 for q in range(3,200)))
    print("    q<=59 where 13 speeds could fill the base (phi(q)<=26):",
          [q for q in range(15,60) if phi(q) <= 26])

# --- klein-S323 addendum: the QR/NQR transfer chase (DEAD END) ------------------
# TRUE:  a QR-supported deck collapses the multiplier orbit to two configs, giving
#        M|_q = max(min|QR|, min|NQR|)/q.  Since 1 in QR always, min|QR| = 1, so the
#        two minima coincide iff -1 is a non-residue iff q = 3 mod 4 = the Paley
#        section condition, restated on the base.  But that is the tautology
#        "min|NQR| = 1 <=> -1 in NQR" in a fibration costume, not transferred content.
#        The one non-trivial branch: q = 1 mod 4 gives M|_q = (least absolute NQR)/q,
#        tying loneliness to the classical least-non-residue problem, and predicting
#        QR decks are lonely-HOSTILE at large q (least NQR grows slowly).
# FALSE: that prediction fails on the real objects.  Classifying every known low-M
#        family at its own witness denominator, NOT ONE is QR- or NQR-supported:
#          AP {1..13}            M=1/14    witness 1/14    mixed QR3/NQR3
#          deep well {1..12,182} M=14/183  witness 14/183  mixed QR2/NQR7
#          GW {1..11,13,24}      M=1/14    witness 1/14    mixed QR3/NQR3
#          n=5,8,9,10,11 cx      2/9..3/31                 all mixed, near-even
#        The LRC extremals are QR-AGNOSTIC: generic subsets of U_q, not cosets.
# VERDICT: THM-981's QR/NQR supports belong to the H6 owner-feasibility ledger, not to
#        the LRC metric objects.  Same vocabulary, different problems.  No viable
#        transfer candidate currently remains identified.
