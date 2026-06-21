#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadB_M2_proof_probe_opus_0621.py  (opus, 2026-06-21, THREAD B)

PROOF PROBE for "consec maximizes M_2" (the strongest single binding even-Krawtchouk moment).

M_2 = 2(E[N]-3)^2 - 3 + 2 Var(N).   K_2(t)=2t^2-12t+15 = 2(t-3)^2-3.
So M_2 = 2 E[(N-3)^2] - 3.  Maximizing M_2 == maximizing E[(N-3)^2] = the 2nd moment of N
about the CENTER 3 (the symmetry point of the 6 inner sectors).  This is a clean object:
  M_2/2 + 3/2 = E[(N-3)^2].   consec maximizes the spread of N about 3.

CANDIDATE PROOF MECHANISMS (test each, exact):

 (M1) POINTWISE-CENTERED coupling.  For E subset {0..maxE} with the SAME k, is there a
      coupling making (N_consec(x)-3)^2 >= (N_E(x)-3)^2 in convex order?  We test convex
      order of N about 3:  for every convex phi, E[phi(N_consec)] >= E[phi(N_E)]?  (Then
      ALL even Krawtchouk moments, being convex-in-(N-3) at the band, follow at once.)
      Convex order <=> same mean AND consec majorizes in the convex-order sense.  But means
      DIFFER (consec has the smallest mean).  So pure convex order is the WRONG tool; we test
      the SHIFTED/centered version and report which order actually holds.

 (M2) THE TWO ENDPOINTS DECOMPOSITION.  E[(N-3)^2] = sum_t (t-3)^2 p_t with weights
      w=(9,4,1,0,1,4,9).  Heaviest at t=0 (w=9) and t=6 (w=9).  kps-S9 PROVED:
      p_0=meas(S7) is consec-max (LEMMA B) and p_6=1/(7(k-1)) is consec-max (LEMMA A).
      Test: does the proved max of the two ENDPOINT cells (p_0,p_6) -- the heaviest weights --
      already DOMINATE the M_2 ordering?  i.e. is M_2 mostly carried by 9 p_0 + 9 p_6, both
      proved consec-max?  Decompose M_2 = 2[9 p_0 + 4 p_1 + p_2 + 0 + p_4 + 4 p_5 + 9 p_6]-3
      and report the endpoint share at consec and at the nearest competitors.

 (M3) THE MIDDLE CELLS.  The non-endpoint weighted mass 2[4p_1+p_2+p_4+4p_5] is what could
      break the endpoint domination.  Test: is 9p_0+9p_6 (consec) - 9p_0(E)-9p_6(E) ALWAYS
      >= the middle-cell deficit?  Report worst case.  (If yes on the finite window, M_2 max
      reduces to the two PROVED endpoint lemmas + a bounded middle correction.)

 (M4) k=8 FULL CLOSURE via endpoints.  L_y(k=8) = 1/16 + (1/40)M_2 + (3/80)M_4.  Express
      L_y directly in p_t (=g(t)): g=(1,0,0,1/10,0,0,1).  So L_y = p_0 + (1/10)p_3 + p_6.
      THE THREE BINDING CELLS ARE p_0, p_3, p_6.  p_0 (LEMMA B) and p_6 (LEMMA A) are PROVED
      consec-max.  So k=8 closure reduces to controlling p_3!  Test: is p_3 consec-max too?
      If p_0,p_3,p_6 are ALL consec-max, then L_y(k=8)=p_0+p_3/10+p_6 is consec-max with a
      THREE-CELL proof (two already PROVED).  THIS IS THE CLEANEST POSSIBLE k=8 STATEMENT.
"""
import sys, itertools
from math import comb, gcd
from functools import reduce
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

def occ(E):
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Enz:
        for m in range(0, 7 * e + 1): bps.add(F(m, 7 * e))
    bps = sorted(bps); p = [F(0)] * 7
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2; L = x1 - x0
        hit = set(int(7 * e * xm) % 7 for e in E)
        free = 6 - len([s for s in hit if s != 0]); p[free] += L
    return p
def prim(E):
    g = reduce(gcd, E); return [e // g for e in E] if g > 1 else E
def banner(t): print("\n" + "=" * 96 + f"\n{t}\n" + "=" * 96)

# ============ (M4) THE THREE-CELL k=8 STATEMENT: g=(1,0,0,1/10,0,0,1), cells p_0,p_3,p_6 ============
banner("(M4) k=8: L_y = p_0 + (1/10)p_3 + p_6.  Are p_0, p_3, p_6 EACH consec-max?")
print("  p_0 (=meas S7): kps-S9 LEMMA B PROVED consec-max.   p_6=1/(7(k-1)): LEMMA A PROVED consec-max.")
print("  Remaining: is p_3 consec-max?  Sweep the exhaustive window.")
for k, maxE in [(8, 16), (9, 14)]:
    pc = occ(list(range(k)))
    b0 = b3 = b6 = bLy = 0; n = 0
    worst_p3 = pc[3]; arg3 = 'consec'
    Lc = pc[0] + F(1,10)*pc[3] + pc[6] if k == 8 else None
    for rest in itertools.combinations(range(1, maxE + 1), k - 1):
        E = prim([0] + list(rest))
        if E[-1] != list(rest)[-1] and reduce(gcd,[0]+list(rest))==1:
            pass
        p = occ([0] + list(rest)); n += 1
        if p[0] > pc[0]: b0 += 1
        if p[3] > pc[3]: b3 += 1; worst_p3 = max(worst_p3, p[3]); arg3 = [0]+list(rest)
        if p[6] > pc[6]: b6 += 1
        if k == 8:
            Ly = p[0] + F(1,10)*p[3] + p[6]
            if Ly > Lc: bLy += 1
    print(f"  k={k} (maxE<={maxE}, {n} sets): p_0-beaters={b0}  p_3-beaters={b3}  p_6-beaters={b6}"
          + (f"  L_y-beaters={bLy}" if k==8 else ""))
    if b3:
        print(f"       p_3(consec)={float(pc[3]):.5f}  max p_3 in window={float(worst_p3):.5f} at {arg3}")

# ============ (M2/M3) endpoint domination of M_2 ============
banner("(M2/M3) M_2 = 2[9p_0+4p_1+p_2+p_4+4p_5+9p_6]-3.  Endpoint (9p_0+9p_6) domination?")
W = {9, 4, 1, 0, 1, 4, 9}
wt = [9, 4, 1, 0, 1, 4, 9]
for k, maxE in [(8, 13), (9, 13)]:
    pc = occ(list(range(k)))
    M2c = 2*sum(wt[t]*pc[t] for t in range(7)) - 3
    end_c = pc[0] + pc[6]  # the proved-max endpoints (sum)
    worst_def = None; nbeat = 0; n = 0
    for rest in itertools.combinations(range(1, maxE + 1), k - 1):
        E = [0] + list(rest); p = occ(E); n += 1
        M2 = 2*sum(wt[t]*p[t] for t in range(7)) - 3
        if M2 > M2c: nbeat += 1
        # decompose: endpoint surplus 9(p0c-p0)+9(p6c-p6) vs middle deficit
        end_surplus = 9*(pc[0]-p[0]) + 9*(pc[6]-p[6])
        mid_change = 4*(pc[1]-p[1]) + (pc[2]-p[2]) + (pc[4]-p[4]) + 4*(pc[5]-p[5])
        total = end_surplus + mid_change  # = (M2c-M2)/2
        if worst_def is None or total < worst_def[0]:
            worst_def = (total, end_surplus, mid_change, E)
    print(f"  k={k} (maxE<={maxE}, {n} sets): M_2-beaters={nbeat}")
    t, es, mc, E = worst_def
    print(f"     tightest: (M2c-M2)/2={float(t):.5f}  endpoint-surplus 9*(p0,p6 gap)={float(es):.5f}  "
          f"middle-change={float(mc):.5f}  at {E}")
    print(f"     => endpoints {'DOMINATE' if es>=-mc else 'do NOT dominate'} (endpoint surplus {'>=' if es>=-mc else '<'} middle deficit)")

# ============ (M1) which stochastic order actually holds ============
banner("(M1) Convex order of N about its mean?  Centered comparison consec vs competitors.")
print("  Pure convex order needs equal means; means differ.  We report mean & E[(N-3)^2] per set.")
for k in [8]:
    pc = occ(list(range(k)))
    muc = sum(t*pc[t] for t in range(7)); cen_c = sum((t-3)**2*pc[t] for t in range(7))
    print(f"  k={k} consec: mean={float(muc):.4f}  E[(N-3)^2]={float(cen_c):.4f}  (M_2=2*that-3={float(2*cen_c-3):.4f})")
    # show a few competitors
    cnt = 0
    for rest in itertools.combinations(range(1, 12), k-1):
        E = [0]+list(rest); p = occ(E)
        mu = sum(t*p[t] for t in range(7)); cen = sum((t-3)**2*p[t] for t in range(7))
        if cen > cen_c - F(1,2):  # near or above
            cnt += 1
            if cnt <= 6:
                print(f"     {E}: mean={float(mu):.4f} E[(N-3)^2]={float(cen):.4f} {'<' if cen<cen_c else '>='}consec")
print("\nDONE.")
