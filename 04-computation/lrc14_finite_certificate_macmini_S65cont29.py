#!/usr/bin/env python3
"""
lrc14_finite_certificate_macmini_S65cont29.py -- THM-702: the EXPLICIT finite certificate
for kps's THM-701 wide-spread recursion closure.

Three bookkeeping items kps left named-not-written, made exact here:
(2) cap-growth table: cap_{k+1} - cap_k >= 2/21 for k=8..12, EXACT rationals;
(3) the explicit far-element threshold via THM-699's PROVEN constant:
    per-peel error <= C(E')/w, C(E') = 672 * (sum E')  [K7 = sum_t C(7,t) t^2 (7-t)/7 = 672],
    budget: allocate half the worst slack s* = min_k (growth_k - 2/21)/1 -> threshold
    W(E') = C(E')/s_half; per-step induction needs NO summability (each peel self-budgets);
(1) the balanced-core margins at the known argmax (consecutive) -- exact Phi values,
    and the HONEST reduction: exhaustive core check is infeasible at B ~ 672 Sum e'/s;
    the core closure = Phi-consec-extremality (the same single extremal lemma as THM-534).
"""
from fractions import Fraction as F
from itertools import combinations

CAP = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91), 11: F(66,91), 12: F(6,7), 13: F(1)}
K7 = sum(F(__import__('math').comb(7,t) * t*t * (7-t), 7) for t in range(1,8))
print(f"K7 (THM-699 IE constant) = {K7} (= 672 exact: {K7 == 672})")
print()
print("(2) CAP-GROWTH TABLE (exact): growth_k = cap_(k+1) - cap_k vs 2/21 =", F(2,21), f"~ {float(F(2,21)):.4f}")
worst_slack = None
for k in range(8, 13):
    g = CAP[k+1] - CAP[k]
    s = g - F(2,21)
    ok = "OK" if s > 0 else "*** FAIL ***"
    if worst_slack is None or s < worst_slack: worst_slack = s
    print(f"  k={k}: growth = {g} ~ {float(g):.4f}, slack = {s} ~ {float(s):.4f}  {ok}")
print(f"  WORST SLACK s* = {worst_slack} ~ {float(worst_slack):.5f} (at k=8) -- ALL POSITIVE: (2) PROVED")
print()
print("(3) THE EXPLICIT THRESHOLD: half-slack budget s*/2 =", worst_slack/2)
print(f"    W(E') = C(E')/(s*/2) = 672*(sum E') / {float(worst_slack/2):.5f} = {float(672/(worst_slack/2)):.0f} * (sum E')")
Wcoef = 672 / (worst_slack/2)
print(f"    => far-element threshold: w >= {float(Wcoef):.0f} * (sum of remaining offsets)")
print("    Per-step induction: Phi(E) <= Phi(E') + 2/21 + C(E')/w <= Phi(E') + 2/21 + s*/2")
print("                        <= cap_k + growth_k  whenever Phi(E') <= cap_k. NO summability needed.")
print()
print("(1) BALANCED-CORE margins at the consec argmax (exact Phi = p0 + p1/3):")
def sector_counts(E, x_den_limit=None):
    # exact p_j for offsets E: measure of {x: exactly j sectors empty}... compute p0,p1 via IE:
    # p0 = meas(no sector empty) -- complement of measS7; p1 = meas(exactly one empty).
    # exact breakpoint sweep over all sector boundaries:
    pts = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for m in range(e):
            for s in range(8):
                pts.add(F(m,e) + F(s, 7*e))
    pts = sorted(p for p in pts if 0 <= p <= 1)
    p = [F(0)]*8
    for a,b in zip(pts, pts[1:]):
        mid = (a+b)/2
        hit = set()
        for e in E:
            if e == 0: hit.add(0); continue
            fr = e*mid - (e*mid).__floor__()
            hit.add(int(fr*7))
        empty = 7 - len(hit)
        p[empty] += b - a
    return p
for k in range(8, 11):
    E = list(range(k))  # consecutive core 0..k-1
    p = sector_counts(E)
    Phi = p[0] + p[1]/3
    m = CAP[k+1] - Phi   # THM-701 induction: Phi(F) <= cap_(|F|+1) (peel target p0(E)<=Phi(E'), |E|=|F|+1)
    print(f"  core |F|={k} consec: p0={float(p[0]):.4f} p1={float(p[1]):.4f} Phi={float(Phi):.4f} cap_(k+1)={float(CAP[k+1]):.4f} margin={float(m):+.4f} {'OK' if m>0 else '*** FAIL ***'}")
print()
print("HONEST SCOPE: (2),(3) PROVED exact. (1) at the known argmax has positive margin;")
print("the exhaustive core sweep at B ~ W*(sum e') is infeasible -- core closure = the")
print("Phi-consec-extremality lemma, THE single remaining extremal statement (as THM-534).")
