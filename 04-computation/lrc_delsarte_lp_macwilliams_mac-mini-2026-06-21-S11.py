#!/usr/bin/env python3
"""
lrc_delsarte_lp_macwilliams — mac-mini-2026-06-21-S11

THE DELSARTE-LP UNIFICATION of the LRC(14) cover bound (unifies kps HYP-2723 relation-code/MDS,
mac-mini HYP-2724 even-Krawtchouk, and THM-534 moment-LP).

THM-534's moment-LP dual g(t) (g(t) >= 1[t=0] on t=0..6 gives measS7(E) <= L_y(E) = sum_r y_r S_r
for every E) expands in the binary KRAWTCHOUK basis K_j(t;6) with ALL-NONNEGATIVE coefficients --
i.e. it is exactly a DELSARTE linear-programming certificate. So the LRC cover bound IS a Delsarte
LP on the 7-cell (Z/7) occupancy scheme; the relation code Lambda(E) (HYP-2723) is the Delsarte
scheme; consec = the LP-tight (anti-MDS) extremal config.

VERIFIED here (exact Fractions): the THM-534 duals for k=8 / k=9,10 / k=11,12,13 are all
Krawtchouk-nonnegative (Delsarte-positive). k=8 is purely EVEN-degree (K0,K2,K4) -> explains
HYP-2724's even-band cleanliness; k>=9 use both parities, all >=0.
"""
from fractions import Fraction as F
from math import comb

def Kraw(j, t, n=6):  # binary Krawtchouk K_j(t) on the 6 inner sectors
    return sum((-1)**i * comb(t, i) * comb(n - t, j - i) for i in range(j + 1))

def kraw_expand(g):   # solve sum_j c_j K_j(t) = g(t), t=0..6, over Q
    n = 7
    Aug = [[F(Kraw(j, t)) for j in range(7)] + [g[t]] for t in range(7)]
    for col in range(n):
        piv = next(r for r in range(col, n) if Aug[r][col] != 0)
        Aug[col], Aug[piv] = Aug[piv], Aug[col]
        pv = Aug[col][col]; Aug[col] = [x / pv for x in Aug[col]]
        for r in range(n):
            if r != col and Aug[r][col] != 0:
                f = Aug[r][col]; Aug[r] = [Aug[r][i] - f * Aug[col][i] for i in range(n + 1)]
    return [Aug[i][n] for i in range(n)]

duals = {
    'k=8       (deg4)': [F((t-1)*(t-2)*(t-4)*(t-5), 40) for t in range(7)],
    'k=9,10    (deg3)': [F(-(t-2)*(t-3)*(t-6), 36) for t in range(7)],
    'k=11,12,13(deg2)': [F((t-3)*(t-4), 12) for t in range(7)],
}
print("THM-534 moment-LP duals as DELSARTE certificates (Krawtchouk-nonnegative):")
for name, g in duals.items():
    c = kraw_expand(g)
    feasible = all(g[t] >= (1 if t == 0 else 0) for t in range(7))
    delsarte = c[0] > 0 and all(c[j] >= 0 for j in range(1, 7))
    even_only = all(c[j] == 0 for j in range(1, 7, 2))
    print(f"  {name}: g={[str(x) for x in g]}")
    print(f"       Krawtchouk coeffs c={[str(x) for x in c]}")
    print(f"       g>=1[t=0]:{feasible}  Delsarte(all c>=0):{delsarte}  even-only:{even_only}")
print("\n=> the LRC cover bound is a DELSARTE LP; the dual is Krawtchouk-nonnegative at every binding k.")
print("   Relation code Lambda(E)=Delsarte scheme (HYP-2723); even-band (HYP-2724) = k=8 even-only case.")
