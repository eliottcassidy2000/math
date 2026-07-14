#!/usr/bin/env python3
"""
lrc14_thm735_bonferroni_multipeel_kps_S128c3.py
===============================================
kind-pasteur-2026-07-13-S128 (cont.3).  THM-735 (i)+(ii): the SIMULTANEOUS (Bonferroni) MULTI-PEEL.

LEMMA (i):  L(E u F) >= (1 - j/7) m_E - Sum_{v in F} |eps_v(E)|,   j=|F| <= 6,
            eps_v(E) = |D_v cap G_E| - m_E/7,   |eps_v(E)|^2 <= (6/49) disc_v(G_E)   [THM-731 on G_E]
CRUDE (ii): disc_v(G_E) <= r_E^2/(3v^2)  [THM-732]  =>  closure whenever
            Sum_{v in F} 1/v < (7-j) m_E / (sqrt2 r_E).
KEY POINT:  every eps is measured against the FIXED body good set G_E -- no peel ever faces a base
carved by the other far arcs. Sequential peeling does (that carving = klein-S289's isolation wall).
So (ii) is CLUSTERING-IMMUNE: F may be consecutive integers.

This script (MISTAKE-136 rule -- verify the chain itself, exactly, in Q):
 (1) for a battery of (E,F): exact m_E, r_E, eps_v, disc_v, L_true; verify EVERY link:
       L_true >= (1-j/7)m_E - Sum|eps_v|          [lemma]
       eps_v^2 <= (6/49) disc_v                    [per-peel covariance bound]
       disc_v <= r_E^2/(3 v^2)                     [crude disc]
     and evaluate the two computable positivity certificates:
       CRUDE:      (1-j/7) m_E - (99/70)(r_E/7) Sum 1/v > 0            [99/70 > sqrt2, safe]
       EXACT-DISC: (1-j/7)^2 m_E^2 > j (6/49) Sum_v disc_v             [Cauchy-Schwarz, safe]
 (2) the sequential-vs-simultaneous DEMO on clustered triples {1..10, c, c+1, c+2}:
     sequential crude (peel b=c+2 from the exact 12-speed base {1..10,c,c+1}) vs Bonferroni (ii).
"""
from fractions import Fraction as F
from functools import lru_cache
import time

ONE = F(1)
SIXTH = F(1, 6)
SQRT2_UP = F(99, 70)          # > sqrt2
assert SQRT2_UP * SQRT2_UP > 2

@lru_cache(maxsize=4096)
def bad_pieces(w):
    r = F(1, 14 * w)
    out = []
    for k in range(w):
        c = F(k, w)
        lo, hi = c - r, c + r
        if lo < 0:
            out.append((F(0), hi)); out.append((lo + 1, ONE))
        elif hi > 1:
            out.append((lo, ONE)); out.append((F(0), hi - 1))
        else:
            out.append((lo, hi))
    return tuple(out)

def good_intervals(speeds):
    pieces = []
    for w in speeds:
        pieces += bad_pieces(w)
    pieces = sorted(pieces)
    comps = []
    for lo, hi in pieces:
        if comps and lo <= comps[-1][1]:
            if hi > comps[-1][1]:
                comps[-1][1] = hi
        else:
            comps.append([lo, hi])
    n = len(comps)
    goods = []
    for i in range(n):
        a = comps[i][1]
        j = (i + 1) % n
        b = comps[j][0] + (ONE if j == 0 else 0)
        if a < b:
            goods.append((a, b))
    return goods, len(goods), sum(b - a for a, b in goods)

def normalize(goods):
    out = []
    for a, b in goods:
        if b <= 1:
            out.append((a, b))
        else:
            out.append((a, ONE)); out.append((F(0), b - 1))
    out.sort()
    return out

def inter_measure(U1, U2):
    tot = F(0); i = j = 0
    while i < len(U1) and j < len(U2):
        a = max(U1[i][0], U2[j][0]); b = min(U1[i][1], U2[j][1])
        if a < b:
            tot += b - a
        if U1[i][1] < U2[j][1]:
            i += 1
        else:
            j += 1
    return tot

def B2(x):
    t = x - (x.numerator // x.denominator)
    return t * t - t + SIXTH

def disc_pair(goods, v):
    E = []
    for a, b in goods:
        E.append((a if a < 1 else a - 1, 1)); E.append((b if b < 1 else b - 1, -1))
    S = F(0)
    for e, s in E:
        for f, u in E:
            S += s * u * B2(v * (e - f))
    return S / (2 * v * v)

def eps_exact(GE_norm, m_E, v):
    return inter_measure(GE_norm, sorted(bad_pieces(v))) - m_E / 7

# ----------------------------------------------------------------------------------------
print("=" * 102)
print("THM-735 (i)+(ii) EXACT VERIFICATION -- simultaneous multi-peel, all links in Q")
print("=" * 102)

BATTERY = [
    (list(range(1, 11)),  [13, 22, 84],            "k=10 extremal body+peels"),
    (list(range(1, 12)),  [13, 84],                "residue, j=2"),
    (list(range(1, 13)),  [182],                   "deep well, j=1"),
    (list(range(1, 11)),  [150, 151, 152],         "clustered j=3, c=150"),
    (list(range(1, 11)),  [200, 201, 202],         "clustered j=3, c=200"),
    (list(range(1, 11)),  [215, 216, 217],         "clustered j=3, c=215"),
    (list(range(1, 11)),  [300, 301, 302],         "clustered j=3, c=300"),
    (list(range(1, 11)),  [500, 501, 502],         "clustered j=3, c=500"),
    (list(range(1, 8)),   [300, 301, 302, 303, 304, 305], "j=6 clustered, body {1..7}"),
]

all_ok = True
for E, Fset, name in BATTERY:
    j = len(Fset)
    assert j <= 6
    goodsE, rE, mE = good_intervals(E)
    GE = normalize(goodsE)
    epss = [eps_exact(GE, mE, v) for v in Fset]
    discs = [disc_pair(goodsE, v) for v in Fset]
    _, _, Ltrue = good_intervals(E + Fset)
    base = (1 - F(j, 7)) * mE
    lemma_ok = Ltrue >= base - sum(abs(e) for e in epss)
    cov_ok = all(e * e <= F(6, 49) * d for e, d in zip(epss, discs))
    crude_ok = all(d <= F(rE * rE, 3) / (v * v) for d, v in zip(discs, Fset))
    cert_crude = base - SQRT2_UP * F(rE, 7) * sum(F(1, v) for v in Fset)
    cert_exact = base * base - j * F(6, 49) * sum(discs)   # >0 => (CS) Sum sqrt((6/49)disc) < base
    ok = lemma_ok and cov_ok and crude_ok
    all_ok &= ok
    print("\n%-34s  E-size=%d j=%d   m_E=%.5f r_E=%d   L_true=%.6f" % (name, len(E), j, float(mE), rE, float(Ltrue)))
    print("   lemma L>=base-Sum|eps|: %s   per-peel eps^2<=(6/49)disc: %s   disc<=r^2/3v^2: %s"
          % (lemma_ok, cov_ok, crude_ok))
    print("   base=(1-j/7)m_E=%.6f   Sum|eps|=%.6f   bound value=%.6f (vs L_true %.6f)"
          % (float(base), float(sum(abs(e) for e in epss)),
             float(base - sum(abs(e) for e in epss)), float(Ltrue)))
    print("   CERT crude: %.6f -> %s     CERT exact-disc (CS): %.3e -> %s"
          % (float(cert_crude), "FIRES (L>0 proved)" if cert_crude > 0 else "no",
             float(cert_exact), "FIRES (L>0 proved)" if cert_exact > 0 else "no"))
assert all_ok
print("\n=> ALL chain links hold EXACTLY on the whole battery.")

# ----------------------------------------------------------------------------------------
print("\n" + "=" * 102)
print("SEQUENTIAL vs SIMULTANEOUS on clustered triples {1..10, c, c+1, c+2}")
print("  sequential crude: peel b=c+2 from exact 12-speed base {1..10,c,c+1}: fires iff b > (99/70)r/(6m)")
print("  simultaneous (ii): fires iff (4/7)m0 - (99/70)(r0/7)(1/c+1/(c+1)+1/(c+2)) > 0")
print("=" * 102)
E0 = list(range(1, 11))
goods0, r0, m0 = good_intervals(E0)
G0 = normalize(goods0)
print("body {1..10}: r=%d m=%s~%.5f" % (r0, m0, float(m0)))
first_sim = None
for c in [120, 150, 180, 200, 210, 215, 220, 250, 300, 400, 600, 1000]:
    Fset = [c, c + 1, c + 2]
    # simultaneous
    cert = (1 - F(3, 7)) * m0 - SQRT2_UP * F(r0, 7) * sum(F(1, v) for v in Fset)
    # sequential: exact base {1..10,c,c+1}
    gseq, rseq, mseq = good_intervals(E0 + [c, c + 1])
    seq_fires = F(c + 2) > SQRT2_UP * F(rseq, 1) / (6 * mseq)
    # exact truth
    _, _, Lt = good_intervals(E0 + Fset)
    if cert > 0 and first_sim is None:
        first_sim = c
    print("  c=%5d: simultaneous %-18s sequential(peel c+2 vs base r=%3d m=%.5f) %-9s  L_true=%.6f"
          % (c, "FIRES(%.5f)" % float(cert) if cert > 0 else "no(%.5f)" % float(cert),
             rseq, float(mseq), "fires" if seq_fires else "NO", float(Lt)))
print("\n  first simultaneous fire among samples: c=%s ; exact threshold V1 = %d (all three >= V1 closed)"
      % (first_sim, int(F(297 * r0, 280 * m0)) + 1))
print("done.")
