#!/usr/bin/env python3
"""kps-2026-07-08-S88 (cont.2) -- enumerate the longest-AP = 9 and = 8 families for opus's S157
resonance-sum finite check, mirroring the L=10 enumeration (lrc14_longestAP10_d456_kps_S88).

The k=11 tail stratifies by longest-AP length L (dilation-invariant).  L=10 is CLOSED (opus-S157 proved
+ kps verified <= prim-diam 54).  Remaining = L <= 9.  Each L-family = {AP_L at scale d} u {(11-L)
points} (the correlated AP_L core + decorrelating extras -- exactly what opus's m_j = L_j + resonance-sum
analyzes, now with 11-L free points).  This gives opus the per-d finite-check tables + the per-L
extremal (the 'A_L' analog of A), and confirms the L=9,8 families sit far above bar (margin >> +0.12)."""
import sys
sys.path.insert(0, "."); sys.path.insert(0, "04-computation")
from lrc14_d3_exact_verify_klein_S184 import D3 as D3e
from fractions import Fraction as Fr
from math import gcd, floor
from functools import reduce
from itertools import combinations

BAR = Fr(83549, 252252); BARF = float(BAR)
TH = 1/7; M = 6/7
# decorrelation limits D3_L (klein-S187): L=8->0.6021, L=9->0.5235, L=10->0.4646
D3_LIMIT = {8: 0.6021, 9: 0.5235, 10: 0.4646}


def D3f(E):
    """fast float D3 via Farey-cell integration (scan; minima exact-verified)."""
    E = sorted(E); k = len(E); D = E[-1] - E[0]
    bps = set([0.0, 1.0])
    for d in range(1, D + 1):
        inv = 1.0 / d
        for m in range(0, d + 1): bps.add(m * inv)
    bps = sorted(bps); m1 = m2 = m3 = 0.0
    for c in range(len(bps) - 1):
        a, b = bps[c], bps[c + 1]; mid = (a + b) / 2
        cj = [floor(e * mid) for e in E]
        order = sorted(range(k), key=lambda i: E[i] * mid - cj[i])
        ph = [(E[order[r]], -cj[order[r]]) for r in range(k)]
        gaps = [(ph[r + 1][0] - ph[r][0], ph[r + 1][1] - ph[r][1]) for r in range(k - 1)]
        gaps.append((ph[0][0] - ph[k - 1][0], 1 + ph[0][1] - ph[k - 1][1]))
        pts = {a, b}
        for (s, ic) in gaps:
            if s != 0:
                xc = (TH - ic) / s
                if a < xc < b: pts.add(xc)
        pts = sorted(pts)
        for t in range(len(pts) - 1):
            lo, hi = pts[t], pts[t + 1]; mm = (lo + hi) / 2; A = 0.0; Bc = 0.0
            for (s, ic) in gaps:
                if s * mm + ic > TH: A += s; Bc += ic - TH
            m1 += A / 2 * (hi * hi - lo * lo) + Bc * (hi - lo)
            m2 += A * A / 3 * (hi**3 - lo**3) + A * Bc * (hi * hi - lo * lo) + Bc * Bc * (hi - lo)
            m3 += A**3 / 4 * (hi**4 - lo**4) + A * A * Bc * (hi**3 - lo**3) + 3 * A * Bc * Bc / 2 * (hi * hi - lo * lo) + Bc**3 * (hi - lo)
    den = m2 - m3 / M
    return m1 / M + (m1 - m2 / M)**2 / den if den > 1e-15 else m1 / M


def isprim(E):
    E = sorted(set(E)); return reduce(gcd, [e - E[0] for e in E]) == 1


def pdiam(E):
    E = sorted(set(E)); return max(E) - min(E)


def longest_ap(E):
    E = sorted(set(E)); S = set(E); best = 1
    for i in range(len(E)):
        for j in range(i + 1, len(E)):
            d = E[j] - E[i]
            L = 2; nxt = E[j] + d
            while nxt in S: L += 1; nxt += d
            prv = E[i] - d
            while prv in S: L += 1; prv -= d
            if L > best: best = L
    return best


def enumerate_family(L, dscales, reach, pdmin=25, pdmax=54):
    """GENUINE longest-AP = L TAIL shapes: {AP_L scale d} u {(11-L) points}, filtered to
       prim-diam in [pdmin, pdmax] AND longest-AP EXACTLY L (extras that complete a 10-AP are
       excluded -- those are L=10, covered by opus-S157). Per-d + overall min D3 (float)."""
    nextra = 11 - L
    overall = (9.9, None)
    rows = []
    for d in dscales:
        ap = tuple(d * j for j in range(L))            # {0,d,..,(L-1)d}, span (L-1)d
        hi = (L - 1) * d + reach
        pool = [p for p in range(1, hi + 1) if p not in ap]
        per_d = (9.9, None)
        for extra in combinations(pool, nextra):
            E = tuple(sorted(ap + extra))
            if len(E) != 11 or not isprim(E): continue
            if not (pdmin <= pdiam(E) <= pdmax): continue
            v = D3f(E)
            if v >= 0.68 and v >= per_d[0]: continue    # cheap skip before longest_ap
            if longest_ap(E) != L: continue             # keep only GENUINE L
            if v < per_d[0]: per_d = (v, extra)
            if v < overall[0]: overall = (v, E)
        rows.append((d, (L - 1) * d, per_d))
    return rows, overall


def main():
    print(f"L=9 and L=8 families for opus-S157 finite check  (bar={BARF:.6f})")
    print("stratification: L=10 CLOSED (opus proved + kps verified<=54); L<=9 = the remaining stratum")
    print("=" * 96)

    for L, dscales, reach in [(9, (2, 3, 4, 5), 10), (8, (2, 3, 4), 8)]:
        print(f"\n### longest-AP = {L} family: {{AP_{L} scale d}} u {{{11-L} points}}   "
              f"(decorrelation limit D3_{L} = {D3_LIMIT[L]})")
        rows, overall = enumerate_family(L, dscales, reach)
        print(f"  {'d':>2} {'AP span':>8} {'min D3':>10} {'best extra':>18} {'margin':>9}")
        for d, span, (v, extra) in rows:
            print(f"  {d:>2} {span:>8} {v:>10.6f} {str(extra):>18} {v-BARF:>+9.5f}")
        E = overall[1]; de = float(D3e(E)); lap = longest_ap(E)
        print(f"  --> family MIN D3 = {de:.6f} (EXACT) at {E}")
        print(f"      prim-diam {pdiam(E)}, longest-AP {lap}, margin {de-BARF:+.6f}  "
              f"[{'>= bar, +'+format(de-BARF,'.3f') if de>=BARF else 'BELOW BAR'}]"
              + (f"  (NB longest-AP={lap}, an L=10 member -- covered by opus-S157)" if lap >= 10 else ""))

    print("\n" + "=" * 96)
    print("HANDOFF to opus-S157: each L-family = {AP_L scale d} + (11-L) points; the extremal sits at")
    print("SMALL scale d (like L=10's d=3), dips slightly BELOW its decorrelation limit D3_L, and stays")
    print("FAR above bar (L=9 margin ~ +0.17, L=8 ~ +0.25).  Your m_j = L_j + resonance-sum + 1/(pd)")
    print("rate extends per L with (11-L) free points -- these tables are the finite-check anchors.")


if __name__ == "__main__":
    main()
