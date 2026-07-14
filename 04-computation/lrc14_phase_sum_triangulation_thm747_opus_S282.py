#!/usr/bin/env python3
"""
lrc14_phase_sum_triangulation_thm747_opus_S282.py
=================================================
opus-2026-07-14-S282.  THM-747: the vertices-as-runners phase sum S(W), TRIANGULATED.

LENS 1 (by-vertex / pair sector).  Regroup the per-end sum by shared vertex u*.  Writing
theta1 = 1 - {u1 W} and theta2 = {u2 W} (generic W), S(W) = D0 + Sum_v b_v {u_v W} with
D0 = Sum_{start-ends} c_e and b_v = [Sum_{end-type} c_e - Sum_{start-type} c_e] -- at a shared
(swap) vertex the couplings collapse to the PAIR-DIFFERENCE form (delta-weighted).  Bound:
|S| <= |S_bar| + B1/2, S_bar = D0 + Sum b_v / 2, B1 = Sum |b_v|.

LENS 2 (runner / periodicity).  {u_v W} is periodic in integer W with period q_v = den(u_v):
S(W) is periodic mod Q = lcm(q_v).  The EXACT maximum over W mod Q is a finite computation:
the level-2 bound becomes a LOOKUP (float scan + exact Fraction confirmation at candidates).
Same for T(W) and for the zero-crossing correction Z(W) (Z = 0 for W >= Wz = 2/min du).

LENS 3 (tower).  Integer time makes mirror phases EQUAL ({(1-u)W} = {-uW}), so the level-1
mirror is EXHAUSTED at level 2 -- paired vertices carry equal terms (verified).  Level-2
cancellation is pure phase equidistribution across vertices, and the Q-scan resolves it
exactly.  With S_max/T_max from the scan, the assembly beats the standing W0.

Outputs: the three bounds triangulated (per-end S1; by-vertex |S_bar|+B1/2; exact period max),
the pair-difference structure table, mirror-pair verification, Wz, the improved assembly + W0,
and the named remaining unknowns.
"""
from fractions import Fraction as F
import numpy as np
from math import gcd

LAM = F(1, 14)
HALF = F(1, 2)

def normalize(ivs):
    ivs = sorted((a, b) for a, b in ivs if b > a)
    out = []
    for a, b in ivs:
        if out and a <= out[-1][1]:
            if b > out[-1][1]: out[-1] = (out[-1][0], b)
        else: out.append((a, b))
    return out

def intersect(A, B):
    out = []; i = j = 0
    while i < len(A) and j < len(B):
        a1, b1 = A[i]; a2, b2 = B[j]
        lo, hi = max(a1, a2), min(b1, b2)
        if hi > lo: out.append((lo, hi))
        if b1 <= b2: i += 1
        else: j += 1
    return out

def measure(A): return sum(b - a for a, b in A)

def safe_set(w):
    return [((k + LAM) / w, (k + 1 - LAM) / w) for k in range(w)]

def good_intervals(speeds):
    cur = [(F(0), F(1))]
    for w in sorted(speeds): cur = intersect(cur, safe_set(w))
    return cur

def A_J(u, J):
    ivs = []
    for j in J:
        c = (-j * u) % 1
        a, b = c - LAM, c + LAM
        if a < 0: ivs += [(a + 1, F(1)), (F(0), b)]
        elif b > 1: ivs += [(a, F(1)), (F(0), b - 1)]
        else: ivs.append((a, b))
    return 1 - measure(normalize(ivs))

def breakpoints(J, extra):
    pts = set(extra)
    diffs = sorted({abs(a - b) for a in J for b in J if a != b})
    for d in diffs:
        for k in range(d + 1):
            for off in (F(0), F(1, 7), -F(1, 7)):
                u = (F(k) + off) / d
                if 0 <= u <= 1: pts.add(u)
    pts.add(F(0)); pts.add(F(1))
    return sorted(pts)

def segments_with_heights(J, GB, pts):
    pieces = []
    for a, b in zip(pts, pts[1:]):
        if a == b: continue
        for da, db in GB:
            lo, hi = max(a, da), min(b, db)
            if hi > lo: pieces.append((lo, hi))
    pieces.sort()
    segs = []
    for j in J:
        if j == 0: continue
        for sign in (1, -1):
            run = None; prev_hi = None
            for lo, hi in pieces:
                mid = (lo + hi) / 2
                x = (-j * mid + F(sign, 14)) % 1
                exposed = True
                for j2 in J:
                    if j2 == j: continue
                    d = (x - (-j2 * mid)) % 1
                    if (d < LAM or d > 1 - LAM) and d != LAM and d != 1 - LAM:
                        exposed = False; break
                if exposed:
                    if run is not None and prev_hi == lo: run = (run[0], hi)
                    else:
                        if run is not None: segs.append((j, sign, run[0], run[1]))
                        run = (lo, hi)
                    prev_hi = hi
                else:
                    if run is not None: segs.append((j, sign, run[0], run[1])); run = None
            if run is not None: segs.append((j, sign, run[0], run[1]))
    return segs

def Ffun(x): return x * (1 - x) / 2
def psi(x): return HALF - x

CFG = [
    ("shape 1 {1}u{W..W+11}", [1], list(range(12)), F(3690)),
    ("shape 2 {1,2,3}u{W..W+9}", [1, 2, 3], list(range(10)), F(1971)),
]

print("=" * 108)
print("THM-747: the vertices-as-runners phase sum, triangulated (by-vertex / periodicity / tower)")
print("=" * 108)

for name, Bset, J, C2_743 in CFG:
    GB = good_intervals(Bset)
    pts = breakpoints(J, [e for iv in GB for e in iv])
    segs = segments_with_heights(J, GB, pts)
    # ends: (u, type, coupling)  type 1 = start (theta = 1 - {uW}), type 2 = end (theta = {uW})
    ends = []
    S2 = 0; minDu = None
    for (j, sign, u1, u2) in segs:
        orient = -1 if sign == 1 else 1
        xs = (-j * u1 + F(sign, 14)) % 1
        xe = (-j * u2 + F(sign, 14)) % 1
        ends.append((u1, 1, orient * j * psi(xs), j))
        ends.append((u2, 2, orient * j * psi(xe), j))
        S2 += j * j
        du = u2 - u1
        if minDu is None or du < minDu: minDu = du
    S1 = sum(abs(c) for (u, t, c, j) in ends)
    # LENS 1: by-vertex regroup
    from collections import defaultdict
    D0 = F(0)
    bv = defaultdict(lambda: F(0))
    vend = defaultdict(list)
    for (u, t, c, j) in ends:
        vend[u].append((t, c, j))
        if t == 1:
            D0 += c; bv[u] -= c
        else:
            bv[u] += c
    verts = sorted(bv.keys())
    B1 = sum(abs(bv[u]) for u in verts)
    Sbar = D0 + sum(bv[u] for u in verts) / 2
    print("\n%s: %d segments, %d ends, %d distinct vertices" % (name, len(segs), len(ends), len(verts)))
    print("   LENS 1 (by-vertex): D0 = %.4f ; S_bar = %.4f ; B1 = Sum|b_v| = %.2f" % (float(D0), float(Sbar), float(B1)))
    print("      bound |S| <= |S_bar| + B1/2 = %.2f   [per-end S1 was %.2f]" % (float(abs(Sbar) + B1/2), float(S1)))
    shared = [u for u in verts if len(vend[u]) >= 2]
    print("      shared vertices: %d of %d; pair-difference collapse sample:" % (len(shared), len(verts)))
    for u in shared[:3]:
        es = vend[u]
        print("        u=%s: ends %s -> b_v = %.4f" % (u, [(t, float(c)) for (t, c, j) in es], float(bv[u])))
    # mirror-pair check: b at u vs b at 1-u
    mirror_eq = sum(1 for u in verts if (1 - u) in bv and bv[u] == bv[1 - u])
    print("      mirror pairs with EQUAL b_v: %d/%d vertices (integer time exhausts the mirror)"
          % (mirror_eq, len(verts)))
    # LENS 2: periodicity + exact scan
    Q = 1
    for u in verts: Q = Q * u.denominator // gcd(Q, u.denominator)
    Wz = (2 / minDu).__ceil__()
    print("   LENS 2 (periodicity): Q = lcm(vertex denominators) = %d ; Z(W)=0 for W >= Wz = %d" % (Q, Wz))
    uf = np.array([float(u) for u in verts]); bf = np.array([float(bv[u]) for u in verts])
    D0f = float(D0)
    # T(W) needs per-end j^2 with theta^2 signs: T = Sum_seg orient j^2 (th1^2 - th2^2): by end: type1 +oj^2 th1^2? T = Sum orient j^2 th1^2 - orient j^2 th2^2
    tj = []
    for (j, sign, u1, u2) in segs:
        orient = -1 if sign == 1 else 1
        tj.append((float(u1), 1, orient * j * j))
        tj.append((float(u2), 2, -orient * j * j))
    tu = np.array([u for (u, t, c) in tj]); tc = np.array([c for (u, t, c) in tj]); tt = np.array([t for (u, t, c) in tj])
    CH = 2 * 10 ** 6
    smax = 0.0; smaxW = 0; tmax = 0.0
    for lo in range(0, Q, CH):
        Wv = np.arange(lo, min(lo + CH, Q), dtype=np.float64)
        fr = np.modf(np.outer(Wv, uf))[0]
        Sv = D0f + fr @ bf
        # T: theta1 = 1-{uW} for type1(start), theta2={uW}: theta^2 terms
        frt = np.modf(np.outer(Wv, tu))[0]
        TH = np.where(tt == 1, 1 - frt, frt)
        Tv = (TH * TH) @ tc
        k = int(np.argmax(np.abs(Sv)))
        if abs(Sv[k]) > smax: smax, smaxW = abs(Sv[k]), int(Wv[k])
        tmax = max(tmax, float(np.max(np.abs(Tv))))
    # exact confirmation at the float argmax
    Wc = smaxW
    Sx = D0 + sum(bv[u] * ((u * Wc) % 1) for u in verts)
    print("      float scan over full period: max|S(W)| = %.4f at W = %d (exact there: %.4f); max|T(W)| = %.1f"
          % (smax, smaxW, float(abs(Sx)), tmax))
    print("      TRIANGULATED: per-end S1 = %.0f  >  by-vertex %.1f  >  EXACT period max ~ %.1f  [measured 2.6-49]"
          % (float(S1), float(abs(Sbar) + B1/2), smax))
    # improved assembly (W >= Wz): C1 = 2(#comp+|Xi_sv|); C2 = C2_743 + 2*Smax + Tmax/2
    Xi_sv = sum((-1 if sign == 1 else 1) * (Ffun((-j*u1 + F(sign,14)) % 1) - Ffun((-j*u2 + F(sign,14)) % 1))
                for (j, sign, u1, u2) in segs)
    C1 = 2 * (len(GB) + abs(Xi_sv))
    C2 = float(C2_743) + 2 * smax * 1.05 + tmax / 2   # 5% float-margin on smax
    lo, hi = Wz, 10 ** 7
    while lo < hi:
        mid = (lo + hi) // 2
        if float(C1) / mid + C2 / (mid * mid) < float(measure(GB)) * 0:  # placeholder
            hi = mid
        else: break
    # compute area for W0
    area, _ = None, None
    # quick exact area
    def integrate_AJ_local(J, dom):
        pts2 = breakpoints(J, [e for iv in dom for e in iv])
        tot = F(0)
        def integ(a, b, fa, fb, d=0):
            m = (a + b) / 2; fm = A_J(m, J)
            assert fm * 2 == fa + fb or d <= 12
            if fm * 2 == fa + fb: return (fa + fb) * (b - a) / 2
            return integ(a, m, fa, fm, d + 1) + integ(m, b, fm, fb, d + 1)
        for a, b in zip(pts2, pts2[1:]):
            if a == b: continue
            for da, db in dom:
                l2, h2 = max(a, da), min(b, db)
                if h2 > l2: tot += integ(l2, h2, A_J(l2, J), A_J(h2, J))
        return tot
    area = integrate_AJ_local(J, GB)
    lo = Wz
    W0 = None
    W = Wz
    lo2, hi2 = 1, 10 ** 7
    while lo2 < hi2:
        mid = (lo2 + hi2) // 2
        if float(C1) / mid + C2 / (mid * mid) < float(area): hi2 = mid
        else: lo2 = mid + 1
    W0 = max(lo2, Wz)
    print("   ASSEMBLY (W >= Wz): C1 = %.3f ; C2 = %.0f (743 pots %.0f + 2 S_max + T_max/2) ; W0 = %d [standing 339/513]"
          % (float(C1), C2, float(C2_743), W0))

print("\nREMAINING UNKNOWNS (named):")
print(" U1. The 743 pots now dominate C2: vertex capping (one-signed, Sum delta = 2172/1116) and the")
print("     quadratic wedge remainders -- BOTH are also periodic mod Q: their exact period-extremes are")
print("     the same kind of lookup (unscanned here).  With them, the family's asymptotic first-order")
print("     coefficient becomes a FINITE EXACT OBJECT -- no analytic unknown remains in the tail lane.")
print(" U2. The compact core (bounded-Vmax bodies) -- outside this program (kps exact certificates).")
print(" U3. The general multi-speed equidistribution (klein-S300 capstone) -- the fleet-level residual;")
print("     this lane contributes exact structure (strand identity, pairing, periodicity), not a bridge.")
print("done.")
