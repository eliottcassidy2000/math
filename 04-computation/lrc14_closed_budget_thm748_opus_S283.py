#!/usr/bin/env python3
"""
lrc14_closed_budget_thm748_opus_S283.py
=======================================
opus-2026-07-14-S283.  THM-748: THE CLOSED BUDGET of the tail lane.

CLAIM.  For W >= Wz (every segment has a full-inside crossing), W (L - Area) decomposes EXACTLY:

  W (L - Area)  =  PHI(W)  +  QPOT(W)  +  KAP(W),

  PHI  = Sum_seg orient (j/W) Sum_{full-inside crossings} psi(h_m)          [linear wedge]
  QPOT = Sum_seg orient (j^2/(W(W+j))) Sum_{full crossings} h_m             [exact quadratic]
         (per-crossing f - g = (j/W) psi(h) + h j^2/(W(W+j)) EXACTLY: the strand crossing is
          sigma* = hW/(W+j), the strip average h - j/2W; no remainder, F'' = -1 not needed here)
  KAP  = Sum over END-EVENTS of exact closed-form kappa(phi), phi = the event's grid phase:
    DEATH (upper path j_up, lower j_lo meet at height x; delta = j_up - j_lo > 0):
        kappa_D = W delta (phi - x)_+ / ((W+j_up)(W+j_lo))  -  delta phi^2 / (2W)
    BIRTH (delta = j_lo - j_up > 0):
        kappa_B = W delta (x - phi)_+ / ((W+j_up)(W+j_lo))  -  delta (1-phi)^2 / (2W)
    SWAP (path kinks from slope -j_A to -j_B at height x, orientation o):
        kappa_S = o [ sigma_valid - x - (j_A phi^2 - j_B (1-phi)^2)/(2W) ],
        sigma_valid = (xW + j_A phi)/(W+j_A) if x <= phi else (xW + j_B phi)/(W+j_B)
    CUT (right G_B endpoint u_G; per exposed run with paths (x_up, j_up), (x_lo, j_lo)):
        g = (x_up - x_lo) phi + (j_up - j_lo) phi^2/(2W)
        f = len( [max(sig_lo, 0), min(sig_up, phi)] ), sig_p = (x_p W + j_p phi)/(W + j_p)
        kappa_C = f - g ;  LEFT G_B endpoint: the mirror formulas.
  All exact rationals given phi = the fractional position of the event in its strip.

REFEREE (the gate): TOTAL_pred(W) vs W (L_engine - Area): EXACT Fraction equality demanded at a
battery of W >= Wz.  On pass: the float scan of the total envelope E(W) over a full period
[Wz, Wz + Q) closes the lane's budget: |L - Area| <= max|E| / W for W in the window, extended
beyond by the explicit 1/W-monotone parts.
"""
from fractions import Fraction as F

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

def integrate_AJ(J, dom):
    pts = breakpoints(J, [e for iv in dom for e in iv])
    tot = F(0)
    def integ(a, b, fa, fb, d=0):
        m = (a + b) / 2; fm = A_J(m, J)
        assert fm * 2 == fa + fb or d <= 12
        if fm * 2 == fa + fb: return (fa + fb) * (b - a) / 2
        return integ(a, m, fa, fm, d + 1) + integ(m, b, fm, fb, d + 1)
    for a, b in zip(pts, pts[1:]):
        if a == b: continue
        for da, db in dom:
            lo, hi = max(a, da), min(b, db)
            if hi > lo: tot += integ(lo, hi, A_J(lo, J), A_J(hi, J))
    return tot, pts

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

def build_events(segs, GB):
    """classify segment-end events.  Returns (events, per-seg data).
    event kinds: ('death', j_up, j_lo, x, u), ('birth', ...), ('swap', o, j_A, j_B, x, u),
    ('cutR', runs, u_G), ('cutL', runs, u_G) with runs = [(x_lo, j_lo, x_up, j_up), ...]."""
    from collections import defaultdict
    GBends = set()
    for a, b in GB: GBends.add(a); GBends.add(b)
    endlist = defaultdict(list)   # (u, x) -> [(side, orient, j)]
    cutlist = defaultdict(list)   # u_G -> [(side, orient, j, x)]
    for (j, sign, u1, u2) in segs:
        orient = -1 if sign == 1 else 1
        xs = (-j * u1 + F(sign, 14)) % 1
        xe = (-j * u2 + F(sign, 14)) % 1
        if u1 in GBends: cutlist[u1].append(('start', orient, j, xs))
        else: endlist[(u1, xs)].append(('start', orient, j))
        if u2 in GBends: cutlist[u2].append(('end', orient, j, xe))
        else: endlist[(u2, xe)].append(('end', orient, j))
    events = []
    unmatched = 0
    for (u, x), mem in endlist.items():
        ends = [(o, j) for (s, o, j) in mem if s == 'end']
        starts = [(o, j) for (s, o, j) in mem if s == 'start']
        # pair swaps: end+start with same orient
        used_e = [False]*len(ends); used_s = [False]*len(starts)
        for ie, (oe, je) in enumerate(ends):
            for isx, (os_, js) in enumerate(starts):
                if used_e[ie] or used_s[isx]: continue
                if oe == os_:
                    events.append(('swap', oe, je, js, x, u))
                    used_e[ie] = True; used_s[isx] = True
        rem_e = [(o, j) for k, (o, j) in enumerate(ends) if not used_e[k]]
        rem_s = [(o, j) for k, (o, j) in enumerate(starts) if not used_s[k]]
        # deaths: pair rem_e with opposite orients
        while len(rem_e) >= 2:
            (o1, j1) = rem_e.pop()
            mate = next((k for k, (o2, j2) in enumerate(rem_e) if o2 == -o1), None)
            if mate is None: rem_e.insert(0, (o1, j1)); break
            (o2, j2) = rem_e.pop(mate)
            jup = j1 if o1 == 1 else j2
            jlo = j2 if o1 == 1 else j1
            events.append(('death', jup, jlo, x, u))
        while len(rem_s) >= 2:
            (o1, j1) = rem_s.pop()
            mate = next((k for k, (o2, j2) in enumerate(rem_s) if o2 == -o1), None)
            if mate is None: rem_s.insert(0, (o1, j1)); break
            (o2, j2) = rem_s.pop(mate)
            jup = j1 if o1 == 1 else j2
            jlo = j2 if o1 == 1 else j1
            events.append(('birth', jup, jlo, x, u))
        # STATIC-PARTNER events (the arc-0 edges at heights 1/14 and 13/14 -- the origin band):
        #   (end, +1): moving upper drops onto the static floor -> death (j_up=j, j_lo=0)
        #   (start,-1): moving lower drops off the static ceiling -> birth (j_up=0, j_lo=j)
        #   (start,+1): moving upper emerges below the ceiling -> swap static->moving (jA=0, jB=j)
        #   (end, -1): moving lower dies into the floor -> swap moving->static (jA=j, jB=0)
        for (o1, j1) in rem_e:
            if o1 == 1: events.append(('death', j1, 0, x, u))
            else: events.append(('swap', -1, j1, 0, x, u))
        for (o1, j1) in rem_s:
            if o1 == -1: events.append(('birth', 0, j1, x, u))
            else: events.append(('swap', 1, 0, j1, x, u))
    for uG, mem in cutlist.items():
        # reconstruct runs at the cut: sort by height; pair (-1 lower, +1 upper) ascending
        side = 'end' if any(s == 'end' for (s, o, j, x) in mem) else 'start'
        pts_ = sorted([(x, o, j) for (s, o, j, x) in mem])
        # pad with the STATIC arc-0 boundaries: a leading +1 needs the static floor below;
        # a trailing -1 needs the static ceiling above
        if pts_ and pts_[0][1] == 1: pts_ = [(F(1, 14), -1, 0)] + pts_
        if pts_ and pts_[-1][1] == -1: pts_ = pts_ + [(F(13, 14), 1, 0)]
        runs = []
        k = 0
        ok = True
        while k < len(pts_):
            if k + 1 < len(pts_) and pts_[k][1] == -1 and pts_[k+1][1] == 1:
                runs.append((pts_[k][0], pts_[k][2], pts_[k+1][0], pts_[k+1][2]))
                k += 2
            else:
                ok = False; k += 1
        if not ok: unmatched += 1
        events.append(('cutR' if side == 'end' else 'cutL', runs, uG))
    return events, unmatched

def total_pred(segs, events, W, exact=True):
    """PHI + QPOT + KAP at integer W (exact Fractions if exact else floats)."""
    conv = (lambda z: z) if exact else float
    PHI = 0; QPOT = 0; KAP = 0
    for (j, sign, u1, u2) in segs:
        orient = -1 if sign == 1 else 1
        m1 = (u1 * W).__ceil__(); m2 = (u2 * W).__floor__() - 1
        if m2 < m1: continue
        K1 = m2 - m1 + 1
        alpha = F(j, W)
        hf = (F(sign, 14) - F(j * m1, W)) % 1
        # no wraps (strong no-wrap lemma): arithmetic sums
        sum_h = K1 * hf - alpha * (K1 - 1) * K1 / 2
        sum_psi = K1 * HALF - sum_h
        PHI += conv(orient * alpha * sum_psi)
        QPOT += conv(orient * F(j * j, W * (W + j)) * sum_h)
    for ev in events:
        if ev[0] in ('death', 'birth'):
            _, jup, jlo, x, u = ev
            phi = (u * W) % 1
            if phi == 0: continue
            if ev[0] == 'death':
                delta = jup - jlo
                fpart = W * delta * max(phi - x, F(0)) / ((W + jup) * (W + jlo))
                gpart = delta * phi * phi / (2 * W)
            else:
                delta = jlo - jup
                fpart = W * delta * max(x - phi, F(0)) / ((W + jup) * (W + jlo))
                gpart = delta * (1 - phi) * (1 - phi) / (2 * W)
            KAP += conv(fpart - gpart)
        elif ev[0] == 'swap':
            _, o, jA, jB, x, u = ev
            phi = (u * W) % 1
            if phi == 0: continue
            if x <= phi: sig = (x * W + jA * phi) / (W + jA)
            else: sig = (x * W + jB * phi) / (W + jB)
            KAP += conv(o * (sig - x - (jA * phi * phi - jB * (1 - phi) * (1 - phi)) / (2 * W)))
        elif ev[0] in ('cutR', 'cutL'):
            _, runs, uG = ev
            phi = (uG * W) % 1
            if phi == 0: continue
            for (xlo, jlo, xup, jup) in runs:
                if ev[0] == 'cutR':
                    g = (xup - xlo) * phi + (jup - jlo) * phi * phi / (2 * W)
                    slo = (xlo * W + jlo * phi) / (W + jlo)
                    sup = (xup * W + jup * phi) / (W + jup)
                    f = max(min(sup, phi) - max(slo, F(0)), F(0))
                else:
                    # left cut: run exists for u > u_G within the strip
                    g = (xup - xlo) * (1 - phi) - (jup - jlo) * (1 - phi) * (1 - phi) / (2 * W)
                    # paths: p(u) = x_p - j_p (u - u_G): strand sigma = x_p - j_p(sigma-phi)/W
                    slo = (xlo * W + jlo * phi) / (W + jlo)
                    sup = (xup * W + jup * phi) / (W + jup)
                    f = max(min(sup, F(1)) - max(slo, phi), F(0))
                KAP += conv(f - g)
    return PHI, QPOT, KAP

# ================= run =================

CFG = [
    ("shape 2 {1,2,3}u{W..W+9}", [1, 2, 3], list(range(10)), 588, [600, 601, 700, 977]),
    ("shape 1 {1}u{W..W+11}", [1], list(range(12)), 924, [930, 1001]),
]

print("=" * 106)
print("THM-748: the closed budget -- exact decomposition W(L-Area) = PHI + QPOT + KAP, referee-gated")
print("=" * 106)

for name, Bset, J, Wz, battery in CFG:
    GB = good_intervals(Bset)
    area, pts = integrate_AJ(J, GB)
    segs = segments_with_heights(J, GB, pts)
    events, unmatched = build_events(segs, GB)
    from collections import Counter
    kinds = Counter(ev[0] for ev in events)
    print("\n%s: %d segments; events: %s; unmatched ends: %d"
          % (name, len(segs), dict(kinds), unmatched))
    print("   REFEREE (exact Fraction equality vs engine):")
    for W in battery:
        body = list(Bset) + [W + j for j in J]
        Lx = measure(good_intervals(body))
        R_true = W * (Lx - area)
        PHI, QPOT, KAP = total_pred(segs, events, W, exact=True)
        pred = PHI + QPOT + KAP
        match = (pred == R_true)
        diff = float(R_true - pred)
        print("     W=%5d: R_true = %+.8f ; pred = %+.8f (PHI %+.6f QPOT %+.6f KAP %+.6f) ; EXACT: %s (diff %.2e)"
              % (W, float(R_true), float(pred), float(PHI), float(QPOT), float(KAP), match, diff))

print("\n(referee results determine whether the scan proceeds -- see .out continuation)")
