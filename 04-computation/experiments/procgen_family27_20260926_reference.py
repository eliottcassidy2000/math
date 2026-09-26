#!/usr/bin/env python3
"""procgen_family27_20260926_reference.py

Independent brute-force reference for procgen_family27_20260926_scan.c at small X.
Every orbit is iterated in full (to 1) with plain Python integers; no memo, no closed forms
for even numbers, no heap for delay records.  It prints the same report format as the C
scanner (REC, RHO2, BLK, HFULL, HWIN, HDT, HDS, HGL, HLAB, HB27, C27, EXL lines), so the two
reports can be compared line by line (EXL float sums approximately).

Usage: python3 -u procgen_family27_20260926_reference.py X
"""
import math
import sys


def T(x):
    return x // 2 if x % 2 == 0 else (3 * x + 1) // 2


def orbit(n):
    o = [n]
    while o[-1] != 1:
        o.append(T(o[-1]))
    return o


def build_roots():
    roots = []
    x = 27
    while True:
        if x not in roots:
            roots.append(x)
        if x == 1:
            break
        x = T(x)
    p = 2051
    for _ in range(6):
        p = 4 * p + 1
        if p < 2 ** 24 and p not in roots:
            roots.append(p)
    tr = 1
    for _ in range(12):
        if tr < 2 ** 24 and tr not in roots:
            roots.append(tr)
        tr = 4 * tr + 1
    p = 27
    for _ in range(9):
        p = 4 * p + 1
        if p < 2 ** 24 and p not in roots:
            roots.append(p)
    return roots


def riser_index(t, n, cap=130):
    i = 0
    while i + 1 < cap:
        j = i + 1
        ok = (t >= n * 2 ** (j // 2)) if j % 2 == 0 else (t * t >= n * n * 2 ** j)
        if not ok:
            break
        i = j
    return i


def main():
    X = int(sys.argv[1])
    roots = build_roots()
    ridx = {r: i for i, r in enumerate(roots)}
    reach = [[0] * len(roots) for _ in roots]
    for i, r in enumerate(roots):
        for v in orbit(r):
            if v in ridx:
                reach[i][ridx[v]] = 1
    inB = [reach[i][ridx[3077]] for i in range(len(roots))]
    out = []
    out.append("ROOTS %d" % len(roots))
    for i, r in enumerate(roots):
        out.append("ROOT %d %d %d" % (i, r, inB[i]))
    p = math.log(2) / math.log(3)
    sig = math.sqrt(p * (1 - p)) * math.log(3)
    NB, NRB, NRI = 34, 512, 130
    H_all = [0] * NB
    HF = [[0] * NRI for _ in range(NB)]
    HW = [[0] * NRI for _ in range(NB)]
    HDT = [[0] * NRB for _ in range(NB)]
    HDS = [[0] * NRB for _ in range(NB)]
    HGL = [[0] * NRB for _ in range(NB)]
    HL = [[0] * len(roots) for _ in range(NB)]
    HB = [[0] * 3 for _ in range(NB)]
    HDTB = [[0] * NRB for _ in range(NB)]
    HGLEN = [dict() for _ in range(NB)]
    C27 = [[0] * 5 for _ in range(NB)]
    EX = {}
    recs = []
    rho2 = []
    best = {"tT": -1, "tS": -1, "dT": -1, "dS": -1, "gT": -1, "gS": -1,
            "gam": -1.0, "glr": -1.0, "rho": -1.0}
    LN27, LN4616 = math.log(27), math.log(4616)
    for n in range(1, X + 1):
        o = orbit(n)
        dT = len(o) - 1
        ones = sum(1 for v in o[:-1] if v % 2 == 1)
        tT = max(o)
        odd_out = [3 * v + 1 for v in o[:-1] if v % 2 == 1]
        tS = max([n] + odd_out)
        w = n.bit_length() - 1
        # window: T^j(n), j <= w, continuing through the 1 -> 2 cycle if needed
        x, wm = n, n
        for _ in range(w):
            x = T(x)
            wm = max(wm, x)
        lab = next(ridx[v] for v in o if v in ridx)
        if n >= 2:
            g = next(j for j in range(1, len(o)) if o[j] < n)
            gC = sum(1 for v in o[:g] if v % 2 == 1)
            gmax = max(o[:g + 1])
        else:
            g, gC, gmax = 0, 0, 1
        # records in increasing n
        for key, val in (("tT", tT), ("tS", tS), ("dT", dT), ("dS", dT + ones)):
            if val > best[key]:
                best[key] = val
                recs.append((key, n, val))
        if n >= 2:
            for key, val in (("gT", g), ("gS", g + gC)):
                if val > best[key]:
                    best[key] = val
                    recs.append((key, n, val))
        if n >= 3 and n % 2 == 1:
            lnm = math.log(n)
            for key, val in (("gam", dT / lnm), ("glr", g / (lnm / math.log(2))), ("rho", math.log(tT) / lnm)):
                if val > best[key]:
                    best[key] = val
                    recs.append((key, n, val))
            if g >= 20:
                u = math.log(gmax / n) / (sig * math.sqrt(g))
                L = min(g, 1200)
                c, s1, s2 = EX.get(L, (0, 0.0, 0.0))
                EX[L] = (c + 1, s1 + u, s2 + u * u)
        b = w
        H_all[b] += 1
        HL[b][lab] += 1
        if inB[lab]:
            HB[b][n % 3] += 1
        HF[b][riser_index(tT, n)] += 1
        HW[b][riser_index(wm, n)] += 1
        if tT * 27 >= 4616 * n:
            C27[b][0] += 1
        if wm * 27 >= 4616 * n:
            C27[b][1] += 1
        if n >= 2:
            lnn = math.log(n)
            HDT[b][min(NRB - 1, int(4.0 * dT / lnn))] += 1
            if inB[lab]:
                HDTB[b][min(NRB - 1, int(4.0 * dT / lnn))] += 1
            HGLEN[b][min(g, 1023)] = HGLEN[b].get(min(g, 1023), 0) + 1
            HDS[b][min(NRB - 1, int(4.0 * (dT + ones) / lnn))] += 1
            HGL[b][min(NRB - 1, int(4.0 * g / (lnn / math.log(2))))] += 1
            if dT * LN27 >= 70.0 * lnn - 1e-9:
                C27[b][2] += 1
            if g * LN27 >= 59.0 * lnn - 1e-9:
                C27[b][3] += 1
            if math.log(tT) * LN27 >= LN4616 * lnn - 1e-12:
                C27[b][4] += 1
            if tT > n * n:
                rho2.append((n, tT, lab))
    for key, n, val in recs:
        if isinstance(val, float):
            out.append("REC %s %d %.9f" % (key, n, val))
        else:
            out.append("REC %s %d %d" % (key, n, val))
    for n, t, lab in rho2:
        out.append("RHO2 %d %d %d" % (n, t, lab))
    out.append("X %d" % X)
    for b in range(NB):
        if not H_all[b]:
            continue
        out.append("BLK %d %d" % (b, H_all[b]))
        out.append("HFULL %d " % b + " ".join(map(str, HF[b])))
        out.append("HWIN %d " % b + " ".join(map(str, HW[b])))
        out.append("HDT %d " % b + " ".join(map(str, HDT[b])))
        out.append("HDS %d " % b + " ".join(map(str, HDS[b])))
        out.append("HGL %d " % b + " ".join(map(str, HGL[b])))
        out.append("HLAB %d " % b + " ".join(map(str, HL[b])))
        out.append("HDTB %d " % b + " ".join(map(str, HDTB[b])))
        for L in sorted(HGLEN[b]):
            out.append("HGLEN %d %d %d" % (b, L, HGLEN[b][L]))
        out.append("HB27 %d %d %d %d" % (b, HB[b][0], HB[b][1], HB[b][2]))
        out.append("C27 %d " % b + " ".join(map(str, C27[b])))
    for L in sorted(EX):
        c, s1, s2 = EX[L]
        out.append("EXL %d %d %.10e %.10e" % (L, c, s1, s2))
    print("\n".join(out))


if __name__ == "__main__":
    main()
