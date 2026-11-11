# Auditor numerics-prereg: independent re-derivation of the A-PRIORI window prediction from the text of
# PREREGISTRATION.md alone (written before reading the predictor of odd_modulus_census.py).
#   W(h1, s) = int_{1/2}^{1} 2 d rho (1/2) max(0, B/f - max(A/f, y1)),  f = x1 rho + 2^s y1
# closed form: W = (1/x1) [ (B-A) ln(min(f1, A/y1)/f0)_+  +  B ln(u1/u0) - y1 (u1 - u0) ],
#   f0 = x1/2 + 2^s y1, f1 = x1 + 2^s y1, u0 = max(f0, A/y1), u1 = min(f1, B/y1).
# F_d = W-share of the smaller states of type d; Phi(g0) = sum_{d | g0} d F_d.
# usage: python3 -B my_predictor.py KHI KLO [registered_json] [stored_window_out]
import sys
import json
import math
import re
from aud_common import psi_table, divisors

MODS = [9, 27, 81, 25, 125, 49, 45, 75]


def W(x1, y1, s, A, B):
    b = (1 << s) * y1
    f0 = x1 / 2.0 + b
    f1 = x1 + b
    w = 0.0
    cut = A / y1
    if f0 < cut:
        w += (B - A) * math.log(min(f1, cut) / f0)
    u0 = max(f0, cut)
    u1 = min(f1, B / y1)
    if u1 > u0:
        w += B * math.log(u1 / u0) - y1 * (u1 - u0)
    return w / x1


def main():
    KHI = int(sys.argv[1])
    KLO = int(sys.argv[2])
    A = float(1 << KLO)
    B = float(1 << KHI)
    Hmax = int(math.isqrt(1 << KHI)) + 10
    psi = psi_table(Hmax)
    tot = [0.0, 0.0]
    byd = [dict(), dict()]  # sheet -> {g: W-mass}
    nstates = [0, 0]
    s = 0
    while (1 << s) < B:
        sh = 0 if s == 0 else 1
        for h1 in range(0, Hmax):
            x1 = 2 * int(psi[h1]) + 1
            y1 = 2 * h1 + 1
            if y1 * (y1 / 2.0 + (1 << s) * y1) > B:
                break
            w = W(x1, y1, s, A, B)
            if w > 0:
                g = math.gcd(x1, y1)
                tot[sh] += w
                byd[sh][g] = byd[sh].get(g, 0.0) + w
                nstates[sh] += 1
        s += 1
    print("a-priori prediction for the window [2^%d, 2^%d): predicted ordered triples sheet 0: %.0f, sheet 1: %.0f (smaller states used: %d, %d)" % (
        KLO, KHI, 2 * tot[0], 2 * tot[1], nstates[0], nstates[1]))
    phi = {}
    F = {}
    for sh in range(2):
        for m in MODS:
            Fd = {d: 0.0 for d in divisors(m)}
            for g, w in byd[sh].items():
                Fd[math.gcd(g, m)] += w / tot[sh]
            F[(sh, m)] = Fd
            for g0 in divisors(m):
                phi["%d,%d,%d" % (sh, m, g0)] = sum(d * Fd[d] for d in divisors(g0))
    if len(sys.argv) > 3 and sys.argv[3] != "-":
        P = json.load(open(sys.argv[3]))
        worst = 0.0
        wk = None
        for k, v in P["phi"].items():
            d = abs(v - phi[k])
            if d > worst:
                worst = d
                wk = k
        print("registered json: %d values of Phi(g0); largest |mine - registered| = %.3e at %s (mine %.10f, registered %.10f)" % (
            len(P["phi"]), worst, wk, phi[wk], P["phi"][wk]))
        worstF = 0.0
        for sh in range(2):
            for i, m in enumerate(P["moduli"]):
                ds = divisors(m)
                reg = P["F"][sh][i]
                for j, d in enumerate(ds):
                    worstF = max(worstF, abs(reg[j] - F[(sh, m)][d]))
        print("registered json: largest |mine - registered| over the type shares F_d = %.3e" % worstF)
    if len(sys.argv) > 4:
        bad = 0
        n = 0
        worst = 0.0
        for ln in open(sys.argv[4]):
            mo = re.match(r"^\s+m =\s+(\d+) sheet (\d): (.*)$", ln)
            if not mo:
                continue
            m = int(mo.group(1))
            sh = int(mo.group(2))
            for part in mo.group(3).split(";"):
                mm = re.match(r"g0=(\d+): ([0-9.]+) \[([0-9.]+), ([0-9.]+)\] \| ([0-9.]+) \| ([0-9.]+) \| ([0-9.]+)", part.strip())
                g0 = int(mm.group(1))
                ap = float(mm.group(5))
                mine = phi["%d,%d,%d" % (sh, m, g0)]
                n += 1
                worst = max(worst, abs(mine - ap))
                if abs(mine - ap) > 0.5000001e-4:
                    bad += 1
        print("stored output %s: a-priori column, %d of %d values differ from mine at 4 decimals (largest difference %.2e)" % (sys.argv[4], bad, n, worst))
    for m in (9, 25, 49):
        print("  m = %d: " % m + "  ".join("sheet %d: " % sh + ", ".join("g0=%d: %.4f" % (g0, phi["%d,%d,%d" % (sh, m, g0)]) for g0 in divisors(m)) for sh in range(2)))


if __name__ == "__main__":
    main()
