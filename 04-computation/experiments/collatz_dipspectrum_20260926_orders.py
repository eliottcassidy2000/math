#!/usr/bin/env python3
"""collatz_dipspectrum_20260926_orders.py -- polynomial orders of the dip spectrum (THM-4498 control;
session collatz-exponent-atlas-20260926, opus, 2026-09-26).

For gamma in (log_4 3, 1) the claim is D_b(X, gamma) = Theta(X^E (log X)^(-1/2)), E = h(gamma/log_2 3),
against Theta(X^(h*) (log X)^(-3/2)) at gamma = 1 (THM-4495). Two controls:
 (1) exact word counts  W_t(gamma) = #{ w in {0,1}^t : o_j log_2 3 - j > (gamma - 1) t for all 1 <= j <= t }
     by a DP on (position, ones), for t up to TMAX: the normalised W_t(gamma) t^(1/2) 2^(-tE) should settle
     in a bounded window (and W_t(1) t^(3/2) 2^(-t h*) likewise, THM-4495);
 (2) the brute-force orbit counts D_+(2^T, gamma) of THM-4487's control (collatz_dipspectrum_20260926.out,
     plus sheet, T = 10..24) normalised by X^E (log_2 X)^(1/2): the same ratios at the sizes where orbits
     were actually counted, next to the alternative normalisation (log_2 X)^(3/2), which drifts.
Also prints the constants of the lower-bound construction (Hoeffding block K, i_0(c)) for the gammas used.
Usage: python3 collatz_dipspectrum_20260926_orders.py [TMAX=600]
"""
import math, os, re, sys

ALPHA = math.log2(3)


def h(p):
    return 0.0 if p <= 0 or p >= 1 else -(p * math.log2(p) + (1 - p) * math.log2(1 - p))


def E(gamma):
    return h(max(0.5, gamma / ALPHA))


def W_dp(t, gamma):
    """words of length t with all prefix sums o_j*ALPHA - j > (gamma-1)*t (float threshold; ties impossible)."""
    thr = (gamma - 1.0) * t
    cur = {0: 1}
    for j in range(1, t + 1):
        nxt = {}
        for o, c in cur.items():
            for step in (0, 1):
                o2 = o + step
                if o2 * ALPHA - j > thr:
                    nxt[o2] = nxt.get(o2, 0) + c
        cur = nxt
    return sum(cur.values())


def main():
    TMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 600
    gammas = [0.82, 0.85, 0.88, 0.91, 0.94, 0.97, 1.0]
    print("gamma:        " + "  ".join("%7.3f" % g for g in gammas))
    print("E(gamma):     " + "  ".join("%7.4f" % E(g) for g in gammas))
    print("(1) exact word counts: W_t(gamma) t^(1/2) 2^(-tE)  [last column gamma = 1 uses t^(3/2)]")
    for t in (25, 50, 100, 150, 200, 300, 400, 500, 600, 800, 1000):
        if t > TMAX:
            break
        row = []
        for g in gammas:
            W = W_dp(t, g)
            pw = 1.5 if g == 1.0 else 0.5
            row.append(2 ** (math.log2(W) + pw * math.log2(t) - t * E(g)) if W > 0 else float('nan'))
        print("   t=%5d:   " % t + "  ".join("%7.4f" % r for r in row))
    # (2) brute-force orbit counts from THM-4487's control
    here = os.path.dirname(os.path.abspath(__file__))
    outp = os.path.join(here, 'collatz_dipspectrum_20260926.out')
    if os.path.exists(outp):
        txt = open(outp, encoding='utf-8', errors='replace').read().split('== sheet 3n-1')[0]
        gam_line = next(l for l in txt.split('\n') if l.strip().startswith('gamma:'))
        gs = [float(x) for x in gam_line.split(':')[1].split()]
        rows = [(int(m.group(1)), [int(x) for x in m.group(2).split()]) for m in re.finditer(r"X=2\^(\d+) D:\s+([0-9 ]+)", txt)]
        print("(2) THM-4487's brute-force D_+(2^T, gamma) normalised: D T^(1/2) / 2^(T E)   (gamma = 0.792 is Korec's endpoint, E = 1)")
        print("    gamma:   " + "  ".join("%7.3f" % g for g in gs))
        for T, D in rows:
            print("    T=%2d:    " % T + "  ".join("%7.4f" % (2 ** (math.log2(d) + 0.5 * math.log2(T) - T * E(g))) for g, d in zip(gs, D)))
        print("    the alternative normalisation D T^(3/2) / 2^(T E) at gamma = 0.91: " + ", ".join("%.2f" % (2 ** (math.log2(D[4]) + 1.5 * math.log2(T) - T * E(gs[4]))) for T, D in rows))
    # constants of the lower-bound construction
    print("(3) lower-bound construction: c = 1 - gamma, Hoeffding tail exp(-c^2 i/(2 alpha^2)); i_0(c) with the tail sum <= 1/2, delta_0 = max(alpha, 0.585 i_0), K_0 = ceil((delta_0 + c + 2)/(c + 0.585))")
    for g in gammas[:-1]:
        c = 1 - g
        q = math.exp(-c * c / (2 * ALPHA * ALPHA))
        i0 = 1
        while q ** i0 / (1 - q) > 0.5:
            i0 += 1
        delta0 = max(ALPHA, math.log2(1.5) * i0)
        K0 = math.ceil((delta0 + c + 2) / (c + math.log2(1.5)))
        print("   gamma=%.2f: c=%.2f, i_0=%d, delta_0=%.2f, K_0=%d, cost 2^(-K_0 E) = %.2e" % (g, c, i0, delta0, K0, 2 ** (-K0 * E(g))))


if __name__ == '__main__':
    main()
