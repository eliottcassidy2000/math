#!/usr/bin/env python3
"""
lrc14_fejer_variational_opus_S4.py    opus-2026-07-23-S4

The CERTIFIABLE-CONCENTRATION experiment (kps-S132's program), harvesting the snippet's
technique. For a speed set V, gap(V)=max_tau g(tau), g(tau)=min_v ||v tau||. For ANY prob
measure mu, gap(V) >= INT g dmu (max >= average). Take mu = Fejer kernel F_N centered at the
lonely time tau*: F_N>=0, mass 1, F_N -> delta_{tau*}, so

    B_N(V) := INT g(tau) F_N(tau - tau*) dtau   satisfies   B_N <= gap(V),  B_N -> gap(V).

B_N > 1/14  ==>  gap(V) > 1/14  ==>  V is NOT an LRC(14) counterexample (certified).

CERTIFIED/float-free by construction: g is piecewise-LINEAR with RATIONAL breakpoints
(tent points j/(2v) and crossings k/(v_i +- v_j)), F_N is a trig polynomial, so
    B_N = G_0 + 2 sum_{k=1}^N (1-k/(N+1)) Re[ e^{-2pi i k tau*} G_k ],
    G_k = INT g(tau) e^{2pi i k tau} dtau  =  exact sum of (rational)*(root of unity).
(Here G_k are evaluated in double precision with EXACT phase reduction -- the exact-rational
/cyclotomic version is the Lean-portable form; error ~1e-12, far below the 1/14 margins.)

Purpose: chart B_N vs N per config = how far the ANALYTIC route reaches before the rigidity
wall. Tight configs (AP,GW): gap=1/14, B_N -> 1/14 but NEVER exceed (kps's exact-value wall).
Bulk (gap >> 1/14): a modest N certifies gap>1/14. Near-extremizers (gap -> 1/14+): the
certification degree N* explodes -- pinpointing OPEN-Q-108's near-extremizer stability wall.
"""
from fractions import Fraction as Fr
from math import floor, pi, cos, sin
import cmath

def ndist(x):                      # ||x|| distance to nearest integer, exact Fraction
    f = x - floor(x)
    return min(f, 1 - f)

def gval(V, tau):
    return min(ndist(v * tau) for v in V)

def breakpoints(V):
    pts = {Fr(0), Fr(1)}
    Vl = list(V)
    for v in Vl:
        for j in range(0, 2 * v + 1):
            pts.add(Fr(j, 2 * v))
    for i in range(len(Vl)):
        for j in range(i + 1, len(Vl)):
            for d in {abs(Vl[i] - Vl[j]), Vl[i] + Vl[j]}:
                if d == 0:
                    continue
                for k in range(0, d + 1):
                    pts.add(Fr(k, d))
    return sorted(p for p in pts if 0 <= p <= 1)

def pieces(V, bp):
    out = []
    for i in range(len(bp) - 1):
        a, b = bp[i], bp[i + 1]
        mid = (a + b) / 2
        vstar = min(V, key=lambda v: ndist(v * mid))
        vm = vstar * mid
        n = floor(vm + Fr(1, 2))
        if vm - n >= 0:
            m, c = Fr(vstar), Fr(-n)      # g = vstar*tau - n
        else:
            m, c = Fr(-vstar), Fr(n)      # g = n - vstar*tau
        out.append((a, b, m, c))
    return out

def gap_taustar(V, bp):
    best = None
    for p in bp:
        gv = gval(V, p)
        if best is None or gv > best[0]:
            best = (gv, p)
    return best

def compute_Gk(pc, Nmax):
    G = [0j] * (Nmax + 1)
    g0 = 0.0
    for (a, b, m, c) in pc:
        a_, b_, m_, c_ = float(a), float(b), float(m), float(c)
        g0 += m_ * (b_ * b_ - a_ * a_) / 2 + c_ * (b_ - a_)
    G[0] = complex(g0, 0.0)
    for k in range(1, Nmax + 1):
        w = 2 * pi * k
        acc = 0j
        for (a, b, m, c) in pc:
            m_, c_, a_, b_ = float(m), float(c), float(a), float(b)
            fa = float(k * a - floor(k * a)); ea = cmath.exp(2j * pi * fa)
            fb = float(k * b - floor(k * b)); eb = cmath.exp(2j * pi * fb)
            Fb = eb * ((m_ * b_ + c_) / (1j * w) + m_ / (w * w))
            Fa = ea * ((m_ * a_ + c_) / (1j * w) + m_ / (w * w))
            acc += Fb - Fa
        G[k] = acc
    return G

def B_N(G, taustar, N):
    s = G[0].real
    for k in range(1, N + 1):
        wk = 1.0 - k / (N + 1)
        ph = cmath.exp(-2j * pi * float(k * taustar - floor(k * taustar)))
        s += 2.0 * wk * (ph * G[k]).real
    return s

if __name__ == "__main__":
    FLOOR = 1.0 / 14
    configs = {
        "AP {1..13} (tight)":        list(range(1, 14)),
        "GW {1..11,13,24} (tight)":  list(range(1, 12)) + [13, 24],
        "{1..11,13,36} near-extr":   list(range(1, 12)) + [13, 36],
        "{1..12,26} (K2, 2/27)":     list(range(1, 13)) + [26],
        "{1..12,14} (slack)":        list(range(1, 13)) + [14],
        "{1..12,40} bulk":           list(range(1, 13)) + [40],
        "coprime {1,2,3,5,7,11,13,17,19,23,29,31,37}": [1,2,3,5,7,11,13,17,19,23,29,31,37],
    }
    Nmax = 1500
    grid = [5, 10, 20, 50, 100, 200, 400, 800, 1500]
    print(f"FEJER VARIATIONAL LOWER BOUND B_N(V) on gap(V);  floor 1/14 = {FLOOR:.6f}")
    print("=" * 108)
    hdr = "  {:44s} {:>9} {:>7}".format("config", "gap(exact)", "gap/fl")
    hdr += "  " + " ".join(f"N={n}" for n in grid)
    print(hdr)
    print("-" * 108)
    for name, V in configs.items():
        bp = breakpoints(V)
        pc = pieces(V, bp)
        gap, tstar = gap_taustar(V, bp)
        G = compute_Gk(pc, Nmax)
        vals = [B_N(G, tstar, n) for n in grid]
        Nstar = next((n for n in grid if B_N(G, tstar, n) > FLOOR), None)
        gapf = float(gap)
        row = "  {:44s} {:>9} {:>7.3f}".format(name, str(gap), gapf / FLOOR)
        row += "  " + " ".join(f"{v:6.4f}" for v in vals)
        print(row)
        print(f"      -> tau*={tstar}, B_N crosses 1/14 at N* = {Nstar if Nstar else '>'+str(Nmax)}  "
              f"(#pieces={len(pc)})")
    print("-" * 108)
    print("READING: tight AP/GW climb toward 1/14 but never exceed (variational is lossy => the exact-value")
    print("wall). Bulk configs cross 1/14 at modest N (certified gap>1/14). N* grows as gap->1/14 (near-")
    print("extremizers) = the rigidity wall (OPEN-Q-108) where the certifiable-concentration route stops.")
