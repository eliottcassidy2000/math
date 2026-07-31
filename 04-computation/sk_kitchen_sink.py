#!/usr/bin/env python3
"""Kitchen-sink PSLQ: pi*S(4), pi*S(5), and the A/B halves of S(4),
against mixed weight 0/1/2 + elliptic (K,E at m=1/2) bases, dependencies removed.
Known internal relations deliberately broken:
  Legendre  2 E K - K^2 = pi/2           -> keep K, E, drop K^2
  Landen    Li2(s2-1)-Li2(1-s2) = (pi^2-4 L2^2)/8  -> keep Li2(s2-1) only
"""
import json, sys
import mpmath as mm
from mpmath import mp

DPS = 380
mp.dps = DPS
V = json.load(open('/tmp/math-wt-coinC2/04-computation/Sk_520.json'))['vals']
Sk = {int(k): mm.mpf(v) for k, v in V.items()}

pi, one = mm.pi, mm.mpf(1)
s2, s5 = mm.sqrt(2), mm.sqrt(5)
l2, L2, a2, G = mm.log(2), mm.log(1 + s2), mm.atan(s2), mm.catalan
K = mm.ellipk(one / 2)
E = mm.ellipe(one / 2)

SINK = {
    "pi": pi, "log2": l2, "L2": L2, "at2": a2,
    "pi^2": pi ** 2, "G": G, "Cl2(pi/4)": mm.clsin(2, pi / 4),
    "log2^2": l2 ** 2, "L2^2": L2 ** 2, "at2^2": a2 ** 2,
    "pi*log2": pi * l2, "pi*L2": pi * L2, "pi*at2": pi * a2,
    "log2*L2": l2 * L2, "at2*log2": a2 * l2, "at2*L2": a2 * L2,
    "Li2(s2-1)": mm.polylog(2, s2 - 1),
    "K": K, "E": E, "pi*K": pi * K, "pi*E": pi * E,
    "K*L2": K * L2, "K*at2": K * a2, "K*log2": K * l2, "E*L2": E * L2,
}
ELL = {"K": K, "E": E, "pi*K": pi * K, "pi*E": pi * E, "K*L2": K * L2,
       "K*at2": K * a2, "K*log2": K * l2, "E*L2": E * L2, "E*at2": E * a2,
       "pi": pi, "pi^2": pi ** 2, "G": G}


def run(T, d, muls, tag, maxc=10 ** 5, tolfrac=0.86):
    vals, labs = [one], ["1"]
    for m, mv in muls:
        if m != "1":
            vals.append(mv); labs.append(m)
        for k, v in d.items():
            vals.append(mv * v); labs.append(k if m == "1" else m + "*" + k)
    rel = mm.pslq([T] + vals, tol=mm.mpf(10) ** (-int(DPS * tolfrac)),
                  maxcoeff=maxc, maxsteps=400000)
    if rel is None:
        print(f"  [{tag}] none  n={len(vals)}"); sys.stdout.flush(); return
    r = rel[0] * T + sum(mm.mpf(c) * v for c, v in zip(rel[1:], vals))
    body = " ".join(("+" if c > 0 else "-") + f" {abs(c)}*{l}"
                    for c, l in zip(rel[1:], labs) if c)
    print(f"  [{tag}] {rel[0]}*X {body}   resid={mm.nstr(r,4)}")
    sys.stdout.flush()


M1 = [("1", one)]
# A-half of S(4):  A = int_0^{pi/2} w dw / sqrt(1+sin^2 w) = (1/sqrt2) int_0^{pi/2} F(u|1/2) du
mp.dps = DPS + 30
A = mm.quad(lambda w: w / mm.sqrt(1 + mm.sin(w) ** 2), [0, mm.pi / 2], maxdegree=12)
mp.dps = DPS
A = +A
B = pi * Sk[4] / 2 - A

print("A =", mm.nstr(A, 40))
print("B =", mm.nstr(B, 40))
run(pi * Sk[4], SINK, M1, "pi*S(4) kitchen-sink")
run(A, ELL, M1, "A elliptic")
run(B, ELL, M1, "B elliptic")
run(pi * Sk[5], SINK, M1, "pi*S(5) kitchen-sink")
