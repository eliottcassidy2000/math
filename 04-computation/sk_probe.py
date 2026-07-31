#!/usr/bin/env python3
"""Fast, unbuffered PSLQ probes.  Usage: sk_probe.py <block>"""
import json, sys
import mpmath as mm
from mpmath import mp

DPS = 400
mp.dps = DPS
V = json.load(open('/tmp/math-wt-coinC2/04-computation/Sk_520.json'))['vals']
Sk = {int(k): mm.mpf(v) for k, v in V.items()}

pi = mm.pi
s2, s3, s5, s6, s10, s15 = (mm.sqrt(n) for n in (2, 3, 5, 6, 10, 15))
one = mm.mpf(1)

A = {
    "pi": pi, "log2": mm.log(2), "log3": mm.log(3), "log5": mm.log(5),
    "L(1+s2)": mm.log(1 + s2), "L(2+s3)": mm.log(2 + s3),
    "L(s3+s2)": mm.log(s3 + s2), "L(s5+2)": mm.log(s5 + 2),
    "L(s5+s2)": mm.log(s5 + s2), "L(phi)": mm.log((1 + s5) / 2),
    "at(s2)": mm.atan(s2), "at(s6)": mm.atan(s6), "at(s5)": mm.atan(s5),
    "at(s10)": mm.atan(s10), "at(s2/3)": mm.atan(s2 / 3),
    "at(s3/5)": mm.atan(s3 / 5), "at(1/2)": mm.atan(one / 2),
    "at(1/3)": mm.atan(one / 3), "at(2s2)": mm.atan(2 * s2),
    # weight 2
    "pi^2": pi ** 2, "G": mm.catalan,
    "log2^2": mm.log(2) ** 2, "pi*log2": pi * mm.log(2),
    "L2^2": mm.log(1 + s2) ** 2, "pi*L2": pi * mm.log(1 + s2),
    "log2*L2": mm.log(2) * mm.log(1 + s2),
    "Cl2(pi/4)": mm.clsin(2, pi / 4), "Cl2(pi/3)": mm.clsin(2, pi / 3),
    "Cl2(pi/6)": mm.clsin(2, pi / 6),
    "at2^2": mm.atan(s2) ** 2, "pi*at2": pi * mm.atan(s2),
    "at2*L2": mm.atan(s2) * mm.log(1 + s2),
    "G14^4/pi^3": mm.gamma(one / 4) ** 4 / pi ** 3,
    "pi^3/G14^4": pi ** 3 / mm.gamma(one / 4) ** 4,
}
M = {"1": one, "s2": s2, "s3": s3, "s5": s5, "s6": s6, "s10": s10, "s15": s15}


def run(T, alph, muls, tag, maxc=10 ** 5, tolfrac=0.90):
    vals, labs = [one], ["1"]
    for m in muls:
        if m != "1":
            vals.append(M[m]); labs.append(m)
        for a in alph:
            vals.append(M[m] * A[a]); labs.append(a if m == "1" else m + "*" + a)
    rel = mm.pslq([T] + vals, tol=mm.mpf(10) ** (-int(DPS * tolfrac)),
                  maxcoeff=maxc, maxsteps=300000)
    if rel is None:
        print(f"  [{tag}] none  n={len(vals)}"); sys.stdout.flush(); return None
    r = rel[0] * T + sum(mm.mpf(c) * v for c, v in zip(rel[1:], vals))
    body = " ".join(("+" if c > 0 else "-") + f" {abs(c)}*{l}"
                    for c, l in zip(rel[1:], labs) if c)
    print(f"  [{tag}] {rel[0]}*X {body}   resid={mm.nstr(r,4)}")
    sys.stdout.flush()
    return rel


W1 = ["pi", "log2", "log3", "L(1+s2)", "L(2+s3)", "L(s3+s2)", "at(s2)"]
W2 = ["pi^2", "G", "log2^2", "pi*log2", "L2^2", "pi*L2", "log2*L2",
      "at2^2", "pi*at2", "at2*L2"]

blk = sys.argv[1]
if blk == "A":                       # k=4 weight-1, growing fields
    T = pi * Sk[4]
    run(T, W1, ["1", "s2", "s3", "s6"], "k4 w1 Q(s2,s3)")
    run(T, W1 + ["at(s6)", "at(1/2)", "at(1/3)", "at(2s2)", "at(s2/3)"],
        ["1", "s2", "s3", "s6"], "k4 w1 wide")
elif blk == "B":                     # k=4 weight-2
    T = pi * Sk[4]
    run(T, W2, ["1", "s2"], "k4 w2 Q(s2)")
    run(T, W2 + ["Cl2(pi/4)", "Cl2(pi/3)", "Cl2(pi/6)"], ["1", "s2"], "k4 w2+Cl")
    run(T, W2, ["1", "s2", "s3", "s6"], "k4 w2 Q(s2,s3)")
elif blk == "C":                     # k=5
    T = pi * Sk[5]
    run(T, ["pi", "log2", "log5", "L(1+s2)", "L(phi)", "L(s5+2)", "L(s5+s2)",
            "at(s2)", "at(s5)"], ["1", "s2", "s5", "s10"], "k5 w1 Q(s2,s5)")
    run(T, W2, ["1", "s2", "s5", "s10"], "k5 w2")
elif blk == "D":                     # even k
    for k in (6, 8, 12):
        T = pi * Sk[k]
        run(T, W1, ["1", "s2", "s3", "s6"], f"k{k} w1")
        run(T, W2, ["1", "s2"], f"k{k} w2")
elif blk == "E":                     # gamma-period mixtures for k=4
    T = pi * Sk[4]
    run(T, ["pi", "log2", "L(1+s2)", "at(s2)", "G14^4/pi^3", "pi^3/G14^4"],
        ["1", "s2", "s3", "s6"], "k4 gamma")
