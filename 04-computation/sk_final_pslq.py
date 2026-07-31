#!/usr/bin/env python3
"""Final broad PSLQ for S(4), S(5) and friends -- weight<=2 at 8th/12th/20th
roots of unity, plus elliptic (Gamma(1/4)) mixtures."""
import json, sys
import mpmath as mm
from mpmath import mp

DPS = 400
mp.dps = DPS
V = json.load(open('/tmp/math-wt-coinC2/04-computation/Sk_520.json'))['vals']
Sk = {int(k): mm.mpf(v) for k, v in V.items()}

pi = mm.pi
one = mm.mpf(1)
s2, s3, s5 = mm.sqrt(2), mm.sqrt(3), mm.sqrt(5)
l2 = mm.log(2)
L2 = mm.log(1 + s2)          # log fundamental unit of Q(sqrt2)
L3 = mm.log(2 + s3)
Lp = mm.log((1 + s5) / 2)
a2 = mm.atan(s2)
G = mm.catalan
K12 = mm.ellipk(one / 2)      # = Gamma(1/4)^2/(4 sqrt(pi))
E12 = mm.ellipe(one / 2)


def base8():
    d = {
        "pi^2": pi ** 2, "G": G, "Cl2(pi/4)": mm.clsin(2, pi / 4),
        "log2^2": l2 ** 2, "pi*log2": pi * l2,
        "L2^2": L2 ** 2, "pi*L2": pi * L2, "log2*L2": l2 * L2,
        "at2^2": a2 ** 2, "pi*at2": pi * a2, "at2*log2": a2 * l2,
        "at2*L2": a2 * L2,
        "Li2(s2-1)": mm.polylog(2, s2 - 1), "Li2(1-s2)": mm.polylog(2, 1 - s2),
    }
    return d


def base12():
    d = base8()
    d.update({"Cl2(pi/3)": mm.clsin(2, pi / 3), "Cl2(pi/6)": mm.clsin(2, pi / 6),
              "L3^2": L3 ** 2, "pi*L3": pi * L3, "log2*L3": l2 * L3,
              "log3^2": mm.log(3) ** 2, "pi*log3": pi * mm.log(3)})
    return d


def base20():
    d = base8()
    d.update({"Cl2(pi/5)": mm.clsin(2, pi / 5), "Cl2(2pi/5)": mm.clsin(2, 2 * pi / 5),
              "Lphi^2": Lp ** 2, "pi*Lphi": pi * Lp, "log2*Lphi": l2 * Lp,
              "log5^2": mm.log(5) ** 2, "pi*log5": pi * mm.log(5),
              "Li2(1/phi)": mm.polylog(2, 2 / (1 + s5))})
    return d


def ell():
    return {"K": K12, "pi*K": pi * K12, "K*L2": K12 * L2, "K*at2": K12 * a2,
            "E": E12, "pi*E": pi * E12, "K^2": K12 ** 2, "pi/K": pi / K12,
            "1/K": 1 / K12, "L2/K": L2 / K12, "at2/K": a2 / K12,
            "pi^2/K": pi ** 2 / K12}


def run(T, d, muls, tag, maxc=10 ** 5, tolfrac=0.88):
    vals, labs = [one], ["1"]
    for m, mv in muls:
        if m != "1":
            vals.append(mv); labs.append(m)
        for k, v in d.items():
            vals.append(mv * v); labs.append(k if m == "1" else m + "*" + k)
    rel = mm.pslq([T] + vals, tol=mm.mpf(10) ** (-int(DPS * tolfrac)),
                  maxcoeff=maxc, maxsteps=400000)
    if rel is None:
        print(f"  [{tag}] none  n={len(vals)}"); sys.stdout.flush(); return None
    r = rel[0] * T + sum(mm.mpf(c) * v for c, v in zip(rel[1:], vals))
    body = " ".join(("+" if c > 0 else "-") + f" {abs(c)}*{l}"
                    for c, l in zip(rel[1:], labs) if c)
    print(f"  [{tag}] {rel[0]}*X {body}   resid={mm.nstr(r,4)}")
    sys.stdout.flush(); return rel


blk = sys.argv[1]
M1 = [("1", one)]
M2 = [("1", one), ("s2", s2)]
if blk == "k4":
    T = pi * Sk[4]
    run(T, base8(), M1, "k4 w2-8th /Q")
    run(T, base8(), M2, "k4 w2-8th /Q(s2)")
    run(T, base12(), M1, "k4 w2-12th /Q")
    run(T, ell(), M2, "k4 elliptic")
elif blk == "k5":
    T = pi * Sk[5]
    run(T, base20(), M1, "k5 w2-20th /Q")
    run(T, base20(), [("1", one), ("s5", s5)], "k5 w2-20th /Q(s5)")
    run(T, ell(), M2, "k5 elliptic")
elif blk == "kk":
    for k in (6, 8, 12):
        T = pi * Sk[k]
        run(T, base8(), M2, f"k{k} w2-8th /Q(s2)")
        run(T, base12(), M1, f"k{k} w2-12th /Q")
