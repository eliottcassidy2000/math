#!/usr/bin/env python3
"""PSLQ hunt for S(4) (and A/B split) at ~500 digits."""
import json, itertools, sys
import mpmath as mm
from mpmath import mp

DPS = 500
mp.dps = DPS
V = json.load(open('/tmp/math-wt-coinC2/04-computation/Sk_520.json'))['vals']
Sk = {int(k): mm.mpf(v) for k, v in V.items()}

pi = mm.pi
s2, s3, s5 = mm.sqrt(2), mm.sqrt(3), mm.sqrt(5)
l2, l3 = mm.log(2), mm.log(3)
L2 = mm.log(1 + s2)
G = mm.catalan

C = {}
C["1"] = mm.mpf(1)
C["pi^2"] = pi ** 2
C["G"] = G
C["log2^2"] = l2 ** 2
C["pi*log2"] = pi * l2
C["L2^2"] = L2 ** 2
C["pi*L2"] = pi * L2
C["log2*L2"] = l2 * L2
C["Cl2(pi/4)"] = mm.clsin(2, pi / 4)
C["pi"] = pi
C["log2"] = l2
C["L2"] = L2
C["sqrt2"] = s2


def run(target, names, tag, maxcoeff=10 ** 6, tolfrac=0.90, sqrt2=False):
    vals, labs = [], []
    for n in names:
        vals.append(C[n]); labs.append(n)
        if sqrt2:
            vals.append(s2 * C[n]); labs.append("sqrt2*" + n)
    tol = mm.mpf(10) ** (-int(DPS * tolfrac))
    rel = mm.pslq([target] + vals, tol=tol, maxcoeff=maxcoeff, maxsteps=10 ** 6)
    if rel is None:
        print(f"  [{tag}] none (n={len(vals)}, maxcoeff={maxcoeff})")
        return None
    resid = rel[0] * target + sum(mm.mpf(c) * v for c, v in zip(rel[1:], vals))
    terms = " ".join(f"{'+' if c>0 else '-'} {abs(c)}*{l}"
                     for c, l in zip(rel[1:], labs) if c)
    print(f"  [{tag}] {rel[0]}*X {terms}    resid={mm.nstr(resid,4)}")
    return rel


W2 = ["1", "pi^2", "G", "log2^2", "pi*log2", "L2^2", "pi*L2", "log2*L2"]

X = pi * Sk[4]
print("pi*S(4)      =", mm.nstr(X, 50))
print("pi*S(4)/sqrt2=", mm.nstr(X / s2, 50))
print()
print("--- pi*S(4) against weight-2 sets")
run(X, W2, "w2")
run(X, W2, "w2+sqrt2", sqrt2=True)
run(X, W2 + ["Cl2(pi/4)"], "w2+Cl2")
run(X, W2 + ["Cl2(pi/4)"], "w2+Cl2+sqrt2", sqrt2=True)
print()
print("--- pi*S(4)/sqrt2 against weight-2 sets")
Y = X / s2
run(Y, W2, "w2")
run(Y, W2 + ["Cl2(pi/4)"], "w2+Cl2")
print()
print("--- mixed weight 1 and 2")
run(X, ["1", "pi", "log2", "L2", "pi^2", "G", "log2^2", "pi*log2", "L2^2",
        "pi*L2", "log2*L2"], "mix", sqrt2=True)
