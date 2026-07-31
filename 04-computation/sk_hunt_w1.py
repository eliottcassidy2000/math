#!/usr/bin/env python3
"""
Weight-1 PSLQ hunt guided by the k<=3 data:

  pi*S(1) = 8 sqrt2/3
  pi*S(2) = 4 log(1+sqrt2)
  pi*S(3) = 2 sqrt3 log(sqrt3+sqrt2) - 2 pi + 6 arctan(sqrt2)      [<= rewritten]

so the natural weight-1 alphabet is
  pi, log2, log3, log(1+sqrt2), log(sqrt3+sqrt2), log(2+sqrt3), arctan(sqrt2)
with coefficients in Q(sqrt2, sqrt3) i.e. multipliers {1, sqrt2, sqrt3, sqrt6}.
"""
import json, sys
import mpmath as mm
from mpmath import mp

DPS = 500
mp.dps = DPS
V = json.load(open('/tmp/math-wt-coinC2/04-computation/Sk_520.json'))['vals']
Sk = {int(k): mm.mpf(v) for k, v in V.items()}

pi = mm.pi
s2, s3, s5, s6 = mm.sqrt(2), mm.sqrt(3), mm.sqrt(5), mm.sqrt(6)
ALPH = {
    "pi": pi,
    "log2": mm.log(2),
    "log3": mm.log(3),
    "log5": mm.log(5),
    "L(1+s2)": mm.log(1 + s2),
    "L(s3+s2)": mm.log(s3 + s2),
    "L(2+s3)": mm.log(2 + s3),
    "L(phi)": mm.log((1 + s5) / 2),
    "L(s5+s2)": mm.log(s5 + s2),
    "at(s2)": mm.atan(s2),
    "at(s3)": mm.atan(s3),          # = pi/3, dependent; used only in probes
    "at(1/2)": mm.atan(mm.mpf(1) / 2),
    "at(1/3)": mm.atan(mm.mpf(1) / 3),
    "at(s2/2)": mm.atan(s2 / 2),
    "at(s5)": mm.atan(s5),
    "at(s3/2)": mm.atan(s3 / 2),
}
MUL = {"1": mm.mpf(1), "s2": s2, "s3": s3, "s6": s6, "s5": s5, "s10": mm.sqrt(10)}


def run(target, alph, muls, tag, maxcoeff=10 ** 6, tolfrac=0.92, quiet=False):
    vals, labs = [mm.mpf(1)], ["1"]
    for m in muls:
        if m != "1":
            vals.append(MUL[m]); labs.append(m)
        for a in alph:
            vals.append(MUL[m] * ALPH[a])
            labs.append(a if m == "1" else f"{m}*{a}")
    tol = mm.mpf(10) ** (-int(DPS * tolfrac))
    rel = mm.pslq([target] + vals, tol=tol, maxcoeff=maxcoeff, maxsteps=10 ** 6)
    if rel is None:
        if not quiet:
            print(f"  [{tag}] none  (n={len(vals)}, maxcoeff={maxcoeff})")
        return None
    resid = rel[0] * target + sum(mm.mpf(c) * v for c, v in zip(rel[1:], vals))
    terms = " ".join(f"{'+' if c>0 else '-'} {abs(c)}*{l}"
                     for c, l in zip(rel[1:], labs) if c)
    print(f"  [{tag}] n={len(vals)}  {rel[0]}*X {terms}")
    print(f"        resid={mm.nstr(resid,4)}  maxcoeff={max(abs(c) for c in rel)}")
    return rel


if __name__ == "__main__":
    # sanity: the engine must rediscover k=3
    print("--- sanity: rediscover pi*S(3)")
    run(pi * Sk[3], ["pi", "L(s3+s2)", "at(s2)"], ["1", "s3"], "k3-check")

    base = ["pi", "log2", "log3", "L(1+s2)", "L(s3+s2)", "L(2+s3)", "at(s2)"]
    for k in [4, 5, 6, 8, 12]:
        print(f"--- pi*S({k})")
        T = pi * Sk[k]
        run(T, base, ["1", "s2"], f"k{k}/Q(s2)")
        run(T, base, ["1", "s2", "s3", "s6"], f"k{k}/Q(s2,s3)")
        run(T, base + ["L(phi)", "L(s5+s2)", "at(s5)"], ["1", "s2", "s5", "s10"],
            f"k{k}/Q(s2,s5)")
