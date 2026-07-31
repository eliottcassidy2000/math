#!/usr/bin/env python3
"""
Lane A1 PSLQ engine for S(k) = sum C(2n,n)C(4n,2n)/((kn+1)64^n).

S(k) = int_0^1 2F1(1/4,3/4;1;t^k) dt, evaluated by AGM (ellipk) + tanh-sinh.
Everything is computed at WORK dps and rounded to TARGET dps for PSLQ so that
quadrature error never masquerades as a relation.
"""
import sys, json, time
import mpmath as mm
from mpmath import mp

WORK_EXTRA = 40


def F(z):
    w = mm.sqrt(1 - z)
    m = (1 - w) / (1 + w)
    return mm.sqrt(2 / (1 + w)) * (2 / mm.pi) * mm.ellipk(m)


def S(k, dps):
    old = mp.dps
    mp.dps = dps + WORK_EXTRA
    kk = mm.mpf(k)
    v = mm.quad(lambda t: F(t ** kk), [0, 1], maxdegree=10 + int(mp.dps / 25))
    mp.dps = dps
    return +v


def compute_table(dps, ks, path):
    t0 = time.time()
    out = {}
    for k in ks:
        out[k] = mm.nstr(S(k, dps), dps - 5)
    json.dump({"dps": dps, "vals": {str(k): v for k, v in out.items()}},
              open(path, "w"))
    print(f"wrote {path} ({time.time()-t0:.1f}s)")
    return out


def load_table(path):
    d = json.load(open(path))
    return d["dps"], {int(k): mm.mpf(v) for k, v in d["vals"].items()}


# ------------------------------------------------------------------ constants
def consts():
    """All constants evaluated at the CURRENT mp.dps."""
    pi = mm.pi
    s2, s3, s5 = mm.sqrt(2), mm.sqrt(3), mm.sqrt(5)
    l2, l3, l5 = mm.log(2), mm.log(3), mm.log(5)
    L2 = mm.log(1 + s2)                       # arcsinh 1
    L3 = mm.log(s3 + s2)
    L23 = mm.log(2 + s3)
    Lphi = mm.log((1 + s5) / 2)
    G = mm.catalan
    C = {
        "1": mm.mpf(1), "pi": pi, "pi^2": pi ** 2,
        "log2": l2, "log3": l3, "log5": l5,
        "L2": L2, "L3": L3, "L23": L23, "Lphi": Lphi,
        "G": G,
        "log2^2": l2 ** 2, "pi*log2": pi * l2,
        "L2^2": L2 ** 2, "pi*L2": pi * L2, "log2*L2": l2 * L2,
        "atanS2": mm.atan(s2), "atan25": mm.atan(s2 / 5),
        "Cl2_pi4": mm.clsin(2, pi / 4),
        "Cl2_pi3": mm.clsin(2, pi / 3),
        "sqrt2": s2, "sqrt3": s3, "sqrt5": s5,
        "G14^4/pi^3": mm.gamma(mm.mpf(1) / 4) ** 4 / pi ** 3,
        "pi^3/G14^4": pi ** 3 / mm.gamma(mm.mpf(1) / 4) ** 4,
    }
    return C


def run(target, names, C, mults=("1",), tag="", tolfrac=0.85, maxcoeff=10 ** 8,
        verbose=True):
    vals, labs = [], []
    for m in mults:
        for n in names:
            vals.append(C[n] if m == "1" else C[m] * C[n])
            labs.append(n if m == "1" else f"{m}*{n}")
    tol = mm.mpf(10) ** (-int(mp.dps * tolfrac))
    rel = mm.pslq([target] + vals, tol=tol, maxcoeff=maxcoeff, maxsteps=2 * 10 ** 6)
    if rel is None:
        if verbose:
            print(f"  [{tag}] none  (n={len(vals)}, dps={mp.dps}, maxcoeff={maxcoeff})")
        return None
    resid = rel[0] * target + sum(mm.mpf(c) * v for c, v in zip(rel[1:], vals))
    terms = " ".join(f"{'+' if c>0 else '-'} {abs(c)}*{l}"
                     for c, l in zip(rel[1:], labs) if c)
    if verbose:
        print(f"  [{tag}] {rel[0]}*X {terms}")
        print(f"       resid={mm.nstr(resid,5)}  maxc={max(abs(c) for c in rel)}")
    return rel, resid


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "table":
        dps = int(sys.argv[2])
        ks = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 20, 24]
        compute_table(dps, ks, f"/tmp/math-wt-coinC2/04-computation/Sk_{dps}.json")
