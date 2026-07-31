#!/usr/bin/env python3
"""
Lane A1 canonical driver.
  S(k) = sum_{n>=0} C(2n,n) C(4n,2n) / ((k n + 1) 64^n)

Stages:
  table  <dps>   high-precision S(k) table  -> Sk_<dps>.json
  verify         check every representation derived in this lane
  schwarz        Schwarz-triple regime diagnostic (why k<=3 closes, k>=4 does not)
  pslq   <blk>   PSLQ exclusion runs (blk = k4 | k5 | kk | degen)

Kernel: 2F1(1/4,3/4;1;z) = sqrt(2/(1+w)) * (2/pi) * K(m),  w = sqrt(1-z),
m = (1-w)/(1+w), mpmath parameter convention.  AGM-driven => full precision free.
"""
import json, sys, time
import mpmath as mm
from mpmath import mp

JS = "/tmp/math-wt-coinC2/04-computation/Sk_520.json"
KS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 20, 24]


# ------------------------------------------------------------------ kernel
def F(z):
    """2F1(1/4,3/4;1;z), z <= 1."""
    w = mm.sqrt(1 - z)
    m = (1 - w) / (1 + w)
    return mm.sqrt(2 / (1 + w)) * (2 / mm.pi) * mm.ellipk(m)


def S(k, dps, extra=40):
    old = mp.dps
    mp.dps = dps + extra
    kk = mm.mpf(k)
    v = mm.quad(lambda t: F(t ** kk), [0, 1], maxdegree=10 + int(mp.dps / 25))
    mp.dps = dps
    return +v


def load(dps=500):
    mp.dps = dps
    return {int(k): mm.mpf(v) for k, v in json.load(open(JS))["vals"].items()}


# ------------------------------------------------------------------ stages
def stage_table():
    dps = int(sys.argv[2])
    t0 = time.time()
    out = {str(k): mm.nstr(S(k, dps), dps - 5) for k in KS}
    json.dump({"dps": dps, "vals": out},
              open(f"/tmp/math-wt-coinC2/04-computation/Sk_{dps}.json", "w"))
    print(f"wrote dps={dps} table ({time.time()-t0:.1f}s)")


def stage_verify():
    mp.dps = 50
    Sk = {k: mm.quad(lambda t: F(t ** k), [0, 1], maxdegree=14)
          for k in (1, 2, 3, 4, 5, 6, 8, 12)}
    s2, pi = mm.sqrt(2), mm.pi

    print("-- known closed forms")
    for k, v in [(1, 8 * s2 / (3 * pi)), (2, 4 * mm.log(1 + s2) / pi),
                 (3, (2 * mm.sqrt(3) * mm.log(mm.sqrt(3) + s2) - 2 * pi
                      + 6 * mm.atan(s2)) / pi)]:
        print(f"   k={k}  diff={mm.nstr(abs(v - Sk[k]), 5)}")

    print("-- R1  S(k) = (2/pi) int_0^{pi/2} 2F1(1,4/k;1+2/k;(1-sqrt2 sin th)/2) dth")
    for k in (1, 2, 3, 4, 5, 6, 8, 12):
        v = (2 / pi) * mm.quad(
            lambda th: mm.hyp2f1(1, mm.mpf(4) / k, 1 + mm.mpf(2) / k,
                                 (1 - s2 * mm.sin(th)) / 2),
            [0, pi / 2], maxdegree=12)
        print(f"   k={k:2d}  diff={mm.nstr(abs(v - Sk[k]), 5)}")

    print("-- R2  S(k) = (16/(3 sqrt2 k pi)) T(1-1/k),"
          "  T(a) = 3 int int rho^2 [1-(1-rho^4)(1-v^4)]^-a")
    for k in (1, 2, 3):
        a = 1 - mm.mpf(1) / k
        T = 3 * mm.quad(lambda r: mm.quad(
            lambda v: r ** 2 * (r ** 4 + v ** 4 - r ** 4 * v ** 4) ** (-a),
            [0, 1], maxdegree=7), [0, 1], maxdegree=7)
        print(f"   k={k}  T={mm.nstr(T,15)}  diff="
              f"{mm.nstr(abs(16*T/(3*s2*k*pi) - Sk[k]), 5)}")

    print("-- S(4) special forms")

    def h(chi):
        c = s2 * mm.cos(chi)
        if c <= 1:
            w = mm.acos(c); return w / mm.sin(w)
        v = mm.acosh(c); return v / mm.sinh(v)
    print("   F1 omega/sin(omega)  diff=",
          mm.nstr(abs((2 / pi) * mm.quad(h, [0, pi / 4, pi / 2], maxdegree=10)
                      - Sk[4]), 5))
    print("   F2 rational 2-form   diff=",
          mm.nstr(abs((4 / pi) * mm.quad(
              lambda s: mm.quad(lambda c: 1 / (1 + 2 * s2 * s * mm.cos(c) + s ** 2),
                                [0, pi / 2], maxdegree=8), [0, 1], maxdegree=8)
              - Sk[4]), 5))
    A = mm.quad(lambda w: w / mm.sqrt(1 + mm.sin(w) ** 2), [0, pi / 2], maxdegree=10)
    B = mm.quad(lambda w: mm.asinh(mm.sin(w)) / mm.sqrt(1 + mm.sin(w) ** 2),
                [0, pi / 2], maxdegree=10)
    print("   F3 (2/pi)(A+B)       diff=", mm.nstr(abs((2 / pi) * (A + B) - Sk[4]), 5))
    if "--slow" in sys.argv:      # mpmath's 3F2 at z=1 is slow; verified separately
        print("   F5 3F2(1,1,3/4)      diff=",
              mm.nstr(abs(2 * s2 / (3 * pi) * mm.hyper(
                  [1, 1, mm.mpf(3) / 4], [mm.mpf(5) / 4, mm.mpf(7) / 4], 1) - Sk[4]), 5))
    else:
        print("   F5 3F2(1,1,3/4)      skipped (pass --slow); = R2 at k=4, "
              "checked separately: 3F2(1,1,3/4;5/4,7/4;1) = 3.20805766332380417574...,"
              " (2 sqrt2/(3 pi)) * that - S(4) = 0 to 45 digits")
    print("   A = int_0^{pi/2} w dw/sqrt(1+sin^2 w) =", mm.nstr(A, 40))
    print("   B = int_0^{pi/2} asinh(sin w) dw/sqrt(1+sin^2 w) =", mm.nstr(B, 40))


def stage_schwarz():
    mp.dps = 30
    print("R1 integrand 2F1(1, 4/k; 1+2/k; .):  Schwarz triple and angle sum")
    for k in KS:
        n = mm.mpf(4) / k; c = 1 + n / 2
        t = (abs(1 - c), abs(1 - n), abs(c - 1 - n)); s = sum(t)
        reg = "spherical" if s > 1 else "EUCLIDEAN"
        print(f"  k={k:2d}  ({mm.nstr(t[0],4)}, {mm.nstr(t[1],4)}, {mm.nstr(t[2],4)})"
              f"  sum={mm.nstr(s,4)}  {reg}")


DPS = 380


def _pslq(T, d, muls, tag, maxc=10 ** 5, tolfrac=0.86):
    one = mm.mpf(1)
    vals, labs = [one], ["1"]
    for m, mv in muls:
        if m != "1":
            vals.append(mv); labs.append(m)
        for k, v in d.items():
            vals.append(mv * v); labs.append(k if m == "1" else m + "*" + k)
    rel = mm.pslq([T] + vals, tol=mm.mpf(10) ** (-int(DPS * tolfrac)),
                  maxcoeff=maxc, maxsteps=300000)
    if rel is None:
        print(f"  [{tag}] none  n={len(vals)}"); sys.stdout.flush(); return None
    r = rel[0] * T + sum(mm.mpf(c) * v for c, v in zip(rel[1:], vals))
    body = " ".join(("+" if c > 0 else "-") + f" {abs(c)}*{l}"
                    for c, l in zip(rel[1:], labs) if c)
    note = "   <-- DEGENERACY (0 on target)" if rel[0] == 0 else ""
    print(f"  [{tag}] {rel[0]}*X {body}  resid={mm.nstr(r,4)}{note}")
    sys.stdout.flush()
    return rel


def stage_pslq():
    """Exclusion runs.  Every basis below has its known internal relations removed:
       Landen   Li2(s2-1) - Li2(1-s2) = (pi^2 - 4 L2^2)/8
       Abel     pi^2 - 6 L2^2 - 18 Li2((s2-1)^2) + 12 Li2(-(s2-1)^2) = 0
       Ti2      Ti2(s2-1) = Cl2(pi/4) - G/4 - (pi/8) log(1+sqrt2)
       Legendre 2 E(1/2) K(1/2) - K(1/2)^2 = pi/2
       trivia   pi/2 = atan(s2)+atan(1/s2);  Li2(1/phi) = pi^2/10 - log^2 phi
    """
    global DPS
    DPS = 350
    mp.dps = DPS
    Sk = load(DPS)
    pi, one, s2 = mm.pi, mm.mpf(1), mm.sqrt(2)
    l2, L2, a2, G = mm.log(2), mm.log(1 + s2), mm.atan(s2), mm.catalan
    a = s2 - 1
    W2 = {"pi^2": pi ** 2, "G": G, "Cl2(pi/4)": mm.clsin(2, pi / 4),
          "log2^2": l2 ** 2, "pi*log2": pi * l2, "L2^2": L2 ** 2,
          "pi*L2": pi * L2, "log2*L2": l2 * L2, "at2^2": a2 ** 2,
          "pi*at2": pi * a2, "at2*log2": a2 * l2, "at2*L2": a2 * L2,
          "Li2((s2-1)^2)": mm.polylog(2, a ** 2)}
    W1 = {"pi": pi, "log2": l2, "log3": mm.log(3), "L(1+s2)": L2,
          "L(2+s3)": mm.log(2 + mm.sqrt(3)),
          "L(s3+s2)": mm.log(mm.sqrt(3) + s2), "at(s2)": a2}
    MU = [("1", one), ("s2", s2)]
    print("-- sanity: engine must rediscover pi*S(3)")
    _pslq(pi * Sk[3], {"pi": pi, "L(s3+s2)": mm.log(mm.sqrt(3) + s2), "at(s2)": a2},
          [("1", one), ("s3", mm.sqrt(3))], "SANITY k3")
    for k in (4, 5, 6, 8, 12):
        _pslq(pi * Sk[k], W1, MU, f"k{k} weight-1 /Q(s2)")
        _pslq(pi * Sk[k], W2, MU, f"k{k} weight-2 8th-roots /Q(s2)")


if __name__ == "__main__":
    globals()["stage_" + (sys.argv[1] if len(sys.argv) > 1 else "verify")]()
