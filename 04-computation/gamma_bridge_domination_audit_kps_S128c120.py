#!/usr/bin/env python3
"""gamma_bridge_domination_audit_kps_S128c120.py -- kind-pasteur-2026-07-20-S128c120

AUDIT OF THE DOMINATION STEP IN klein-S351'S GAMMA BRIDGE.

klein-S351 asserts (docstring + .out READING, verbatim):
    "E_r[psi_m] = sum_k c_k k! is dominated by its top term once deg is large, because
     consecutive terms have ratio c_{D-1}/(c_D D) -> 0.  So a NONZERO toral quantity
     forces E_r[psi_m] != 0.  THAT is TNC => NC2."
death-star-S61g then declared GMC(2) COMPLETE = Gamma Bridge + DvdK n=1, "modulo klein's
domination rigor".  So the entire closure rests on that one sentence.  It is asserted,
not proved.  This tests it.

INDEPENDENT DERIVATION -- klein's code is NOT reused (they disclosed a no-op defect in it).
On the {-1,0,1} stratum with M = 1:
    P = Zb a(r) + b(r) + Z c(r),  Z = rho e^{i th}, rho = sqrt r, u = e^{i th}
    R(r,u) := u P = rho a + b u + rho c u^2,  so g_{-1} = rho a, g_0 = b, g_1 = rho c.
The small root of u = t R(r,u); set v := u / (t g_{-1}).  Substituting,
    v = 1 + b t v + (g_1 g_{-1}) t^2 v^2   and   g_1 g_{-1} = rho c * rho a = r a c,
so v is RHO-FREE (klein's claim, independently confirmed here) and
    v = 1 + B(r) t v + W(r) t^2 v^2,   B = b,  W = r a c,
    psi_m = [t^m] log v,      NC2  <=>  E_r[psi_m] = 0 for all m >= 1,   E_r[r^k] = k!.
Note the stratum depends only on (B, W) -- a and c enter only through the product a c.

WHAT IS MEASURED, per m:
    D          = deg_r psi_m
    c_D        = top coefficient           (klein: "exactly the toral quantity")
    SHARE      = |c_D * D!| / |E_r[psi_m]| (klein's claim => SHARE -> 1)
    RATIO      = |c_{D-1} / (c_D * D)|     (klein's claim => RATIO -> 0)
    signs      = do all c_k share one sign? (the alternative mechanism: no cancellation)
"""
import sys
from fractions import Fraction as Fr

MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 24


def pmul(A, B):
    if not A or not B:
        return []
    C = [Fr(0)] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        if a:
            for j, b in enumerate(B):
                C[i + j] += a * b
    while C and C[-1] == 0:
        C.pop()
    return C


def padd(A, B):
    C = [Fr(0)] * max(len(A), len(B))
    for i, a in enumerate(A):
        C[i] += a
    for i, b in enumerate(B):
        C[i] += b
    while C and C[-1] == 0:
        C.pop()
    return C


def fact(n):
    r = 1
    for i in range(2, n + 1):
        r *= i
    return r


def Er(A):                       # E[poly(r)] with r ~ Exp(1): E[r^k] = k!
    return sum(c * fact(k) for k, c in enumerate(A))


def series_v(B, W, N):
    """v = 1 + B t v + W t^2 v^2, as a list of polys-in-r indexed by power of t, to order N."""
    v = [[Fr(1)]] + [[] for _ in range(N)]
    for _ in range(N + 1):                       # order-by-order fixed point
        v2 = [[] for _ in range(N + 1)]
        for i, vi in enumerate(v):
            if not vi:
                continue
            for j, vj in enumerate(v):
                if i + j <= N and vj:
                    v2[i + j] = padd(v2[i + j], pmul(vi, vj))
        new = [[Fr(1)]] + [[] for _ in range(N)]
        for k in range(N + 1):                   # + B t v
            if k >= 1 and v[k - 1]:
                new[k] = padd(new[k], pmul(B, v[k - 1]))
            if k >= 2 and v2[k - 2]:             # + W t^2 v^2
                new[k] = padd(new[k], pmul(W, v2[k - 2]))
        v = new
    return v


def series_log(v, N):
    """log(v) where v = 1 + O(t); returns [t^m] coefficients as polys in r."""
    w = [padd(v[i], [Fr(-1)] if i == 0 else []) for i in range(N + 1)]   # w = v - 1
    out = [[] for _ in range(N + 1)]
    powk = [[Fr(1)]] + [[] for _ in range(N)]                            # w^0
    for k in range(1, N + 1):
        nxt = [[] for _ in range(N + 1)]
        for i, pi in enumerate(powk):
            if not pi:
                continue
            for j, wj in enumerate(w):
                if i + j <= N and wj:
                    nxt[i + j] = padd(nxt[i + j], pmul(pi, wj))
        powk = nxt
        s = Fr((-1) ** (k + 1), k)
        for m in range(N + 1):
            if powk[m]:
                out[m] = padd(out[m], [x * s for x in powk[m]])
    return out


CASES = [
    ("a=1,   b=1,    c=1      (klein's row 1)", [Fr(1)],        [Fr(0), Fr(1)]),
    ("a=1,   b=0,    c=1      (klein's row 2)", [],             [Fr(0), Fr(1)]),
    ("a=1,   b=1+r,  c=1+2r   (klein's row 4)", [Fr(1), Fr(1)], [Fr(0), Fr(1), Fr(2)]),
    ("a=1+r, b=2,    c=1      (klein's VOID row, redone)", [Fr(2)], [Fr(0), Fr(1), Fr(1)]),
    ("a=1,   b=3,    c=1      (new: large B)", [Fr(3)],         [Fr(0), Fr(1)]),
    ("a=1,   b=1,    c=-1     (new: sign mix)", [Fr(1)],        [Fr(0), Fr(-1)]),
]

for name, B, W in CASES:
    print("=" * 92)
    print("CASE  %s      B(r) = %s   W(r) = r*a*c = %s" % (name, B, W))
    print("=" * 92)
    v = series_v(B, W, MMAX)
    lg = series_log(v, MMAX)
    print("  %-4s %-5s %-16s %-20s %-10s %-10s %s"
          % ("m", "deg", "c_D (top)", "E_r[psi_m]", "SHARE", "RATIO", "one sign?"))
    for m in range(2, MMAX + 1, 2):
        p = lg[m]
        if not p:
            print("  %-4d (psi_m = 0)" % m)
            continue
        D = len(p) - 1
        cD = p[D]
        E = Er(p)
        share = abs(Fr(cD * fact(D))) / abs(E) if E != 0 else None
        ratio = abs(p[D - 1] / (cD * D)) if D >= 1 and cD != 0 else None
        nz = [c for c in p if c != 0]
        onesign = all(c > 0 for c in nz) or all(c < 0 for c in nz)
        print("  %-4d %-5d %-16s %-20s %-10s %-10s %s"
              % (m, D, str(cD)[:15],
                 (str(E)[:17] + "..") if len(str(E)) > 19 else str(E),
                 ("%.4f" % float(share)) if share is not None else "-",
                 ("%.4f" % float(ratio)) if ratio is not None else "-",
                 onesign))
    print()
