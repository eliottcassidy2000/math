#!/usr/bin/env python3
"""
mm94 lane: symmetrized profiles of genuine tensor characters on polynomial multiplication C(a,b).

Quantum functionals F_theta (Christandl-Vrana-Zuiddam, JAMS 2023) are spectral points; on a free tensor
they equal Strassen's support functional, and for the tight tensor C(a,b) (support {(i,j,i+j)})
    F_theta(C(a,b)) = max_{P on support} 2^{theta_1 H(P_X) + theta_2 H(P_Y) + theta_3 H(P_Z)}.
Every F_theta has dot-product exponents p_X = 1 - theta_1 etc., so t = 2/3 and the paper's profile is
    P_theta(a,b) = ( prod_{sigma in S_3} F_{sigma theta}(C(a,b)) )^(1/4).
Checks (NUMERICAL): P_theta >= P_min (least element), P_theta satisfies (C), (T) and the k = 5 gadget;
diagonal exponent of P_theta is 3/2 (vs 4/3 for P_min).  Collatz slice a = 2: C(2,b) is the carry-free
tensor of n -> 3n = n + 2n; profiles are asymptotically s*(b + 1/2), the coordinate that linearises 3x+1.
"""
import numpy as np, itertools, math, sys

LN2 = math.log(2.0)

def entropy_bits(p):
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())

def F_theta(a, b, th, iters=4000, tol=1e-13):
    """max over distributions on {(i,j)} of sum_k th_k H(marginal_k), in bits; returns 2^value."""
    I = np.repeat(np.arange(a), b); J = np.tile(np.arange(b), a); K = I + J
    P = np.full(a * b, 1.0 / (a * b))
    def val(P):
        return (th[0] * entropy_bits(np.bincount(I, P, a)) + th[1] * entropy_bits(np.bincount(J, P, b))
                + th[2] * entropy_bits(np.bincount(K, P, a + b - 1)))
    v = val(P)
    for it in range(iters):
        PX = np.bincount(I, P, a); PY = np.bincount(J, P, b); PZ = np.bincount(K, P, a + b - 1)
        g = -(th[0] * np.log2(PX[I] + 1e-300) + th[1] * np.log2(PY[J] + 1e-300) + th[2] * np.log2(PZ[K] + 1e-300))
        eta = 1.0
        while True:                                 # mirror ascent with backtracking
            Q = P * np.exp(eta * LN2 * (g - g.max()))
            Q /= Q.sum()
            vq = val(Q)
            if vq >= v - 1e-15 or eta < 1e-6:
                break
            eta *= 0.5
        gap = float((P * g).sum())                  # = v ; KKT gap: max g - v
        P, vold, v = Q, v, vq
        if g.max() - gap < tol:
            break
    return 2.0 ** v

def P_theta(a, b, th, cache={}):
    key = (a, b, th)
    if key not in cache:
        prod = 1.0
        for s in itertools.permutations(th):
            prod *= F_theta(a, b, s)
        cache[key] = prod ** 0.25
    return cache[key]

def Pmin(a, b):
    m, M = min(a, b), max(a, b)
    s = 1.0
    for j in range(1, m):
        s *= 1 + 1 / (3 * j)
    return s * (2 * M + m - 1) / 2

thetas = [(1/3, 1/3, 1/3), (0.5, 0.25, 0.25), (0.6, 0.3, 0.1), (1.0, 0.0, 0.0)]
A = 7
print("Check P_theta(1,b) = b, P_theta >= P_min, concavity, tripling, k=5 gadget (tolerance 1e-7 relative)")
for th in thetas:
    ok_bd, ok_dom, ok_conc, ok_tri, ok_k5 = True, True, True, True, True
    worst_ratio = 1e9
    for a in range(1, A + 1):
        for b in range(1, 3 * A + 4):
            v = P_theta(a, b, th)
            if a == 1: ok_bd &= abs(v - b) < 1e-6 * b
            ok_dom &= v >= Pmin(a, b) * (1 - 1e-7)
            worst_ratio = min(worst_ratio, v / Pmin(a, b))
        for b in range(2, 3 * A + 3):
            ok_conc &= 2 * P_theta(a, b, th) >= (P_theta(a, b + 1, th) + P_theta(a, b - 1, th)) * (1 - 1e-7)
        for h in range(1, 4):
            if 3 * h + a - 1 <= 3 * A + 3:
                ok_tri &= P_theta(a, 3 * h + a - 1, th) >= 3 * P_theta(a, h, th) * (1 - 1e-7)
            if 5 * h + 2 * (a - 1) <= 3 * A + 3:
                ok_k5 &= P_theta(a, 5 * h + 2 * (a - 1), th) >= 5 * P_theta(a, h, th) * (1 - 1e-7)
    print(f"  theta={tuple(round(x,3) for x in th)}: P(1,b)=b {ok_bd}; >=P_min {ok_dom} (min ratio {worst_ratio:.4f}); concave {ok_conc}; tripling {ok_tri}; k=5 gadget {ok_k5}")

print("Diagonal growth: P(a,a) for theta = (1/3,1/3,1/3) and flattening vs P_min")
prev = None
for a in (2, 4, 8, 16, 24, 32):
    v = P_theta(a, a, thetas[0]); f = math.sqrt(a * a * (2 * a - 1)); m = Pmin(a, a)
    line = f"  a={a:3d}: sym {v:10.4f}  flat {f:10.4f}  P_min {m:10.4f}"
    if prev:
        line += f"   local exponents: sym {math.log(v/prev[0])/math.log(a/prev[3]):.4f}  flat {math.log(f/prev[1])/math.log(a/prev[3]):.4f}  P_min {math.log(m/prev[2])/math.log(a/prev[3]):.4f}"
    prev = (v, f, m, a)
    print(line)

print("Collatz slice a = 2 (carry-free tensor of n -> 3n): P(2,b) - s*(b+1/2), s = asymptotic slope")
for th, name in ((thetas[0], "sym"), (thetas[2], "(.6,.3,.1)")):
    vals = {b: P_theta(2, b, th) for b in (50, 100, 200, 400)}
    s = (vals[400] - vals[200]) / 200
    print(f"  {name}: slope ~ {s:.6f} (sqrt 2 = {math.sqrt(2):.6f}); offsets P(2,b)/s - b: " +
          ", ".join(f"b={b}: {vals[b]/s - b:.4f}" for b in vals))
print(f"  least element: P_min(2,b) = (4/3)(b + 1/2) exactly for b >= 2; flattening: sqrt(2b(b+1)) = sqrt2*sqrt((b+1/2)^2 - 1/4)")
