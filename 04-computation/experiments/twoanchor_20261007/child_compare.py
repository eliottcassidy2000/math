#!/usr/bin/env python3
"""Haar comparison of deletion children on residual sources with prescribed two-run length J
(mac-mini-2026-10-07-twoanchor). Sources n = 2^K t - 1 with 400-bit random t (2-adic regime); for each D the
equal-time merge of the D-chain (source run end x vs child run end y_D) is tested within post-run depth s."""
import random, sys
def v2(x): return (x & -x).bit_length() - 1
def T(x): return (3*x + 1) >> 1 if x & 1 else x >> 1
def absorb_time(x, y, k, smax):
    for s in range(smax):
        if x == y and k == 0: return s
        k += (x & 1) - (y & 1); x, y = T(x), T(y)
    return None
rnd = random.Random(42)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
S = 400
print(f"post-run depth budget {S}; {N} sources per J; entries = fraction certified")
print("   J   r  | D=1    D=2    D=3    D=4    D=5    D=6   any D<=8 | median post-run depth (D=1, D=3)")
for J in (1, 2, 3, 4, 6, 10, 16, 24, 32):
    for r in (1, 2):
        j = 2*J + r
        if j < 3: continue
        cert = {D: 0 for D in range(1, 9)}; anyc = 0; d1 = []; d3 = []
        for _ in range(N):
            K = rnd.randint(12, 40)
            M = 1 << (j - 1)
            t0 = pow(3, -(K-1), M)
            while True:
                t = t0 + M * rnd.getrandbits(400)
                if t & 1 and v2(3**(K-1)*t - 1) == j - 1: break
            x = 2*3**(K-1)*t - 1
            got = False
            for D in range(1, 9):
                y = (x + 1)//3**D - 1
                a = absorb_time(x, y, D, 2*J + S)
                if a is not None:
                    cert[D] += 1; got = True
                    if D == 1: d1.append(a - 2*J)
                    if D == 3: d3.append(a - 2*J)
            anyc += got
        med = lambda L: sorted(L)[len(L)//2] if L else None
        print(f"  {J:3d}  {r}  | " + "  ".join(f"{cert[D]/N:.3f}" for D in range(1, 7)) + f"   {anyc/N:.3f}  | {med(d1)}, {med(d3)}", flush=True)
