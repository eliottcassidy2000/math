#!/usr/bin/env python3
"""Audit A: the results note's '(Haar-weighted) deletion children D <= 8 certify about 66% of the whole residual branch
within 400 post-run steps' -- own Monte Carlo: residual sources with j - 3 ~ Geom(1/2) (Haar given j >= 3), K uniform in
[12, 40], 400-bit random t; certified if some D <= 8 chain absorbs within 2J + 400 Terras steps after x.
Also reports the 3-adic depth-0 certificate m = (2n-1)/3 (available iff 3 | t), to type the 'any K-uniform rule' claim."""
import random
from a2_barrier import make_source, T

def run(N=1200, S=400, seed=8):
    rnd = random.Random(seed)
    cert = 0; three = 0; byj = {}
    for _ in range(N):
        j = 3
        while rnd.getrandbits(1): j += 1
        J = (j - 1)//2
        K = rnd.randint(12, 40)
        n, t, x = make_source(rnd, K, j, 400)
        ok = False
        for D in range(1, 9):
            y = (x + 1)//3**D - 1
            u, v, k = x, y, D
            for s in range(2*J + S):
                if u == v and k == 0:
                    ok = True; break
                k += (u & 1) - (v & 1); u, v = T(u), T(v)
            if ok: break
        cert += ok
        key = min(j, 7)
        a, b = byj.get(key, (0, 0)); byj[key] = (a + ok, b + 1)
        if t % 3 == 0:
            m = (2*n - 1)//3
            assert 3*m == 2*n - 1 and m < n and ((3*m + 1) // 2 == n)   # U(m) = n with letter 1
            three += 1
    print(f"Haar-weighted residual branch, {N} sources: certified by some D <= 8 within 400 post-run steps: {cert/N:.3f} "
          f"(+-{(cert/N*(1-cert/N)/N)**0.5:.3f})")
    print("   by type: " + ", ".join(f"j={'>=7' if k == 7 else k}: {a}/{b}={a/b:.2f}" for k, (a, b) in sorted(byj.items())))
    print(f"   3-adic depth-0 certificate m = (2n-1)/3 < n with U(m) = n available on {three}/{N} sources (3 | t)")

if __name__ == "__main__":
    run()
