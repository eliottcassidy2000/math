#!/usr/bin/env python3
"""THM-790 referee (opus-S305): verify every step of the blue parity law proof
exactly, on ALL grid-symmetric tilings at n = 4..7.

The proof chain (to be canonized):
  L1 (reflection identity): for grid-sym t, the tile out-degree satisfies
      e_v + e_{rho(v)} = (U+L)_v,  rho(v) = n+1-v  =>  a_{rho(v)} = -a_v,
      where a_v = (U+L)_v - 2 e_v is the half score-shift (d_v(tbar) = d_v + 2a_v).
  L2 (antisymmetric scores): d_{rho(v)} = -d_v for grid-sym t
      (=> blue classes have complement-symmetric score multisets).
  L3 (drop decomposition): Dx = 8 * Sum_{v < rho(v)} a_v (d_v + a_v).
  L4 (parity): odd n:  a(d+a) = a mod 2 and Sum_half (U+L)_v is odd => Dx = 8 mod 16;
              even n:  d odd => a(d+a) = a(a+1) = 0 mod 2 => Dx = 0 mod 16.
  L5 (transitive legs): for the transitive tiling, a_v(d_v+a_v) = 0 for all
      interior v; only v = 1 contributes -(n-2): Dx = -8(n-2).
  L6 (half-tiling count): #blue tilings = 2^{(m+f)/2}, f = floor((n-1)/2)
      (fixed anti-diagonal tiles x+y = n+1).
"""
from collections import defaultdict

def check(name, cond, fails=[0]):
    print(("  [PASS] " if cond else "  [FAIL] ") + name)
    if not cond: fails[0] += 1
    return fails[0]

for n in range(4, 8):
    V = list(range(1, n+1))
    tiles = [(x, y) for y in range(1, n-1) for x in range(n, y+1, -1) if x-y >= 2]
    m = len(tiles)
    ti = {t: i for i, t in enumerate(tiles)}
    refl = [ti[(n-y+1, n-x+1)] for (x, y) in tiles]
    f_count = sum(1 for i, t in enumerate(tiles) if refl[i] == i)
    U = {v: max(0, v-2) for v in V}
    L = {v: max(0, n-v-1) for v in V}
    def sco(bits):
        s = {v: (1 if v >= 2 else 0) for v in V}
        for i, (x, y) in enumerate(tiles):
            if bits >> i & 1: s[x] += 1
            else: s[y] += 1
        return s
    def etile(bits):
        e = {v: 0 for v in V}
        for i, (x, y) in enumerate(tiles):
            if bits >> i & 1: e[x] += 1
            else: e[y] += 1
        return e
    full = (1 << m) - 1
    blues = [b for b in range(1 << m)
             if all((b >> i & 1) == (b >> refl[i] & 1) for i in range(m))]
    ok1 = ok2 = ok3 = ok4 = True
    for b in blues:
        e = etile(b); s = sco(b)
        d = {v: 2*s[v] - (n-1) for v in V}
        a = {v: U[v] + L[v] - 2*e[v] for v in V}
        for v in V:
            r = n+1-v
            if e[v] + e[r] != U[v] + L[v]: ok1 = False
            if d[r] != -d[v]: ok2 = False
        sb = sco(b ^ full); db = {v: 2*sb[v] - (n-1) for v in V}
        dx = sum(db[v]**2 for v in V) - sum(d[v]**2 for v in V)
        half = [v for v in V if v < n+1-v]
        if dx != 8 * sum(a[v]*(d[v]+a[v]) for v in half): ok3 = False
        if n % 2 == 1 and dx % 16 != 8: ok4 = False
        if n % 2 == 0 and dx % 16 != 0: ok4 = False
    check(f"n={n} L1 reflection identity e_v + e_rho = U+L on all {len(blues)} blue tilings", ok1)
    check(f"n={n} L2 antisymmetric centered scores d_rho = -d", ok2)
    check(f"n={n} L3 drop decomposition Dx = 8 Sum_half a(d+a)", ok3)
    check(f"n={n} L4 parity law (odd: 8 mod 16 / even: 0 mod 16)", ok4)
    # L4 arithmetic ingredient: Sum_half (U+L)_v odd at odd n
    if n % 2 == 1:
        sh = sum(U[v]+L[v] for v in V if v < n+1-v)
        check(f"n={n} L4' Sum_half (U+L) = {sh} odd", sh % 2 == 1)
    # L5 transitive legs
    e = etile(full); s = sco(full)
    d = {v: 2*s[v]-(n-1) for v in V}; a = {v: U[v]+L[v]-2*e[v] for v in V}
    contrib = {v: a[v]*(d[v]+a[v]) for v in V if v < n+1-v}
    check(f"n={n} L5 transitive drop concentrated at v=1: {contrib} sum=-(n-2)",
          contrib[1] == -(n-2) and all(c == 0 for v, c in contrib.items() if v != 1))
    # L6 half-tiling count
    check(f"n={n} L6 #blue = 2^((m+f)/2): {len(blues)} = 2^{(m+f_count)//2}",
          len(blues) == 2 ** ((m + f_count)//2))
print("done")
