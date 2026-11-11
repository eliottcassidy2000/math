#!/usr/bin/env python3
"""Auditor numbers-vs-outputs: S_n, D_n of THM-5128 (E) by an own transfer (state = two lowest digits read so far,
digits read from the top), n <= 48, against brute force for n <= 22; zero pattern of D_n modulo 12."""
import numpy as np

def brute(nb):
    size = 1 << nb
    q = np.zeros(size, dtype=np.int64)
    for h in range(1, size):
        half = h >> 1
        q[h] = q[half] + (((half + 1) >> 1) if (h & 1) else 0)
    x = 2 * np.arange(size, dtype=np.int64) + 1 - 2 * q
    out = {}
    for n in range(1, nb + 1):
        s = int(np.where(q[: 1 << n] % 2 == 0, 1, -1).sum())
        d = int((x[: 1 << n] % 4 == 1).sum() - (x[: 1 << n] % 4 == 3).sum())
        out[n] = (s, d)
    return out

def transfer(nmax):
    # q(2h + 1) = q(h) + ceil(h/2); ceil(h/2) mod 2 = (h0 + h1) mod 2 for h = ...h1 h0; x = 2h + 1 - 2q, so
    # x mod 4 = 1 iff (h - q) even iff (new lowest digit + q) even
    cur = {(0, 0, 0): 1}     # (lowest digit, second lowest, q mod 2) -> count
    out = {}
    for n in range(0, nmax + 1):
        s = sum(c * (1 if st[2] == 0 else -1) for st, c in cur.items())
        d = sum(c * (1 if (st[0] + st[2]) % 2 == 0 else -1) for st, c in cur.items())
        out[n] = (s, d)
        nxt = {}
        for (a, b, qq), c in cur.items():
            for dig in (0, 1):
                q2 = (qq + (dig * ((a + b) % 2))) % 2
                key = (dig, a, q2)
                nxt[key] = nxt.get(key, 0) + c
        cur = nxt
    return out

b = brute(22)
t = transfer(48)
print("transfer = brute force for n <= 22:", all(b[n] == t[n] for n in range(1, 23)))
print("D_18, D_20, D_27 =", t[18][1], t[20][1], t[27][1])
print("n <= 48 with D_n = 0:", [n for n in range(1, 49) if t[n][1] == 0])
print("  residues mod 12:", sorted({n % 12 for n in range(1, 49) if t[n][1] == 0}), " and D_n != 0 at residues", sorted({n % 12 for n in range(1, 49) if t[n][1] != 0}))
print("n <= 48 with S_n = 0:", [n for n in range(1, 49) if t[n][0] == 0])
corank = {1: 1, 2: 0, 3: 1, 4: 2, 5: 1, 0: 0}
ok = True
for n in range(1, 49):
    r = n - corank[n % 6]
    for v in t[n]:
        ok &= v == 0 or abs(v) == 2 ** (n - r // 2)
print("S_n, D_n in {0, +-2^(n - rank/2)} with corank 1,0,1,2,1,0 for n = 1..5,0 mod 6, n <= 48:", ok)
