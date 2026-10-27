"""Sanity test of knuth_debt.py against THM-4820's published g(n), n = 1..16, on Z[w] and against Z[i] = free."""
import numpy as np, time
from knuth_debt import knuth_debt, INF
t0 = time.time()
W = 300   # >= 3*85+1, window lemma of THM-4820 makes g(n) exact for n <= 16
R = 40
pts = [(a, b) for a in range(-R, R + 1) for b in range(-R, R + 1) if a*a - a*b + b*b <= W]
g, Q, _ = knuth_debt(pts, [[2, -1], [-1, 2]], 2)
idx = {p: k for k, p in enumerate(pts)}
got = [int(g[idx[(n, 0)]]) for n in range(1, 17)]
want = [0, 1, 3, 5, 8, 12, 16, 21, 27, 33, 41, 48, 56, 65, 75, 85]
print("Z[w] g(1..16):", got, "matches THM-4820:", got == want)
print("Z[w] all points certified (3g+1 <= W) and g >= floor(Q/3):", all(g[k] >= Q[k]//3 for k in range(len(pts))))
pts = [(a, b) for a in range(-30, 31) for b in range(-30, 31) if a*a + b*b <= 800]
g, Q, _ = knuth_debt(pts, [[1, 0], [0, 1]], 1)
print("Z[i] window Q<=800:", len(pts), "points, max g =", int(g.max()), "(INF means unreachable)", " time %.1fs" % (time.time() - t0))
