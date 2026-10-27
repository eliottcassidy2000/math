"""Audit THM-4824 (D) and HYP-9428 (3): L_2 = {u = v mod 2}, Q = u^2 + 2 v^2, units (+-1, +-1).
(1) exhaustive box checks of the potential Phi = 2 floor(Q/8); (2) exceptional norms <= 891 (W = 900, window lemma
Q <= 4C + 7); (3) on-axis points (2k, 0) for k <= 29 (W = 3500) and, as an extension, k <= KMAX."""
import sys, time
import numpy as np
from collections import defaultdict
from knuth_debt import knuth_debt, INF

Qf = lambda u, v: u*u + 2*v*v
Phi = lambda q: 2 * (q // 8)
R = 26
box = [(u, v) for u in range(-R, R + 1) for v in range(-R, R + 1) if (u - v) % 2 == 0]
res = {Qf(*p) % 8 for p in box}
bad_pot = bad_even = bad_perp = 0
for (u, v) in box:
    q1 = Qf(u, v)
    for (s, t) in box:
        q2 = Qf(s, t)
        b = 2 * (u*s + 2*v*t)
        if b % 2: bad_even += 1
        if q1 % 8 == 4 and q2 % 8 == 4 and b == 0: bad_perp += 1
        if Phi(q1 + q2 + b) - Phi(q1) - Phi(q2) > max(b, 0): bad_pot += 1
print(f"(1) box |u|,|v| <= {R} ({len(box)} pts, {len(box)**2} ordered pairs): Q mod 8 residues {sorted(res)}; "
      f"odd beta {bad_even}; perpendicular 4-mod-8 pairs {bad_perp}; potential violations {bad_pot}")

def solve(W):
    Rr = int(W ** 0.5) + 2
    pts = [(u, v) for u in range(-Rr, Rr + 1) for v in range(-Rr, Rr + 1) if (u - v) % 2 == 0 and Qf(u, v) <= W]
    g, Q, _ = knuth_debt(pts, [[1, 0], [0, 2]], 1)
    return pts, g, Q

t0 = time.time()
W = 900
pts, g, Q = solve(W)
cert = (4 * g + 7 <= W)
allle = [k for k in range(len(pts)) if Q[k] <= 891]
print(f"(2) W = {W}: {len(pts)} points; points with Q <= 891: {len(allle)}, all certified exact: {all(cert[k] for k in allle)}; "
      f"certified points: {int(cert.sum())}, max certified norm {int(Q[cert].max())}  [{time.time()-t0:.1f}s]")
assert all(g[k] >= Phi(int(Q[k])) for k in range(len(pts)))
exc = defaultdict(set); tot = defaultdict(int); excpts = defaultdict(list)
for k in range(len(pts)):
    if Q[k] <= 891:
        q = int(Q[k]); tot[q] += 1
        e = int(g[k]) - Phi(q)
        if e:
            exc[q].add(e); excpts[q].append(pts[k])
claimed = [83, 100, 107, 147, 211, 227, 251, 307, 347, 419, 475, 507, 547, 587, 643, 659, 676, 747, 787, 811, 827]
mine = sorted(exc)
print("    exceptional norms:", [(q, sorted(exc[q])) for q in mine])
print("    matches the THM-4824 (D) list:", mine == claimed, "; excess-4 norms:", [q for q in mine if 4 in exc[q]],
      "; excess values used:", sorted(set().union(*exc.values())))
partial = [(q, len(excpts[q]), tot[q]) for q in mine if len(excpts[q]) != tot[q]]
print("    norms where only SOME points are exceptional (norm, #exceptional, #points):", partial)
print("    first-quadrant exceptional representatives:", {q: sorted(p for p in excpts[q] if p[0] >= 0 and p[1] >= 0) for q in mine})

KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 29
W = max(3500, 4 * (2 * (4 * KMAX * KMAX // 8) + 6) + 7)
t0 = time.time()
pts, g, Q = solve(W)
idx = {p: k for k, p in enumerate(pts)}
print(f"(3) on-axis (2k, 0), W = {W}, {len(pts)} points  [{time.time()-t0:.1f}s]")
isprime = lambda n: n > 1 and all(n % d for d in range(2, int(n**0.5) + 1))
rows = []
for k in range(1, KMAX + 1):
    j = idx[(2*k, 0)]
    q = int(Q[j]); gg = int(g[j]); ok = 4*gg + 7 <= W
    e = gg - Phi(q)
    rows.append((k, q, gg, e, ok))
print("    k: excess (cert)  ->", " ".join(f"{k}:{e}{'' if ok else '?'}" for k, q, gg, e, ok in rows))
print("    excess 4 at k =", [k for k, q, gg, e, ok in rows if e == 4], "; excess 2 at k =", [k for k, q, gg, e, ok in rows if e == 2],
      "; other nonzero:", [(k, e) for k, q, gg, e, ok in rows if e not in (0, 2, 4)], "; all certified:", all(r[4] for r in rows))
print("    primes p <= KMAX by class mod 8 -> excess:", {p: (p % 8, next(e for k, q, gg, e, ok in rows if k == p)) for p in range(2, KMAX + 1) if isprime(p)})
print("    THM table k <= 29 (g):", [(k, gg) for k, q, gg, e, ok in rows if k <= 29])
