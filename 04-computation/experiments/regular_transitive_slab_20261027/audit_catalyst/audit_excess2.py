"""Window-free lower bounds g(x) >= Phi(x) + 1 and g(x) >= Phi(x) + 2 by closure (THM-4820 / HYP-9426 (3) audit).

Independent Python/numpy re-implementation (the C program audit_lattice.c is only compared against, not used):

  tight merge:   beta = 1 with r = s = 1, or -(r+s) <= beta <= 0          (r = Q(y) mod 3, s = Q(z) mod 3)
  G  = least set containing 0 and the six units, closed under tight splits with both parts nonzero in G.
       g(x) = Phi(x)  iff  x in G;  all-tight trees have node norms <= Q(x) + 1  (confinement lemma).
  Slack-1 lemma: a tree of total slack <= 1 has node norms <= Q(x) + 6, and (following zero-child merges) its first
       nontrivial split is either tight with one part in G1 and the other in G, or of slack exactly 1 with both parts
       in G.  So for a point x with NO tight split, g(x) >= Phi(x) + 2 iff every slack-1 split of x has a part outside G.

Usage: python3 audit_excess2.py N1 xa xb [xa xb ...]    (closure on all points of norm <= N1)
Optionally compares G with column 6 of audit_lat_W<N1>.txt if present.
"""
import os
import sys
import time

import numpy as np

N1 = int(sys.argv[1])
targets = [(int(sys.argv[i]), int(sys.argv[i + 1])) for i in range(2, len(sys.argv), 2)]
t0 = time.time()
R = int((4 * N1 / 3) ** 0.5) + 3
aa, bb = np.meshgrid(np.arange(-R, R + 1), np.arange(-R, R + 1), indexing="ij")
aa, bb = aa.ravel().astype(np.int64), bb.ravel().astype(np.int64)
qq = aa * aa - aa * bb + bb * bb
keep = qq <= N1
PA, PB, PQ = aa[keep], bb[keep], qq[keep]
order = np.argsort(PQ, kind="stable")
PA, PB, PQ = PA[order], PB[order], PQ[order]
NP = len(PA)
grid = -np.ones((2 * R + 1, 2 * R + 1), dtype=np.int64)
grid[PA + R, PB + R] = np.arange(NP)
print(f"{NP} points of norm <= {N1}")


def tight_mask(qx, qy, qz):
    beta = qx - qy - qz
    r, s = qy % 3, qz % 3
    return ((beta == 1) & (r == 1) & (s == 1)) | ((beta <= 0) & (beta >= -(r + s)))


# enumerate tight splits with both parts nonzero (parts of a tight split have norm <= Q(x) + 1)
TX, TY, TZ = [], [], []
cum = np.searchsorted(PQ, np.arange(N1 + 3), side="right")
for i in range(NP):
    xa, xb, qx = int(PA[i]), int(PB[i]), int(PQ[i])
    m = cum[min(qx + 1, N1)]
    ya, yb, qy = PA[:m], PB[:m], PQ[:m]
    za, zb = xa - ya, xb - yb
    qz = za * za - za * zb + zb * zb
    ok = (qz <= N1) & ~((ya == 0) & (yb == 0)) & ~((za == 0) & (zb == 0))
    ok &= tight_mask(qx, qy, qz)
    idx = np.nonzero(ok)[0]
    if len(idx):
        zi = grid[za[idx] + R, zb[idx] + R]
        TX.append(np.full(len(idx), i))
        TY.append(idx)
        TZ.append(zi)
TX, TY, TZ = np.concatenate(TX), np.concatenate(TY), np.concatenate(TZ)
print(f"{len(TX)} ordered tight splits ({time.time() - t0:.1f}s)")
G = np.zeros(NP, dtype=bool)
for u in [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)]:
    G[grid[u[0] + R, u[1] + R]] = True
rounds = 0
while True:
    rounds += 1
    new = G[TY] & G[TZ] & ~G[TX]
    if not new.any():
        break
    G[TX[new]] = True
exact = PQ + (PQ % 3 == 0) <= N1
print(f"closure converged in {rounds} rounds; points with exact status: {int(exact.sum())}; "
      f"not in G (exceptional, window-free): {int((~G & exact).sum())}")

fnc = os.environ.get("AUDIT_CMP", f"audit_lat_W{N1}.txt")
if os.path.exists(fnc):
    dis = 0
    n = 0
    for line in open(fnc):
        t = line.split()
        if t[0] != "P":
            continue
        a, b, q, g, ing = map(int, t[1:])
        if q > N1 or ing == -1:
            continue
        j = grid[a + R, b + R]
        n += 1
        if bool(G[j]) != (ing == 1):
            dis += 1
    print(f"comparison with the C closure in {fnc}: {n} points, {dis} disagreements")


def phi(q):
    return q // 3


for (xa, xb) in targets:
    qx = xa * xa - xa * xb + xb * xb
    assert qx + 6 <= N1, "window too small for the slack-1 lemma"
    m = cum[min(qx + 4, N1)]
    ya, yb, qy = PA[:m], PB[:m], PQ[:m]
    za, zb = xa - ya, xb - yb
    qz = za * za - za * zb + zb * zb
    ok = (qz <= N1) & ~((ya == 0) & (yb == 0)) & ~((za == 0) & (zb == 0))
    beta = qx - qy - qz
    slack = np.maximum(beta, 0) + qy // 3 + qz // 3 - qx // 3
    tight = ok & (slack == 0)
    s1 = ok & (slack == 1)
    idx = np.nonzero(s1)[0]
    zi = grid[za[idx] + R, zb[idx] + R]
    both_ex = exact[idx] & exact[zi]
    bothG = G[idx] & G[zi]
    print(f"x = {xa} + {xb}w (norm {qx}, Phi {phi(qx)}): in G: {bool(G[grid[xa + R, xb + R]])}; "
          f"ordered tight splits {int(tight.sum())}; ordered slack-1 splits {len(idx)} "
          f"(beta values {sorted(set(beta[idx].tolist()))}); all parts have exact G-status: {bool(both_ex.all())}; "
          f"slack-1 splits with BOTH parts in G: {int(bothG.sum())}")
    if int(tight.sum()) == 0 and int(bothG.sum()) == 0 and bool(both_ex.all()):
        print(f"   => g({xa} + {xb}w) >= Phi + 2 = {phi(qx) + 2} (window-free)")
    bad_parts = sorted({(int(PA[j]), int(PB[j])) for j in list(idx) + list(zi) if not G[j]},
                       key=lambda v: v[0] * v[0] - v[0] * v[1] + v[1] * v[1])
    print(f"   exceptional parts used by its slack-1 splits (norms): "
          f"{sorted({a * a - a * b + b * b for a, b in bad_parts})}")
print(f"done ({time.time() - t0:.1f}s)")
