"""a3_offset_transport.py (AUDIT) - are the plus carriers of Proposition 2.5 literally the carriers 3, 21, 35 of
multiplicand_saturation_chanpak_20261111/chanpak/OFFSET_TRANSPORT.md?  Build base (rows = out-neighbour masks),
selected set O, close with a one-vertex tail and the closer (B => t => v, v -> O), and test isomorphism with the
packet's kernels (codes 172278145, 71647617 at order 8; 17332994050 at order 9) or their converses."""
import sys, os
from aud_common import *
out = []
def P(s):
    print(s); out.append(s)
def from_rows(rows):
    n = len(rows)
    T = [set(j for j in range(n) if r >> j & 1) for r in rows]
    for i in range(n):
        for j in range(n):
            if i != j: assert (j in T[i]) != (i in T[j]), "not a tournament"
    return T
def closing(B, O, nx=1):
    nb = len(B); N = nb + nx + 1; v = nb + nx
    o = [set() for _ in range(N)]
    for a in range(nb):
        o[a] |= set(B[a]); o[a] |= set(range(nb, nb + nx))
        if a in O: o[v].add(a)
        else: o[a].add(v)
    for a in range(nx):
        o[nb + a] |= set(range(nb + a + 1, nb + nx)); o[nb + a].add(v)
    return o
def conv(A):
    m = len(A)
    return [set(j for j in range(m) if j != i and i in A[j]) for i in range(m)]
def iso(A, B):
    m = len(A)
    if sorted(len(s) for s in A) != sorted(len(s) for s in B): return False
    img = [-1] * m; used = [False] * m
    def rec(k):
        if k == m: return True
        for v in range(m):
            if used[v] or len(B[v]) != len(A[k]): continue
            if all((k in A[i]) == (v in B[img[i]]) for i in range(k)):
                img[k] = v; used[v] = True
                if rec(k + 1): return True
                used[v] = False
        return False
    return rec(0)
data = [("3", (40, 33, 3, 6, 15, 28), {3}, 8, 172278145), ("21", (42, 28, 57, 48, 33, 2), {3, 4}, 8, 71647617), ("35", (86, 124, 40, 113, 100, 65, 4), {6}, 9, 17332994050)]
for name, rows, O, m, code in data:
    B = from_rows(rows)
    K = closing(B, O)
    Kp = decode(code, m)
    P("carrier %s of OFFSET_TRANSPORT.md closed with a one-vertex tail (order %d): H = %d; packet kernel code %d: H = %d; isomorphic: %s; to the converse: %s"
      % (name, len(K), ham_py(K), code, ham_py(Kp), iso(K, Kp), iso(K, conv(Kp))))
open(os.path.join(HERE, "a3_offset_transport.out"), "w").write("\n".join(out) + "\n")
