"""Independent exhaustive census for THM-4819 (E) at n = 4, 6 (no pruning at all) and n = 8 (pruned only by the
trivial score formula s'(v) = n-1-s(v) + [v first] - [v last]).

For every balanced class (gentourng -q -d(m-1)Dm n), every vertex v and BOTH roles (source / sink), decide whether some
all-black native word of length exactly K = m-1 ends with v a source (resp. sink) and T - v strong.  A native move along a
Hamilton path Q reverses every arc not on Q; Q is black unless its positional reversal is an anti-automorphism.
Prints per order: instances (vertex, role) that succeed, split into the dipole-forced role and the other role.
Usage: python3 census_small_even.py n
"""
import subprocess, sys
from itertools import permutations

n = int(sys.argv[1]); m = n // 2; K = m - 1

def read_classes():
    out = subprocess.run(["gentourng", "-q", f"-d{m-1}D{m}", str(n)], capture_output=True, text=True).stdout.split()
    res = []
    for s in out:
        A = [[0] * n for _ in range(n)]; k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if s[k] == '1': A[i][j] = 1
                else: A[j][i] = 1
                k += 1
        res.append((s, A))
    return res

def ham_paths(A):
    res = []
    def rec(p, used):
        if len(p) == n: res.append(tuple(p)); return
        for w in range(n):
            if not used >> w & 1 and A[p[-1]][w]: p.append(w); rec(p, used | 1 << w); p.pop()
    for s in range(n): rec([s], 1 << s)
    return res

def blue(A, P):
    sig = {P[i]: P[n - 1 - i] for i in range(n)}
    return all(A[x][y] == A[sig[y]][sig[x]] for x in range(n) for y in range(n) if x != y)

def native(A, Q):
    B = [[(1 - A[i][j]) if i != j else 0 for j in range(n)] for i in range(n)]
    for i in range(n - 1):
        a, b = Q[i], Q[i + 1]; B[a][b] = 1; B[b][a] = 0
    return B

def strong_without(A, v):
    W = [u for u in range(n) if u != v]
    def reach(s, fwd):
        seen = {s}; st = [s]
        while st:
            x = st.pop()
            for y in W:
                if y not in seen and (A[x][y] if fwd else A[y][x]): seen.add(y); st.append(y)
        return seen
    return len(reach(W[0], True)) == len(W) and len(reach(W[0], False)) == len(W)

def search(A, v, role, depth, prune):
    if depth == K:
        s = sum(A[v]); ok_role = (s == n - 1) if role == 'source' else (s == 0)
        return ok_role and strong_without(A, v)
    for Q in ham_paths(A):
        if prune:
            # score after the move; reaching |s - (n-1)/2| = m - 1/2 at depth K needs +1 of |.| per move
            s2 = n - 1 - sum(A[v]) + (Q[0] == v) - (Q[-1] == v)
            if abs(2 * s2 - (n - 1)) != 2 * (depth + 1) + 1: continue
        if blue(A, Q): continue
        if search(native(A, Q), v, role, depth + 1, prune): return True
    return False

classes = read_classes()
forced_ok = forced_tot = other_ok = other_tot = 0; fails = []
for idx, (s, A) in enumerate(classes):
    for v in range(n):
        sv = sum(A[v])
        forced = ('source' if m % 2 == 1 else 'sink') if sv == m else ('sink' if m % 2 == 1 else 'source')
        for role in ('source', 'sink'):
            if n >= 8 and role != forced: continue   # n >= 8: forced role only (the other role is excluded by the
                                                     # dipole parity; checked without pruning at n = 4, 6)
            ok = search(A, v, role, 0, prune=(n >= 8))
            if role == forced:
                forced_tot += 1; forced_ok += ok
                if not ok: fails.append((s, v, sv, role))
            else:
                other_tot += 1; other_ok += ok
print(f"n={n}: balanced classes {len(classes)}; forced-role instances {forced_tot}, succeed {forced_ok}; "
      f"other-role instances {other_tot}, succeed {other_ok}")
for f in fails: print("   forced-role failure: class", f[0], "vertex", f[1], "score", f[2], "role", f[3])
