"""Wave 58.  Conserved-charge digit test for two-hub spine blocks (wave 45 charge principle applied to the two-hub transfer).
State per hub order: (s0,s1,s2) in {E,G,B}^3 (27); two orders -> 54 dims.  Letter = (x_i, y_i).
Appending block beta: coordinates c -> c N_beta on the reachable row space R; H = c F_e.
A digit block needs N F = lam F + d Phi, N Phi = Phi  <=>  u = (N - I) F is an eigenvector of N (eigenvalue lam)."""
from fractions import Fraction as Fr
import itertools
def idx(s0, s1, s2): return s0 + 3 * s1 + 9 * s2
def letter_matrix(A, B):
    """27x27 transition for one spine vertex with first-hub bit A, second-hub bit B (A=1: first hub beats vertex)."""
    M = [[0] * 27 for _ in range(27)]
    for s0, s1, s2 in itertools.product(range(3), repeat=3):
        i = idx(s0, s1, s2)
        M[i][idx(1 if A == 0 else 2, s1, s2)] += 1
        if s1 != 0 or A == 1: M[i][idx(s0, 1 if B == 0 else 2, s2)] += 1
        if s2 != 0 or B == 1: M[i][idx(s0, s1, 1)] += 1
    return M
def final(ab):
    return [1 if (2 not in (s0, s1, s2) and (s1 != 0 or ab)) else 0 for s2 in range(3) for s1 in range(3) for s0 in range(3)]
# fix index order: idx = s0 + 3 s1 + 9 s2  -> iterate s2 outer
def final_vec(ab):
    v = [0] * 27
    for s0, s1, s2 in itertools.product(range(3), repeat=3):
        if 2 not in (s0, s1, s2) and (s1 != 0 or ab): v[idx(s0, s1, s2)] = 1
    return v
letters = [(x, y) for x in (0, 1) for y in (0, 1)]
A54 = {}
for (x, y) in letters:
    M1 = letter_matrix(x, y); M2 = letter_matrix(y, x)   # order 1: u1 first (A=x,B=y); order 2: u2 first (A=y,B=x)
    M = [[0] * 54 for _ in range(54)]
    for i in range(27):
        for j in range(27):
            M[i][j] = M1[i][j]; M[27 + i][27 + j] = M2[i][j]
    A54[(x, y)] = M
def fvec(e):
    return final_vec(e) + final_vec(1 - e)    # e = 1: u1 -> u2
init = [0] * 54; init[0] = 1; init[27] = 1
def vecmat(v, M): return [sum(v[i] * M[i][j] for i in range(len(v)) if v[i]) for j in range(len(M[0]))]
# reachable row space by Gaussian elimination over Fractions
basis = []; piv = []
def reduce(v):
    v = [Fr(x) for x in v]
    for b, p in zip(basis, piv):
        if v[p] != 0:
            f = v[p] / b[p]; v = [vi - f * bi for vi, bi in zip(v, b)]
    return v
def add(v):
    r = reduce(v)
    for p, x in enumerate(r):
        if x != 0:
            basis.append(r); piv.append(p); return True
    return False
queue = [init]; add(init); raw = [init]
while queue:
    v = queue.pop()
    for l in letters:
        w = vecmat(v, A54[l])
        if add(w): queue.append(w); raw.append(w)
d = len(basis); print("reachable row space dimension:", d)
# coordinates: express any reachable row r as c.Rm with Rm = raw basis vectors (use 'raw' as basis: independent set)
Rm = raw[:d]
# solve c from r: build matrix and invert via Fractions
import copy
def solve_coords(r):
    # solve c Rm = r ; use least-squares-free exact elimination on transposed system
    n = d; m = 54
    Mt = [[Fr(Rm[i][j]) for i in range(n)] + [Fr(r[j])] for j in range(m)]   # m equations, n unknowns
    rows = Mt; prow = 0; pivcols = []
    for col in range(n):
        sel = None
        for rr in range(prow, m):
            if rows[rr][col] != 0: sel = rr; break
        if sel is None: continue
        rows[prow], rows[sel] = rows[sel], rows[prow]
        pv = rows[prow][col]; rows[prow] = [x / pv for x in rows[prow]]
        for rr in range(m):
            if rr != prow and rows[rr][col] != 0:
                f = rows[rr][col]; rows[rr] = [a - f * b for a, b in zip(rows[rr], rows[prow])]
        pivcols.append(col); prow += 1
    c = [Fr(0)] * n
    for i, col in enumerate(pivcols): c[col] = rows[i][n]
    for rr in range(prow, m): assert rows[rr][n] == 0, "not in span"
    return c
N = {}
for l in letters:
    N[l] = [solve_coords(vecmat(Rm[i], A54[l])) for i in range(d)]   # row i: coords of Rm[i] A_l
F = {e: [sum(Rm[i][j] * fvec(e)[j] for j in range(54)) for i in range(d)] for e in (0, 1)}
def matmul(P, Q): return [[sum(P[i][k] * Q[k][j] for k in range(len(Q)) if P[i][k]) for j in range(len(Q[0]))] for i in range(len(P))]
def matvec(P, v): return [sum(P[i][k] * v[k] for k in range(len(v)) if P[i][k]) for i in range(len(P))]
print("checking blocks...")
found = []
for L in range(1, 6):
    for blk in itertools.product(letters, repeat=L):
        Nb = N[blk[0]]
        for l in blk[1:]: Nb = matmul(Nb, N[l])
        for e in (0, 1):
            Fe = F[e]; NF = matvec(Nb, Fe); u = [a - b for a, b in zip(NF, Fe)]
            if all(x == 0 for x in u):
                continue
            Nu = matvec(Nb, u)
            # check Nu = lam u
            lam = None; ok = True
            for a, b in zip(Nu, u):
                if b == 0:
                    if a != 0: ok = False; break
                else:
                    r = Fr(a) / b
                    if lam is None: lam = r
                    elif r != lam: ok = False; break
            if ok and lam is not None:
                found.append((L, blk, e, lam))
print("blocks (length <= 5) with a conserved-charge digit structure:", len(found))
from collections import Counter
print(Counter((f[0], f[3]) for f in found).most_common(20))
for f in found[:20]: print(f)

# ---- part 2: middle blocks with fixed reader suffixes ----
import itertools
from collections import Counter, defaultdict
blocks = {}
for L in range(1, 5):
    for blk in itertools.product(letters, repeat=L):
        Nb = N[blk[0]]
        for l in blk[1:]: Nb = matmul(Nb, N[l])
        blocks[blk] = Nb
found = defaultdict(list)
for Ls in range(0, 4):
    for suf in itertools.product(letters, repeat=Ls):
        Ns = None
        for l in suf: Ns = N[l] if Ns is None else matmul(Ns, N[l])
        for e in (0, 1):
            Fs = F[e] if Ns is None else matvec(Ns, F[e])
            for blk, Nb in blocks.items():
                NF = matvec(Nb, Fs); u = [a - b for a, b in zip(NF, Fs)]
                if all(x == 0 for x in u): continue
                Nu = matvec(Nb, u); lam = None; ok = True
                for a, b in zip(Nu, u):
                    if b == 0:
                        if a != 0: ok = False; break
                    else:
                        r = Fr(a) / b
                        if lam is None: lam = r
                        elif r != lam: ok = False; break
                if ok and lam is not None:
                    Phi = tuple(a - lam * b for a, b in zip(NF, Fs))   # = d * Phi_unit
                    # normalise direction
                    nz = next((x for x in Phi if x != 0), None)
                    key = (suf, e, lam, tuple(x / nz for x in Phi) if nz else None)
                    found[key].append((blk, nz))
print("reader/charge structures (suffix, e, lambda, Phi-direction) with >= 2 digit blocks:")
cnt = 0
for key, lst in found.items():
    if len(lst) >= 2:
        cnt += 1
        if cnt <= 25: print(f"  suffix={key[0]} e={key[1]} lambda={key[2]} digits(d)={sorted(set(str(d) for _, d in lst))[:12]} blocks={len(lst)}")
print("total keys:", len(found), " with >=2 blocks:", cnt)
print("lambda values seen:", Counter(k[2] for k in found).most_common(10))
