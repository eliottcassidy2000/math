#!/usr/bin/env python3
# audit J1 (THM-4911 (6)): exact eight-point sphere of the signature simplex (S2 from corner powers), Pow(O_T) = n(n-1/s) with s from its definition, Pow(p_v) = 4(n-1+Delta_v) with Delta_v by brute-force 2-paths; regime scan over all classes n = 4..7.
import itertools, subprocess
from fractions import Fraction as Fr
import numpy as np

GT = "/opt/homebrew/bin/gentourng"
def classes(n):
    m = n * (n - 1) // 2; iu = np.triu_indices(n, 1)
    out = subprocess.run([GT, "-q", str(n)], capture_output=True).stdout
    a = np.frombuffer(out, dtype=np.uint8).reshape(-1, m + 1)[:, :m] - 48
    A = np.zeros((len(a), n, n), dtype=np.int64); A[:, iu[0], iu[1]] = a; A[:, iu[1], iu[0]] = 1 - a
    return A
def sub(a, b): return [x - y for x, y in zip(a, b)]
def add(a, b): return [x + y for x, y in zip(a, b)]
def mul(c, a): return [c * x for x in a]
def dt(a, b): return sum(x * y for x, y in zip(a, b))
def solve(M, b):
    m = len(M); A = [list(r) + [bb] for r, bb in zip(M, b)]
    for c in range(m):
        p = next((r for r in range(c, m) if A[r][c] != 0), None)
        if p is None: raise ZeroDivisionError
        A[c], A[p] = A[p], A[c]
        for r in range(m):
            if r != c and A[r][c] != 0:
                f = A[r][c] / A[c][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [A[i][m] / A[i][i] for i in range(m)]
def cc(P):
    p0 = P[0]; D = [sub(p, p0) for p in P[1:]]
    lam = solve([[dt(x, y) for y in D] for x in D], [dt(x, x) / 2 for x in D])
    return add(p0, [sum(l * d[j] for l, d in zip(lam, D)) for j in range(len(p0))])
def det3(G):
    return (G[0][0]*(G[1][1]*G[2][2]-G[1][2]*G[2][1]) - G[0][1]*(G[1][0]*G[2][2]-G[1][2]*G[2][0]) + G[0][2]*(G[1][0]*G[2][1]-G[1][1]*G[2][0]))

def setup(A, v):
    n = len(A); M = A - A.T
    P = [[Fr(int(k == i) + int(M[k, i])) for k in range(n)] for i in range(n)]      # p_i = column i of I + M
    order = [v] + [j for j in range(n) if j != v]; V = [P[i] for i in order]
    m = n - 1
    O = cc(V); R2 = dt(sub(V[0], O), sub(V[0], O)); G = mul(Fr(1, n), [sum(c) for c in zip(*V)])
    Fc = V[1:]; Op = cc(Fc); Gp = mul(Fr(1, m), [sum(c) for c in zip(*Fc)])
    vv, w = sub(V[0], O), sub(Op, O); s = sub(mul(m, Gp), mul(m, O))
    basis = [vv, s, w]; Gb = [[dt(x, y) for y in basis] for x in basis]
    return dict(n=n, m=m, M=M, V=V, O=O, R2=R2, G=G, Op=Op, Gp=Gp, Fc=Fc, basis=basis, Gb=Gb, regime=det3(Gb) != 0)

def full_check(A, v):
    S = setup(A, v); n, m, V, O, R2, G, Op, Gp, Fc, basis = (S[k] for k in ("n", "m", "V", "O", "R2", "G", "Op", "Gp", "Fc", "basis"))
    rp2 = dt(sub(Fc[0], Op), sub(Fc[0], Op))
    dV = sub(G, V[0]); mu = -2 * dt(sub(V[0], O), dV) / dt(dV, dV); Apt = add(V[0], mul(mu, dV))
    w = sub(Op, O); B = add(O, mul(R2 / dt(w, w), w)); Cs = add(mul(m + 1, G), mul(-m, V[0]))
    g = sub(G, O); D = add(O, mul(m + 1, g)); E = add(O, mul(R2 / dt(g, g), g)); F = add(O, mul(m + 1, w))
    H = add(Op, mul(m, sub(Gp, Op))); gg = sub(Gp, Op); X = add(Op, mul(rp2 / dt(gg, gg), gg))
    rows, rhs = [], []
    for P_, c in ((O, m + 1), (G, m), (Gp, m - 1), (Op, m)):
        y = sub(P_, O); rows.append([-2 * dt(b, y) for b in basis] + [Fr(1)]); rhs.append(c * (R2 - dt(y, y)) - dt(y, y))
    sol = solve(rows, rhs); al, kap = sol[:3], sol[3]
    def pw(P_):
        y = sub(P_, O); return dt(y, y) - 2 * sum(a * dt(b, y) for a, b in zip(al, basis)) + kap
    on8 = all(pw(P_) == 0 for P_ in (Apt, B, Cs, D, E, F, H, X))
    # s from its definition
    Mm = S["M"]; Gm = [[Fr(int(i == j) + sum(int(Mm[k, i]) * int(Mm[k, j]) for k in range(n))) for j in range(n)] for i in range(n)]
    s = sum(solve(Gm, [Fr(1)] * n))
    okO = pw(O) == n * R2 == n * (n - 1 / s)
    # Delta_v by brute-force enumeration of directed 2-paths j->k->v and v->k->j
    Delta = sum(1 for j, k in itertools.permutations([t for t in range(n) if t != v], 2) if A[j, k] and A[k, v]) + \
            sum(1 for k, j in itertools.permutations([t for t in range(n) if t != v], 2) if A[v, k] and A[k, j])
    okv = pw(V[0]) == 4 * (n - 1 + Delta) == sum(dt(sub(V[0], q), sub(V[0], q)) for q in Fc)
    return on8, okO, okv

tot_checked = tot_ok = 0
for n in (4, 5, 6, 7):
    A_all = classes(n); noneq = 0; noreg = []; inreg_pairs = 0
    for b, A in enumerate(A_all):
        M = A - A.T; x = M.sum(1); y = M @ x
        if (y == y[0]).all(): continue                       # equicentric: G = O, outside the regime
        noneq += 1
        regs = [setup(A, v)["regime"] for v in range(n)]
        inreg_pairs += sum(regs)
        if not any(regs): noreg.append((b, tuple(sorted(A.sum(1).tolist()))))
        # exact transfer check on every in-regime vertex for n <= 5, first in-regime vertex for n = 6, 7 (every 4th class)
        if n <= 5 or b % 4 == 0:
            for v in range(n):
                if regs[v]:
                    r = full_check(A, v); tot_checked += 1; tot_ok += all(r)
                    if n > 5: break
    print(f"n={n}: {noneq} non-equicentric classes; in-regime (T,v) pairs {inreg_pairs}; non-equicentric classes with NO vertex in the regime: {len(noreg)}"
          + (f" e.g. score sequences {[s for _, s in noreg[:6]]}" if noreg else ""))
print(f"(6) exact transfer: {tot_ok}/{tot_checked} in-regime (T, v) pairs with all 8 labels on S2, Pow(O_T) = n r^2 = n(n - 1/s), Pow(p_v) = 4(n - 1 + Delta_v) = sum |p_v - p_j|^2")
