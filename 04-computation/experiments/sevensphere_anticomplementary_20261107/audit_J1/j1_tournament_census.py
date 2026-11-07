#!/usr/bin/env python3
# audit J1 (THM-4911 (1)-(4), HYP-9500 table): independent vectorised census over all gentourng classes n <= 10 (equicentric, regular, strong via Landau, non-strong/almost-regular structure, lambda mod 4, equilateral c, delta-edges and s = 1/dist^2 for n <= 9).
import itertools, subprocess, sys, time
from fractions import Fraction as Fr
import numpy as np

GT = "/opt/homebrew/bin/gentourng"
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10
CHUNK = 200000

def chunks(n, strong_only=False):
    m = n * (n - 1) // 2
    args = [GT, "-q"] + (["-c"] if strong_only else []) + [str(n)]
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    iu = np.triu_indices(n, 1)
    while True:
        buf = p.stdout.read((m + 1) * CHUNK)
        if not buf: break
        a = np.frombuffer(buf, dtype=np.uint8).reshape(-1, m + 1)
        assert (a[:, m] == 10).all()
        bits = (a[:, :m] - 48).astype(np.int8)
        A = np.zeros((len(bits), n, n), dtype=np.int8)
        A[:, iu[0], iu[1]] = bits; A[:, iu[1], iu[0]] = 1 - bits
        yield A
    p.wait()

def strong_landau(A):
    d = np.sort(A.sum(2), axis=1).astype(np.int64); cs = np.cumsum(d, axis=1)
    n = A.shape[1]; k = np.arange(1, n)
    return (cs[:, :-1] > k * (k - 1) // 2).all(1)

def components(A):
    n = len(A); R = (A.astype(bool) | np.eye(n, dtype=bool))
    for k in range(n): R = R | (R[:, [k]] & R[[k], :])
    comps = []; seen = set()
    for i in range(n):
        if i in seen: continue
        c = [j for j in range(n) if R[i, j] and R[j, i]]; comps.append(c); seen |= set(c)
    return comps

def exact_s_and_dist(A):
    # exact: s = u^T (I + M^T M)^{-1} u ; O_T = projection of 0 onto aff{p_i} computed independently by Gram-Schmidt-free normal equations
    n = len(A); M = [[int(A[i][j]) - int(A[j][i]) for j in range(n)] for i in range(n)]
    P = [[Fr(int(i == k) + M[k][i]) for k in range(n)] for i in range(n)]       # P[i] = p_i (column i of I+M)
    Gm = [[sum(P[i][k] * P[j][k] for k in range(n)) for j in range(n)] for i in range(n)]
    def solve(Mx, b):
        m = len(Mx); Ag = [list(r) + [bb] for r, bb in zip(Mx, b)]
        for c in range(m):
            piv = next(r for r in range(c, m) if Ag[r][c] != 0); Ag[c], Ag[piv] = Ag[piv], Ag[c]
            for r in range(m):
                if r != c and Ag[r][c] != 0:
                    f = Ag[r][c] / Ag[c][c]; Ag[r] = [x - f * y for x, y in zip(Ag[r], Ag[c])]
        return [Ag[i][m] / Ag[i][i] for i in range(m)]
    s = sum(solve(Gm, [Fr(1)] * n))
    Q = [[P[j][k] - P[0][k] for k in range(n)] for j in range(1, n)]
    QQ = [[sum(a * b for a, b in zip(q1, q2)) for q2 in Q] for q1 in Q]
    coef = solve(QQ, [-sum(a * b for a, b in zip(q, P[0])) for q in Q])
    Ot = [P[0][k] + sum(c * q[k] for c, q in zip(coef, Q)) for k in range(n)]
    d2 = sum(x * x for x in Ot)
    r2 = [sum((P[i][k] - Ot[k])**2 for k in range(n)) for i in range(n)]
    x = [sum(M[i]) for i in range(n)]
    return s, d2, r2, sum(t * t for t in x)

t0 = time.time()
viol = {key: 0 for key in ("1", "2", "3a", "3b", "3c", "3d", "lam=|x|^2/n", "4", "bounds", "strong-cross")}
print("n | classes | equicentric | non-strong | almost-reg | strong | lambda (all) | lambda (strong) | equilateral c=-1, c=0")
for n in range(3, NMAX + 1):
    tot = eq = ns = ar = st = 0; lam_all = set(); lam_strong = set(); c_m1 = c_0 = 0; eq_reps = []
    for A in chunks(n):
        tot += len(A)
        M = (A - A.transpose(0, 2, 1)).astype(np.int16)
        x = M.sum(2); y = np.einsum('bij,bj->bi', M, x)
        regular = (x == 0).all(1); equic = (y == y[:, :1]).all(1)
        if n % 2 == 1: viol["3a"] += int((equic != regular).sum())
        # (3c) almost regular: x in {+-1}
        if n % 2 == 0:
            alm = (np.abs(x) == 1).all(1)
            if alm.any():
                Aa = A[alm].astype(np.int16); Sm = (x[alm] == 1).astype(np.int16); k = n // 2
                wS = np.einsum('bij,bj->bi', Aa, Sm); wSb = np.einsum('bij,bj->bi', Aa, 1 - Sm)
                inS = Sm.astype(bool)
                okS = np.where(inS, (wS == (k - 1) // 2) & (wSb == (k + 1) // 2), (wSb == (k - 1) // 2) & (wS == (k - 1) // 2))
                struct = okS.all(1) & (k % 2 == 1) & (Sm.sum(1) == k)
                viol["3c"] += int((struct != equic[alm]).sum())
        # (1) and (2) for n <= 9: delta-edges (combinatorial) and float s/dist
        if n <= 9:
            P = (np.eye(n, dtype=np.int16)[None] + M)
            ham = (P[:, :, :, None] != P[:, :, None, :]).sum(1)              # Hamming distance between columns i, j
            sep = (A[:, :, :, None] != A[:, :, None, :])                      # [b, k, i, j]: A[k,i] != A[k,j]
            idx = np.arange(n)
            mask = (idx[:, None, None] != idx[None, :, None]) & (idx[:, None, None] != idx[None, None, :])   # k not in {i,j}
            delta = (sep & mask[None]).sum(1)
            off = ~np.eye(n, dtype=bool)
            viol["1"] += int(((ham != 1 + delta) & off[None]).any((1, 2)).sum())
            viol["1"] += int(((P.astype(np.int32)**2).sum(1) != n).any(1).sum())
            Gm = np.einsum('bki,bkj->bij', P.astype(np.float64), P.astype(np.float64))
            s = np.linalg.solve(Gm, np.ones((len(A), n, 1)))[:, :, 0].sum(1)
            Pf = P.astype(np.float64); Q = Pf[:, :, 1:] - Pf[:, :, :1]
            QQ = np.einsum('bki,bkj->bij', Q, Q); rhs = -np.einsum('bki,bk->bi', Q, Pf[:, :, 0])
            coef = np.linalg.solve(QQ, rhs[:, :, None])[:, :, 0]
            Ot = Pf[:, :, 0] + np.einsum('bki,bi->bk', Q, coef)
            d2 = (Ot**2).sum(1); r2 = ((Pf - Ot[:, :, None])**2).sum(1)
            viol["2"] += int((np.abs(d2 * s - 1) > 1e-9).sum()) + int((np.abs(r2 - (n - 1 / s)[:, None]) > 1e-8).any(1).sum())
            lower = n * n / (n + (x.astype(np.float64)**2).sum(1))
            viol["bounds"] += int((s < lower - 1e-9).sum() + (s > n + 1e-9).sum())
            viol["bounds"] += int(((np.abs(s - lower) < 1e-9) != equic).sum()) + int(((np.abs(s - n) < 1e-9) != regular).sum())
        # equicentric details
        sl = strong_landau(A)
        for b in np.nonzero(equic)[0]:
            eq += 1; lam = int(-y[b, 0])
            if lam * n != int((x[b].astype(np.int64)**2).sum()): viol["lam=|x|^2/n"] += 1
            if n % 2 == 0:
                lam_all.add(lam)
                if (lam - (n - 1)) % 4: viol["3d"] += 1
            comps = components(A[b])
            if (len(comps) == 1) != bool(sl[b]): viol["strong-cross"] += 1
            if len(comps) > 1:
                ns += 1
                regs = [((A[b][np.ix_(c, c)].sum(1) * 2) == len(c) - 1).all() for c in comps]
                if not (len(comps) == 2 and all(regs) and lam == len(comps[0]) * len(comps[1])): viol["3b"] += 1
            else:
                st += 1
                if n % 2 == 0: lam_strong.add(lam)
            if n % 2 == 0 and (np.abs(x[b]) == 1).all(): ar += 1
            # (4) equilateral (equilateral implies equicentric, so it suffices to look here)
            Pm = np.eye(n, dtype=np.int64) + M[b].astype(np.int64); Gm1 = Pm.T @ Pm
            offv = {int(Gm1[i, j]) for i, j in itertools.combinations(range(n), 2)}
            if len(offv) == 1:
                c = offv.pop(); AA = A[b].astype(np.int64) @ A[b].astype(np.int64).T
                common = {int(AA[i, j]) for i, j in itertools.combinations(range(n), 2)}
                dr = (x[b] == 0).all() and len(common) == 1
                if c == -1 and dr: c_m1 += 1
                elif c == 0 and (Gm1 == n * np.eye(n, dtype=np.int64)).all() and not (x[b] == 0).all(): c_0 += 1
                else: viol["4"] += 1
            elif (x[b] == 0).all():
                AA = A[b].astype(np.int64) @ A[b].astype(np.int64).T
                if len({int(AA[i, j]) for i, j in itertools.combinations(range(n), 2)}) == 1: viol["4"] += 1   # doubly regular but not equilateral
            if n <= 7: eq_reps.append(A[b].copy())
        # (4) direct check over ALL classes for n <= 8 (not only equicentric ones)
        if n <= 8:
            P = (np.eye(n, dtype=np.int16)[None] + M).astype(np.int32)
            Gm2 = np.einsum('bki,bkj->bij', P, P); off = ~np.eye(n, dtype=bool)
            offvals = Gm2[:, off]; equil = (offvals == offvals[:, :1]).all(1)
            viol["4"] += int((equil & ~equic).sum())
            viol["4"] += int((equil & ~np.isin(offvals[:, 0], [-1, 0])).sum())
    print(f"{n} | {tot} | {eq} | {ns} | {ar} | {st} | {sorted(lam_all)} | {sorted(lam_strong)} | {c_m1}, {c_0}   [{time.time()-t0:.0f}s]", flush=True)
    if n <= 7:
        # exact (2) on every class for n <= 7
        bad = 0
        for A in chunks(n):
            for b in range(len(A)):
                s, d2, r2, xx = exact_s_and_dist(A[b])
                ok = (d2 * s == 1) and all(r == n - 1 / s for r in r2) and (Fr(n * n, n + xx) <= s <= n)
                bad += not ok
        viol["2"] += bad
        print(f"   exact (2) on all {tot} classes: {'OK' if bad == 0 else 'VIOLATIONS %d' % bad}")
    if n in (8, 9, 10) and n <= NMAX:
        # cross-check strong equicentric count with gentourng -c (strong classes only)
        cnt = 0
        for A in chunks(n, strong_only=True):
            M = (A - A.transpose(0, 2, 1)).astype(np.int16); x = M.sum(2); y = np.einsum('bij,bj->bi', M, x)
            cnt += int((y == y[:, :1]).all(1).sum())
        print(f"   gentourng -c {n}: strong equicentric classes = {cnt} (Landau/components count above: {st})  [{time.time()-t0:.0f}s]", flush=True)
        if cnt != st: viol["strong-cross"] += 1
print("violations:", viol)
