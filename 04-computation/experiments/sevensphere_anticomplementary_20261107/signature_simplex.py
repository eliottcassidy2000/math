#!/usr/bin/env python3
# THM-4911 (mac-mini-2026-11-07-sevensphere): the signature simplex of a tournament, p_i = (I+M)e_i in {+-1}^n.
# Checks, over all isomorphism classes from nauty gentourng (GENTOURNG env var or /opt/homebrew/bin/gentourng):
#  (1) |p_i|^2 = n and |p_i - p_j|^2 = 4(1 + delta_ij), delta_ij = #vertices beating exactly one of i, j   (n <= 8)
#  (2) s(T) = u^T (I+M^TM)^(-1) u = 1/dist(0, aff{p_i})^2, circumradius^2 = n - 1/s                           (n <= 8)
#  (3a) odd n: equicentric (M^2 u parallel to u) <=> regular                                                   (n <= 9)
#  (3b) non-strong & equicentric <=> exactly two strong components, both regular                               (n <= 10)
#  (3c) almost-regular (scores n/2-1, n/2) & equicentric <=> halves regular + biregular cross; iff n = 2 mod 4  (n <= 10)
#  (3d) even n: -M^2 u = lambda u with lambda = n - 1 (mod 4)                                                   (n <= 10)
#  (4) equilateral <=> I+M^TM = (n-c)I + cJ with c in {-1 (doubly regular), 0 (I+M skew-Hadamard)}             (n <= 8)
#  (5) QR7: signature vectors are -1 exactly on the Fano lines i + {1,2,4}; pairwise Hamming distance 4
#  (6) the eight-point sphere (THM-4910) of the signature simplex at a vertex v, in exact rationals:
#      Pow(O_T) = n r_T^2 = n(n - 1/s), Pow(p_v) = 4(n - 1 + Delta_v), Delta_v = (#2-paths into v) + (#2-paths out of v)
# Output: signature_simplex.out
import itertools, os, subprocess, sys, time
from fractions import Fraction as Fr
import numpy as np

GENTOURNG = os.environ.get("GENTOURNG", "/opt/homebrew/bin/gentourng")
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10

def classes(n):
    p = subprocess.Popen([GENTOURNG, "-q", str(n)], stdout=subprocess.PIPE, text=True)
    iu = np.triu_indices(n, 1)
    for line in p.stdout:
        b = np.frombuffer(line.strip().encode(), dtype=np.uint8) - 48
        A = np.zeros((n, n), dtype=np.int64); A[iu] = b; A[(iu[1], iu[0])] = 1 - b
        yield A
    p.wait()

def strong_parts(A):
    n = len(A); R = (A | np.eye(n, dtype=np.int64)).astype(bool)
    for k in range(n): R = R | (R[:, [k]] & R[[k], :])
    comps, seen = [], set()
    for i in range(n):
        if i in seen: continue
        c = [j for j in range(n) if R[i, j] and R[j, i]]; comps.append(c); seen |= set(c)
    return comps

t0 = time.time()
viol = {k: 0 for k in ["(1)", "(2)", "(3a)", "(3b)", "(3c)", "(3d)", "(4)"]}
census = {}
for n in range(3, NMAX + 1):
    eqc = nonstrong = almost = 0; lam_set = set(); equil = []
    for A in classes(n):
        M = A - A.T; u = np.ones(n, dtype=np.int64); x = M @ u; y = M @ x
        P = np.eye(n, dtype=np.int64) + M
        if n <= 8:
            Gm = P.T @ P
            if not (np.diag(Gm) == n).all(): viol["(1)"] += 1
            for i, j in itertools.combinations(range(n), 2):
                delta = int(np.sum((A[:, i] != A[:, j])[[k for k in range(n) if k not in (i, j)]]))
                if int(((P[:, i] - P[:, j])**2).sum()) != 4 * (1 + delta): viol["(1)"] += 1
            Gf = Gm.astype(float); w = np.linalg.solve(Gf, np.ones(n)); s = w.sum()
            Oc = P @ (w / s)
            if abs((Oc @ Oc) * s - 1) > 1e-9 or any(abs(((P[:, i] - Oc)**2).sum() - (n - 1/s)) > 1e-8 for i in range(n)):
                viol["(2)"] += 1
            off = {int(Gm[i, j]) for i, j in itertools.combinations(range(n), 2)}
            if len(off) == 1:
                c = off.pop(); regular = (x == 0).all()
                common_out = {int((A[i] & A[j]).sum()) for i, j in itertools.combinations(range(n), 2)}
                ok = (c == -1 and regular and len(common_out) == 1) or (c == 0 and (Gm == n * np.eye(n, dtype=np.int64)).all())
                if not ok: viol["(4)"] += 1
                equil.append(c)
        equic = bool((y == y[0]).all())
        regular = bool((x == 0).all())
        if n % 2 == 1 and equic != regular: viol["(3a)"] += 1
        if equic:
            eqc += 1
            lam = -int(y[0])
            if n % 2 == 0:
                lam_set.add(lam)
                if (lam - (n - 1)) % 4 != 0: viol["(3d)"] += 1
            comps = strong_parts(A)
            if len(comps) > 1:
                nonstrong += 1
                if len(comps) != 2 or any(not (x_sub == 0).all() for x_sub in
                                          [(A[np.ix_(c, c)] - A[np.ix_(c, c)].T) @ np.ones(len(c), dtype=np.int64) for c in comps]):
                    viol["(3b)"] += 1
            sc = A.sum(1)
            if n % 2 == 0 and set(sc.tolist()) == {n // 2 - 1, n // 2}:
                almost += 1
                S = [i for i in range(n) if sc[i] == n // 2]; Sb = [i for i in range(n) if sc[i] == n // 2 - 1]
                k = n // 2
                okS = len(S) == k and (A[np.ix_(S, S)].sum(1) == (k - 1) // 2).all() and (A[np.ix_(Sb, Sb)].sum(1) == (k - 1) // 2).all()
                okX = (A[np.ix_(S, Sb)].sum(1) == (k + 1) // 2).all() and (A[np.ix_(Sb, S)].sum(1) == (k - 1) // 2).all()
                if not (okS and okX and n % 4 == 2): viol["(3c)"] += 1
        elif n % 2 == 0 and set(A.sum(1).tolist()) == {n // 2 - 1, n // 2}:
            # almost-regular but not equicentric: check the structure is NOT present (characterisation)
            sc = A.sum(1); S = [i for i in range(n) if sc[i] == n // 2]; Sb = [i for i in range(n) if sc[i] == n // 2 - 1]; k = n // 2
            if len(S) == k and k % 2 == 1 and (A[np.ix_(S, S)].sum(1) == (k - 1) // 2).all() and (A[np.ix_(Sb, Sb)].sum(1) == (k - 1) // 2).all() \
               and (A[np.ix_(S, Sb)].sum(1) == (k + 1) // 2).all():
                viol["(3c)"] += 1
    census[n] = (eqc, nonstrong, almost, sorted(lam_set), equil)
    print(f"n={n}: equicentric {eqc} (non-strong {nonstrong}, almost-regular {almost}); even-n lambdas {sorted(lam_set)}; "
          f"equilateral c-values {equil if n <= 8 else 'n/a'}  [{time.time() - t0:.0f}s]", flush=True)
print("violations:", viol)

# (5) QR7 and the Fano lines
n = 7; QR = {1, 2, 4}
A = np.array([[1 if (j - i) % 7 in QR else 0 for j in range(7)] for i in range(7)])
M = A - A.T; P = np.eye(7, dtype=int) + M
lines = [{(i + q) % 7 for q in QR} for i in range(7)]
ok5 = all({k for k in range(7) if P[k, i] == -1} == lines[i] for i in range(7))
ham = {int((P[:, i] != P[:, j]).sum()) for i, j in itertools.combinations(range(7), 2)}
print(f"(5) QR7 signature vectors are -1 exactly on the Fano lines i+{{1,2,4}}: {ok5}; pairwise Hamming distances {ham}")

# (6) eight-point sphere of the signature simplex, exact rationals
def sub(a, b): return [p - q for p, q in zip(a, b)]
def add(a, b): return [p + q for p, q in zip(a, b)]
def mul(c, a): return [c * p for p in a]
def dt(a, b): return sum(p * q for p, q in zip(a, b))
def solve(Mx, b):
    m = len(Mx); Aug = [list(map(Fr, r)) + [Fr(bb)] for r, bb in zip(Mx, b)]
    for c in range(m):
        piv = next(r for r in range(c, m) if Aug[r][c] != 0); Aug[c], Aug[piv] = Aug[piv], Aug[c]
        for r in range(m):
            if r != c and Aug[r][c] != 0:
                f = Aug[r][c] / Aug[c][c]; Aug[r] = [p - f * q for p, q in zip(Aug[r], Aug[c])]
    return [Aug[i][m] / Aug[i][i] for i in range(m)]
def circumcenter(Pts):
    p0 = Pts[0]; Dm = [sub(p, p0) for p in Pts[1:]]
    lam = solve([[dt(a, b) for b in Dm] for a in Dm], [dt(a, a) / 2 for a in Dm])
    return add(p0, [sum(l * d[j] for l, d in zip(lam, Dm)) for j in range(len(p0))])
import random
random.seed(7)
checked = good = 0
for n in (5, 6, 7, 8):
    tries = 0
    while tries < 6:
        A = np.zeros((n, n), dtype=int)
        for i, j in itertools.combinations(range(n), 2):
            if random.random() < 0.5: A[i, j] = 1
            else: A[j, i] = 1
        M = A - A.T
        if ((M @ (M @ np.ones(n, dtype=int))) == (M @ (M @ np.ones(n, dtype=int)))[0]).all():
            continue                                          # equicentric: outside the regime
        V = [[Fr(int(e)) for e in (np.eye(n, dtype=int) + M)[:, i]] for i in range(n)]
        v = 0; O = circumcenter(V); R2 = dt(sub(V[0], O), sub(V[0], O)); m = n - 1
        G = mul(Fr(1, n), [sum(c) for c in zip(*V)]); Fc = V[1:]
        Op = circumcenter(Fc); Gp = mul(Fr(1, m), [sum(c) for c in zip(*Fc)])
        basis = [sub(V[0], O), sub(G, O), sub(Op, O)]
        Gb = [[dt(a, b) for b in basis] for a in basis]
        det = (Gb[0][0] * (Gb[1][1] * Gb[2][2] - Gb[1][2] * Gb[2][1]) - Gb[0][1] * (Gb[1][0] * Gb[2][2] - Gb[1][2] * Gb[2][0])
               + Gb[0][2] * (Gb[1][0] * Gb[2][1] - Gb[1][1] * Gb[2][0]))
        if det == 0 or G == O or Gp == Op:
            continue                                          # outside the regime
        tries += 1; checked += 1
        d = sub(G, V[0]); mu = -2 * dt(sub(V[0], O), d) / dt(d, d); Apt = add(V[0], mul(mu, d))
        w = sub(Op, O); B = add(O, mul(R2 / dt(w, w), w)); Cs = add(V[0], mul(m + 1, sub(G, V[0])))
        D = add(O, mul(m + 1, sub(G, O))); g = sub(G, O); E = add(O, mul(R2 / dt(g, g), g))
        F = add(O, mul(m + 1, sub(Op, O))); H = add(Op, mul(m, sub(Gp, Op)))
        rp2 = dt(sub(Fc[0], Op), sub(Fc[0], Op)); gg = sub(Gp, Op); X = add(Op, mul(rp2 / dt(gg, gg), gg))
        def coords(Pt):
            return solve(Gb, [dt(sub(Pt, O), b) for b in basis])
        rows, rhs = [], []
        for Pt in (D, F, Cs, H):
            c = coords(Pt); x2 = sum(c[i] * c[j] * Gb[i][j] for i in range(3) for j in range(3))
            rows.append([-2 * c[0], -2 * c[1], -2 * c[2], 1]); rhs.append(-x2)
        sol = solve(rows, rhs); uu, kap = sol[:3], sol[3]
        def pw(Pt):
            c = coords(Pt); return sum(c[i] * c[j] * Gb[i][j] for i in range(3) for j in range(3)) - 2 * sum(c[i] * uu[i] for i in range(3)) + kap
        on8 = all(pw(Pt) == 0 for Pt in (Apt, B, Cs, D, E, F, H, X))
        s = 1 / dt(sub([Fr(0)] * n, O), sub([Fr(0)] * n, O))    # s = 1/|O_T|^2
        Delta = int(sum(A[k, v] * A[:, k].sum() for k in range(n)) + sum(A[v, k] * A[k, :].sum() for k in range(n)))
        okO = pw(O) == n * R2 and R2 == n - 1 / s
        okv = pw(V[0]) == 4 * (n - 1 + Delta)
        good += on8 and okO and okv
print(f"(6) eight-point sphere of the signature simplex (exact): {good}/{checked} tournaments (n = 5..8) with all 8 labels on S2, "
      f"Pow(O_T) = n r^2 = n(n - 1/s), Pow(p_v) = 4(n - 1 + Delta_v)")
print("all done")
