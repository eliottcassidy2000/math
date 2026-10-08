#!/usr/bin/env python3
"""Audit G (independent reimplementation; stdlib only, exact rationals).

THM-4609 statements 1 and 4:
 (a) For every prime d in a range and EVERY unit set Z (|Z| >= 1) the lag covariances of the independent configuration
     (v_j = 0 on Z, v_j = e_k on the k-th non-unit position) equal (1/d)(2I - A_b), A_b = adjacency of the lag-b cycle
     restricted to the non-unit positions; A_b is a disjoint union of paths.  Checked by building C_b from the roots.
 (b) The explicit rank-3 forms Q_AP, Q_nonAP are positive definite and balanced (Sylvester, exact) on their families.
 (c) EVERY 3-subset of non-unit positions for d = 5..31 (prime) is balanced by Q_AP (ordered end, middle, end) or
     Q_nonAP: checked directly on the actual C_b's (no appeal to the family description).
 (d) rank >= 4: Q = I balanced, checked exactly on every unit set (up to AGL(1,d)) for d = 5, 7, 11, 13.
 (e) Z empty (all d multipliers independent): shifted configuration, Q = I balanced (d = 5, 7, 11, 13).
 (f) alpha_max for the forms, and the claim 'Q = I gives equality' (is it equality or strict failure?).
"""
from fractions import Fraction as Fr
import itertools, math

def det(M):
    n = len(M); A = [row[:] for row in M]; dv = Fr(1)
    for c in range(n):
        pr = next((i for i in range(c, n) if A[i][c] != 0), None)
        if pr is None: return Fr(0)
        if pr != c: A[c], A[pr] = A[pr], A[c]; dv = -dv
        dv *= A[c][c]
        for i in range(c + 1, n):
            f = A[i][c] / A[c][c]
            if f: A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return dv

def posdef(S):
    # Sylvester: all leading principal minors > 0 (S symmetric)
    for i in range(len(S)):
        for j in range(len(S)):
            assert S[i][j] == S[j][i]
    return all(det([row[:k] for row in S[:k]]) > 0 for k in range(1, len(S) + 1))

def inv(M):
    n = len(M); A = [[Fr(x) for x in row] + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        pr = next(i for i in range(c, n) if A[i][c] != 0); A[c], A[pr] = A[pr], A[c]
        pv = A[c][c]; A[c] = [x / pv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]; A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return [row[n:] for row in A]

def trace_prod(C, Qi):
    n = len(C); return sum(C[a][b] * Qi[b][a] for a in range(n) for b in range(n))

def balanced_on(Q, Cs):
    """S_b = (1/2) tr(C Q^-1) Q - C positive definite for every nonzero C."""
    Q = [[Fr(x) for x in row] for row in Q]
    if not posdef(Q): return False
    Qi = inv(Q)
    for C in Cs:
        if all(x == 0 for row in C for x in row): continue
        t = trace_prod(C, Qi)
        S = [[t / 2 * Q[a][b] - C[a][b] for b in range(len(Q))] for a in range(len(Q))]
        if not posdef(S): return False
    return True

def config(d, nonunits):
    rho = len(nonunits); v = []
    for j in range(d):
        if j in nonunits:
            e = [Fr(0)] * rho; e[nonunits.index(j)] = Fr(1); v.append(e)
        else:
            v.append([Fr(0)] * rho)
    return v

def lagcov(v, d, b):
    rho = len(v[0]); C = [[Fr(0)] * rho for _ in range(rho)]
    for j in range(d):
        w = [x - y for x, y in zip(v[(j + b) % d], v[j])]
        for p in range(rho):
            if w[p]:
                for q in range(rho):
                    if w[q]: C[p][q] += w[p] * w[q] / d
    return C

def lag_graph(d, nonunits, b):
    rho = len(nonunits); A = [[0] * rho for _ in range(rho)]
    for a in range(rho):
        for c in range(rho):
            if a != c and ((nonunits[c] - nonunits[a]) % d in (b % d, (-b) % d)): A[a][c] = 1
    return A

def is_union_of_paths(A):
    n = len(A); deg = [sum(r) for r in A]
    if any(x > 2 for x in deg): return False
    # acyclic: edges = vertices - components
    seen = [False] * n; comps = 0
    for s in range(n):
        if seen[s]: continue
        comps += 1; st = [s]; seen[s] = True
        while st:
            x = st.pop()
            for y in range(n):
                if A[x][y] and not seen[y]: seen[y] = True; st.append(y)
    edges = sum(deg) // 2
    return edges == n - comps

def K(n, edges):
    M = [[Fr(2) if i == j else Fr(0) for j in range(n)] for i in range(n)]
    for i, j in edges: M[i][j] -= 1; M[j][i] -= 1
    return M

def eig_sym_float(M):
    # small symmetric matrices: use numpy if present, else Jacobi
    try:
        import numpy as np
        return sorted(np.linalg.eigvalsh(np.array([[float(x) for x in r] for r in M])).tolist())
    except Exception:
        raise

def alpha_max(Q, Cs):
    import numpy as np
    Qf = np.array([[float(x) for x in r] for r in Q]); w, U = np.linalg.eigh(Qf)
    Qmh = U @ np.diag(w ** -0.5) @ U.T
    best = 1e9
    for C in Cs:
        Cf = np.array([[float(x) for x in r] for r in C])
        if not Cf.any(): continue
        Sg = Qmh @ Cf @ Qmh; ev = np.linalg.eigvalsh(Sg)
        best = min(best, ev.sum() / ev.max() - 2)
    return best

def primes_upto(n):
    return [p for p in range(5, n + 1) if all(p % q for q in range(2, int(p ** .5) + 1))]

if __name__ == '__main__':
    Q_AP = [[6, -2, 0], [-2, 6, -1], [0, -1, 5]]          # basis (end, middle, end)
    Q_NAP = [[4, -1, -1], [-1, 5, 0], [-1, 0, 5]]
    I3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    AP = [K(3, [(0, 1), (1, 2)]), K(3, [(0, 2)]), K(3, [])]
    NAP = [K(3, [(0, 1)]), K(3, [(1, 2)]), K(3, [(0, 2)]), K(3, [])]
    print("== (b) explicit forms on the universal families ==")
    print("Q_AP posdef:", posdef([[Fr(x) for x in r] for r in Q_AP]), " det =", det([[Fr(x) for x in r] for r in Q_AP]))
    print("Q_nonAP posdef:", posdef([[Fr(x) for x in r] for r in Q_NAP]), " det =", det([[Fr(x) for x in r] for r in Q_NAP]))
    for name, Q, fam in (('Q_AP on AP family', Q_AP, AP), ('Q_AP on AP minus 2I', Q_AP, AP[:2]),
                         ('Q_nonAP on nonAP family', Q_NAP, NAP), ('Q_nonAP on nonAP minus 2I', Q_NAP, NAP[:3])):
        print(f"  {name}: balanced = {balanced_on(Q, fam)};  alpha_max = {alpha_max(Q, fam):.4f}")
    # the minors of S explicitly, for the record
    for name, Q, fam in (('AP', Q_AP, AP), ('nonAP', Q_NAP, NAP)):
        Qf = [[Fr(x) for x in r] for r in Q]; Qi = inv(Qf)
        for k_, C in enumerate(fam):
            t = trace_prod(C, Qi); S = [[t / 2 * Qf[a][b] - C[a][b] for b in range(3)] for a in range(3)]
            mins = [det([row[:k] for row in S[:k]]) for k in (1, 2, 3)]
            print(f"     {name} member {k_}: tr(KQ^-1) = {t}, leading minors of S = {[str(m) for m in mins]}")
    print("== (f) Q = I on the families: equality or strict failure? ==")
    for name, fam in (('AP', AP), ('nonAP', NAP)):
        for k_, C in enumerate(fam):
            ev = eig_sym_float(C); print(f"  {name} member {k_}: eigenvalues {[round(x, 4) for x in ev]}, tr/2 = {sum(ev)/2:.4f}, "
                                         f"lambda_max - tr/2 = {max(ev) - sum(ev)/2:+.4f}")
    print("== (a),(c) every 3-subset of non-unit positions, primes d <= 31: covariance formula, paths, balance ==")
    for d in primes_upto(31):
        n_ap = n_nap = 0; bad = 0; formula_bad = 0; path_bad = 0
        for P in itertools.combinations(range(d), 3):
            P = list(P); v = config(d, P)
            Cs = [lagcov(v, d, b) for b in range(1, d)]
            for b, C in zip(range(1, d), Cs):
                A = lag_graph(d, P, b)
                pred = [[(Fr(2) * (a == c) - A[a][c]) / d for c in range(3)] for a in range(3)]
                if pred != C: formula_bad += 1
                if not is_union_of_paths(A): path_bad += 1
            # AP?  find middle
            mid = None
            for k in range(3):
                o = [P[x] for x in range(3) if x != k]
                if (o[0] + o[1] - 2 * P[k]) % d == 0: mid = k
            if mid is not None:
                n_ap += 1; ends = [x for x in range(3) if x != mid]; order = [ends[0], mid, ends[1]]
                Q = Q_AP
            else:
                n_nap += 1; order = [0, 1, 2]; Q = Q_NAP
            # express Q in the configuration's basis (permute)
            Qp = [[0] * 3 for _ in range(3)]
            for a in range(3):
                for c in range(3): Qp[order[a]][order[c]] = Q[a][c]
            if not balanced_on(Qp, Cs): bad += 1
        print(f"  d = {d:2d}: 3-sets AP {n_ap}, nonAP {n_nap};  covariance-formula mismatches {formula_bad}; "
              f"non-path lag graphs {path_bad}; NOT balanced by the explicit form: {bad}", flush=True)
    print("== (d) rank >= 4 with Q = I, every unit set up to AGL(1,d) ==")
    for d in (5, 7, 11, 13):
        reps = set()
        for zs in range(1, d - 3):
            for Z in itertools.combinations(range(d), zs):
                best = None
                for a in range(1, d):
                    for c in range(d):
                        img = tuple(sorted((a * z + c) % d for z in Z))
                        if best is None or img < best: best = img
                reps.add(best)
        nbad = 0; worst = 0
        for Z in reps:
            P = [j for j in range(d) if j not in Z]; rho = len(P); v = config(d, P)
            Cs = [lagcov(v, d, b) for b in range(1, d)]
            I = [[int(i == j) for j in range(rho)] for i in range(rho)]
            ok = balanced_on(I, Cs)
            nbad += (not ok)
            for C in Cs:
                ev = eig_sym_float(C); worst = max(worst, max(ev) * d)
        print(f"  d = {d}: {len(reps)} AGL-orbits with rank >= 4; Q = I fails on {nbad}; max lambda_max(2I - A) = {worst:.4f} (< 4)", flush=True)
    print("== (e) Z empty: d independent multipliers, shifted by -v_0 ==")
    for d in (5, 7, 11, 13):
        rho = d - 1
        v = [[Fr(0)] * rho] + [[Fr(int(k == j - 1)) for k in range(rho)] for j in range(1, d)]   # v_j - v_0 in basis
        Cs = [lagcov(v, d, b) for b in range(1, d)]
        I = [[int(i == j) for j in range(rho)] for i in range(rho)]
        print(f"  d = {d}: rank {rho}, Q = I balanced: {balanced_on(I, Cs)}")
