#!/usr/bin/env python3
"""audit H, task F: THM-4609 statement 5 (standard-form hierarchy) and the non-translation rows of statement 6.

Independent configuration on Z/d (d prime): unit positions Z (vector 0), the other rho positions carry independent unit
vectors.  For a coupling pi(j) = a j + b, D = sum_j (v_pi(j) - v_j)(v_pi(j) - v_j)^T.  Q = I balances D iff
lambda_max(D) < tr(D)/2.  Checked for EVERY unit set Z (no translation WLOG), every rank 3..d-1, every coupling of three
classes: translations (a = 1), the largest odd-order subgroup, all units.  Exact decision at ties (integer Bareiss
leading minors of (tr/2) I - D).  Also checks the structural facts used in the proof: lambda_max <= 4, = 4 iff an even
cycle lies inside the non-unit positions, tr D = 2 #(moved non-unit positions), at most one fixed point for a != 1.
Then: squares-group forms (exact), the d = 7 squares-group orbit {0,1,3} (min coupling rank), the statement-6 rows.
Usage: python3 g_hierarchy.py [d ...]"""
import itertools, math, sys, time
from fractions import Fraction as Fr
import numpy as np

def bareiss_pd(S):
    """exact: all leading principal minors of the integer matrix S positive (Sylvester)."""
    n = len(S)
    M = [row[:] for row in S]
    prev = 1
    for k in range(n):
        if M[k][k] <= 0:
            # leading minor k+1 equals M[k][k] in Bareiss (fraction-free) elimination
            return False
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return True

def Dmat(d, Z, a, b):
    nonunit = [p for p in range(d) if p not in Z]
    idx = {p: t for t, p in enumerate(nonunit)}
    n = len(nonunit)
    D = np.zeros((n, n), dtype=np.int64)
    for j in range(d):
        w = np.zeros(n, dtype=np.int64)
        pj = (a * j + b) % d
        if pj in idx: w[idx[pj]] += 1
        if j in idx: w[idx[j]] -= 1
        D += np.outer(w, w)
    return D, nonunit

def cycles(d, a, b):
    seen = set(); cyc = []
    for s in range(d):
        if s in seen: continue
        c = [s]; seen.add(s); x = (a * s + b) % d
        while x != s:
            c.append(x); seen.add(x); x = (a * x + b) % d
        cyc.append(c)
    return cyc

def odd_subgroup(d):
    units = list(range(1, d))
    best = {1}
    for g in units:
        H = {1}; x = g
        while x != 1: H.add(x); x = x * g % d
        if len(H) % 2 == 1 and len(H) > len(best): best = H
    return sorted(best)

def hierarchy(d):
    t0 = time.time()
    classes = {'translations': [1], 'odd-order': odd_subgroup(d), 'all units': list(range(1, d))}
    struct_bad = 0
    res = {}
    for rho in range(3, d):
        res[rho] = {}
        for cname, G in classes.items():
            allok = True; fails = 0; ties = 0
            for Z in itertools.combinations(range(d), d - rho):
                Zs = set(Z)
                for a in G:
                    for b in range(d):
                        if (a, b) == (1, 0): continue
                        D, nonunit = Dmat(d, Zs, a, b)
                        tr = int(np.trace(D))
                        ev = np.linalg.eigvalsh(D.astype(float))
                        lmax = ev[-1]
                        # structural facts
                        cyc = cycles(d, a, b)
                        fixed = [c[0] for c in cyc if len(c) == 1]
                        moved_nonunit = sum(1 for p in nonunit if (a * p + b) % d != p)
                        even_inside = any(len(c) % 2 == 0 and all(p not in Zs for p in c) for c in cyc)
                        if not (len(fixed) <= (1 if a != 1 else 0) and tr == 2 * moved_nonunit and lmax <= 4 + 1e-9
                                and (abs(lmax - 4) < 1e-9) == even_inside):
                            struct_bad += 1
                        half = tr // 2
                        if abs(lmax - half) > 1e-7:
                            ok = lmax < half
                        else:
                            ties += 1
                            S = [[(half if i == j else 0) - int(D[i][j]) for j in range(rho)] for i in range(rho)]
                            ok = bareiss_pd(S)
                        if not ok:
                            allok = False; fails += 1
            res[rho][cname] = (allok, fails, ties)
        print(f"d = {d}, rank {rho}: Q = I balances every coupling for every unit set? " +
              "; ".join(f"{c}: {res[rho][c][0]} ({res[rho][c][1]} failing (Z, coupling) pairs, {res[rho][c][2]} exact ties)" for c in classes),
              flush=True)
    print(f"d = {d}: structural-fact violations {struct_bad}  [{time.time() - t0:.1f}s]", flush=True)

def balanced_by(Q, D):
    """exact: tr(D Q^-1) Q - 2 D  positive definite  (Q, D integer)."""
    Qf = [[Fr(x) for x in row] for row in Q]
    n = len(Q)
    # inverse
    A = [Qf[i][:] + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(i for i in range(c, n) if A[i][c] != 0); A[c], A[p] = A[p], A[c]
        pv = A[c][c]; A[c] = [x / pv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]; A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    Qi = [row[n:] for row in A]
    t = sum(Fr(int(D[i][j])) * Qi[j][i] for i in range(n) for j in range(n))
    S = [[t * Qf[i][j] - 2 * int(D[i][j]) for j in range(n)] for i in range(n)]
    den = 1
    for row in S:
        for x in row: den = den * x.denominator // math.gcd(den, x.denominator)
    return bareiss_pd([[int(x * den) for x in row] for row in S])

def squares_forms():
    tests = [(7, {0, 1, 2}, [[7, -1, -1, -1], [-1, 7, -1, -1], [-1, -1, 7, -1], [-1, -1, -1, 7]]),
             (11, set(range(7)), [[13, -1, -1, -2], [-1, 12, -1, -1], [-1, -1, 12, -1], [-2, -1, -1, 12]]),
             (11, {0, 1, 2, 3, 4, 5, 7}, [[13, -1, -1, -2], [-1, 12, -1, -1], [-1, -1, 12, -2], [-2, -1, -2, 12]])]
    for d, Z, Q in tests:
        G = sorted({x * x % d for x in range(1, d)})
        ok = all(balanced_by(Q, Dmat(d, Z, a, b)[0]) for a in G for b in range(d) if (a, b) != (1, 0) and np.any(Dmat(d, Z, a, b)[0]))
        print(f"squares group mod {d} = {G}, units {sorted(Z)}: form {Q} balances every non-identity coupling exactly: {ok}", flush=True)
    # the other AGL-orbit of rank-4 unit sets for d = 7 and the squares group
    for d in (7, 11):
        G = sorted({x * x % d for x in range(1, d)})
        seen = set()
        for Z in itertools.combinations(range(d), d - 4):
            # AGL(1,d) orbit representative
            orb = min(tuple(sorted((u * z + t) % d for z in Z)) for u in range(1, d) for t in range(d))
            if orb in seen: continue
            seen.add(orb)
            ranks = [int(np.linalg.matrix_rank(Dmat(d, set(orb), a, b)[0].astype(float))) for a in G for b in range(d) if (a, b) != (1, 0)]
            print(f"   d = {d}, squares group, rank-4 unit-set orbit {orb}: min coupling covariance rank {min(ranks)}"
                  f"{' -> no one-step form exists' if min(ranks) <= 2 else ''}", flush=True)

def row6():
    from hcore import Map
    for d, m, cls, Q in [(7, [1, 1, 1, 2, 11, 23, 29], 'squares', [[7, -1, -1, -1], [-1, 7, -1, -1], [-1, -1, 7, -1], [-1, -1, -1, 7]]),
                         (7, [1, 1, 2, 11, 23, 29, 37], 'squares', None), (7, [1, 2, 3, 5, 11, 13, 17], 'all', None),
                         (7, [1, 1, 1, 8, 15, 22, 29], 'translations', None), (7, [1, 1, 1, 1, 8, 15, 22], 'translations', 'AP'),
                         (5, [1, 1, 6, 11, 16], 'translations', 'AP')]:
        mp = Map(d, m)
        G = sorted({(x * pow(y, -1, d)) % d for x in m for y in m})
        Gset = set(G)
        # generated group
        H = {1}; todo = [1]
        while todo:
            h = todo.pop()
            for g in G:
                y = h * g % d
                if y not in H: H.add(y); todo.append(y)
        Z = {j for j in range(d) if m[j] == 1}
        others = [x for x in m if x != 1]
        from hcore import lattice_rank, primes_of, expvec
        pr = primes_of(m)
        indep = lattice_rank([expvec(x, pr) for x in others]) == len(others)
        if Q == 'AP':
            Q = None; QAP = True
        else:
            QAP = False
        rho = d - len(Z)
        if QAP:
            nonunit = [p for p in range(d) if p not in Z]
            # order (end, middle, end) for an AP p, p+c, p+2c
            for c in range(1, d):
                for p0 in nonunit:
                    if sorted([p0, (p0 + c) % d, (p0 + 2 * c) % d]) == sorted(nonunit):
                        order = [p0, (p0 + c) % d, (p0 + 2 * c) % d]
            QAPm = [[6, -2, 0], [-2, 6, -1], [0, -1, 5]]
            # build D in the (end, middle, end) order
            def Dord(a, b):
                Dm, nu = Dmat(d, Z, a, b)
                perm = [nu.index(p) for p in order]
                return Dm[np.ix_(perm, perm)]
            ok = all(balanced_by(QAPm, Dord(a, b)) for a in sorted(H) for b in range(d) if (a, b) != (1, 0))
            form = 'Q_AP (end, middle, end) ' + str(order)
        else:
            Qm = Q if Q is not None else np.eye(rho, dtype=np.int64).tolist()
            ok = all(balanced_by(Qm, Dmat(d, Z, a, b)[0]) for a in sorted(H) for b in range(d) if (a, b) != (1, 0))
            form = str(Qm)
        Lam = sum(math.log(x / d) for x in m) / d
        print(f"THM-4609 (6) row Z_{d} m = {m}: residues {[x % d for x in m]}, generated coupling group {sorted(H)}, "
              f"rank {mp.rank}, independent {indep}, prod {math.prod(m)} vs d^d {d ** d}, Lambda {Lam:+.4f}, r = {mp.r}; "
              f"form {form} balances every non-identity coupling a in G exactly: {ok}", flush=True)

if __name__ == '__main__':
    ds = [int(x) for x in sys.argv[1:]] or [7, 11]
    for d in ds:
        hierarchy(d)
    squares_forms()
    row6()
