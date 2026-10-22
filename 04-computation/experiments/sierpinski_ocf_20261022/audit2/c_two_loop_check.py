# Audit2 (B, C): exact comparison of the brute-force excess-0 and excess-1 cluster sums with
#   X_0 = (1/2) log det(I - z^2 Shat^2) - (1/2) log det(I + z^2 Shat^2)                    (THM-4757 (3))
#   X_1 = figure-eight + theta - d^T G d                                                    (THM-4757 (4))
# coefficient by coefficient in z, as exact rationals.  Two independent brute-force methods:
#   P = set-partition enumeration (b_partitions.py),  M = blow-up scaling m^(-e) (b_mscaling.py).
# Test matrices:
#   (i)   regular tournaments: n = 7 (all 3 classes), n = 9 (all 15), n = 11 (R_11, P_11 and 6 non-vertex-transitive
#         classes from gentourng);
#   (ii)  the dumbbell: for EVERY regular tournament d_2 = 2 diag(Shat^2) is constant and G 1 = 0, so -d^T G d has no
#         z^L term for L < 14 (first term: loops 6, 6 and bridge 2); tested at L = 14, 16 by method M on the
#         non-circulant regular 7-tournament in c_dumbbell_n7_L16.py (a regular 9-tournament at z^14 exceeds the memory
#         budget of this audit);
#   (iii) generic "Eulerian" skew matrices (entries in {0, +-1}, zero row sums, non-constant degrees): the derivation
#         of THM-4757 (4) only uses M^T = -M and M 1 = 0, and here the dumbbell already enters at z^6.
import sys, time, random
from fractions import Fraction
import numpy as np
sys.path.insert(0, '.')
from b_partitions import cluster_sums
from b_mscaling import mscaling
from common import rot_A, paley_A, skew, code_to_A, regular_classes, X0_series, X1_series, int_powers

def compare(S, label, Lp, Lm, out):
    t0 = time.time()
    K = max(Lp, Lm)
    x0 = X0_series(S, K); x1 = X1_series(S, K, parts=True)
    resP = cluster_sums(S, Lp, regular=True, emax=1, odd=False, Lmin=2, verbose=False) if Lp >= 2 else {}
    resM = mscaling(S, Lm, verbose=False) if Lm >= 2 else {}
    rows = []; ok = True
    for L in range(2, K + 1, 2):
        f8, th, db, tot = x1[L]
        cells = []
        if L in resP:
            e0 = resP[L].get(0, Fraction(0)); e1 = resP[L].get(1, Fraction(0))
            okP = (e0 == x0[L]) and (e1 == tot); ok &= okP
            cells.append('P:' + ('ok' if okP else f'FAIL({float(e0):.6e},{float(e1):.6e})'))
        if L in resM:
            e0 = resM[L].get(0, Fraction(0)); e1 = resM[L].get(1, Fraction(0))
            okM = (e0 == x0[L]) and (e1 == tot) and all(e >= 0 for e in resM[L])
            ok &= okM
            cells.append('M:' + ('ok' if okM else f'FAIL({float(e0):.6e},{float(e1):.6e})'))
        rows.append(f"    z^{L:<2d} X0 = {float(x0[L]):+.10e}  X1 = {float(tot):+.10e}  [F8 {float(f8):+.3e}, theta {float(th):+.3e}, "
                    f"dumbbell {float(db):+.3e}]  " + ' '.join(cells))
    msg = f"  {label}: " + ('ALL EXACT' if ok else 'MISMATCH') + f"  ({time.time() - t0:.1f}s)"
    print(msg); out.append(msg)
    for r in rows: print(r); out.append(r)
    sys.stdout.flush()
    return ok

def eulerian_skew(n, seed, ncycles):
    """sum of random directed cycles (vertex-disjoint within a cycle), skew part, entries clipped to {0, +-1};
    retried until every vertex has in-degree = out-degree and degrees are not all equal"""
    rng = random.Random(seed)
    while True:
        M = np.zeros((n, n), dtype=np.int64)
        good = True
        for _ in range(ncycles):
            k = rng.randint(3, n)
            cyc = rng.sample(range(n), k)
            for a, b in zip(cyc, cyc[1:] + cyc[:1]):
                if M[a, b] != 0 or M[b, a] != 0: good = False
                M[a, b] = 1; M[b, a] = -1
        if not good: continue
        if (M.sum(axis=1) == 0).all():
            deg = (M != 0).sum(axis=1)
            if len(set(deg.tolist())) > 1: return M

if __name__ == '__main__':
    out = []; allok = True
    print('(i) regular tournaments'); out.append('(i) regular tournaments')
    c7 = regular_classes(7)
    for i, c in enumerate(c7):
        allok &= compare(skew(code_to_A(c)), f'n=7 class {i} ({c})', 10, 12, out)
    c9 = regular_classes(9)
    for i, c in enumerate(c9):
        allok &= compare(skew(code_to_A(c)), f'n=9 class {i} ({c})', 8, 10, out)
    c11 = regular_classes(11)
    for i in (0, 1, 2, 100, 500, 1000):
        allok &= compare(skew(code_to_A(c11[i])), f'n=11 gentourng class #{i}', 8, 8, out)
    allok &= compare(skew(rot_A(11)), 'R_11', 8, 8, out)
    allok &= compare(skew(paley_A(11)), 'P_11', 8, 8, out)
    # (ii) the first dumbbell terms of a regular tournament (z^14, z^16): see c_dumbbell_n7_L16.py (+ .out); a regular
    # 9-tournament at z^14 needs more than the 1.5 GB budget of this audit.
    print('(iii) generic Eulerian skew matrices (dumbbell from z^6)'); out.append('(iii) Eulerian skew matrices')
    for n, seed, nc, Lm in ((7, 1, 3, 12), (8, 2, 4, 12), (9, 3, 4, 10), (9, 4, 6, 10)):
        M = eulerian_skew(n, seed, nc)
        deg = (M != 0).sum(axis=1).tolist()
        allok &= compare(M, f'Eulerian skew n={n} seed={seed} degrees {deg}', 8, Lm, out)
    print('OVERALL: ' + ('ALL EXACT' if allok else 'MISMATCHES FOUND'))
