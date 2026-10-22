"""Auditor: HYP-9370 small cases with independent code (forest_audit.exe: permutation-product DP, cross-checked by the
linear-forest set-partition DP on every small case).  Exact Phi for C_3, the almost-regular 4-tournament, C_5;
all labelled regular 7-tournaments and almost-regular 6-tournaments (enumerated, grouped into classes with a mass
check); the 85 almost-regular 8-classes; regular 9-tournaments as my own completions of the 85 classes.
Usage: python zeros_small.py > zeros_small.out"""
import os, sys, math, itertools
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tourn

HERE = os.path.dirname(os.path.abspath(__file__))
UP = os.path.join(HERE, '..')
w = sp.symbols('w')


def exact_phi(N, n):
    return sum(sp.Rational(N[2 * k], n ** (2 * k)) * w ** k for k in range((n - 1) // 2 + 1))


def report_exact(name, A):
    n = A.shape[0]
    code = tourn.to_code(A)
    (H, c, N2), = tourn.forest_numbers([code], n, m2=True)
    N = tourn.N_from_c(c, n)
    assert N == N2, (N, N2)
    P = exact_phi(N, n)
    roots = sp.Poly(P, w).all_roots() if sp.degree(P, w) > 0 else []
    tz = min((sp.sqrt(abs(r)) for r in roots), key=lambda x: float(x)) if roots else None
    print(f'{name}: n={n}, H={H}, c_L={c}, N_L={N}, P(w)={sp.expand(P)}, roots in w {roots}, '
          f't0 = {sp.nsimplify(tz)} = {float(tz):.6f}')


def enumerate_scores(n, target):
    """all codes of n-tournaments whose sorted score vector equals target."""
    E = n * (n - 1) // 2
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    codes = np.arange(2 ** E, dtype=np.int64)
    sc = np.zeros((n, 2 ** E), dtype=np.int16)
    for e, (i, j) in enumerate(pairs):
        bit = ((codes >> (E - 1 - e)) & 1).astype(np.int16)     # code string position e = bit E-1-e
        sc[i] += bit
        sc[j] += 1 - bit
    ss = np.sort(sc, axis=0)
    ok = np.all(ss == np.array(sorted(target), dtype=np.int16)[:, None], axis=0)
    return [format(int(x), f'0{E}b') for x in codes[ok]]


def zero_stats(name, n, codes, group=True):
    res = tourn.forest_numbers(codes, n)
    rows = []
    for code, (H, c, _) in zip(codes, res):
        N = tourn.N_from_c(c, n)
        t0, wr, kind = tourn.smallest_zero(N, n)
        rows.append((code, H, tuple(N), t0, kind))
    t = np.array([r[3] for r in rows])
    kinds = [r[4] for r in rows]
    print(f'{name}: {len(rows)} tournaments; t0 min {t.min():.4f} median {np.median(t):.4f} max {t.max():.4f}; '
          f'kinds {dict((k, kinds.count(k)) for k in set(kinds))}; any |z|<=1.1: {(t <= 1.1).any()}')
    return rows


def main():
    # exact small cases
    report_exact('C_3', tourn.circulant(3, [1]))
    ar4 = tourn.from_code(open(os.path.join(UP, 'ar4_classes.txt')).read().split()[0])
    report_exact('almost-regular 4 (file)', ar4)
    report_exact('C_5 (rotational)', tourn.circulant(5, [1, 2]))
    report_exact('P_7 (Paley)', tourn.circulant(7, [1, 2, 4]))
    report_exact('C_7(1,2,3)', tourn.circulant(7, [1, 2, 3]))
    # all labelled almost-regular 4-tournaments
    for n, target, nm in ((4, [1, 1, 2, 2], 'almost-regular 4'), (6, [2, 2, 2, 3, 3, 3], 'almost-regular 6'),
                          (5, [2] * 5, 'regular 5'), (7, [3] * 7, 'regular 7')):
        codes = enumerate_scores(n, target)
        rows = zero_stats(f'all labelled {nm} ({len(codes)})', n, codes)
        # exact classes by brute-force canonical form (canon_brute.exe), mass check per class
        import subprocess
        out = subprocess.run([os.path.join(HERE, 'canon_brute.exe'), str(n)], input=chr(10).join(codes) + chr(10),
                             capture_output=True, text=True, check=True).stdout.split(chr(10))
        can = [l.split() for l in out if l.strip()]
        cls = {}
        for (cc, aut, sc), r in zip(can, rows):
            cls.setdefault(cc, []).append((int(aut), int(sc), r))
        for cc, mem in sorted(cls.items(), key=lambda kv: kv[1][0][2][3]):
            aut, sc, r = mem[0]
            ok = len(mem) == math.factorial(n) // aut
            tv = [m[2][3] for m in mem]
            print(f'   class {cc}: H={r[1]} N={list(r[2])} |Aut|={aut} labelled {len(mem)} = {n}!/|Aut|: {ok}; '
                  f'self-converse {bool(sc)}; t0 = {r[3]:.6f} ({r[4]}); constant on class {max(tv) - min(tv) < 1e-12}')
        tcls = np.array([mem[0][2][3] for mem in cls.values()])
        print(f'   {len(cls)} classes; per-class t0 min {tcls.min():.4f} median {np.median(tcls):.4f} max {tcls.max():.4f}; '
              f'self-converse {sum(mem[0][1] for mem in cls.values())}')
    # 85 almost-regular 8-classes (file) and my completions to regular 9
    codes8 = [l.split()[0] for l in open(os.path.join(UP, 'ar8_classes.txt'))]
    H8 = [int(l.split()[1]) for l in open(os.path.join(UP, 'ar8_classes.txt'))]
    rows8 = zero_stats('almost-regular 8, all 85 classes (file)', 8, codes8)
    print(f'   H agrees with file: {sum(r[1] == h for r, h in zip(rows8, H8))}/85')
    i = int(np.argmin([r[3] for r in rows8]))
    print(f'   minimising class {rows8[i][0]}: N = {list(rows8[i][2])}, t0 = {rows8[i][3]:.6f}')
    r0 = rows8[0]
    co = [float(x) for x in tourn.phi_w_coeffs(list(r0[2]), 8)]
    print(f'   first class P_T(w) coefficients {[round(x, 5) for x in co]}  (HYP-9370 quotes 1 + 0.75w + 0.0703w^2 + 0.00049w^3)')
    reg9 = [tourn.to_code(tourn.complete_regular(tourn.from_code(c))) for c in codes8]
    rows9 = zero_stats('regular 9, my completions of the 85 classes (with repeats)', 9, reg9)
    # dedupe regular 9 by brute-force canonical form
    import subprocess
    out = subprocess.run([os.path.join(HERE, 'canon_brute.exe'), '9'], input='\n'.join(reg9) + '\n',
                         capture_output=True, text=True, check=True).stdout.split('\n')
    can = [l.split() for l in out if l.strip()]
    classes = {}
    for (cc, aut, sc), r in zip(can, rows9):
        classes.setdefault(cc, (int(aut), int(sc), r[3], r[4]))
    mass = sum(math.factorial(9) // v[0] for v in classes.values())
    t = np.array([v[2] for v in classes.values()])
    print(f'   distinct regular 9-classes {len(classes)}, mass {mass}, self-converse {sum(v[1] for v in classes.values())}; '
          f'per-class t0 min {t.min():.4f} median {np.median(t):.4f} max {t.max():.4f}; '
          f'kinds {[v[3] for v in classes.values()].count("neg")}/{len(classes)} negative real')


if __name__ == '__main__':
    main()
