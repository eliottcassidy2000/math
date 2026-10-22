"""Auditor helpers (independent of the S54 lane's code).  Tournaments as 0/1 numpy matrices A (A[i,j] = 1 iff i -> j);
code = upper-triangle bits row-major, '1' iff i -> j (i < j) -- the repository's convention."""
import itertools, os, subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, 'forest_audit.exe')


def from_code(code):
    m = len(code)
    n = int(round((1 + (1 + 8 * m) ** 0.5) / 2))
    A = np.zeros((n, n), dtype=int)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code[c] == '1':
                A[i, j] = 1
            else:
                A[j, i] = 1
            c += 1
    return A


def to_code(A):
    n = A.shape[0]
    return ''.join('1' if A[i, j] else '0' for i in range(n) for j in range(i + 1, n))


def circulant(m, conn):
    A = np.zeros((m, m), dtype=int)
    for i in range(m):
        for s in conn:
            A[i, (i + s) % m] = 1
    assert np.all(A + A.T + np.eye(m, dtype=int) == 1), 'not a tournament'
    return A


def lex(Aout, Ain):
    """lexicographic product Aout[Ain]: each vertex of Aout replaced by a copy of Ain."""
    a, b = Aout.shape[0], Ain.shape[0]
    A = np.zeros((a * b, a * b), dtype=int)
    for x in range(a):
        for y in range(a):
            for i in range(b):
                for j in range(b):
                    if x == y:
                        A[x * b + i, y * b + j] = Ain[i, j]
                    else:
                        A[x * b + i, y * b + j] = Aout[x, y]
    return A


def delete(A, v):
    keep = [i for i in range(A.shape[0]) if i != v]
    return A[np.ix_(keep, keep)]


def complete_regular(A):
    """unique regular completion of an almost-regular tournament of even order n (scores n/2, n/2 - 1)."""
    n = A.shape[0]
    s = A.sum(1)
    assert sorted(set(s)) == [n // 2 - 1, n // 2], s
    B = np.zeros((n + 1, n + 1), dtype=int)
    B[:n, :n] = A
    for v in range(n):
        if s[v] == n // 2 - 1:
            B[v, n] = 1          # low-score vertices beat the new vertex
        else:
            B[n, v] = 1
    assert np.all(B.sum(1) == n // 2)
    return B


def forest_numbers(codes, n, m2=False):
    """run forest_audit.exe; returns list of (H, c_list, N_list or None)."""
    inp = '\n'.join(codes) + '\n'
    args = [EXE, str(n)] + (['m2'] if m2 else [])
    out = subprocess.run(args, input=inp, capture_output=True, text=True, check=True).stdout
    res = []
    for line in out.strip().splitlines():
        parts = line.split('|')
        a = list(map(int, parts[0].split()))
        H, c = a[0], a[1:]
        N2 = list(map(int, parts[1].split())) if len(parts) > 1 else None
        res.append((H, c, N2))
    return res


def N_from_c(c, n):
    from math import factorial
    N = []
    for L, x in enumerate(c):
        f = factorial(n - L)
        assert x % f == 0, (L, x, f)
        N.append(x // f)
    return N


def phi_w_coeffs(N, n):
    """P_T(w) coefficients (ascending) with w = z^2: [w^k] = N_(2k) / n^(2k)."""
    from fractions import Fraction
    return [Fraction(N[2 * k], n ** (2 * k)) for k in range((n - 1) // 2 + 1)]


def smallest_zero(N, n):
    """returns (|z| of the smallest zero, w-root, kind) with kind in {'neg', 'pos', 'complex'}."""
    co = [float(x) for x in phi_w_coeffs(N, n)]
    while len(co) > 1 and co[-1] == 0:
        co.pop()
    if len(co) == 1:
        return float('inf'), None, 'none'
    r = np.roots(co[::-1])
    k = np.argmin(np.abs(r))
    w = r[k]
    if abs(w.imag) <= 1e-9 * max(1.0, abs(w)):
        kind = 'neg' if w.real < 0 else 'pos'
    else:
        kind = 'complex'
    return float(np.sqrt(abs(w))), w, kind


def all_automorphisms(A, anti=False):
    """count automorphisms (anti=False) or anti-automorphisms (anti=True: A[x,y] = A[p(y), p(x)]) by backtracking."""
    n = A.shape[0]
    s = A.sum(1)
    target = A.T if anti else A
    ts = target.sum(1)
    cnt = 0
    p = [-1] * n
    used = [False] * n

    def bt(i):
        nonlocal cnt
        if i == n:
            cnt += 1
            return
        for y in range(n):
            if used[y] or ts[y] != s[i]:
                continue
            ok = True
            for j in range(i):
                if A[i, j] != target[y, p[j]]:
                    ok = False
                    break
            if ok:
                p[i] = y; used[y] = True
                bt(i + 1)
                used[y] = False; p[i] = -1
    bt(0)
    return cnt
