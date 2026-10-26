# Quadratic reciprocity and the 4-vertex census of tournaments (owner prompt 2026-10-25: "quadratic reciprocity just
# describes that a tournament on 4 vertices has 4 isomorphism classes").  Precise statements checked here:
#  (L) For a regular tournament on n = 3 (mod 4) vertices, TT4 (transitive 4-sets) >= n(n-1)(n-3)(n-7)/64, with
#      equality iff the tournament is doubly regular.  (Each TT4 has a unique source v and the other three form a
#      transitive triple in N+(v); t(N+(v)) = sum_w C(d_w, 2) is minimised exactly by regular out-neighbourhoods.)
#  (G) For p = 3 (mod 4) the Paley skew matrix has eigenvalues chi(k) g (k != 0), so S^2 = J - pI  <=>  g^2 = -p:
#      the doubly regular 4-census of P_p is the Gauss-sum identity g^2 = (-1/p) p behind the Gauss-sum proof of QR.
#  (Q) With g^2 = p* and the Frobenius congruence g^q = chi_p(q) g (mod q), (p*/q) = (q/p): checked for all odd primes
#      p, q < 200, and the sign rule (p/q)(q/p) = -1 iff p = q = 3 (mod 4), i.e. iff both carry Paley tournaments.
# mac-mini-2026-10-25-pyth.
import itertools, math, random, subprocess, sys
import numpy as np
from math import comb

def tt4(A):
    n = A.shape[0]; tot = 0
    for v in range(n):
        out = [w for w in range(n) if A[v, w]]
        sub = A[np.ix_(out, out)]
        d = sub.sum(axis=1)
        tot += int(sum(comb(int(k), 2) for k in d))
    return tot
def bound(n): return n * (n - 1) * (n - 3) * (n - 7) // 64
def doubly_regular(A):
    n = A.shape[0]; M = A @ A.T
    off = M[~np.eye(n, dtype=bool)]
    return bool(np.all(off == (n - 3) // 4))
out = []
# (L) on every regular tournament of order 7 and 11 (gentourng), and on random regular 15-, 19-tournaments
for n in (7, 11):
    m = (n - 1) // 2
    codes = subprocess.run(['gentourng', '-q', f'-d{m}', f'-D{m}', str(n)], capture_output=True, text=True).stdout.split()
    iu = np.triu_indices(n, 1); viol = 0; eq_drt = 0; eq_non = 0; drt_non_eq = 0; vals = []
    for c in codes:
        A = np.zeros((n, n), dtype=int); b = np.array([int(ch) for ch in c]); A[iu[0], iu[1]] = b; A[iu[1], iu[0]] = 1 - b
        t = tt4(A); dr = doubly_regular(A); vals.append(t)
        if t < bound(n): viol += 1
        if t == bound(n) and dr: eq_drt += 1
        if t == bound(n) and not dr: eq_non += 1
        if dr and t != bound(n): drt_non_eq += 1
    out.append(f"(L) n={n}: {len(codes)} regular classes; bound {bound(n)}; TT4 range [{min(vals)}, {max(vals)}]; below bound {viol}; "
               f"equality & doubly regular {eq_drt}; equality & not DR {eq_non}; DR but not equal {drt_non_eq}")
rng = random.Random(7)
for n in (15, 19, 23):
    d = (n - 1) // 2
    A = np.array([[1 if 0 < (j - i) % n <= d else 0 for j in range(n)] for i in range(n)])
    viol = 0; mn = None
    for s in range(300):
        for _ in range(20 * n * n):
            u, v, w = rng.randrange(n), rng.randrange(n), rng.randrange(n)
            if len({u, v, w}) == 3 and A[u, v] and A[v, w] and A[w, u]: A[u, v] = A[v, w] = A[w, u] = 0; A[v, u] = A[w, v] = A[u, w] = 1
        t = tt4(A); mn = t if mn is None else min(mn, t)
        if t < bound(n) or (t == bound(n)) != doubly_regular(A): viol += 1
    out.append(f"(L) n={n}: 300 random regular tournaments; min TT4 {mn} vs bound {bound(n)}; violations of the equivalence {viol}")
# Paley and the doubly regular codes of the repo: TT4 equals the bound
for q in (7, 11, 19, 23, 31, 43):
    sq = {(x * x) % q for x in range(1, q)}
    A = np.array([[1 if (j - i) % q in sq else 0 for j in range(q)] for i in range(q)])
    out.append(f"(L) Paley P_{q}: TT4 = {tt4(A)}, bound {bound(q)}, doubly regular {doubly_regular(A)}")
# (G) Gauss sums: g^2 = -p for p = 3 (mod 4), and the eigenvalues of S
for p in (7, 11, 19, 23, 31, 43, 47):
    chi = lambda x: 0 if x % p == 0 else (1 if pow(x, (p - 1) // 2, p) == 1 else -1)
    g = sum(chi(x) * np.exp(2j * np.pi * x / p) for x in range(p))
    S = np.array([[chi(j - i) for j in range(p)] for i in range(p)], dtype=float)
    ev = np.linalg.eigvals(S)
    out.append(f"(G) p={p}: g^2 = {g * g:.6f} (-p = {-p}); S^2 = J - pI: {np.allclose(S @ S, np.ones((p, p)) - p * np.eye(p))}; "
               f"|eigenvalues| = sqrt(p): {np.allclose(np.sort(np.abs(ev))[1:], np.sqrt(p))}")
# (Q) reciprocity via Gauss sums and the Paley sign rule
from sympy import primerange, legendre_symbol
bad = 0; bad_sign = 0; cnt = 0
for p in primerange(3, 200):
    ps = p if p % 4 == 1 else -p
    for q in primerange(3, 200):
        if q == p: continue
        cnt += 1
        if legendre_symbol(ps % q, q) != legendre_symbol(q % p, p): bad += 1
        prod = legendre_symbol(p % q, q) * legendre_symbol(q % p, p)
        if (prod == -1) != (p % 4 == 3 and q % 4 == 3): bad_sign += 1
out.append(f"(Q) odd primes p != q < 200: (p*/q) = (q/p) fails {bad} of {cnt}; '(p/q)(q/p) = -1 iff both = 3 mod 4 (both Paley tournaments)' fails {bad_sign}")
print('\n'.join(out))
