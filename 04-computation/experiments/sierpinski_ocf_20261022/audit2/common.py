# Audit2 common helpers (independent of the session's scripts): tournaments, skew matrices, exact series of the
# THM-4757 one- and two-loop formulas.
from fractions import Fraction
import numpy as np

def circulant_A(n, conn):
    conn = {c % n for c in conn}
    return np.array([[1 if (j - i) % n in conn else 0 for j in range(n)] for i in range(n)], dtype=np.int64)

def paley_A(q):
    sq = {(x * x) % q for x in range(1, q)}
    return circulant_A(q, sq)

def rot_A(n):
    return circulant_A(n, range(1, (n - 1) // 2 + 1))

def code_to_A(code):
    """gentourng / drt_results upper-triangle code: '1' at (i,j), i<j, iff i -> j"""
    m = len(code); n = int(round((1 + (1 + 8 * m) ** 0.5) / 2))
    assert n * (n - 1) // 2 == m
    A = np.zeros((n, n), dtype=np.int64); k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code[k] == '1': A[i, j] = 1
            else: A[j, i] = 1
            k += 1
    return A

def skew(A):
    return (A - A.T).astype(np.int64)

def is_regular(A):
    n = A.shape[0]
    return all(int(A[i].sum()) == (n - 1) // 2 for i in range(n)) and n % 2 == 1

def int_powers(S, K):
    """exact integer powers S^0..S^K as numpy object arrays"""
    S = np.array(S, dtype=object)
    n = S.shape[0]
    P = [np.array(np.eye(n, dtype=np.int64), dtype=object)]
    for _ in range(K): P.append(P[-1].dot(S))
    return P

def X0_series(S, K):
    """exact [z^L] of (1/2) log det(I - z^2 Shat^2) - (1/2) log det(I + z^2 Shat^2) = -sum_(m odd) tr(Shat^(2m)) z^(2m)/m"""
    n = S.shape[0]
    P = int_powers(S, K)
    out = {}
    for L in range(1, K + 1):
        if L % 2 == 0 and (L // 2) % 2 == 1:
            m = L // 2
            out[L] = Fraction(-int(np.trace(P[L])), m * n ** L)
        else:
            out[L] = Fraction(0)
    return out

def X1_series(S, K, parts=False):
    """exact [z^L] (L <= K) of the THM-4757 (4) two-loop formula with Shat -> z Shat:
       X_1 = sum_v [2 G_vv H_vv - (3/2) G_vv^2] + sum_(u,v) G_uv (H_uv^2 - O_uv^2 - (2/3) G_uv^2) - d^T G d,
       G = -Shat^2 (I + Shat^2)^-1, H = Shat^2 (I - Shat^2)^-1, O = Shat (I - Shat^2)^-1, d = diag(H - G).
       Built from the power series of G, H, O directly (not from the kernel table)."""
    n = S.shape[0]
    P = int_powers(S, K)
    Z = np.zeros((n, n), dtype=object)
    G = [Z.copy() for _ in range(K + 1)]; H = [Z.copy() for _ in range(K + 1)]; O = [Z.copy() for _ in range(K + 1)]
    # numerators: coefficient of z^k is (integer matrix)/n^k
    for k in range(1, K + 1):
        if k % 2 == 0:
            G[k] = (-1) ** (k // 2) * P[k]          # -X/(1+X) = sum_(j>=1) (-1)^j X^j, X = z^2 Shat^2
            H[k] = P[k]                              # X/(1-X)
        else:
            O[k] = P[k]                              # Shat/(1 - X)
    D = [np.array([H[k][v, v] - G[k][v, v] for v in range(n)], dtype=object) for k in range(K + 1)]
    res = {}
    for L in range(1, K + 1):
        f8 = Fraction(0); th = Fraction(0); db = Fraction(0)
        for a in range(1, L):
            b = L - a
            ga = np.array([G[a][v, v] for v in range(n)], dtype=object)
            hb = np.array([H[b][v, v] for v in range(n)], dtype=object)
            gb = np.array([G[b][v, v] for v in range(n)], dtype=object)
            f8 += 2 * int((ga * hb).sum()) - Fraction(3, 2) * int((ga * gb).sum())
        for a in range(1, L - 1):
            for b in range(1, L - a):
                c = L - a - b
                if c < 1: continue
                th += int((G[a] * (H[b] * H[c] - O[b] * O[c])).sum()) - Fraction(2, 3) * int((G[a] * G[b] * G[c]).sum())
                db -= int(D[a].dot(G[b]).dot(D[c]))
        den = n ** L
        if parts: res[L] = (f8 / den, th / den, db / den, (f8 + th + db) / den)
        else: res[L] = (f8 + th + db) / den
    return res

def regular_classes(n):
    import subprocess
    m = (n - 1) // 2
    out = subprocess.run(['gentourng', '-q', f'-d{m}', f'-D{m}', str(n)], capture_output=True, text=True).stdout
    return [c for c in out.split() if c]
