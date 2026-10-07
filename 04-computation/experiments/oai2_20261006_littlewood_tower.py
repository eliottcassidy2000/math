#!/usr/bin/env python3
"""The Sierpinski skew-Hadamard tower (HYP-9162 / THM-4557) against the Littlewood flatness problem of openai/math #076
(and its companions), plus the cyclotomic stacking partition for #155 (mac-mini, 2026-10-06; independent re-check of the
session reader's claims).

Tower: H_2 = [[1, 1], [-1, 1]], H_(2N) = [[H, H], [-H^T, H^T]] (order 2^k, skew Hadamard: H + H^T = 2I).
Row/column polynomials r_j(z) = sum_y H[j, y] z^y, c_j(z) = sum_y H[y, j] z^y, w = z^N.
Merit factor F(x) = N^2 / (2 sum_(u>=1) C_u^2), C_u the aperiodic autocorrelations; R(x) = ||x||_4^4 / N^2 = 1 + 1/F.

  1. doubling: rows of H_(2N) are r_j (1 + w) and (r_j - 2 z^j)(1 - w); columns are c_j -+ w r_j (exact, k <= 9).
  2. closed form (proved by induction on k; checked k <= 11):
       r_i(y) = w_i(y) - 2 sum_(l : i_l = 1) [y = i mod 2^l] (-1)^(sum_(m >= l) i_m y_m),   w_i(y) = (-1)^(popcount(i & y)).
  3. Morse step x -> x (1 + eps z^m) (m = len x): R' = (3/2 + 2 eps s) R and s' = f(eps s), f(t) = (1 + 4t)/(6 + 8t),
     s = sum_(u=1)^(m-1) C_u C_(m-u) / ||x||_4^4 (exact integer identities, random and exhaustive tests).
     Consequences: F(x (1 + eps z^m)) <= 2F/(F + 1) < 2 for len x >= 2 (Cauchy-Schwarz |s| <= (1 - 1/R)/2); every ROW of
     H_(2^k), k >= 2, has F < 2 (columns are Golay-type steps c_j -+ w r_j and can reach F = 4); Walsh rows have s in [0, 1/4] and F <= 1/((3/2)^floor(k/2) - 1).
  4. tower rows: max F decreases (Thue-Morse regime); sup/sqrt N grows (NUMERICAL, k <= 13).
  5. Rudin-Shapiro switching D = diag((-1)^(sum y_m y_(m+1))): every row of D H D has sup <= (3 + 2 sqrt 2) sqrt(2N)
     (proved: Golay pieces on dyadic progressions) and hence F >= 1/(33 + 24 sqrt 2); numerically sup <= 3.0 sqrt N.
  6. no two rows of H_(2^k) (k >= 2) are Golay complementary (proved); rows + columns, 3 <= k <= 7 (exhaustive).
  7. H_16 (tower) has no +-1 vector d with |H d| = 4 entrywise (exhaustive, 2^15 vectors with d_0 = 1); Sylvester H_16 has 448.
  8. Paley/Legendre: a rotation of a Legendre sequence (any value at 0) is Barker exactly for p in {3, 7, 11} among
     p = 3 mod 4, p < 200 (all larger lengths are excluded by Turyn-Storer 1961).
  9. cyclotomic stacking partition: for prime q = 1 mod s, the order-s cyclotomic classes (0 added to C_0) each satisfy
     E - E = Z/q whenever q - 3s + 1 > (s - 1)(s - 2) sqrt q (Jacobi-sum bound on the cyclotomic numbers (t, t));
     all failures for s <= 8, q < 6000 lie inside the bound.
Run: python3 oai2_20261006_littlewood_tower.py   (about 1-2 min)
"""
import itertools, math, random
import numpy as np

OK = True


def check(cond, msg):
    global OK
    print(("  ok   " if cond else "  FAIL ") + msg, flush=True)
    OK &= bool(cond)


def tower(k):
    H = np.array([[1, 1], [-1, 1]], dtype=np.int64)
    for _ in range(k - 1):
        H = np.block([[H, H], [-H.T, H.T]])
    return H


def sylvester(k):
    H = np.array([[1]], dtype=np.int64)
    for _ in range(k):
        H = np.block([[H, H], [H, -H]])
    return H


def acorr(X):
    """exact aperiodic autocorrelations C_0..C_(N-1) of the rows of X (int64)."""
    X = np.atleast_2d(X)
    m, N = X.shape
    L = 1 << int(math.ceil(math.log2(2 * N)))
    F = np.fft.rfft(X.astype(np.float64), n=L, axis=1)
    C = np.fft.irfft(np.abs(F) ** 2, n=L, axis=1)[:, :N]
    Ci = np.rint(C).astype(np.int64)
    assert np.max(np.abs(C - Ci)) < 1e-3
    return Ci


def merit(X):
    C = acorr(X)
    N = C.shape[1]
    return N * N / (2.0 * (C[:, 1:] ** 2).sum(axis=1))


def supnorm_ub(X, over=8):
    """upper bound for max_|z|=1 |P| from M = over*N samples: ||P|| <= max_grid / cos(pi N / (2M))."""
    X = np.atleast_2d(X)
    m, N = X.shape
    M = over * N
    a = np.abs(np.fft.fft(X.astype(np.float64), n=M, axis=1))
    return a.max(axis=1) / math.cos(math.pi * N / (2 * M))


def poly_mul(a, b):
    return np.convolve(a, b)


print("1. the doubling, as polynomials")
good = True
for k in range(1, 10):
    H, H2 = tower(k), tower(k + 1)
    N = H.shape[0]
    assert np.array_equal(H + H.T, 2 * np.eye(N, dtype=np.int64)) and np.array_equal(H @ H.T, N * np.eye(N, dtype=np.int64))
    onepw = np.zeros(N + 1, dtype=np.int64); onepw[0] = 1; onepw[N] = 1
    onemw = np.zeros(N + 1, dtype=np.int64); onemw[0] = 1; onemw[N] = -1
    for j in range(N):
        r, c = H[j], H[:, j]
        e = np.zeros(N, dtype=np.int64); e[j] = 1
        good &= np.array_equal(poly_mul(r, onepw)[:2 * N], H2[j]) and np.array_equal(poly_mul(r - 2 * e, onemw)[:2 * N], H2[N + j])
        good &= np.array_equal(c, 2 * e - r)
        good &= np.array_equal(np.concatenate([c, -r]), H2[:, j]) and np.array_equal(np.concatenate([c, r]), H2[:, N + j])
check(good, "rows of H_(2N): r_j (1 + w), (r_j - 2 z^j)(1 - w); columns c_j -+ w r_j; c_j = 2 z^j - r_j (k = 1..9)")

print("2. closed form of the tower rows")


def closed_form(k):
    N = 1 << k
    y = np.arange(N)
    pc = np.vectorize(lambda v: bin(v).count("1"))
    R = np.zeros((N, N), dtype=np.int64)
    for i in range(N):
        row = (-1) ** pc(i & y)
        for l in range(k):
            if (i >> l) & 1:
                mask = (y % (1 << l)) == (i % (1 << l))
                hi = (i >> l) << l
                row = row - 2 * mask * (-1) ** pc(hi & y)
        R[i] = row
    return R


check(all(np.array_equal(closed_form(k), tower(k)) for k in range(1, 12)),
      "r_i(y) = w_i(y) - 2 sum_(l: i_l = 1) [y = i mod 2^l] (-1)^(sum_(m>=l) i_m y_m) for every row, k = 1..11")

print("3. the Morse-step transfer identity")


def Rs(x):
    C = acorr(np.array(x))[0]
    m = len(x)
    n4 = int(C[0] ** 2 + 2 * (C[1:] ** 2).sum())
    S = int(sum(C[u] * C[m - u] for u in range(1, m)))
    return n4, S


rng = random.Random(7)
good = True
tests = [list(x) for m in range(1, 9) for x in itertools.product((1, -1), repeat=m)] + \
        [[rng.choice((1, -1)) for _ in range(rng.randint(9, 300))] for _ in range(400)]
for x in tests:
    m = len(x)
    n4, S = Rs(x)
    for eps in (1, -1):
        y = x + [eps * v for v in x]
        n4y, Sy = Rs(y)
        good &= (n4y == 6 * n4 + 8 * eps * S) and (Sy == n4 + 4 * eps * S)
check(good, f"||y||_4^4 = 6||x||_4^4 + 8 eps S and S_y = ||x||_4^4 + 4 eps S for y = (x, eps x): {len(tests)} sequences "
            "(all of length <= 8, 400 random up to 300) [i.e. R' = (3/2 + 2 eps s) R, s' = (1 + 4 eps s)/(6 + 8 eps s)]")
good = True
for x in tests[:600]:
    if len(x) < 2:      # length 1: F = inf, and (x, eps x) has F = 2 exactly
        continue
    F = merit(np.array(x))[0]
    for eps in (1, -1):
        Fy = merit(np.array(x + [eps * v for v in x]))[0]
        good &= Fy < 2 and Fy <= 2 * F / (F + 1) + 1e-12
check(good, "F(x (1 + eps z^m)) <= 2F/(F+1) < 2 for every x of length >= 2 among the same sequences (Cauchy-Schwarz)")
check(all(merit(tower(k)).max() < 2 for k in range(2, 12)), "every row of H_(2^k) has F < 2, k = 2..11 (rows only: "
      "columns are Golay-type steps, and columns of H_4, H_8 reach F = 4)")
good = True
for k in range(2, 13):
    W = sylvester(k)
    Fm = merit(W[1:]).max()       # row 0 is all ones (F tiny); bound holds for all rows anyway
    good &= merit(W).max() <= 1.0 / (1.5 ** (k // 2) - 1) + 1e-12
check(good, "Walsh (Sylvester) rows: F <= 1/((3/2)^floor(k/2) - 1), k = 2..12")

print("4. tower rows: the Thue-Morse regime (NUMERICAL)")
rows = []
for k in range(8, 14):
    H = tower(k)
    N = H.shape[0]
    Fr = np.concatenate([merit(H[s:s + 512]) for s in range(0, N, 512)])
    sup = np.concatenate([supnorm_ub(H[s:s + 256], 4) for s in range(0, N, 256)])
    rows.append((k, N, Fr.max(), Fr.mean(), sup.min() / math.sqrt(N)))
    print(f"   k={k:2d} N={N:5d}: max F {Fr.max():.4f}, mean F {Fr.mean():.4f}, min over rows of sup/sqrtN <= {sup.min()/math.sqrt(N):.2f}")
check(all(rows[i + 1][2] < rows[i][2] for i in range(len(rows) - 1)), "max row F strictly decreasing for k = 8..13")

print("5. Rudin-Shapiro switching flattens every row")
bound = (3 + 2 * math.sqrt(2)) * math.sqrt(2)
worst = []
good = True
for k in range(1, 13):
    N = 1 << k
    y = np.arange(N)
    rs = np.array([(-1) ** sum(((v >> m) & 1) * ((v >> (m + 1)) & 1) for m in range(k)) for v in range(N)], dtype=np.int64)
    H = tower(k)
    G = rs[:, None] * H * rs[None, :]
    assert np.array_equal(G + G.T, 2 * np.eye(N, dtype=np.int64))
    sup = np.concatenate([supnorm_ub(G[s:s + 256], 8) for s in range(0, N, 256)]) / math.sqrt(N)
    F = np.concatenate([merit(G[s:s + 512]) for s in range(0, N, 512)]) if N > 1 else np.array([math.inf])
    worst.append((k, sup.max(), F.min(), F.mean(), F.max()))
    good &= sup.max() <= bound and F.min() >= 1 / (33 + 24 * math.sqrt(2))
for w in worst[2:]:
    print(f"   k={w[0]:2d}: max over rows sup/sqrtN <= {w[1]:.3f}; F in [{w[2]:.4f}, {w[4]:.4f}], mean {w[3]:.4f}")
check(good, f"every row of D H D has sup/sqrtN <= (3 + 2 sqrt2) sqrt2 = {bound:.4f} and F >= 1/(33 + 24 sqrt2) = {1/(33+24*math.sqrt(2)):.4f}, k <= 12")

print("6. Golay pairs")


def is_golay(a, b):
    Ca, Cb = acorr(a)[0], acorr(b)[0]
    return np.all(Ca[1:] + Cb[1:] == 0)


good = True
for k in range(3, 8):
    H = tower(k)
    V = [H[i] for i in range(H.shape[0])] + [H[:, j].copy() for j in range(H.shape[0])]
    A = acorr(np.array(V))[:, 1:]
    # pairs with C_a(u) + C_b(u) = 0 for all u >= 1
    for i in range(len(V)):
        hit = np.all(A[i][None, :] + A == 0, axis=1)
        good &= not hit.any()
check(good, "no Golay complementary pair among the rows and columns of H_(2^k), 3 <= k <= 7")

print("7. bent switchings of order 16")
D = np.array(list(itertools.product((1, -1), repeat=15)), dtype=np.int64)
D = np.hstack([np.ones((D.shape[0], 1), dtype=np.int64), D])
nt = int(np.all(np.abs(D @ tower(4).T) == 4, axis=1).sum())
ns = int(np.all(np.abs(D @ sylvester(4).T) == 4, axis=1).sum())
mx = int(np.abs(D @ tower(4).T).min(axis=1).max())
check(nt == 0 and ns == 448, f"#d (d_0 = 1) with |H_16 d| = 4 entrywise: tower {nt}, Sylvester {ns}; best tower min |(Hd)_i| = {mx}")

print("8. Legendre rotations that are Barker sequences")


def legendre(p, v0):
    return [v0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1) for a in range(p)]


def barker(x):
    C = acorr(np.array(x))[0]
    return np.all(np.abs(C[1:]) <= 1)


hits = []
for p in [q for q in range(3, 200) if all(q % d for d in range(2, int(q ** 0.5) + 1)) and q % 4 == 3]:
    for v0 in (1, -1):
        L = legendre(p, v0)
        if any(barker(L[r:] + L[:r]) for r in range(p)):
            hits.append(p)
            break
check(hits == [3, 7, 11], f"p = 3 mod 4 < 200 with a Barker rotation of the Legendre sequence: {hits}")

print("9. cyclotomic stacking partition")


def prime_factors(n):
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    return out + ([n] if n > 1 else [])


good = True
for s in range(2, 9):
    fails = []
    for q in range(s + 1, 6000, s):
        if not all(q % d for d in range(2, int(q ** 0.5) + 1)):
            continue
        pf = prime_factors(q - 1)
        g = next(g for g in range(2, q) if all(pow(g, (q - 1) // r, q) != 1 for r in pf))
        ind = np.zeros(q, dtype=np.int64)
        v = 1
        for e in range(q - 1):
            ind[v] = e % s
            v = v * g % q
        okq = True
        for t in range(s):
            E = (ind == t).astype(np.float64)
            E[0] = 1.0 if t == 0 else 0.0
            fE = np.fft.rfft(E)
            corr = np.rint(np.fft.irfft(fE * np.conj(fE), n=q)).astype(np.int64)   # corr[x] = #{(a, b) in E^2 : a - b = x}
            okq &= bool((corr > 0).all())
        if not okq:
            fails.append(q)
            good &= not (q - 3 * s + 1 > (s - 1) * (s - 2) * math.sqrt(q))
    print(f"   s={s}: failures {fails}")
check(good, "every failure (s <= 8, q < 6000) violates q - 3s + 1 > (s-1)(s-2) sqrt q (the bound is never contradicted)")
print("ALL CHECKS PASSED" if OK else "SOME CHECK FAILED")
