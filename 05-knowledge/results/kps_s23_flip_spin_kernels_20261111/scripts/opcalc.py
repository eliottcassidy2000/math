"""opcalc.py - the closing operator on cover polynomials (derivation check, brute force).

Cover polynomial C_T(k) = number of ordered k-tuples of disjoint possibly-empty paths covering T.
Claim: for a pair (B, O), X any tournament, T = C(B, O; X) = (B => X => v, v -> O):
   C_T(k) = k [ u_k C_X(k+1) + r_k C_X(k) + s_k C_X(k-1) ],
   u_k = a'_k,  r_k = a_k - 2 a'_k + b'_k,  s_k = a'_k - a_k + b_k - b'_k,
where over (k+1)-tuples (L_1..L_{k+1}) of disjoint possibly-empty paths covering B:
   a_k : L_2 empty or starting in O;                a'_k: L_2 nonempty starting in O;
   b_k : (L_1 empty or ending outside O) and (L_2 empty or starting in O);
   b'_k: (L_1 empty or ending outside O) and L_2 nonempty starting in O.
3-cycle: B = one vertex, O = that vertex:  (Phi C)(k) = k [ C(k+1) + k C(k) - C(k-1) ]."""
import random, itertools
from towers import tables, level_up, rand_t, TT


def cover_tables(A, K):
    n = len(A)
    full = (1 << n) - 1
    en, st, hs = tables(A) if n else ([[]], [[]], [1])
    f = [[0] * (full + 1) for _ in range(K + 1)]
    f[0][0] = 1
    for k in range(1, K + 1):
        for S in range(full + 1):
            t = 0
            sub = S
            while True:
                t += hs[sub] * f[k - 1][S ^ sub]
                if sub == 0:
                    break
                sub = (sub - 1) & S
            f[k][S] = t
    return en, st, hs, f


def cover_poly_values(A, K):
    n = len(A)
    en, st, hs, f = cover_tables(A, K)
    return [f[k][(1 << n) - 1] for k in range(K + 1)]


def closing_coeffs(A, O, K):
    """u_k, r_k, s_k for k = 1..K-? ; needs (k+1)-tuples -> uses f up to K-1 for the other k-1 slots."""
    n = len(A)
    full = (1 << n) - 1
    en, st, hs, f = cover_tables(A, K)
    Oset = set(O)
    endI = [0] * (full + 1)
    startO = [0] * (full + 1)
    for S in range(1, full + 1):
        endI[S] = sum(en[S][v] for v in range(n) if (S >> v & 1) and v not in Oset)
        startO[S] = sum(st[S][v] for v in range(n) if (S >> v & 1) and v in Oset)
    out = {}
    for k in range(1, K + 1):
        a = ap = b = bp = 0
        for S1 in range(full + 1):
            g1_any = hs[S1]
            g1_I = 1 if S1 == 0 else endI[S1]
            rest = full ^ S1
            S2 = rest
            while True:
                rem = f[k - 1][rest ^ S2]
                if rem:
                    g2_ne = startO[S2] if S2 else 0
                    g2 = 1 if S2 == 0 else startO[S2]
                    a += g1_any * g2 * rem
                    ap += g1_any * g2_ne * rem
                    b += g1_I * g2 * rem
                    bp += g1_I * g2_ne * rem
                if S2 == 0:
                    break
                S2 = (S2 - 1) & rest
        out[k] = (ap, a - 2 * ap + bp, ap - a + b - bp)
    return out


if __name__ == "__main__":
    rng = random.Random(23)
    bad = 0
    cnt = 0
    for _ in range(120):
        n = rng.randint(1, 4)
        B = rand_t(n, rng)
        O = [i for i in range(n) if rng.random() < 0.5]
        x = rng.randint(1, 3)
        X = rand_t(x, rng)
        K = 3
        co = closing_coeffs(B, O, K)
        CX = cover_poly_values(X, K + 1)
        T = level_up(B, O, X)
        CT = cover_poly_values(T, K)
        for k in range(1, K + 1):
            u, r, s = co[k]
            cnt += 1
            bad += (CT[k] != k * (u * CX[k + 1] + r * CX[k] + s * CX[k - 1]))
    print("closing operator C_T(k) = k[u_k C_X(k+1) + r_k C_X(k) + s_k C_X(k-1)]: cases", cnt, "failures", bad)
    # the 3-cycle operator and the coefficients of small pairs
    for name, B, O in [("point, O = {it}", [[0]], [0]), ("point, O empty", [[0]], []), ("TT_2, O = {source}", TT(2), [0]),
                       ("TT_2, O = {sink}", TT(2), [1]), ("TT_2, O = both", TT(2), [0, 1]),
                       ("C_3, O = one vertex", [[0, 1, 0], [0, 0, 1], [1, 0, 0]], [0]),
                       ("C_3, O = two vertices", [[0, 1, 0], [0, 0, 1], [1, 0, 0]], [0, 1])]:
        co = closing_coeffs(B, O, 5)
        print("  %-24s (u_k, r_k, s_k) for k = 1..5: %s   u_k s_(k+1): %s" % (name, [co[k] for k in range(1, 6)], [co[k][0] * co[k + 1][2] for k in range(1, 5)]))
    # nested 3-cycle tower: a(j) = (Phi^j k)(1); polynomials as coefficient lists
    def shift(p, c):   # p(k + c)
        from math import comb
        q = [0] * len(p)
        for d, a in enumerate(p):
            for i in range(d + 1):
                q[i] += a * comb(d, i) * c ** (d - i)
        return q
    def phi(p):
        a = shift(p, 1)
        b = shift(p, -1)
        inner = [0] * (len(p) + 1)
        for i in range(len(p)):
            inner[i] += a[i] - b[i]
            inner[i + 1] += p[i]
        return [0] + inner
    p = [0, 1]
    seq = [1]
    polys = []
    for j in range(1, 8):
        p = phi(p)
        seq.append(sum(p))
        polys.append(p)
    print("nested 3-cycle tower H(Phi^j(pt)), j = 0..7:", seq)
    print("cover polynomial of Phi^3(pt) (coefficients of k^0..k^7):", polys[2])
    # literal check for j <= 3: transitive order s_j..s_1, c, v_1..v_j with the arcs s_i -> v_i reversed
    from towers import H
    for j in range(1, 5):
        m = 2 * j + 1
        A = [[1 if a < b else 0 for b in range(m)] for a in range(m)]
        for i in range(j):
            a, b = i, m - 1 - i
            A[a][b], A[b][a] = 0, 1
        print("   literal H for j = %d (order %d): %d" % (j, m, H(A)))
