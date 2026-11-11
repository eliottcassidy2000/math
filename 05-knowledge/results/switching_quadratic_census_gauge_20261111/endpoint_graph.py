#!/usr/bin/env python3
"""Random controls of THM-5139 (the endpoint graph of a switching class), orders 3-10.
  (1) e_{T^X} = e_T + Lambda 1_X, Lambda symmetric, zero diagonal, even rows
  (2) hc(T^X) = hc(T) + sum_{a in X}(1 + e_T(a)) + sum_{pairs in X} Lambda(a,b)            (mod 2)
  (3) Lambda(u,w) = 1 + <a+b, e_Q> + a^T Lambda_Q b,  Q = T-u-w                               (n >= 4)
  (4) (H(T^X) - H(T))/2 = sum_{a in X} lambda_a + sum_{pairs in X} M(a,b)  with lambda, M built from the
      subtournaments on odd vertex sets                                                       (n <= 8)
  (5) M(a,b) = sum_{v in N+(a) xor N+(b)} E_R(v) + iota_a^T M_R iota_b,  R = T-a-b            (n >= 5)
  (6) order 6: M(a,b) = 1 + d_a + d_b
usage: endpoint_graph.py [seed]"""
import random, sys, itertools

def rand_t(n, rng):
    out = [0]*n
    for i in range(n):
        for j in range(i+1, n):
            if rng.getrandbits(1): out[i] |= 1 << j
            else: out[j] |= 1 << i
    return out

def switch(out, X):
    n = len(out); full = (1 << n) - 1; o2 = [0]*n
    for i in range(n):
        if (X >> i) & 1: o2[i] = ((out[i] & X) | (~out[i] & ~X & full)) & ~(1 << i)
        else:            o2[i] = ((out[i] & ~X) | (~out[i] & X & full)) & ~(1 << i)
    return o2

def conv(out):
    n = len(out); full = (1 << n) - 1
    return [(~out[i]) & full & ~(1 << i) for i in range(n)]

def endtab(out):
    """dp[S][v] = number of paths with vertex set S ending at v"""
    n = len(out); dp = [[0]*n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        row = dp[S]
        for v in range(n):
            c = row[v]
            if not c: continue
            t = out[v] & ~S
            while t:
                b = t & -t; dp[S | b][b.bit_length() - 1] += c; t ^= b
    return dp

def H(out): return sum(endtab(out)[(1 << len(out)) - 1])

def evec(out):
    n = len(out); full = (1 << n) - 1
    t = endtab(out)[full]; s = endtab(conv(out))[full]
    return [(s[i] + t[i]) % 2 for i in range(n)]

def hc(out):
    n = len(out)
    if n < 3: return 0
    dp = [[0]*n for _ in range(1 << n)]; dp[1][0] = 1
    for S in range(1, 1 << n, 2):
        row = dp[S]
        for v in range(n):
            c = row[v]
            if not c: continue
            t = out[v] & ~S
            while t:
                b = t & -t; dp[S | b][b.bit_length() - 1] += c; t ^= b
    return sum(dp[(1 << n) - 1][v] for v in range(1, n) if out[v] & 1)

def sub(out, U):
    idx = {v: i for i, v in enumerate(U)}; o = [0]*len(U)
    for v in U:
        for w in U:
            if v != w and (out[v] >> w) & 1: o[idx[v]] |= 1 << idx[w]
    return o

def Lam_of(out):
    n = len(out); e0 = evec(out); L = [[0]*n for _ in range(n)]
    for w in range(n):
        ew = evec(switch(out, 1 << w))
        for u in range(n): L[u][w] = (e0[u] + ew[u]) % 2
    return e0, L

def Mmat(o):
    n = len(o); H0 = H(o); Ha = [H(switch(o, 1 << a)) for a in range(n)]
    M = [[0]*n for _ in range(n)]
    for a, b in itertools.combinations(range(n), 2):
        M[a][b] = M[b][a] = ((H0 + Ha[a] + Ha[b] + H(switch(o, (1 << a) | (1 << b)))) // 2) % 2
    return M

def Eodd(o):
    n = len(o); e = endtab(o); s = endtab(conv(o)); res = [0]*n
    for S in range(1, 1 << n):
        k = bin(S).count('1')
        if k % 2 == 1 and k >= 3:
            for v in range(n): res[v] ^= (e[S][v] + s[S][v]) & 1
    return res

if __name__ == "__main__":
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 20261111)
    for n in range(3, 11):
        trials = 30 if n <= 7 else (8 if n <= 9 else 3)
        c = [0]*6; t = [0]*6
        for _ in range(trials):
            o = rand_t(n, rng); e0, L = Lam_of(o); h0 = hc(o); H0 = H(o)
            ok = all(L[u][w] == L[w][u] for u in range(n) for w in range(n)) and all(L[u][u] == 0 for u in range(n)) and all(sum(r) % 2 == 0 for r in L)
            ok2 = True
            for _2 in range(6):
                X = rng.getrandbits(n); xs = [i for i in range(n) if (X >> i) & 1]; oX = switch(o, X)
                ok &= (evec(oX) == [(e0[u] + sum(L[u][w] for w in xs)) % 2 for u in range(n)])
                ok2 &= (hc(oX) % 2 == (h0 + sum(1 + e0[a] for a in xs) + sum(L[a][b] for a, b in itertools.combinations(xs, 2))) % 2)
            t[0] += 1; c[0] += ok; t[1] += 1; c[1] += ok2
            if n >= 4:
                ok3 = True
                for u, w in itertools.combinations(range(n), 2):
                    Q = [v for v in range(n) if v not in (u, w)]; eq, Lq = Lam_of(sub(o, Q))
                    a = [(o[u] >> v) & 1 for v in Q]; b = [(o[w] >> v) & 1 for v in Q]
                    val = (1 + sum((a[i] + b[i]) * eq[i] for i in range(len(Q))) + sum(a[i] * Lq[i][j] * b[j] for i in range(len(Q)) for j in range(len(Q)))) % 2
                    ok3 &= (val == L[u][w])
                t[2] += 1; c[2] += ok3
            if n <= 8:
                M = [[0]*n for _ in range(n)]; lam = [0]*n
                for k in range(3, n + 1, 2):
                    for U in itertools.combinations(range(n), k):
                        eU, LU = Lam_of(sub(o, list(U)))
                        for i, a in enumerate(U):
                            lam[a] ^= (1 + eU[i]) % 2
                            for j, b in enumerate(U): M[a][b] ^= LU[i][j]
                ok4 = (M == Mmat(o))
                for _2 in range(6):
                    X = rng.getrandbits(n); xs = [i for i in range(n) if (X >> i) & 1]
                    ok4 &= (((H(switch(o, X)) - H0) // 2) % 2 == (sum(lam[a] for a in xs) + sum(M[a][b] for a, b in itertools.combinations(xs, 2))) % 2)
                t[3] += 1; c[3] += ok4
            if 5 <= n <= 9:
                M = Mmat(o); ok5 = True; sc = [bin(x).count('1') for x in o]
                for a, b in itertools.combinations(range(n), 2):
                    R = [v for v in range(n) if v not in (a, b)]; oR = sub(o, R); MR = Mmat(oR); E = Eodd(oR)
                    ia = [(o[a] >> v) & 1 for v in R]; ib = [(o[b] >> v) & 1 for v in R]
                    val = (sum(E[i] for i in range(len(R)) if ia[i] != ib[i]) + sum(ia[i] * MR[i][j] * ib[j] for i in range(len(R)) for j in range(len(R)))) % 2
                    ok5 &= (val == M[a][b])
                t[4] += 1; c[4] += ok5
                if n == 6:
                    t[5] += 1; c[5] += all(M[a][b] == (1 + sc[a] + sc[b]) % 2 for a, b in itertools.combinations(range(n), 2))
        print("n=%d  (1) affine+even graph %d/%d  (2) hc quadratic %d/%d  (3) recursion (**) %d/%d  (4) H mod 4 formula %d/%d  (5) Hessian recursion %d/%d%s"
              % (n, c[0], t[0], c[1], t[1], c[2], t[2], c[3], t[3], c[4], t[4], ("  (6) order-6 closed form %d/%d" % (c[5], t[5])) if n == 6 else ""))
