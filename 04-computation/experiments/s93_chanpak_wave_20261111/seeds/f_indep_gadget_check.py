#!/usr/bin/env python3
"""S93 lane seeds, FINISHER: independent recount (pure python, no numba, no code of the lane) of the 14 unit neutral
gadgets of order 8 over non-transitive X (c_neutral_unit_8.txt), i.e. of the refutation of HYP-9662 (D), first sentence.

Definitions (THM-5130): Pi_2 = ordered pairs (A, B) of disjoint possibly empty paths covering X;
Q(e, f) = ordered pairs (D1, D2) of nonempty paths covering X, D1 ending at e, D2 starting at f;
Q3(e, f) = the same with a third NONEMPTY path D3;  W = sum over ordered pairs e != f of P of Q,  W3 likewise;
kappa = Pi_2 - 2 W - 2 W3,  r_X = H(X) - W.
Then the sign-preserving gadget (colcode 19771611, P mask 50) is put on the seed (C_3, {0, 1}) and the closed
tournaments T_1, T_2 of the new pair are built from the definitions and counted: expected 1 + 2^t * 132458."""


def decode(n, code):
    out = [0] * n
    b = 0
    for j in range(1, n):
        for i in range(j):
            if (code >> b) & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            b += 1
    return out


def tables(n, out):
    full = (1 << n) - 1
    end = [[0] * n for _ in range(full + 1)]
    start = [[0] * n for _ in range(full + 1)]
    h = [0] * (full + 1)
    h[0] = 1
    for S in range(1, full + 1):
        if S & (S - 1) == 0:
            v = S.bit_length() - 1
            end[S][v] = start[S][v] = 1
            h[S] = 1
            continue
        t = 0
        for v in range(n):
            if (S >> v) & 1:
                R = S ^ (1 << v)
                a = b = 0
                for u in range(n):
                    if (R >> u) & 1:
                        if (out[u] >> v) & 1:
                            a += end[R][u]
                        else:
                            b += start[R][u]
                end[S][v] = a
                start[S][v] = b
                t += a
        h[S] = t
    return end, start, h


def hcount(n, out):
    full = (1 << n) - 1
    dp = [dict() for _ in range(full + 1)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, full + 1):
        for v, c in dp[S].items():
            rest = out[v] & ~S & full
            while rest:
                low = rest & -rest
                u = low.bit_length() - 1
                T = S | low
                dp[T][u] = dp[T].get(u, 0) + c
                rest ^= low
    return sum(dp[full].values())


ok = True
rows = [tuple(int(x) for x in line.split()) for line in open("c_neutral_unit_8.txt") if not line.startswith("#")]
for (n, code, mask, szP, H0, Pi0, W0, r0, k0, u0, s0, u30, e0, nu0) in rows:
    out = decode(n, code)
    end, start, h = tables(n, out)
    full = (1 << n) - 1
    H = h[full]
    Pi2 = sum(h[S] * h[full ^ S] for S in range(full + 1))
    P = [v for v in range(n) if (mask >> v) & 1]
    W = W3 = 0
    for e in P:
        for f in P:
            if e == f:
                continue
            for S3 in range(full + 1):
                if (S3 >> e) & 1 or (S3 >> f) & 1:
                    continue
                rest = full ^ S3
                tot = 0
                S = rest
                while S:
                    if (S >> e) & 1 and not (S >> f) & 1:
                        tot += end[S][e] * start[rest ^ S][f]
                    S = (S - 1) & rest
                if S3 == 0:
                    W += tot
                else:
                    W3 += h[S3] * tot
    kappa = Pi2 - 2 * W - 2 * W3
    good = (H == H0 and Pi2 == Pi0 and W == W0 and kappa == 0 and H - W == r0 and len(P) == 3 and H > 1)
    ok &= good
    print(f"X = {code}, P mask {mask}: H = {H}, Pi_2 = {Pi2}, W = {W}, W3 = {W3}, kappa = {kappa}, r_X = {H - W}, |P| = {len(P)}: {'agrees' if good else 'DISAGREES'}")
print("all 14 gadgets: kappa = 0, |P| = 3, r_X = +-1 over a non-transitive X:", ok)

# the sign-preserving gadget on the seed (C_3, {0, 1}): B_new = C_3 => X => v, v -> {0, 1}, vertex 2 -> v, X -> v
X = decode(8, 19771611)
Pm = 50
nb = 3 + 8 + 1
base = [0] * nb
base[0] |= 1 << 1
base[1] |= 1 << 2
base[2] |= 1 << 0
for b in range(3):
    for x in range(8):
        base[b] |= 1 << (3 + x)
for x in range(8):
    base[3 + x] |= X[x] << 3
    base[3 + x] |= 1 << 11
base[11] |= (1 << 0) | (1 << 1)
base[2] |= 1 << 11
Pnew = Pm << 3
res = []
for t in (1, 2):
    N = nb + t + 1
    T = base[:] + [0] * (t + 1)
    w = nb + t
    for b in range(nb):
        for k in range(t):
            T[b] |= 1 << (nb + k)
        if (Pnew >> b) & 1:
            T[w] |= 1 << b
        else:
            T[b] |= 1 << w
    for k in range(t):
        for l in range(k + 1, t):
            T[nb + k] |= 1 << (nb + l)
        T[nb + k] |= 1 << w
    assert all(((T[a] >> b) & 1) + ((T[b] >> a) & 1) == 1 for a in range(N) for b in range(a + 1, N))
    res.append(hcount(N, T))
print(f"sign-preserving gadget on (C_3, two vertices): H(T_1), H(T_2) = {res}; expected {[1 + 2 * 132458, 1 + 4 * 132458]}: {res == [264917, 529833]}")
