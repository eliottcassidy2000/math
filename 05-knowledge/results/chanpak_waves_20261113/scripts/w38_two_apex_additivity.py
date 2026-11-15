#!/usr/bin/env python3
"""Wave 38: are fibre swaps additive inside two-apex sources?

Two-apex source (chanpak_pair_line_census_20261115/TWO_APEX_COARSE_PLACEMENT.md): a transitive spine 0 < 1 < ... < m-1,
apex A and apex B with A -> B; A beats spine vertex i iff a_i = 1, B beats it iff b_i = 1; a_0 = 1, a_{m-1} = 0, b != 0.
The triangle reader is F = H + Q, Q = Pi_2 = P2 + 2H.

Inflation: replace the spine vertices at positions p1 < p2 by modules X1, X2 (each vertex of a module inherits the
spine vertex's arcs).  Every count is then multilinear in the port vectors v(X) = (P_1(X), P_2(X), ...), where P_k(X)
counts ordered k-tuples of disjoint nonempty paths covering X.

A fibre swap replaces X by X' with H(X') = H(X).  Swaps in the two modules are additive iff the mixed difference
    D = F(X1', X2') - F(X1', X2) - F(X1, X2') + F(X1, X2)
vanishes.  By multilinearity D = sum_{k, l >= 2} c_kl dP_k(X1) dP_l(X2), where c_kl counts the quotient covers that
split module 1 into k pieces and module 2 into l pieces.  The script measures D on random two-apex sources and
extracts c_kl exactly with transitive probe modules (P_k(TT_s) = k! S(s, k)).
Usage: python3 w38_two_apex_additivity.py HQ_BINARY
"""
import itertools
import random
import subprocess
import sys
from fractions import Fraction
from math import factorial

HQ = sys.argv[1]
rng = random.Random(20261115)


def hq_batch(tours):
    """tours: list of out-mask lists. Returns list of (H, Q)."""
    inp = "".join(f"{len(t)} " + " ".join(map(str, t)) + "\n" for t in tours)
    out = subprocess.run([HQ], input=inp.encode(), capture_output=True, check=True).stdout.decode().split("\n")
    return [tuple(map(int, l.split())) for l in out if l.strip()]


def gentourng(n):
    """All n-tournaments (isomorphism classes) as out-mask lists."""
    txt = subprocess.run(["gentourng", "-q", str(n)], capture_output=True, check=True).stdout.decode().split()
    res = []
    for line in txt:
        out = [0] * n
        k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if line[k] == "1":
                    out[i] |= 1 << j
                else:
                    out[j] |= 1 << i
                k += 1
        res.append(out)
    return res


def ham_table(out):
    n = len(out)
    N = 1 << n
    dp = [[0] * n for _ in range(N)]
    for v in range(n):
        dp[1 << v][v] = 1
    H = [0] * N
    H[0] = 1
    for S in range(1, N):
        for v in range(n):
            c = dp[S][v]
            if c:
                H[S] += c
                nb = out[v] & ~S
                while nb:
                    u = (nb & -nb).bit_length() - 1
                    nb &= nb - 1
                    dp[S | 1 << u][u] += c
    return H


def ports(out):
    """P_k for k = 1..n: ordered k-tuples of disjoint nonempty paths covering the tournament."""
    n = len(out)
    N = 1 << n
    H = ham_table(out)
    cur = {0: 1}
    # f_k[S] = number of ordered k-tuples covering S
    fk = [0] * N
    fk[0] = 1
    res = []
    for k in range(1, n + 1):
        g = [0] * N
        for S in range(1, N):
            # choose the last block B as a nonempty subset of S
            B = S
            tot = 0
            while B:
                tot += fk[S ^ B] * H[B]
                B = (B - 1) & S
            g[S] = tot
        fk = g
        res.append(fk[N - 1])
    return res


def two_apex(m, a, b, mods):
    """mods: dict position -> module out-mask list. Returns out-mask list of the inflated source."""
    verts = []  # (kind, pos, local index)
    for i in range(m):
        if i in mods:
            for t in range(len(mods[i])):
                verts.append(("s", i, t))
        else:
            verts.append(("s", i, 0))
    verts.append(("A", -1, 0))
    verts.append(("B", -1, 0))
    n = len(verts)
    out = [0] * n
    for x in range(n):
        for y in range(x + 1, n):
            kx, px, tx = verts[x]
            ky, py, ty = verts[y]
            if kx == "s" and ky == "s":
                if px == py:
                    xwins = (mods[px][tx] >> ty) & 1
                else:
                    xwins = px < py
            elif kx == "s" and ky == "A":
                xwins = not ((a >> px) & 1)
            elif kx == "s" and ky == "B":
                xwins = not ((b >> px) & 1)
            else:  # A vs B
                xwins = True
            if xwins:
                out[x] |= 1 << y
            else:
                out[y] |= 1 << x
    return out


def random_words(m):
    while True:
        a = rng.randrange(1 << m) | 1
        a &= ~(1 << (m - 1))
        b = rng.randrange(1, 1 << m)
        return a, b


def stirling2(n, k):
    return sum((-1) ** j * factorial(k) // (factorial(j) * factorial(k - j)) * (k - j) ** n for j in range(k + 1)) // factorial(k)


def TT(s):
    return [sum(1 << j for j in range(i + 1, s)) for i in range(s)]


def main():
    # fibre pairs of modules at orders 5 and 6: equal H, different P2
    pairs = []
    for k in (5, 6):
        allk = gentourng(k)
        byH = {}
        for out in allk:
            pv = ports(out)
            byH.setdefault(pv[0], []).append((pv, out))
        for h, lst in sorted(byH.items()):
            lst = sorted(lst)
            for (p1, o1), (p2, o2) in itertools.combinations(lst, 2):
                if p1[1] != p2[1]:
                    pairs.append((k, h, p1, o1, p2, o2))
    print(f"module fibre pairs (orders 5, 6, equal H, different P2): {len(pairs)}")
    for (k, h, p1, o1, p2, o2) in pairs[:6]:
        print(f"  order {k}, H = {h}: ports {p1} vs {p2}")

    # (1) mixed differences on random two-apex sources
    stats = {}
    tours = []
    meta = []
    for trial in range(240):
        m = rng.choice([6, 7, 8])
        a, b = random_words(m)
        p1, p2 = sorted(rng.sample(range(m), 2))
        (k1, h1, v1, X1, v1p, X1p) = rng.choice(pairs)
        (k2, h2, v2, X2, v2p, X2p) = rng.choice(pairs)
        if m + 2 + (k1 - 1) + (k2 - 1) > 20:
            continue
        for M1, M2 in ((X1, X2), (X1p, X2), (X1, X2p), (X1p, X2p)):
            tours.append(two_apex(m, a, b, {p1: M1, p2: M2}))
        region = ("A" if (a >> p1) & 1 else "a") + ("B" if (b >> p1) & 1 else "b") + "|" + \
                 ("A" if (a >> p2) & 1 else "a") + ("B" if (b >> p2) & 1 else "b")
        meta.append((m, a, b, p1, p2, region))
    res = hq_batch(tours)
    zero = 0
    ratios = []
    rows = []
    for t, (m, a, b, p1, p2, region) in enumerate(meta):
        F = [h + q for (h, q) in res[4 * t: 4 * t + 4]]
        D = F[3] - F[1] - F[2] + F[0]
        d1, d2 = F[1] - F[0], F[2] - F[0]
        if D == 0:
            zero += 1
        if d1 and d2:
            ratios.append(abs(D) / min(abs(d1), abs(d2)))
        rows.append((m, p1, p2, region, d1, d2, D))
    ratios.sort()
    print(f"\n(1) random inflated two-apex sources: {len(meta)} tests; mixed difference D = 0 in {zero}")
    if ratios:
        print(f"    |D| / min(|dF1|, |dF2|): median {ratios[len(ratios) // 2]:.3f}, "
              f"10% {ratios[len(ratios) // 10]:.3f}, 90% {ratios[9 * len(ratios) // 10]:.3f}")
    for r in rows[:8]:
        print("    m=%d p=(%d,%d) apex-bits %s: dF1=%d dF2=%d D=%d" % r)

    # (2) exact multilinear coefficients c_kl with transitive probes, and the additivity criterion
    print("\n(2) exact coefficients c_kl (pieces of module 1, module 2) of F = H + Q, from transitive probes")
    K = 4
    for trial in range(6):
        m = rng.choice([6, 7])
        a, b = random_words(m)
        p1, p2 = sorted(rng.sample(range(m), 2))
        tours = [two_apex(m, a, b, {p1: TT(s), p2: TT(t)}) for s in range(1, K + 1) for t in range(1, K + 1)]
        res = hq_batch(tours)
        Fv = [[res[(s - 1) * K + (t - 1)][0] + res[(s - 1) * K + (t - 1)][1] for t in range(1, K + 1)] for s in range(1, K + 1)]
        # P_k(TT_s) = k! S(s, k), lower triangular in (s, k): solve Fv = P C P^T
        P = [[Fraction(factorial(k) * stirling2(s, k)) for k in range(1, K + 1)] for s in range(1, K + 1)]
        # invert lower-triangular P
        Pinv = [[Fraction(0)] * K for _ in range(K)]
        for i in range(K):
            Pinv[i][i] = 1 / P[i][i]
            for j in range(i):
                Pinv[i][j] = -sum(P[i][l] * Pinv[l][j] for l in range(j, i)) / P[i][i]
        C = [[sum(Pinv[i][s] * Fv[s][t] * Pinv[j][t] for s in range(K) for t in range(K)) for j in range(K)] for i in range(K)]
        region = ("A" if (a >> p1) & 1 else "a") + ("B" if (b >> p1) & 1 else "b") + "|" + \
                 ("A" if (a >> p2) & 1 else "a") + ("B" if (b >> p2) & 1 else "b")
        print(f"  m={m} a={a:0{m}b} b={b:0{m}b} p=({p1},{p2}) bits {region}:")
        for i in range(K):
            print("     " + " ".join(f"{str(C[i][j]):>8}" for j in range(K)))
        # verify multilinearity on a real fibre pair
        (k1, h1, v1, X1, v1p, X1p) = rng.choice([p for p in pairs if p[0] == 5])
        (k2, h2, v2, X2, v2p, X2p) = rng.choice([p for p in pairs if p[0] == 5])
        res2 = hq_batch([two_apex(m, a, b, {p1: M1, p2: M2}) for M1, M2 in ((X1, X2), (X1p, X2), (X1, X2p), (X1p, X2p))])
        F = [h + q for (h, q) in res2]
        Dobs = F[3] - F[1] - F[2] + F[0]
        dv1 = [x - y for x, y in zip(v1p, v1)] + [0] * K
        dv2 = [x - y for x, y in zip(v2p, v2)] + [0] * K
        Dpred = sum(C[i][j] * dv1[i] * dv2[j] for i in range(K) for j in range(K))
        print(f"     multilinear check on order-5 fibre pairs: D observed {Dobs}, predicted {Dpred}")


main()
