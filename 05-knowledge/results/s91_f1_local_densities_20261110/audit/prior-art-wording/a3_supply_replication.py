#!/usr/bin/env python3
"""Auditor prior-art-wording, check a3: independent recomputation of THM-5129 (B), (C), (D) and of the two claims
that no committed output of the session supports: "the poorest offsets (35, 39, 47, 63) are those with only two
strong realisers" and "correlation 0.68 between the covered share and the number of strong classes (23 <= h <= 73)".

Input data (read only, no code imported): the session's list of class representatives reps_le8.npz
(code: for i < j, bit j(j-1)/2 + i is 1 iff i -> j).  Everything else is this script's own code:
H, U_C(y), Pi_2 by subset dynamic programming; strong connectivity; the product law of THM-5129 (A).
"""
import sys
import numpy as np
from numba import njit

REPS = sys.argv[1]


@njit(cache=False)
def analyse(code, n, res_u):
    # returns H, Pi_2, strong flag; res_u[y] = U(y)
    out = np.zeros(n, dtype=np.int64)
    for j in range(n):
        base = j * (j - 1) // 2
        for i in range(j):
            if (code >> (base + i)) & 1:
                out[i] |= (1 << j)
            else:
                out[j] |= (1 << i)
    full = (1 << n) - 1
    # f[S, v] = number of Hamiltonian paths of T[S] ending at v ; g[S, v] = number starting at v
    f = np.zeros((full + 1, n), dtype=np.int64)
    g = np.zeros((full + 1, n), dtype=np.int64)
    hS = np.zeros(full + 1, dtype=np.int64)
    hS[0] = 1
    for S in range(1, full + 1):
        tot = 0
        for v in range(n):
            if not (S >> v) & 1:
                continue
            R = S ^ (1 << v)
            if R == 0:
                f[S, v] = 1
                g[S, v] = 1
            else:
                a = 0
                b = 0
                for u in range(n):
                    if (R >> u) & 1:
                        if (out[u] >> v) & 1:      # u -> v : path of R ending at u, then v
                            a += f[R, u]
                        if (out[v] >> u) & 1:      # v -> u : v, then path of R starting at u
                            b += g[R, u]
                f[S, v] = a
                g[S, v] = b
            tot += f[S, v]
        hS[S] = tot
    H = hS[full]
    pi2 = 0
    for S in range(full + 1):
        pi2 += hS[S] * hS[full ^ S]
    for yv in range(n):
        u = 0
        for S in range(1, full + 1):
            if (S >> yv) & 1:
                u += g[S, yv] * hS[full ^ S]
        res_u[yv] = u
    # strong: reachability from 0 forward and backward
    reach = 1
    changed = True
    while changed:
        changed = False
        for v in range(n):
            if (reach >> v) & 1:
                nr = reach | out[v]
                if nr != reach:
                    reach = nr
                    changed = True
    fwd = reach == full
    reach = 1
    changed = True
    while changed:
        changed = False
        for v in range(n):
            if not (reach >> v) & 1:
                if out[v] & reach:
                    reach |= (1 << v)
                    changed = True
    strong = fwd and (reach == full)
    return H, pi2, strong


def main():
    z = np.load(REPS, allow_pickle=False)
    LIM = 73
    root = {1: {1}}            # Root(h): slopes U_C(y) of strong C with H(C) = h
    pis = {}                   # h -> set of Pi_2(C), C strong of order >= 3 with H(C) = h
    nstrong = {}               # h -> number of strong classes with H = h (orders >= 3)
    least = {}
    counts = {}
    for n in range(1, 9):
        codes = z[f"n{n}"]
        ns = 0
        res_u = np.zeros(n, dtype=np.int64)
        for code in codes:
            H, pi2, strong = analyse(int(code), n, res_u)
            if n == 1:
                assert H == 1 and pi2 == 2 and res_u[0] == 1
            if strong and n >= 3:
                ns += 1
                least[n] = min(least.get(n, 10**9), H)
                if H <= LIM:
                    nstrong[H] = nstrong.get(H, 0) + 1
                    root.setdefault(H, set()).update(int(u) for u in res_u)
                    pis.setdefault(H, set()).add(int(pi2))
        counts[n] = (len(codes), ns)
    print("classes / strong classes by order:", counts)
    print("least H of a strong tournament, orders 3..8:", [least[n] for n in range(3, 9)])
    ok = [counts[n][0] for n in range(1, 9)] == [1, 1, 2, 4, 12, 56, 456, 6880]
    ok &= [counts[n][1] for n in range(3, 9)] == [1, 1, 6, 35, 353, 6008]
    ok &= [least[n] for n in range(3, 9)] == [3, 5, 9, 15, 25, 45]
    # attained numbers <= LIM: products of H of strong tournaments
    strongH = sorted(nstrong)
    att = {1}
    grew = True
    while grew:
        grew = False
        for a in list(att):
            for h in strongH:
                if a * h <= LIM and a * h not in att:
                    att.add(a * h)
                    grew = True
    missing = [m for m in range(1, LIM + 1, 2) if m not in att]
    print("odd numbers <= 73 that are not attained:", missing)
    ok &= missing == [7, 21]
    suf = {1: {1}}

    def Suf(d):
        if d in suf:
            return suf[d]
        res = set()
        for h in strongH:
            if h > 1 and d % h == 0:
                rest = Suf(d // h)
                for p in pis[h]:
                    for w in rest:
                        res.add(p * w)
        suf[d] = res
        return res

    def slopes0(D):
        res = set()
        for d1 in att:
            if D % d1:
                continue
            m = D // d1
            for h in [1] + strongH:
                if m % h:
                    continue
                d3 = m // h
                ws = Suf(d3)
                if not ws:
                    continue
                for u in root[h]:
                    for w in ws:
                        res.add(d1 * u * w)
        return res

    def oddparts(D):
        best = {}
        for U in slopes0(D):
            j = 0
            d = U
            while d % 2 == 0:
                d //= 2
                j += 1
            best[d] = min(best.get(d, 99), j)
        return best

    table = {D: oddparts(D) for D in sorted(att)}
    exp_counts = {3: 1, 5: 5, 9: 8, 17: 11, 33: 34, 65: 66}
    for D, c in exp_counts.items():
        print(f"offset {D}: {len(table[D])} odd parts (THM-5129: {c});",
              " ".join(f"{d}:{j}" for d, j in sorted(table[D].items())[:14]), "...")
        ok &= len(table[D]) == c
    ok &= table[5] == {1: 5, 3: 1, 5: 0, 7: 0, 9: 0}
    ok &= table[9] == {1: 4, 5: 1, 9: 0, 11: 0, 15: 0, 21: 0, 23: 0, 25: 0}
    ok &= table[17] == {5: 2, 7: 2, 9: 1, 11: 2, 17: 0, 21: 0, 31: 0, 33: 0, 55: 0, 61: 0, 69: 0}
    print("least exponent of d = 9 at offsets 5, 9, 17, 33, 65:", [table[D].get(9) for D in (5, 9, 17, 33, 65)],
          "  bound from U >= D alone (least j with 9 * 2^j >= D):",
          [next(j for j in range(20) if 9 * 2**j >= D) for D in (5, 9, 17, 33, 65)])
    ok &= [table[D].get(9) for D in (5, 9, 17, 33, 65)] == [0, 0, 1, 2, 5]

    def is_slope(D, U):
        if D not in table:
            return False
        j = 0
        d = U
        while d % 2 == 0:
            d //= 2
            j += 1
        return d in table[D] and j >= table[D][d]

    # (C)
    zero = []
    for N in range(3, 225, 2):
        r = sum(1 for D in range(1, N // 3 + 1, 2) if is_slope(D, (N - D) // 2))
        if r == 0:
            zero.append(N)
    print("odd N in [3, 225) with no representation D + 2U:", zero)
    ok &= zero == [7, 11, 13, 21, 47, 67]
    # (D)
    share = {}
    line = []
    for h in sorted(att):
        if h == 7 or h == 21:
            continue
        ints = list(range(h, (3 * h) // 2 + 1))
        c = sum(1 for U in ints if is_slope(h, U))
        share[h] = c / len(ints)
        line.append(f"{h}:{c}/{len(ints)}")
    print("slopes in [h, 3h/2]:", " ".join(line))
    hs = [h for h in sorted(share) if 23 <= h <= 73]
    xs = np.array([share[h] for h in hs])
    ys = np.array([nstrong.get(h, 0) for h in hs], dtype=float)
    corr = float(np.corrcoef(xs, ys)[0, 1])
    print("number of strong classes with H = h, 23 <= h <= 73:", " ".join(f"{h}:{nstrong.get(h, 0)}" for h in hs))
    print(f"Pearson correlation of the covered share with the number of strong classes, 23 <= h <= 73: {corr:.3f}")
    ranks_x = xs.argsort().argsort()
    ranks_y = ys.argsort().argsort()
    print(f"(correlation with log of the number of strong classes: "
          f"{float(np.corrcoef(xs, np.log(np.maximum(ys, 1)))[0, 1]):.3f})")
    two = [h for h in sorted(share) if 15 <= h <= 73 and nstrong.get(h, 0) == 2]
    print("offsets 15 <= h <= 73 with exactly two strong classes:", two)
    poor = sorted(share, key=lambda h: share[h])[:8]
    print("the eight poorest offsets (share, strong classes):",
          [(h, round(share[h], 3), nstrong.get(h, 0)) for h in poor])
    m = [share[h] for h in sorted(share) if 37 <= h <= 73]
    print(f"mean share on 37 <= h <= 73: {np.mean(m):.3f} (min {min(m):.3f}, max {max(m):.3f})")
    print("ALL REPLICATION CHECKS PASSED" if ok else "REPLICATION CHECK FAILED")


if __name__ == "__main__":
    main()
