#!/usr/bin/env python3
"""THM-5125 (cycle-gas specialisations) and the tested-and-refuted bridges of HYP-9657.  Exact checks, all labelled
tournaments of orders 3..6.

For a tournament T with adjacency matrix A (a_ij = 1 iff i -> j), a *cycle packing* is a set of pairwise
vertex-disjoint directed cycles (lengths >= 3); an *odd packing* uses odd cycles only.

(G1) per(I + A) = number of cycle packings;  det(I + A) = sum over packings of (-1)^(number of even cycles);
     H(T) = sum over odd packings of 2^(number of cycles)  (the odd-cycle formula, checked against path counting);
     H(T) = sum over odd packings of (-1)^(number of cycles)  (mod 3).
(G2) Cubic criterion, pure adjacency: nullity_F2(A^2 + A + I) > 0  iff  P_0 = P_1 = P_2 (mod 2), where P_j is the
     number of cycle packings (the empty one included) covering a number of vertices congruent to j mod 3.
(G3) Cubic criterion, Monsky convention M = A + diag(row parities): nullity_F2(M^2 + M + I) > 0 iff
     sum over cycle packings of omega^(e - o) = 0 in F_4, e / o = numbers of uncovered vertices of even / odd
     out-degree; equivalently Q_0 = Q_1 = Q_2 (mod 2) with Q_j = #packings with e - o = j mod 3.
(B)  Bridges tested and REFUTED: H mod 4 or mod 8 is not a function of the F_2 data
     (corank of the F_2-Laplacian, nullity(M^2+M+I), nullity(I+A)); "nullity(M^2+M+I) > 0 implies 3 | H" fails
     first at order 5 (H = 13), and with the pure adjacency matrix first at order 6 (H = 17).
Run:  python3 cycle_gas.py > cycle_gas.out
"""
import itertools
import sys
from collections import Counter, defaultdict


def rank_f2(rows, n):
    rows = list(rows)
    r = 0
    for col in range(n):
        piv = None
        for i in range(r, len(rows)):
            if rows[i] >> col & 1:
                piv = i
                break
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i] >> col & 1:
                rows[i] ^= rows[r]
        r += 1
    return r


def mul_f2(A, B, n):
    C = []
    for i in range(n):
        row = 0
        a = A[i]
        k = 0
        while a:
            if a & 1:
                row ^= B[k]
            a >>= 1
            k += 1
        C.append(row)
    return C


def phi3_nullity(A, n):
    A2 = mul_f2(A, A, n)
    return n - rank_f2([A2[i] ^ A[i] ^ (1 << i) for i in range(n)], n)


def ham_paths(n, beats):
    size = 1 << n
    cnt = [[0] * n for _ in range(size)]
    for v in range(n):
        cnt[1 << v][v] = 1
    for mask in range(1, size):
        row = cnt[mask]
        for v in range(n):
            c = row[v]
            if c:
                nxt = beats[v] & ~mask
                while nxt:
                    w = nxt & -nxt
                    cnt[mask | w][w.bit_length() - 1] += c
                    nxt ^= w
    return sum(cnt[size - 1])


def directed_cycles(n, beats):
    """all directed cycles as (vertex mask, length); each cycle once (rooted at its least vertex)."""
    out = []
    for start in range(n):
        stack = [(start, 1 << start, 1)]
        while stack:
            v, mask, length = stack.pop()
            nxt = beats[v]
            if length >= 3 and (nxt >> start) & 1:
                out.append((mask, length))
            cand = nxt & ~mask
            while cand:
                w = cand & -cand
                cand ^= w
                wi = w.bit_length() - 1
                if wi > start:
                    stack.append((wi, mask | w, length + 1))
    return out


def packings(cycles):
    """yield (covered mask, number of cycles, number of even cycles) over all sets of pairwise disjoint cycles."""
    cycles = sorted(cycles)
    out = [(0, 0, 0)]

    def rec(idx, mask, k, ev):
        for j in range(idx, len(cycles)):
            m, length = cycles[j]
            if m & mask:
                continue
            item = (mask | m, k + 1, ev + (1 if length % 2 == 0 else 0))
            out.append(item)
            rec(j + 1, item[0], item[1], item[2])

    rec(0, 0, 0, 0)
    return out


def det_int(M):
    n = len(M)
    total = 0
    for perm in itertools.permutations(range(n)):
        prod = 1
        for i in range(n):
            prod *= M[i][perm[i]]
            if prod == 0:
                break
        if prod:
            inv = sum(1 for i in range(n) for j in range(i + 1, n) if perm[i] > perm[j])
            total += -prod if inv & 1 else prod
    return total


def per_int(M):
    n = len(M)
    total = 0
    for perm in itertools.permutations(range(n)):
        prod = 1
        for i in range(n):
            prod *= M[i][perm[i]]
            if prod == 0:
                break
        total += prod
    return total


def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        beats = [0] * n
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1:
                beats[i] |= 1 << j
            else:
                beats[j] |= 1 << i
        yield beats


def main():
    ok = True
    for n in range(3, 7):
        bad = Counter()
        tab8 = defaultdict(Counter)
        ex_mon = None
        ex_adj = None
        count = 0
        for beats in all_tournaments(n):
            count += 1
            cyc = directed_cycles(n, beats)
            pk = packings(cyc)
            H = ham_paths(n, beats)
            # (G1)
            ocf = sum(2 ** k for (m, k, ev) in pk if ev == 0)          # odd packings only
            if ocf != H:
                bad["ocf"] += 1
            eul = sum((-1) ** k for (m, k, ev) in pk if ev == 0)
            if (H - eul) % 3 != 0:
                bad["mod3"] += 1
            if n <= 5:
                IA = [[(1 if i == j else 0) + (beats[i] >> j & 1) for j in range(n)] for i in range(n)]
                if per_int(IA) != len(pk):
                    bad["per"] += 1
                if det_int(IA) != sum((-1) ** ev for (m, k, ev) in pk):
                    bad["det"] += 1
            # (G2) pure adjacency
            P = [0, 0, 0]
            for (m, k, ev) in pk:
                P[bin(m).count("1") % 3] ^= 1
            nul_adj = phi3_nullity(beats, n)
            if (nul_adj > 0) != (P[0] == P[1] == P[2]):
                bad["G2"] += 1
            # (G3) Monsky convention
            par = [bin(beats[i]).count("1") & 1 for i in range(n)]
            M = [beats[i] ^ (par[i] << i) for i in range(n)]
            Q = [0, 0, 0]
            for (m, k, ev) in pk:
                e = sum(1 for v in range(n) if not (m >> v) & 1 and par[v] == 0)
                o = sum(1 for v in range(n) if not (m >> v) & 1 and par[v] == 1)
                Q[(e - o) % 3] ^= 1
            nul_mon = phi3_nullity(M, n)
            if (nul_mon > 0) != (Q[0] == Q[1] == Q[2]):
                bad["G3"] += 1
            # (B) bridges
            L = M                                                   # F_2 Laplacian D + A
            corank = n - 1 - rank_f2(L, n)
            nulIA = n - rank_f2([beats[i] ^ (1 << i) for i in range(n)], n)
            tab8[(corank, nul_mon, nulIA)][H % 8] += 1
            if nul_mon > 0 and H % 3 != 0 and ex_mon is None:
                ex_mon = (list(beats), H)
            if nul_adj > 0 and H % 3 != 0 and ex_adj is None:
                ex_adj = (list(beats), H)
        good = (sum(bad.values()) == 0)
        ok &= good
        print(f"n={n}: {count} labelled tournaments; identity failures {dict(bad)}", "OK" if good else "FAIL")
        print(f"   (B) first tournament with nullity(M^2+M+I) > 0 and 3 not dividing H: {ex_mon}; pure adjacency: {ex_adj}")
        nonfunc = sum(1 for key in tab8 if len({h % 4 for h in tab8[key]}) > 1)
        print(f"   (B) classes of (Laplacian corank - 1, nullity(M^2+M+I), nullity(I+A)): {len(tab8)}; "
              f"classes on which H mod 4 is not constant: {nonfunc}")
        for key in sorted(tab8):
            print("        ", key, dict(sorted(tab8[key].items())))
    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
