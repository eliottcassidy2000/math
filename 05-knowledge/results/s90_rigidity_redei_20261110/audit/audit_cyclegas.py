#!/usr/bin/env python3
"""Independent audit of THM-5125 statement 4 (cycle-gas identities) and HYP-9657 (C) (refuted bridges).
Own code.  Differences from the packet's cycle_gas.py:
  * cycle packings are enumerated as permutations sigma with sigma(i) = i or i -> sigma(i) (the packet enumerates
    cycles and then disjoint families);
  * per(I + A) is computed by Ryser's formula and det(I + A) by fraction-free elimination (the packet expands both
    over permutations);
  * the polynomial identities over F_2 are checked coefficientwise (det(xI + A), det(xI + M) by the
    Faddeev-LeVerrier recursion over Z, reduced mod 2), not only at x = omega;
  * the F_4 statements are checked by Gaussian elimination over F_4 (dim ker(A + omega I)), and
    nullity_F2(A^2 + A + I) = 2 dim_F4 ker(A + omega I);
  * the orientation independence of the Monsky nullity is checked on every tournament (converse orientation);
  * Hamiltonian paths by depth-first enumeration.
Universe: every labelled tournament of orders 3..6.
"""
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction

MUL4 = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]]     # F_4 = {0, 1, w = 2, w^2 = 3}, addition = xor
INV4 = [0, 1, 3, 2]


def tournaments(n):
    prs = [(i, j) for j in range(1, n) for i in range(j)]
    for code in range(1 << len(prs)):
        out = [0] * n
        for t, (i, j) in enumerate(prs):
            if code >> t & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
        yield out


def hp_dfs(out):
    n = len(out)
    full = (1 << n) - 1
    total = 0
    stack = [(v, 1 << v) for v in range(n)]
    while stack:
        v, used = stack.pop()
        if used == full:
            total += 1
            continue
        c = out[v] & ~used
        while c:
            low = c & -c
            c ^= low
            stack.append((low.bit_length() - 1, used | low))
    return total


def packings(out):
    """all sigma with sigma(i) = i or i -> sigma(i); returns list of (n_cycles, n_even_cycles, covered mask)."""
    n = len(out)
    res = []
    sigma = [0] * n

    def rec(i, used):
        if i == n:
            seen = 0
            ncyc = 0
            nev = 0
            cov = 0
            for s in range(n):
                if seen >> s & 1 or sigma[s] == s:
                    continue
                length = 0
                t = s
                while not seen >> t & 1:
                    seen |= 1 << t
                    cov |= 1 << t
                    t = sigma[t]
                    length += 1
                ncyc += 1
                if length % 2 == 0:
                    nev += 1
            res.append((ncyc, nev, cov))
            return
        if not used >> i & 1:
            sigma[i] = i
            rec(i + 1, used | (1 << i))
        c = out[i] & ~used
        while c:
            low = c & -c
            c ^= low
            sigma[i] = low.bit_length() - 1
            rec(i + 1, used | low)

    rec(0, 0)
    return res


def ryser_permanent(M):
    n = len(M)
    total = 0
    for S in range(1, 1 << n):
        prod = 1
        for i in range(n):
            s = 0
            row = M[i]
            for j in range(n):
                if S >> j & 1:
                    s += row[j]
            prod *= s
            if prod == 0:
                break
        total += (-1) ** (n - bin(S).count("1")) * prod
    return total


def det_exact(M):
    n = len(M)
    A = [[Fraction(x) for x in row] for row in M]
    det = Fraction(1)
    for c in range(n):
        piv = None
        for r in range(c, n):
            if A[r][c] != 0:
                piv = r
                break
        if piv is None:
            return 0
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            det = -det
        det *= A[c][c]
        for r in range(c + 1, n):
            if A[r][c] != 0:
                f = A[r][c] / A[c][c]
                for k in range(c, n):
                    A[r][k] -= f * A[c][k]
    return int(det)


def charpoly_xI_plus(M):
    """coefficients c_0..c_n (c_k = coefficient of x^(n-k)) of det(xI + M) over Z, Faddeev-LeVerrier on -M."""
    n = len(M)
    N = [[-x for x in row] for row in M]
    coeffs = [1]
    Mk = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    for k in range(1, n + 1):
        # NM = N * Mk
        NM = [[sum(N[i][t] * Mk[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
        tr = sum(NM[i][i] for i in range(n))
        assert tr % k == 0
        c = -tr // k
        coeffs.append(c)
        Mk = [[NM[i][j] + (c if i == j else 0) for j in range(n)] for i in range(n)]
    return coeffs


def rank_gf2(rows):
    rows = list(rows)
    rank = 0
    while rows:
        piv = rows.pop()
        if piv == 0:
            continue
        rank += 1
        low = piv & -piv
        rows = [r ^ piv if r & low else r for r in rows]
    return rank


def mul_gf2(A, B):
    res = []
    for a in A:
        r = 0
        k = 0
        while a:
            if a & 1:
                r ^= B[k]
            a >>= 1
            k += 1
        res.append(r)
    return res


def phi3_nullity(M):
    n = len(M)
    M2 = mul_gf2(M, M)
    return n - rank_gf2([M2[i] ^ M[i] ^ (1 << i) for i in range(n)])


def rank_f4(M):
    n = len(M)
    A = [row[:] for row in M]
    rank = 0
    for c in range(n):
        piv = None
        for r in range(rank, n):
            if A[r][c]:
                piv = r
                break
        if piv is None:
            continue
        A[rank], A[piv] = A[piv], A[rank]
        inv = INV4[A[rank][c]]
        A[rank] = [MUL4[inv][x] for x in A[rank]]
        for r in range(n):
            if r != rank and A[r][c]:
                f = A[r][c]
                A[r] = [x ^ MUL4[f][y] for x, y in zip(A[r], A[rank])]
        rank += 1
    return rank


def polymul_mod2(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                if b:
                    r[i + j] ^= 1
    return r


def main():
    t0 = time.time()
    ok = True
    first_mon = None
    first_adj = None
    for n in range(3, 7):
        bad = Counter()
        tab = defaultdict(Counter)
        cnt = 0
        cex_mon = Counter()
        cex_adj = Counter()
        ex_mon = None
        ex_adj = None
        for out in tournaments(n):
            cnt += 1
            A = [[(out[i] >> j) & 1 for j in range(n)] for i in range(n)]
            IA = [[A[i][j] + (1 if i == j else 0) for j in range(n)] for i in range(n)]
            pk = packings(out)
            H = hp_dfs(out)
            # (a)
            if ryser_permanent(IA) != len(pk):
                bad["per"] += 1
            if det_exact(IA) != sum((-1) ** ev for (k, ev, cov) in pk):
                bad["det"] += 1
            if sum(2 ** k for (k, ev, cov) in pk if ev == 0) != H:
                bad["ocf"] += 1
            if (H - sum((-1) ** k for (k, ev, cov) in pk if ev == 0)) % 3:
                bad["mod3"] += 1
            # (b) coefficientwise over F_2, then the cubic criterion
            cp = charpoly_xI_plus(A)
            pi = [0] * (n + 1)
            for (k, ev, cov) in pk:
                pi[bin(cov).count("1")] ^= 1
            if [c % 2 for c in cp] != pi:
                bad["charpoly_A"] += 1
            P = [0, 0, 0]
            for (k, ev, cov) in pk:
                P[bin(cov).count("1") % 3] ^= 1
            nulA = phi3_nullity(out)
            Aw = [[A[i][j] ^ (2 if i == j else 0) for j in range(n)] for i in range(n)]
            kerw = n - rank_f4(Aw)
            if nulA != 2 * kerw:
                bad["F4_A"] += 1
            if (nulA > 0) != (P[0] == P[1] == P[2]):
                bad["crit_A"] += 1
            # (c) Monsky convention
            par = [bin(out[i]).count("1") & 1 for i in range(n)]
            Mrows = [out[i] ^ (par[i] << i) for i in range(n)]
            Mint = [[A[i][j] + (par[i] if i == j else 0) for j in range(n)] for i in range(n)]
            cpM = [c % 2 for c in charpoly_xI_plus(Mint)]          # c_k = coeff of x^(n-k)
            rhs = [0] * (n + 1)                                    # indexed by degree
            Q = [0, 0, 0]
            for (k, ev, cov) in pk:
                e = sum(1 for v in range(n) if not cov >> v & 1 and par[v] == 0)
                o = sum(1 for v in range(n) if not cov >> v & 1 and par[v] == 1)
                Q[(e - o) % 3] ^= 1
                poly = [0] * e + [1]                               # x^e
                for _ in range(o):
                    poly = polymul_mod2(poly, [1, 1])              # (1 + x)
                for d, c in enumerate(poly):
                    rhs[d] ^= c
            if [cpM[n - d] for d in range(n + 1)] != rhs:
                bad["charpoly_M"] += 1
            nulM = phi3_nullity(Mrows)
            Mw = [row[:] for row in A]
            for i in range(n):
                Mw[i][i] = par[i] ^ 2                              # par_i + omega in F_4
            if nulM != 2 * (n - rank_f4(Mw)):
                bad["F4_M"] += 1
            if (nulM > 0) != (Q[0] == Q[1] == Q[2]):
                bad["crit_M"] += 1
            # orientation independence: converse tournament
            full = (1 << n) - 1
            conv = [full & ~out[i] & ~(1 << i) for i in range(n)]
            parc = [bin(conv[i]).count("1") & 1 for i in range(n)]
            if phi3_nullity([conv[i] ^ (parc[i] << i) for i in range(n)]) != nulM:
                bad["orientation_M"] += 1
            if phi3_nullity(conv) != nulA:
                bad["orientation_A"] += 1
            # bridges
            corank = n - 1 - rank_gf2(Mrows)
            nulIA = n - rank_gf2([out[i] ^ (1 << i) for i in range(n)])
            tab[(corank, nulM, nulIA)][H % 8] += 1
            if nulM > 0 and H % 3:
                cex_mon[H] += 1
                if ex_mon is None:
                    ex_mon = (list(out), H)
            if nulA > 0 and H % 3:
                cex_adj[H] += 1
                if ex_adj is None:
                    ex_adj = (list(out), H)
        good = sum(bad.values()) == 0
        ok &= good
        print(f"n={n}: {cnt} labelled tournaments; identity failures {dict(bad)}", "OK" if good else "FAIL",
              f"[{time.time()-t0:.0f}s]", flush=True)
        non4 = sum(1 for key in tab if len({h % 4 for h in tab[key]}) > 1)
        non8 = sum(1 for key in tab if len(tab[key]) > 1)
        print(f"   classes of (Laplacian corank - 1, nullity(M^2+M+I), nullity(I+A)): {len(tab)}; H mod 4 not constant on "
              f"{non4}; H mod 8 not constant on {non8}")
        print(f"   'nullity(M^2+M+I) > 0 => 3 | H': counterexamples by H: {dict(sorted(cex_mon.items()))}; first in this "
              f"enumeration: {ex_mon}")
        print(f"   'nullity(A^2+A+I) > 0 => 3 | H': counterexamples by H: {dict(sorted(cex_adj.items()))}; first in this "
              f"enumeration: {ex_adj}")
        if cex_mon and first_mon is None:
            first_mon = (n, min(cex_mon))
        if cex_adj and first_adj is None:
            first_adj = (n, min(cex_adj))
    print("first order with a Monsky counterexample, least H:", first_mon, "; pure adjacency:", first_adj)
    ok &= first_mon == (5, 13) and first_adj is not None and first_adj[0] == 6

    # the packet's two named tournaments
    t5 = [18, 8, 3, 5, 14]
    t6 = [32, 1, 3, 7, 15, 30]
    for t in (t5, t6):
        n = len(t)
        assert all(((t[i] >> j) & 1) + ((t[j] >> i) & 1) == 1 for i in range(n) for j in range(i + 1, n))
        par = [bin(t[i]).count("1") & 1 for i in range(n)]
        print(f"   named tournament {t}: H = {hp_dfs(t)}, scores {[bin(x).count('1') for x in t]}, nullity(M^2+M+I) = "
              f"{phi3_nullity([t[i] ^ (par[i] << i) for i in range(n)])}, nullity(A^2+A+I) = {phi3_nullity(t)}")
    ok &= hp_dfs(t5) == 13 and hp_dfs(t6) == 17
    print("AUDIT CYCLE GAS:", "ALL CHECKS PASSED" if ok else "SOME CHECK FAILED", f"[{time.time()-t0:.0f}s]")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
