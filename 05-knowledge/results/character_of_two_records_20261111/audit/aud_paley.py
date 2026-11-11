"""AUDIT: the Paley theorem of section 4 (THM-5136 (P)), recomputed independently.

usage: python3 aud_paley.py /path/to/aud_ham   (compiled aud_ham.c)
"""
import sys
from aud_common import (gf_field, paley, fold, paths_ending_in, ham_exact, ham_python, check_tournament)

HAM = sys.argv[1] if len(sys.argv) > 1 else None


def cycles_by_definition(A, psi):
    """Enumerate every directed odd cycle; return ({length: count}, {length: invariant count})."""
    n = len(A)
    tot, inv = {}, {}
    path = []
    on = [False] * n

    def test():
        L = len(path)
        tot[L] = tot.get(L, 0) + 1
        succ = {path[i]: path[(i + 1) % L] for i in range(L)}
        for u, w in succ.items():
            if succ.get(psi[w]) != psi[u]:
                return
        inv[L] = inv.get(L, 0) + 1

    def dfs(root):
        u = path[-1]
        if len(path) >= 3 and A[u][root]:
            test()
        for w in range(root + 1, n):
            if not on[w] and A[u][w]:
                on[w] = True
                path.append(w)
                dfs(root)
                path.pop()
                on[w] = False

    for root in range(n):
        path[:] = [root]
        on[root] = True
        dfs(root)
        on[root] = False
    return tot, inv


def halfpath_counts(A, psi, v0):
    """F_k = number of (x_1..x_k): v0 -> x_1, directed path, at most one vertex per pair, x_k -> psi(x_k)."""
    n = len(A)
    pair = {}
    idx = 0
    for a in range(n):
        if a != v0 and a not in pair:
            pair[a] = pair[psi[a]] = idx
            idx += 1
    m = idx
    layer = {}
    for a in range(n):
        if a != v0 and A[v0][a]:
            layer[(1 << pair[a], a)] = 1
    F = [0] * (m + 1)
    for k in range(1, m + 1):
        nxt = {}
        for (mask, u), x in layer.items():
            if A[u][psi[u]]:
                F[k] += x
            for w in range(n):
                if w != v0 and A[u][w] and not (mask >> pair[w]) & 1:
                    key = (mask | (1 << pair[w]), w)
                    nxt[key] = nxt.get(key, 0) + x
        layer = nxt
    return F


def binom(n, k):
    from math import comb
    return comb(n, k)


known_H = {3: 3, 7: 189, 11: 95095, 19: 1172695746915, 23: 15760206976379349}

for q in (3, 7, 11, 19, 23, 27, 31):
    els, add, mul, neg = gf_field(q) if q != 31 else gf_field(31)
    A = paley(q)
    n = q
    psi = tuple(neg(a) for a in els)
    Z, D1, plus, minus, v0 = fold(A, psi)
    m = (q - 1) // 2
    squares = sorted(set(mul(a, a) for a in els if a != 0))
    assert v0 == 0 and sorted(plus) == squares
    # (1) D_1 = everything iff q = 3 mod 8, empty iff q = 7 mod 8
    d1_all, d1_none = all(D1), not any(D1)
    assert (q % 8 == 3 and d1_all) or (q % 8 == 7 and d1_none), (q, D1)
    two = add(1, 1)
    chi2 = 1 if two in squares else -1
    assert (chi2 == -1) == (q % 8 == 3)
    # (2) a -> a^2 is an isomorphism from Z onto P_q[squares]
    iso = all(Z[i][j] == A[mul(plus[i], plus[i])][mul(plus[j], plus[j])] for i in range(m) for j in range(m) if i != j)
    sq_bij = sorted(mul(a, a) for a in plus) == squares
    # H(Z)
    HZ = ham_python(Z) if m <= 15 else None
    zp = paths_ending_in(Z, D1)
    # F_k
    F = halfpath_counts(A, psi, v0)
    if q <= 11:
        tot, inv = cycles_by_definition(A, psi)
        Fdef = [0] + [inv.get(2 * k + 1, 0) for k in range(1, m + 1)]
        assert Fdef == F, (q, Fdef, F)
        pair_ok = all((tot.get(2 * k + 1, 0) - Fdef[k]) % 2 == 0 for k in range(1, m + 1))
        cyc_txt = "; cycles by definition: c_odd lengths %s, invariant %s (half-path DP agrees: True; c = F mod 2: %s)" % (
            [tot.get(2 * k + 1, 0) for k in range(1, m + 1)], Fdef[1:], pair_ok)
    else:
        cyc_txt = ""
    par = "".join(str(F[k] % 2) for k in range(1, m + 1))
    pred = "".join(str((binom(m, k) if q % 8 == 3 else 0) % 2) for k in range(1, m + 1))
    zpar = "".join(str(zp[k] % 2) for k in range(1, m + 1))
    Ftot = sum(F)
    line = "q = %2d (q mod 8 = %d, chi(2) = %+d): |D_1| = %d of %d; Z = P_q[squares] via a -> a^2: %s (squaring bijective on squares: %s); H(Z) = %s" % (
        q, q % 8, chi2, sum(D1), m, iso, sq_bij, HZ)
    print(line)
    print("   F_k mod 2 = %s ; predicted [q=3 mod 8] C(m,k) = %s ; paths of Z ending in D_1 mod 2 = %s ; agree: %s" % (
        par, pred, zpar, par == pred == zpar))
    print("   F = %d ; F mod 4 = %d ; (q+1)/4 mod 4 = %d ; HYP (c) predicts %d%s" % (
        Ftot, Ftot % 4, ((q + 1) // 4) % 4, (((q + 1) // 4) % 4 if q % 8 == 3 else 0), cyc_txt))
    if HAM and q <= 23:
        H = ham_exact(A, HAM)
        e = q * (q - 1) // 2
        K = H // e
        print("   H(P_q) = %d (repository value %s: %s); H mod 4 = %d, mod 8 = %d; q(q-1)/2 mod 8 = %d; 2 - chi(2) = %d; H = 1 + 2F mod 4: %s; K_q = H/(q(q-1)/2) = %d (exact division: %s), K_q mod 4 = %d, mod 8 = %d" % (
            H, known_H.get(q), H == known_H.get(q), H % 4, H % 8, e % 8, 2 - chi2, H % 4 == (1 + 2 * Ftot) % 4, K, H % e == 0, K % 4, K % 8))
    if q <= 11:
        # paths starting with the arc 0 -> 1, directly
        nn = len(A)
        cnt = 0
        def ext(last, used):
            global cnt
            if used == (1 << nn) - 1:
                cnt += 1
                return
            for w in range(nn):
                if not (used >> w) & 1 and A[last][w]:
                    ext(w, used | (1 << w))
        ext(1, 0b11)
        print("   paths starting 0 -> 1 (direct DFS): %d" % cnt)
