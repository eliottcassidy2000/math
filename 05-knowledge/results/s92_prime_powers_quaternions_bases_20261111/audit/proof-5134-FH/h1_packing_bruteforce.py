#!/usr/bin/env python3
"""Auditor proof-5134-FH.  THM-5134 (H)(i): brute-force odd-cycle packing polynomials.

Independent code, written from the definitions (THM-4940 Setting, THM-5070 Setting and (2)) before
reading the session's every_base.py.

I(Omega(T), z) = sum over sets of pairwise vertex-disjoint directed odd cycles of T of z^(number of cycles).

Two independent evaluations:
  (E) explicit: list every directed odd cycle by depth-first search (each cycle once, rooted at its least
      vertex), then enumerate all sets of pairwise disjoint cycles by backtracking;
  (D) subset DP: cyc[S] = number of directed Hamiltonian cycles of T[S] (|S| odd), then
      P[mask] = P[mask - v] + z * sum_{S in mask, v in S} cyc[S] P[mask - S], v the least vertex of mask.
Sanity: I(Omega(T), 2) = number of Hamiltonian paths of T (OCF), computed by a third DP.

Constructions (my reading of THM-5070 Setting and (2)):
  T(W) = TT_k + u : chain c_1 -> c_2 -> ... -> c_k (c_i -> c_j for i < j), k = bit length of W;
         u -> c_i iff the i-th bit of W, most significant first, is 1; otherwise c_i -> u.
  B(W') = TT_l + w : chain d_1, ..., d_l (d_i -> d_j for i < j), l = bit length of W';
         d_i -> w iff bit i-1 of W' (least significant first) is 1; otherwise w -> d_i.
  T (.) B : every vertex of T beats every vertex of B, except w -> u.
F1 coordinates: W = 2^s1 (2 h1 + 1) - 1, W' = 2^s2 (2 h2 + 1) - 1, s = s1 + s2  (so W = 2 h1 when s1 = 0).
  R(h) = T(2h);  J(h1, h2, s) = T(W) (.) B(W').

Claims tested:
  (a) I(Omega(T(W)), z) = 1 + a(W) z,  a(0) = 0, a(2W+1) = a(W), a(2W) = a(W) + ceil(W/2);  a(2h) = psi(h);
      same for B(W').
  (b) THM-5070 (2):  I(T(W) (.) B(W')) = (1 + z a(W))(1 + z a(W')) + z [ceil(W/2)(1 + floor(W'/2)) + (1 + floor(W/2)) ceil(W'/2)].
  (c) THM-5134 (H):  I(J(h1,h2,s)) = (1 + psi1 z)(1 + psi2 z) + z (2^s y1 y2 - [s = 0])/2.
"""
import sys
import random
import numpy as np
from numba import njit

sys.setrecursionlimit(100000)


# ---------- arithmetic ----------
def psi(h):
    if h == 0:
        return 0
    if h % 2 == 0:
        return psi(h // 2) + h // 2
    k = h // 2
    return psi(k) + k // 2 + 1


def a_word(W):
    if W == 0:
        return 0
    if W % 2 == 1:
        return a_word(W // 2)
    return a_word(W // 2) + (W // 2 + 1) // 2


# ---------- constructions ----------
def build_T(W):
    """returns (adj, port) ; vertices 0..k-1 = chain c_1..c_k, vertex k = u."""
    k = W.bit_length()
    n = k + 1
    adj = [[0] * n for _ in range(n)]
    for i in range(k):
        for j in range(i + 1, k):
            adj[i][j] = 1
    u = k
    for i in range(k):
        bit = (W >> (k - 1 - i)) & 1          # most significant first
        if bit:
            adj[u][i] = 1
        else:
            adj[i][u] = 1
    return adj, u


def build_B(Wp):
    """returns (adj, port) ; vertices 0..l-1 = chain d_1..d_l, vertex l = w."""
    l = Wp.bit_length()
    n = l + 1
    adj = [[0] * n for _ in range(n)]
    for i in range(l):
        for j in range(i + 1, l):
            adj[i][j] = 1
    w = l
    for i in range(l):
        bit = (Wp >> i) & 1                    # least significant first
        if bit:
            adj[i][w] = 1
        else:
            adj[w][i] = 1
    return adj, w


def glue(adjT, u, adjB, w):
    n1, n2 = len(adjT), len(adjB)
    n = n1 + n2
    adj = [[0] * n for _ in range(n)]
    for i in range(n1):
        for j in range(n1):
            adj[i][j] = adjT[i][j]
    for i in range(n2):
        for j in range(n2):
            adj[n1 + i][n1 + j] = adjB[i][j]
    for i in range(n1):
        for j in range(n2):
            adj[i][n1 + j] = 1
    adj[u][n1 + w] = 0
    adj[n1 + w][u] = 1
    return adj


def is_tournament(adj):
    n = len(adj)
    for i in range(n):
        if adj[i][i]:
            return False
        for j in range(i + 1, n):
            if adj[i][j] + adj[j][i] != 1:
                return False
    return True


# ---------- (E) explicit enumeration ----------
def odd_cycles(adj):
    """all directed odd cycles, each once, as vertex bitmasks (with repetition for different cycles on one set)."""
    n = len(adj)
    out = []
    nb = [[j for j in range(n) if adj[i][j]] for i in range(n)]

    def dfs(start, v, mask, length):
        for x in nb[v]:
            if x == start:
                if length >= 3 and length % 2 == 1:
                    out.append(mask)
            elif x > start and not (mask >> x) & 1:
                dfs(start, x, mask | (1 << x), length + 1)

    for s in range(n):
        dfs(s, s, 1 << s, 1)
    return out


def packing_poly_explicit(adj):
    cyc = odd_cycles(adj)
    cyc.sort()
    m = len(cyc)
    coef = [0] * (len(adj) // 3 + 2)

    def rec(i, used, k):
        coef[k] += 1
        for j in range(i, m):
            if not (cyc[j] & used):
                rec(j + 1, used | cyc[j], k + 1)

    rec(0, 0, 0)
    while len(coef) > 1 and coef[-1] == 0:
        coef.pop()
    return coef, m


# ---------- (D) subset DP ----------
@njit(cache=True)
def ham_cycle_counts(adj, n):
    """cyc[S] = number of directed Hamiltonian cycles of T[S]."""
    full = 1 << n
    cyc = np.zeros(full, dtype=np.int64)
    for s in range(n):
        # cycles whose least vertex is s : paths from s inside vertices > s
        size = 1 << n
        dp = np.zeros((size, n), dtype=np.int64)
        dp[1 << s, s] = 1
        for mask in range(size):
            if not (mask >> s) & 1:
                continue
            if mask & ((1 << s) - 1):
                continue
            for v in range(s, n):
                c = dp[mask, v]
                if c == 0:
                    continue
                if adj[v, s] and mask != (1 << s):
                    cyc[mask] += c
                for x in range(s + 1, n):
                    if adj[v, x] and not (mask >> x) & 1:
                        dp[mask | (1 << x), x] += c
    return cyc


@njit(cache=True)
def packing_poly_dp(adj, n):
    cyc = ham_cycle_counts(adj, n)
    full = 1 << n
    deg = n // 3 + 1
    P = np.zeros((full, deg + 1), dtype=np.int64)
    P[0, 0] = 1
    pop = np.zeros(full, dtype=np.int64)
    for mask in range(1, full):
        pop[mask] = pop[mask >> 1] + (mask & 1)
    for mask in range(1, full):
        low = mask & (-mask)
        rest = mask ^ low
        for d in range(deg + 1):
            P[mask, d] = P[rest, d]
        # submasks S of mask containing low, |S| odd >= 3
        sub = rest
        while True:
            S = sub | low
            if pop[S] >= 3 and pop[S] % 2 == 1 and cyc[S] != 0:
                c = cyc[S]
                rem = mask ^ S
                for d in range(deg):
                    P[mask, d + 1] += c * P[rem, d]
            if sub == 0:
                break
            sub = (sub - 1) & rest
    return P[full - 1]


@njit(cache=True)
def ham_paths(adj, n):
    full = 1 << n
    dp = np.zeros((full, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for mask in range(full):
        for v in range(n):
            c = dp[mask, v]
            if c == 0:
                continue
            for x in range(n):
                if adj[v, x] and not (mask >> x) & 1:
                    dp[mask | (1 << x), x] += c
    tot = 0
    for v in range(n):
        tot += dp[full - 1, v]
    return tot


def poly_dp(adj):
    n = len(adj)
    A = np.array(adj, dtype=np.int64)
    c = [int(x) for x in packing_poly_dp(A, n)]
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    return c


def poly_mul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            r[i + j] += x * y
    return r


def trim(c):
    c = list(c)
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    return c


def poly_add(p, q):
    r = [0] * max(len(p), len(q))
    for i, x in enumerate(p):
        r[i] += x
    for i, x in enumerate(q):
        r[i] += x
    return trim(r)


def evalp(c, z):
    return sum(x * z ** i for i, x in enumerate(c))


def main():
    out = []
    P = out.append
    random.seed(92092)
    P("THM-5134 (H)(i): brute-force packing polynomials (auditor proof-5134-FH)")
    P("=" * 78)

    # ---- 0. the two evaluations agree, and OCF, on random tournaments ----
    nfail = 0
    ntest = 0
    for n in range(3, 10):
        for _ in range(40 if n <= 8 else 15):
            adj = [[0] * n for _ in range(n)]
            for i in range(n):
                for j in range(i + 1, n):
                    if random.random() < 0.5:
                        adj[i][j] = 1
                    else:
                        adj[j][i] = 1
            ce, m = packing_poly_explicit(adj)
            cd = poly_dp(adj)
            hp = int(ham_paths(np.array(adj, dtype=np.int64), n))
            ntest += 1
            if ce != cd or evalp(cd, 2) != hp:
                nfail += 1
    P(f"0. random tournaments of orders 3..9: explicit enumeration = subset DP, and I(Omega, 2) = #Hamiltonian paths (OCF): "
      f"{ntest} tested, {nfail} failures")

    # ---- (a) one-apex ----
    fail_a = 0
    cnt_a = 0
    cnt_a_explicit = 0
    for W in range(0, 2 ** 11):                      # T(W): up to 12 vertices
        adjT, u = build_T(W)
        assert is_tournament(adjT)
        c = poly_dp(adjT)
        want = trim([1, a_word(W)])
        cnt_a += 1
        if c != want:
            fail_a += 1
        if len(adjT) <= 9:
            ce, _ = packing_poly_explicit(adjT)
            cnt_a_explicit += 1
            if ce != want:
                fail_a += 1
        adjB, w = build_B(W)
        cb = poly_dp(adjB)
        if cb != want:
            fail_a += 1
    P(f"(a) I(Omega(T(W))) = I(Omega(B(W))) = 1 + a(W) z for every W < 2^11 (orders <= 12): {cnt_a} words, "
      f"{cnt_a_explicit} also by explicit enumeration; failures {fail_a}")
    bad = [h for h in range(0, 2 ** 14) if a_word(2 * h) != psi(h)]
    P(f"    a(2h) = psi(h) for every h < 2^14: {len(bad) == 0};   R(h) = T(2h): I = 1 + psi(h) z  (h < 2^10 covered above)")
    # trailing ones do not change a
    bad2 = [(h, s) for h in range(0, 512) for s in range(0, 6) if a_word(2 ** s * (2 * h + 1) - 1) != psi(h)]
    P(f"    a(2^s (2h+1) - 1) = psi(h) for h < 512, s < 6: {len(bad2) == 0}")

    # ---- (b), (c) joins ----
    def test_join(W, Wp, explicit=False):
        adjT, u = build_T(W)
        adjB, w = build_B(Wp)
        adj = glue(adjT, u, adjB, w)
        assert is_tournament(adj)
        n = len(adj)
        c = poly_dp(adj)
        cW = (W + 1) // 2 * (1 + Wp // 2) + (1 + W // 2) * ((Wp + 1) // 2)
        want = poly_add(poly_mul([1, a_word(W)], [1, a_word(Wp)]), [0, cW])
        ok = (c == want)
        if explicit:
            ce, _ = packing_poly_explicit(adj)
            ok = ok and (ce == want)
        hp = int(ham_paths(np.array(adj, dtype=np.int64), n))
        ok = ok and (evalp(c, 2) == hp)
        return ok, n, c

    # (b) all word pairs with total order <= 10, explicit too for order <= 9
    cnt_b = 0
    fail_b = 0
    by_order = {}
    for W in range(0, 2 ** 8):
        for Wp in range(0, 2 ** 8):
            n = W.bit_length() + Wp.bit_length() + 2
            if n > 10:
                continue
            ok, n, c = test_join(W, Wp, explicit=(n <= 9))
            cnt_b += 1
            by_order[n] = by_order.get(n, 0) + 1
            if not ok:
                fail_b += 1
    P(f"(b) THM-5070 (2) on ALL word pairs (W, W') with total order <= 10: {cnt_b} joins, failures {fail_b}; "
      f"by order {dict(sorted(by_order.items()))}")
    # random pairs of orders 11, 12, 13
    cnt_b2 = 0
    fail_b2 = 0
    by_order2 = {}
    tries = 0
    while cnt_b2 < 600 and tries < 100000:
        tries += 1
        n = random.choice([11, 12, 12, 13])
        k = random.randint(0, n - 2)
        l = n - 2 - k
        W = random.randint(2 ** (k - 1), 2 ** k - 1) if k >= 1 else 0
        Wp = random.randint(2 ** (l - 1), 2 ** l - 1) if l >= 1 else 0
        ok, nn, c = test_join(W, Wp)
        assert nn == n
        cnt_b2 += 1
        by_order2[n] = by_order2.get(n, 0) + 1
        if not ok:
            fail_b2 += 1
    P(f"    random word pairs of orders 11-13: {cnt_b2} joins, failures {fail_b2}; by order {dict(sorted(by_order2.items()))}")

    # (c) F1 coordinates
    cnt_c = 0
    fail_c = 0
    by_s = {}
    by_order3 = {}
    examples = []
    for h1 in range(0, 40):
        for h2 in range(0, 40):
            for s in range(0, 5):
                for s1 in range(0, s + 1):
                    s2 = s - s1
                    W = 2 ** s1 * (2 * h1 + 1) - 1
                    Wp = 2 ** s2 * (2 * h2 + 1) - 1
                    n = W.bit_length() + Wp.bit_length() + 2
                    if n > 12:
                        continue
                    ok, n, c = test_join(W, Wp, explicit=(n <= 9))
                    y1, y2 = 2 * h1 + 1, 2 * h2 + 1
                    p1, p2 = psi(h1), psi(h2)
                    seam = (2 ** s * y1 * y2 - (1 if s == 0 else 0))
                    assert seam % 2 == 0
                    want = poly_add(poly_mul([1, p1], [1, p2]), [0, seam // 2])
                    ok = ok and (c == want)
                    # N_b at b = 2: x1 x2 + 2^s y1 y2 - [s=0]
                    x1, x2 = 2 * p1 + 1, 2 * p2 + 1
                    ok = ok and (evalp(c, 2) == x1 * x2 + 2 ** s * y1 * y2 - (1 if s == 0 else 0))
                    for b in (3, 4, 5, 10):
                        X1, X2 = 1 + b * p1, 1 + b * p2
                        ok = ok and (evalp(c, b) == X1 * X2 + b * seam // 2)
                    cnt_c += 1
                    by_s[s] = by_s.get(s, 0) + 1
                    by_order3[n] = by_order3.get(n, 0) + 1
                    if not ok:
                        fail_c += 1
                    if len(examples) < 6 and h1 >= 3 and h2 >= 2 and s == len(examples) % 3:
                        examples.append((h1, h2, s, s1, n, c))
    P(f"(c) THM-5134 (H) formula on every (h1, h2, s; split s1 + s2 = s) with h1, h2 < 40, s <= 4, order <= 12: "
      f"{cnt_c} joins, failures {fail_c}")
    P(f"    by s: {dict(sorted(by_s.items()))}; by order: {dict(sorted(by_order3.items()))}")
    P("    (each also checked: I(2) = x1 x2 + 2^s y1 y2 - [s=0] = #Hamiltonian paths, and I(b) = X1 X2 + b (2^s y1 y2 - [s=0])/2 for b = 3, 4, 5, 10)")
    for (h1, h2, s, s1, n, c) in examples:
        P(f"    example: h1 = {h1}, h2 = {h2}, s = {s} (s1 = {s1}), order {n}: I = {c}  "
          f"[psi = {psi(h1)}, {psi(h2)}; (2^s y1 y2 - [s=0])/2 = {(2 ** s * (2 * h1 + 1) * (2 * h2 + 1) - (1 if s == 0 else 0)) // 2}]")
    # distinct tournaments at orders 11, 12 with s = 0, 1, 2
    c1112 = sum(v for k, v in by_order3.items() if k >= 11)
    P(f"    joins of order 11 or 12 among them: {c1112}")
    tot_fail = nfail + fail_a + fail_b + fail_b2 + fail_c + len(bad) + len(bad2)
    P("")
    P("ALL CHECKS PASS" if tot_fail == 0 else f"FAILURES: {tot_fail}")
    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
