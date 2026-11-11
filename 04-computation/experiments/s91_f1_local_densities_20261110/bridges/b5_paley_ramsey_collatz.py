#!/usr/bin/env python3
"""S91 bridges (5): the Paley / Ramsey / Collatz identities of understand_sweep_numbers.txt items (f)1-(f)5, own code.

All checks are exact integer computations.  CITED (not re-proved here): R(5) = 14 with the unique extremal ST_13
(Reid-Parker 1970) and R(6) = 28 with the unique extremal ST_27 (Sanchez-Flores 1994); the script verifies the
lower-bound side (ST_13 has no transitive 5-subtournament, P_27 has no transitive 6-subtournament) and, exhaustively,
R(3) = 4 and R(4) = 8 with the uniqueness of P_7.

Conventions: Paley tournament P_q (q = 3 mod 4): i -> j iff j - i is a nonzero square of F_q.  F_27 = F_3[t]/(t^3 - t - 1).
Skew matrix S: S[i][j] = +1 if i -> j, -1 if j -> i, 0 on the diagonal.  Syracuse map U(n) = odd part of 3n + 1.
Legendre tournament on primes = 3 mod 4: p -> q iff p is a quadratic residue modulo q (THM-5125).
"""
import itertools
import sys

import numpy as np
from numba import njit

RESULTS = []


def check(name, cond, detail=""):
    RESULTS.append(bool(cond))
    print(("PASS  " if cond else "FAIL  ") + name + (("   [" + detail + "]") if detail else ""), flush=True)


def order(a, m):
    k, v = 1, a % m
    while v != 1:
        v = v * a % m
        k += 1
    return k


def factor(n):
    out = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


# ---------------------------------------------------------------- 1. orders and factorisations
check("1a ord_27(2) = ord_19(2) = 18; 2 generates both unit groups (phi(27) = phi(19) = 18)", order(2, 27) == 18 and order(2, 19) == 18)
check("1b 2^9 + 1 = 513 = 27 * 19, 2^18 - 1 = 27 * 7 * 19 * 73, 2^6 - 1 = 9 * 7, 19 = 3^3 - 2^3",
      2 ** 9 + 1 == 27 * 19 and 2 ** 18 - 1 == 27 * 7 * 19 * 73 and 2 ** 6 - 1 == 63 and 19 == 27 - 8)
check("1c the primes p with ord_p(2) | 18 are exactly 3, 7, 19, 73", sorted(factor(2 ** 18 - 1)) == [3, 7, 19, 73])
iso = {pow(2, a, 27): pow(2, a, 19) for a in range(18)}
check("1d 2^a mod 27 -> 2^a mod 19 is a well defined group isomorphism (Z/27)^* -> (Z/19)^*",
      len(iso) == 18 and len(set(iso.values())) == 18 and all(iso[(u * v) % 27] == (iso[u] * iso[v]) % 19 for u in iso for v in iso))
check("1e primitive cube roots of unity: 2 mod 7, 3 mod 13, 3/2 = 11 mod 19, 22/15 = 140 mod 1039",
      pow(2, 3, 7) == 1 and pow(3, 3, 13) == 1 and (3 * pow(2, -1, 19)) % 19 == 11 and pow(11, 3, 19) == 1
      and (22 * pow(15, -1, 1039)) % 1039 == 140 and pow(140, 3, 1039) == 1 and 140 != 1)


# ---------------------------------------------------------------- 2. Collatz
def U(n):
    n = 3 * n + 1
    while n % 2 == 0:
        n //= 2
    return n


def syracuse(n):
    out = [n]
    while n != 1:
        n = U(n)
        out.append(n)
    return out


def collatz(n):
    out = [n]
    while n != 1:
        n = n // 2 if n % 2 == 0 else 3 * n + 1
        out.append(n)
    return out


check("2a Syracuse orbit of 7: 7, 11, 17, 13, 5, 1; so U^3(7) = 13", syracuse(7) == [7, 11, 17, 13, 5, 1])
ok = True
for a in range(1, 41):
    v = 2 ** a - 1
    for _ in range(a):
        v = U(v)
    t = (3 ** a - 1) // 2
    while t % 2 == 0:
        t //= 2
    ok &= (v == t)
check("2b THM-4556: U^a(2^a - 1) = odd part of (3^a - 1)/2 for a = 1..40", ok)
s1039 = syracuse(1039)
check("2c Syracuse orbit of 1039 = 1039, 1559, 2339, 3509, 329, 247, 371, 557, 209, 157, 59, 89, 67, 101, 19, 29, 11, 17, 13, 5, 1",
      s1039 == [1039, 1559, 2339, 3509, 329, 247, 371, 557, 209, 157, 59, 89, 67, 101, 19, 29, 11, 17, 13, 5, 1])
check("2d it passes through 67, 19, 13, 5 and merges with the orbit of 7 at 11 (it contains 11, 17, 13, 5, 1 but not 7)",
      all(v in s1039 for v in (67, 19, 13, 5)) and s1039[-5:] == syracuse(7)[1:] and 7 not in s1039)
o27 = collatz(27)
o1039 = collatz(1039)
meet = [v for v in o1039 if v in set(o27)][0]
s27 = syracuse(27)
check("2e the 3x+1 orbits of 27 (111 steps, maximum 9232) and of 1039 first meet at 40; the Syracuse orbit of 27 avoids 7, 13, 19 and enters 5 from 53",
      len(o27) - 1 == 111 and max(o27) == 9232 and meet == 40 and all(v not in s27 for v in (7, 13, 19)) and s27[-3:] == [53, 5, 1],
      "first common value %d" % meet)
LIM = 200000
targets = [5, 7, 11, 13, 17, 19, 27, 67, 1039]
bit = {t: 1 << i for i, t in enumerate(targets)}
memo = {1: 0}
cnt = {t: 0 for t in targets}
for n in range(1, LIM + 1):
    stack = []
    v = n
    while v not in memo:
        stack.append(v)
        v = v // 2 if v % 2 == 0 else 3 * v + 1
    acc = memo[v]
    for w in reversed(stack):
        if w in bit:
            acc |= bit[w]
        memo[w] = acc
    for t in targets:
        if memo[n] & bit[t]:
            cnt[t] += 1
rates = {t: cnt[t] / LIM for t in targets}
print("      fraction of n <= 200000 whose 3x+1 orbit passes through t:", {t: round(v, 4) for t, v in rates.items()})
check("2f base rates: 17.7 % of n <= 200000 pass through 19 (helper's figure); 93-95 % through 5; about half through 13",
      abs(rates[19] - 0.177) < 0.001 and 0.93 < rates[5] < 0.95 and 0.4 < rates[13] < 0.7)


# ---------------------------------------------------------------- 3. Paley tournaments
def paley_prime(qq):
    sq = {(v * v) % qq for v in range(1, qq)}
    return [[1 if ((j - i) % qq) in sq else 0 for j in range(qq)] for i in range(qq)]


def f27():
    """elements of F_27 as triples (a0, a1, a2) = a0 + a1 t + a2 t^2 with t^3 = t + 1."""
    els = list(itertools.product(range(3), repeat=3))

    def mul(u, v):
        c = [0] * 5
        for i in range(3):
            for j in range(3):
                c[i + j] += u[i] * v[j]
        # t^3 = t + 1, t^4 = t^2 + t
        c[1] += c[3]
        c[0] += c[3]
        c[2] += c[4]
        c[1] += c[4]
        return (c[0] % 3, c[1] % 3, c[2] % 3)

    def sub(u, v):
        return tuple((a - b) % 3 for a, b in zip(u, v))

    return els, mul, sub


els, mul27, sub27 = f27()
zero = (0, 0, 0)
check("3a t^3 - t - 1 is irreducible over F_3 (no root), so the 27 triples form the field F_27",
      all((v ** 3 - v - 1) % 3 != 0 for v in range(3)))
t = (0, 1, 0)
pw = (1, 0, 0)
ords = 0
for k in range(1, 27):
    pw = mul27(pw, t)
    if pw == (1, 0, 0):
        ords = k
        break
squares27 = {mul27(u, u) for u in els if u != zero}
check("3b F_27^* has 13 nonzero squares and -1 is not a square (27 = 3 mod 4); ord(t) = %d divides 26" % ords,
      len(squares27) == 13 and (2, 0, 0) not in squares27 and 26 % ords == 0)
P27 = [[1 if sub27(v, u) in squares27 else 0 for v in els] for u in els]
PAL = {3: paley_prime(3), 7: paley_prime(7), 11: paley_prime(11), 19: paley_prime(19), 27: P27}


def is_tournament(adj):
    n = len(adj)
    return all(adj[i][i] == 0 for i in range(n)) and all(adj[i][j] + adj[j][i] == 1 for i in range(n) for j in range(i + 1, n))


def skew_identity(adj):
    """S S^T = q I - J and S 1 = 0 (exact): the eigenvalues of S are 0 and +- i sqrt(q), each (q-1)/2 times."""
    n = len(adj)
    S = np.array([[adj[i][j] - adj[j][i] for j in range(n)] for i in range(n)], dtype=np.int64)
    return bool(np.array_equal(S @ S.T, n * np.eye(n, dtype=np.int64) - np.ones((n, n), dtype=np.int64))) and bool(np.all(S.sum(axis=1) == 0))


check("3c P_3, P_7, P_11, P_19, P_27 are tournaments", all(is_tournament(PAL[k]) for k in PAL))
check("3d S S^T = q I - J and S 1 = 0 for q = 3, 7, 11, 19, 27: skew eigenvalue moduli sqrt 3, sqrt 7, sqrt 11, sqrt 19, sqrt 27",
      all(skew_identity(PAL[k]) for k in PAL))


def reduced_forms(D):
    out = []
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a) == 0:
                c = (b * b - D) // (4 * a)
                if c >= a and not (b < 0 and a == c):
                    out.append((a, b, c))
        a += 1
    return out


def gcd3(f):
    from math import gcd
    return gcd(gcd(f[0], abs(f[1])), f[2])


h = {D: len([f for f in reduced_forms(D) if gcd3(f) == 1]) for D in (-3, -7, -12, -19, -27, -52, -75, -28, -8, -4)}
check("3e class numbers h(D) = 1 for D = -3, -7, -12, -19, -27 (and -4, -8, -28); h(-52) = 2, h(-75) = 2",
      all(h[D] == 1 for D in (-3, -7, -12, -19, -27, -4, -8, -28)) and h[-52] == 2 and h[-75] == 2, str(h))
check("3f (q + 1)/4 = 1, 2, 5, 7 for q = 3, 7, 19, 27; 13 = 1 mod 4 has no Paley tournament",
      [(qq + 1) // 4 for qq in (3, 7, 19, 27)] == [1, 2, 5, 7] and 13 % 4 == 1)


# ---------------------------------------------------------------- 4. Ramsey numbers for transitive subtournaments
def max_transitive(adj):
    """size of the largest transitive subtournament (exhaustive search by extension of acyclic vertex sets)."""
    n = len(adj)
    out = [sum(1 << j for j in range(n) if adj[i][j]) for i in range(n)]
    best = 0
    # a transitive subtournament is a chain v1 -> v2 -> ... with every later vertex dominated by all earlier ones
    stack = [(1 << i, out[i], 1) for i in range(n)]
    while stack:
        used, cand, size = stack.pop()
        if size > best:
            best = size
        c = cand
        while c:
            low = c & (-c)
            j = low.bit_length() - 1
            stack.append((used | low, cand & out[j], size + 1))
            c ^= low
    return best


C3 = paley_prime(3)
ST13 = [[1 if ((j - i) % 13) in {1, 2, 3, 5, 6, 9} else 0 for j in range(13)] for i in range(13)]
check("4a ST_13 = circulant {1,2,3,5,6,9} is a tournament", is_tournament(ST13))
mt = [max_transitive(C3), max_transitive(PAL[7]), max_transitive(ST13), max_transitive(P27)]
check("4b largest transitive subtournaments of C_3, P_7, ST_13, P_27 have 2, 3, 4, 5 vertices: R(k) - 1 >= 3, 7, 13, 27 for k = 3..6",
      mt == [2, 3, 4, 5], str(mt))
check("4c P_19 and P_11 contain transitive subtournaments of size %d and %d (19 is not in the sequence 3, 7, 13, 27)"
      % (max_transitive(PAL[19]), max_transitive(PAL[11])), max_transitive(PAL[19]) >= 5 and max_transitive(PAL[11]) >= 4)


@njit(cache=False)
def scan_fixed_path(n, k):
    """tournaments on n vertices containing 0 -> 1 -> ... -> n-1; returns the codes with no transitive k-subtournament
    (at most 4096 of them) and their number."""
    pi = np.zeros(64, np.int64)
    pj = np.zeros(64, np.int64)
    m = 0
    for i in range(n):
        for j in range(i + 2, n):
            pi[m] = i
            pj[m] = j
            m += 1
    total = 1 << m
    found = np.zeros(4096, np.int64)
    nf = 0
    out = np.zeros(n, np.int64)
    for code in range(total):
        for i in range(n):
            out[i] = 0
        for i in range(n - 1):
            out[i] |= 1 << (i + 1)
        for b in range(m):
            if (code >> b) & 1:
                out[pi[b]] |= 1 << pj[b]
            else:
                out[pj[b]] |= 1 << pi[b]
        # search for a transitive k-set: chains v1 -> ... -> vk with v_j dominated by all earlier
        has = False
        if k == 3:
            for a in range(n):
                ca = out[a]
                for b in range(n):
                    if (ca >> b) & 1 and (ca & out[b]) != 0:
                        has = True
                        break
                if has:
                    break
        else:                                   # k == 4
            for a in range(n):
                ca = out[a]
                for b in range(n):
                    if (ca >> b) & 1:
                        cab = ca & out[b]
                        if cab == 0:
                            continue
                        for c in range(n):
                            if (cab >> c) & 1 and (cab & out[c]) != 0:
                                has = True
                                break
                        if has:
                            break
                if has:
                    break
        if not has:
            if nf < 4096:
                found[nf] = code
            nf += 1
    return found, nf


_, n3_3 = scan_fixed_path(3, 3)
_, n3_4 = scan_fixed_path(4, 3)
check("4d R(3) = 4: a TT_3-free tournament on 3 vertices exists (C_3), none on 4 vertices (exhaustive)", n3_3 == 1 and n3_4 == 0)
found7, n4_7 = scan_fixed_path(7, 4)
_, n4_8 = scan_fixed_path(8, 4)
check("4e R(4) = 8: no TT_4-free tournament on 8 vertices (all 2^21 completions of a Hamiltonian path)", n4_8 == 0)


def decode(n, code):
    prs = [(i, j) for i in range(n) for j in range(i + 2, n)]
    adj = [[0] * n for _ in range(n)]
    for i in range(n - 1):
        adj[i][i + 1] = 1
    for b, (i, j) in enumerate(prs):
        if (code >> b) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj


def isomorphic_to(adj, ref):
    n = len(adj)
    for perm in itertools.permutations(range(n)):
        if all(adj[i][j] == ref[perm[i]][perm[j]] for i in range(n) for j in range(n)):
            return True
    return False


iso7 = all(isomorphic_to(decode(7, int(found7[i])), PAL[7]) for i in range(n4_7))
check("4f the TT_4-free tournaments on 7 vertices containing 0 -> ... -> 6 number 9 = H(P_7) * (7!/21) / 7!, and all are copies of P_7 (uniqueness)",
      n4_7 == 9 and iso7, "%d found" % n4_7)


# ---------------------------------------------------------------- 5. Hamiltonian path counts
@njit(cache=False)
def ham(adj, n):
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), np.int64)
    for i in range(n):
        dp[1 << i, i] = 1
    for mask in range(1, full + 1):
        for last in range(n):
            c = dp[mask, last]
            if c == 0:
                continue
            for nxt in range(n):
                if (mask >> nxt) & 1 == 0 and adj[last, nxt]:
                    dp[mask | (1 << nxt), nxt] += c
    tot = 0
    for last in range(n):
        tot += dp[full, last]
    return tot


def H(adj, verts=None):
    if verts is None:
        verts = list(range(len(adj)))
    a = np.array([[adj[i][j] for j in verts] for i in verts], dtype=np.uint8)
    return int(ham(a, len(verts)))


h3, h7, h11 = H(PAL[3]), H(PAL[7]), H(PAL[11])
check("5a H(P_3) = 3, H(P_7) = 189 = 7 * 27, H(P_11) = 95095 = 5 * 7 * 11 * 13 * 19, H(P_11)/55 = 1729 = 7 * 13 * 19",
      (h3, h7, h11) == (3, 189, 95095) and 189 == 7 * 27 and 95095 == 5 * 7 * 11 * 13 * 19 and 95095 // 55 == 1729 == 7 * 13 * 19)
sub6 = sorted({H(PAL[7], [v for v in range(7) if v != d]) for d in range(7)})
sub5 = sorted({H(PAL[7], [v for v in range(7) if v not in pr]) for pr in itertools.combinations(range(7), 2)})
sub4 = sorted({H(PAL[7], [v for v in range(7) if v not in tr]) for tr in itertools.combinations(range(7), 3)})
check("5b THM-258: all 21 five-vertex subtournaments of P_7 have H = 13; the six-vertex ones have H = 45; four-vertex ones 3 or 5",
      sub5 == [13] and sub6 == [45] and sub4 == [3, 5], "%s %s %s" % (sub6, sub5, sub4))
h13 = H(ST13)
z13 = [[1 if ((j - i) % 13) in {1, 3, 5, 7, 9, 11} else 0 for j in range(13)] for i in range(13)]
h13max = H(z13)
check("5c THM-128: H(ST_13) = 3703011 and H(Z_13{1,3,5,7,9,11}) = 3711175", h13 == 3703011 and h13max == 3711175, "%d %d" % (h13, h13max))
h19 = H(PAL[19])
check("5d H(P_19) = 1172695746915 = 171 * 5 * 7 * 11 * 23 * 774463 (171 = 9 * 19 = |Aut P_19|): the factor 13 of H(P_11) does not persist",
      h19 == 1172695746915 and h19 == 171 * 5 * 7 * 11 * 23 * 774463 and h19 % 13 != 0, str(factor(h19)))
# regular 7-tournaments: H in {171, 175, 189}
found_reg = set()
prs7 = [(i, j) for i in range(7) for j in range(i + 2, 7)]
for code in range(1 << 15):
    sc = [0] * 7
    for i in range(6):
        sc[i] += 1
    okc = True
    for b, (i, j) in enumerate(prs7):
        if (code >> b) & 1:
            sc[i] += 1
        else:
            sc[j] += 1
    if sc == [3] * 7:
        found_reg.add(H(decode(7, code)))
check("5e THM-053: regular 7-tournaments have H in {171, 175, 189} = {9 * 19, 25 * 7, 27 * 7}", sorted(found_reg) == [171, 175, 189])


# ---------------------------------------------------------------- 6. Legendre tournament on the primes = 3 mod 4 of the list
def legendre(a, p):
    v = pow(a % p, (p - 1) // 2, p)
    return 1 if v == 1 else -1


pr3 = [3, 7, 19, 1039]
arcs = {(p, qq): legendre(p, qq) == 1 for p in pr3 for qq in pr3 if p != qq}
check("6a reciprocity: for distinct primes = 3 mod 4 exactly one of p -> q, q -> p holds", all(arcs[(p, qq)] != arcs[(qq, p)] for p in pr3 for qq in pr3 if p < qq))
order_t = [7, 19, 1039, 3]
check("6b the Legendre tournament on {3, 7, 19, 1039} is transitive with order 7 -> 19 -> 1039 -> 3",
      all(arcs[(order_t[i], order_t[j])] for i in range(4) for j in range(i + 1, 4)))

if all(RESULTS):
    print("ALL CHECKS PASSED (%d checks)" % len(RESULTS))
    sys.exit(0)
print("CHECK FAILED (%d of %d)" % (RESULTS.count(False), len(RESULTS)))
sys.exit(1)
