#!/usr/bin/env python3
"""Tournament characterizations of Fermat and Wieferich primes (THM-4644); Paley Hamiltonian path counts and the
orbit-count hypothesis (HYP-9273) with controls; Paley graphs on Fermat primes and Ramsey; the (2,3,k) chain, perfect 6,
primary pseudoperfect 42; the tau(n) = 4 identity sigma(n)/n - 2 = triangle excess.

Session opus-2026-10-14-S33.  Note: 05-knowledge/results/number_tournament_bridges_20261014.md
Reproduce: python3 04-computation/experiments/number_tournament_bridges_20261014.py   (about 3 minutes, 4 GB)
"""
import math, itertools, time
from fractions import Fraction as Fr
import numpy as np
import numba as nb
import sympy as sp

FAIL = []


def check(label, ok):
    if not ok:
        FAIL.append(label)
    return ok


def odd_part(n):
    while n % 2 == 0:
        n //= 2
    return n


# ---------- 1. Fermat primes = primes all of whose circulant tournaments are rigid ----------
def circulant_tournaments(p):
    """all connection sets S with S + (-S) = Z_p^*: one choice per pair {x, -x}"""
    pairs = [(x, p - x) for x in range(1, (p - 1) // 2 + 1)]
    for bits in range(1 << len(pairs)):
        yield frozenset(pairs[i][(bits >> i) & 1] for i in range(len(pairs)))


def aut_order_bruteforce(p, S):
    """|Aut| of the circulant tournament (Z_p, S), by brute force over all permutations fixing 0 composed with translations
    (vertex-transitive: |Aut| = p * |Stab(0)|); feasible for p <= 11"""
    adj = [[(j - i) % p in S for j in range(p)] for i in range(p)]
    cnt = 0
    for perm in itertools.permutations(range(1, p)):
        f = (0,) + perm
        if all(adj[f[i]][f[j]] == adj[i][j] for i in range(p) for j in range(p) if i != j):
            cnt += 1
    return p * cnt


def multiplier_stabilizer(p, S):
    return [a for a in range(1, p) if frozenset(a * x % p for x in S) == S]


print('== 1. Fermat primes and rigid circulant tournaments')
print('   Claim: p is a Fermat prime  <=>  every circulant tournament on Z_p has Aut = Z_p.')
for p in (3, 5, 7, 11):
    sizes = sorted({aut_order_bruteforce(p, S) for S in circulant_tournaments(p)})
    print(f'   p = {p:2d}: Aut orders of all {2 ** ((p - 1) // 2)} circulant tournaments (brute force): {sizes}')
    check(f'brute {p}', (sizes == [p]) == (p in (3, 5)))
rows = []
for p in sp.primerange(3, 400):
    stabs = {len(multiplier_stabilizer(p, S)) for S in circulant_tournaments(p)} if p <= 31 else None
    rigid_by_theory = (odd_part(p - 1) == 1)
    if stabs is not None:
        check(f'stab {p}', (stabs == {1}) == rigid_by_theory)
    rows.append((p, rigid_by_theory))
print('   primes p < 400 with p - 1 a power of 2 (all circulant tournaments rigid):', [p for p, r in rows if r])
print('   multiplier stabilizers checked exhaustively for p <= 31: rigid exactly at p = 3, 5, 17:', not [x for x in FAIL if x.startswith('stab')])
print('   number of rotational tournaments up to isomorphism at a Fermat prime p (free Z_(p-1) action): 2^((p-1)/2)/(p-1):',
      {p: (f'2^{(p - 1) // 2 - int(math.log2(p - 1))}' if p > 17 else 2 ** ((p - 1) // 2) // (p - 1)) for p in (3, 5, 17, 257, 65537)})

# ---------- 2. Wieferich primes and doubling-invariant circulant tournaments on Z_(p^2) ----------
print('\n== 2. Wieferich primes: circulant tournaments on Z_(p^2) invariant under g = the odd-order part of <2>')
print('   (g = 2^(ord_p(2)/d), d = odd part of ord_p(2); invariant tournaments exist because <g> has odd order, so -1 is not in <g>)')


def invariant_tournament_count_exponent(p):
    d_full = sp.n_order(2, p)
    d = odd_part(d_full)
    q = p * p
    g = pow(2, d_full // d, q)
    og = sp.n_order(g, q)
    unit_orbits = (p * (p - 1)) // og
    nonunit_orbits = (p - 1) // d          # p*(Z/p)^*, on which g acts with order d
    return d_full, d, og, unit_orbits, nonunit_orbits, (unit_orbits + nonunit_orbits) // 2


wief = []
for p in list(sp.primerange(3, 4000)):
    d_full, d, og, uo, nuo, e = invariant_tournament_count_exponent(p)
    is_w = pow(2, p - 1, p * p) == 1
    generic = (p - 1) // d                 # exponent when ord_(p^2)(g) = p d
    if is_w or p in (3, 7, 31, 127):
        wief.append((p, d_full, d, og, e, generic, is_w))
    check(f'wief {p}', (e != generic) == is_w)
for p, d_full, d, og, e, generic, is_w in wief:
    print(f'   p = {p:5d}: ord_p(2) = {d_full:5d}, odd part d = {d:5d}, ord_(p^2)(g) = {og:8d};  invariant tournaments on Z_(p^2): 2^{e}'
          f'  (2^{generic} unless Wieferich);  Wieferich: {is_w}')
print('   below 4000 the exponent jumps exactly at the Wieferich primes 1093 and 3511:', check('wieferich jump', not [x for x in FAIL if x.startswith('wief')]))
p = 3511
print(f'   3511: ord_3511(2) = {sp.n_order(2, 3511)} = (p-1)/2, p = {3511 % 8} mod 8, p = {3511 % 4} mod 4: <2> = the squares, so the doubling-invariant',
      'circulant tournaments on Z_3511 are the Paley tournament and its reverse (as <2> = {1,2,4} gives P_7 at p = 7)')
check('3511 paley', sp.n_order(2, 3511) == 1755 and 3511 % 4 == 3)
print(f'   1093 = (3^7 - 1)/2: {(3 ** 7 - 1) // 2 == 1093};  ord_1093(3) = {sp.n_order(3, 1093)};  ord_1093(2) = {sp.n_order(2, 1093)} = 4 * 91, 91 = 7 * 13;'
      f'  1093 mod 21 = {1093 % 21}, 3511 mod 21 = {3511 % 21} (both in <2> mod 21 = {sorted({pow(2, k, 21) for k in range(6)})}, i.e. both split in Q(sqrt(-7)))')

# ---------- 3. Hamiltonian paths of Paley tournaments ----------
print('\n== 3. Hamiltonian paths of the Paley tournaments P_p (p = 3 mod 4)')


@nb.njit(cache=True)
def ham_paths_from0(adjmask, n):
    """number of Hamiltonian paths starting at vertex 0 (bitmask DP)"""
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), dtype=np.int64)
    dp[1, 0] = 1
    for mask in range(1, 1 << n):
        if (mask & 1) == 0:
            continue
        for v in range(n):
            c = dp[mask, v]
            if c == 0:
                continue
            nxt = adjmask[v] & ~mask
            while nxt:
                w = 0
                b = nxt & (-nxt)
                while (b >> w) != 1:
                    w += 1
                dp[mask | b, w] += c
                nxt ^= b
    tot = 0
    for v in range(n):
        tot += dp[full, v]
    return tot


hp = {}
for p in (3, 7, 11, 19, 23):
    qr = {x * x % p for x in range(1, p)}
    adj = np.array([sum(1 << j for j in range(p) if (j - i) % p in qr) for i in range(p)], dtype=np.int64)
    t0 = time.time()
    h0 = ham_paths_from0(adj, p)
    H = p * h0
    aut = p * (p - 1) // 2
    hp[p] = H
    print(f'   p = {p:2d}: H(P_p) = {H}  (= {aut} * {H // aut}; |Aut| = p(p-1)/2 divides H: {H % aut == 0});  H/p mod p = {(H // p) % p}, '
          f'H/|Aut| mod p = {(H // aut) % p};  {time.time() - t0:.1f} s')
    check(f'aut divides H {p}', H % aut == 0)

print('   orbit counts H/|Aut| mod 8 for the Paley primes:', {p: (hp[p] // (p * (p - 1) // 2)) % 8 for p in hp})
check('paley mod 8', all((hp[p] // (p * (p - 1) // 2)) % 8 == 1 for p in hp))

print('\n   control: every circulant tournament on Z_7 and Z_11 (|Aut| = p |M|, M the multiplier stabilizer; Alspach)')
for p in (7, 11):
    seen_cls = {}
    for S in circulant_tournaments(p):
        M = len(multiplier_stabilizer(p, S))
        adj = np.array([sum(1 << j for j in range(p) if (j - i) % p in S) for i in range(p)], dtype=np.int64)
        H = p * ham_paths_from0(adj, p)
        key = (M, H)
        seen_cls[key] = seen_cls.get(key, 0) + 1
    for (M, H), cnt in sorted(seen_cls.items()):
        aut = p * M
        print(f'     p = {p}: {cnt:2d} connection sets with |Aut| = {aut:3d}, H = {H:7d}, H/|Aut| = {H // aut} = {(H // aut) % 8} mod 8')


@nb.njit(cache=True)
def ham_paths_from0_mod(adjmask, n, M):
    """Hamiltonian paths starting at vertex 0, modulo M (uint16 table, masks containing vertex 0 only)"""
    half = 1 << (n - 1)
    dp = np.zeros((half, n), dtype=np.uint16)
    dp[0, 0] = 1                      # index = mask >> 1 for masks with bit 0 set
    for idx in range(half):
        mask = (idx << 1) | 1
        for v in range(n):
            c = np.int64(dp[idx, v])
            if c == 0:
                continue
            nxt = adjmask[v] & ~mask
            while nxt:
                b = nxt & (-nxt)
                w = 0
                while (b >> w) != 1:
                    w += 1
                j = (mask | b) >> 1
                t = np.int64(dp[j, w]) + c
                if t >= M:
                    t -= M
                dp[j, w] = t
                nxt ^= b
    tot = 0
    for v in range(n):
        tot += np.int64(dp[half - 1, v])
    return tot % M


print('\n   all circulant tournaments on Z_p: residues of the orbit count H/|Aut| (|Aut| = p |M|) mod 8, by p mod 4')
from collections import Counter as _C
mod8_by_p = {}
for p in (3, 5, 7, 11, 13, 17, 19):
    res = _C()
    for S in circulant_tournaments(p):
        M = len(multiplier_stabilizer(p, S))
        adj = np.array([sum(1 << j for j in range(p) if (j - i) % p in S) for i in range(p)], dtype=np.int64)
        H = p * ham_paths_from0(adj, p)
        aut = p * M
        assert H % aut == 0
        res[(H // aut) % 8] += 1
    mod8_by_p[p] = dict(res)
    print(f'     p = {p:2d} (p = {p % 4} mod 4): {2 ** ((p - 1) // 2)} connection sets, residues of H/|Aut| mod 8: {dict(sorted(res.items()))}')
print('     -> the guess "every circulant on Z_p with p = 3 mod 4 has orbit count 1 mod 8" holds at p = 3, 7, 11 and is REFUTED at p = 19')
check('mod8 refuted at 19', set(mod8_by_p[19]) != {1} and all(set(mod8_by_p[p]) == {1} for p in (3, 7, 11)))

print('\n   control: the doubly regular tournaments of the skew-Hadamard doubling tower (gilbreath_fermat_tower_20260926_tournaments.py;')
print('   |Aut(T_4)| = 21 there, not vertex-transitive).  Aut acts freely on Hamiltonian paths of any tournament.')
for kk, nn_expect in ((3, 7), (4, 15)):
    Hm = np.array([[1, 1], [-1, 1]], dtype=np.int64)
    for _ in range(kk - 1):
        Hm = np.vstack([np.hstack([Hm, Hm]), np.hstack([-Hm.T, Hm.T])])
    A = (Hm[1:, 1:] == 1)
    np.fill_diagonal(A, False)
    nn = A.shape[0]
    adj = np.array([sum(1 << j for j in range(nn) if A[i, j]) for i in range(nn)], dtype=np.int64)
    # total Hamiltonian paths: sum over start vertices (relabel so the start is vertex 0)
    tot = 0
    for s in range(nn):
        perm = [s] + [v for v in range(nn) if v != s]
        pos = {v: i for i, v in enumerate(perm)}
        adj_s = np.array([sum(1 << pos[j] for j in range(nn) if A[perm[i], j]) for i in range(nn)], dtype=np.int64)
        tot += ham_paths_from0(adj_s, nn)
    aut = 21
    print(f'     T_{kk} (n = {nn}): H = {tot}, |Aut| = {aut}, H/|Aut| = {tot // aut} = {(tot // aut) % 8} mod 8')
    if nn == 15:
        check('DR15 not 1 mod 8', tot % aut == 0 and (tot // aut) % 8 == 5)

print('\n   the Paley tournament on GF(27) = F_3[x]/(x^3 - x - 1): H mod 8424 = 8 * 1053 (|Aut| = 27 * 13 * 3 = 1053)')
def f27_mul(a, b):
    A = [(a // 3 ** i) % 3 for i in range(3)]
    B = [(b // 3 ** i) % 3 for i in range(3)]
    prod = [0] * 5
    for i in range(3):
        for j in range(3):
            prod[i + j] += A[i] * B[j]
    for k in (4, 3):                  # x^3 = x + 1
        c = prod[k]
        prod[k] = 0
        prod[k - 3] += c
        prod[k - 2] += c
    return sum((prod[i] % 3) * 3 ** i for i in range(3))


def f27_sub(a, b):
    return sum((((a // 3 ** i) - (b // 3 ** i)) % 3) * 3 ** i for i in range(3))


sq27 = {f27_mul(x, x) for x in range(1, 27)}
check('27 squares', len(sq27) == 13 and f27_sub(0, 1) not in sq27)
adj27 = np.array([sum(1 << j for j in range(27) if j != i and f27_sub(j, i) in sq27) for i in range(27)], dtype=np.int64)
t0 = time.time()
h0 = ham_paths_from0_mod(adj27, 27, 8424)
H27 = (27 * h0) % 8424
print(f'     H(P_27) mod 8424 = {H27};  divisible by 1053: {H27 % 1053 == 0};  orbit count under the full Aut (with Frobenius): H/1053 = {H27 // 1053} mod 8;'
      f'  under the affine group x -> a x + b (a a square, order 351): H/351 = {(3 * (H27 // 1053)) % 8} mod 8;  {time.time() - t0:.0f} s')
check('P27 mod 8', H27 % 1053 == 0 and (3 * (H27 // 1053)) % 8 == 1)

# ---------- 4. Paley graphs at Fermat primes and Ramsey numbers ----------
print('\n== 4. Paley graphs on Fermat primes: clique numbers')


def clique_number_paley(q):
    qr = {x * x % q for x in range(1, q)}
    nbrs = sorted(qr)                 # neighbourhood of 0; the graph is vertex-transitive, so a maximum clique may contain 0
    best = [1]

    def extend(clique, cand):
        if len(clique) > best[0]:
            best[0] = len(clique)
        if len(clique) + len(cand) <= best[0]:
            return
        for i, v in enumerate(cand):
            extend(clique + [v], [w for w in cand[i + 1:] if (w - v) % q in qr])
    extend([0], nbrs)
    return best[0]


for q in (5, 13, 17, 29, 37, 41, 257):
    t0 = time.time()
    w = clique_number_paley(q)
    print(f'   P_{q}: clique number {w} (= independence number, self-complementary), so R({w + 1},{w + 1}) > {q};  {time.time() - t0:.1f} s')
    if q == 5:
        check('P5', w == 2)
    if q == 17:
        check('P17', w == 3)
print('   R(3,3) = 6 and R(4,4) = 18: the Paley graphs on the Fermat primes 5 and 17 are the extremal graphs (classical; uniqueness cited, not re-checked)')

# ---------- 5. the (2,3,k) chain, perfect and primary pseudoperfect numbers ----------
print('\n== 5. Triangle groups (2,3,k), perfect 6 and primary pseudoperfect 42')
for k in range(3, 9):
    chi = Fr(1, 2) + Fr(1, 3) + Fr(1, k) - 1
    kind = 'spherical' if chi > 0 else ('Euclidean' if chi == 0 else 'hyperbolic')
    extra = f', |group| = 2/chi = {2 / chi}' if chi > 0 else (f', Hurwitz: |Aut| <= (2g-2)/|chi| = {1 / (-chi) * 2}(g-1)' if chi < 0 else '')
    print(f'   (2,3,{k}): orbifold Euler characteristic {chi} -> {kind}{extra}')
pseudo = [n for n in range(2, 60000) if sp.factorint(n) and all(e == 1 for e in sp.factorint(n).values())
          and sum(Fr(1, pp) for pp in sp.factorint(n)) + Fr(1, n) == 1]
print('   primary pseudoperfect numbers n < 60000 (sum_{p|n} 1/p + 1/n = 1):', pseudo)
check('pseudo', pseudo[:4] == [2, 6, 42, 1806])
print('   chi(2,3,7) = -1/42 with 42 = 2*3*7 primary pseudoperfect;  chi(2,3,6) = 0 with 6 perfect;  Sylvester 2, 3, 7, 43, 1807 = 13 * 139')
print('   sigma(n)/n = sum_{d | n} 1/d for every n, so sigma(n)/n - 2 = sum_{d | n, d > 1} 1/d - 1 is a tautology;  the triangle reading')
print('   needs exactly three divisors above 1, i.e. tau(n) = 4 (n = pq or p^3):')
for n in (6, 15, 21, 35, 77, 4, 8, 27):
    divs = [d for d in sp.divisors(n) if d > 1]
    exc = sum(Fr(1, d) for d in divs) - 1
    ab = Fr(sp.divisor_sigma(n), n) - 2
    if sp.divisor_count(n) == 4:
        print(f'   n = {n:2d} (tau = 4): divisors > 1 {divs}, triangle excess/pi = {exc}, sigma(n)/n - 2 = {ab}, equal: {check("tau4 " + str(n), exc == ab)}')
    else:
        print(f'   n = {n:2d} (tau = {sp.divisor_count(n)}): divisors > 1 {divs}, not a triangle (excluded)')
abund4 = [n for n in range(2, 100000) if sp.divisor_count(n) == 4 and sp.divisor_sigma(n) > 2 * n]
perf4 = [n for n in range(2, 100000) if sp.divisor_count(n) == 4 and sp.divisor_sigma(n) == 2 * n]
print(f'   tau(n) = 4, n < 10^5: abundant {abund4}, perfect {perf4} (only 6 = 2*3 <-> (2,3,6); no spherical case)')
check('tau4 abundant none', abund4 == [] and perf4 == [6])

print('\nFAILED CHECKS:', FAIL if FAIL else 'none')
print('ALL CHECKS PASSED' if not FAIL else 'SOME CHECKS FAILED')
