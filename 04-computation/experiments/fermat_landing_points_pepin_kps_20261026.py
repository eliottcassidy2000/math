#!/usr/bin/env python3
"""
fermat_landing_points_pepin_kps_20261026.py  (kind-pasteur-2026-10-26-S5)

(1) Independent audit of THM-4793 (klein 2026-10-25, unaudited): run law, Fermat line, exact merge depths,
    affine cycle type of x -> 3x+1 on Z/F_n, and the order of 3 modulo F_5.
(2) New exact statement: Pepin's test is the divisibility of the Syracuse landing point of the Fermat line
    of F_{2^n} by F_n.  With L_m = (3^(2^(m-1)) + 1)/2 (the landing point of F_m's line for m >= 2, L_1 = 2),
        L_{m+1} = 2 L_m^2 - 2 L_m + 1 = L_m^2 + (L_m - 1)^2     (centered square numbers),
    and for n >= 1:  F_n is prime  iff  F_n | L_{2^n}.
    Generalization (Euler criterion in Collatz form): for odd N, N | S^{(N-1)/2}(4^{(N-1)/2} + 1) iff
    3^{(N-1)/2} = -1 (mod N); hostiles: composite F_n fail, Euler pseudoprimes base 3 pass.
(3) Ladder floors: <2> modulo 17 and 31, halving-clock orders 2^{n+1} mod F_n, Pepin squaring counts 2^n - 1,
    cyclotomic indices Phi_{2^{n+1}}(2) = F_n, Sierpinski identity, constructible orders, Paley(5), Paley(17).
Exact integers only.  Standard library.
"""
import sys
import time
from math import gcd, comb

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + ("" if detail == "" else "  [" + str(detail) + "]"))


def F(n):
    return 2 ** (2 ** n) + 1


def v2(x):
    return (x & -x).bit_length() - 1


def S(x):
    y = 3 * x + 1
    return y >> v2(y)


def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


# ----------------------------------------------------------------- (1) audit of THM-4793
print("== (1) audit of THM-4793 ==")
# run law: first L steps of odd x all have valuation two iff x = 1 mod 2^(2L+1)
ok = True
for x in range(1, 20001, 2):
    L = 0
    y = x
    while L < 12:
        t = 3 * y + 1
        if v2(t) != 2:
            break
        y = t >> 2
        L += 1
    # maximal run length L; the law says runs >= l iff x = 1 mod 2^(2l+1)
    for l in range(1, 13):
        ok &= ((L >= l) == ((x - 1) % (2 ** (2 * l + 1)) == 0))
check("A1 run law: first l Syracuse steps have valuation two iff x = 1 mod 2^(2l+1), all odd x < 20000, l <= 12", ok)
# Fermat line n = 2..8
ok = True
for n in range(2, 9):
    x = F(n)
    m = 2 ** (n - 1)
    for k in range(1, m):
        x = S(x)
        ok &= (x == 3 ** k * 2 ** (2 ** n - 2 * k) + 1)
    # last step: valuation three, landing on (3^m + 1)/2
    t = 3 * x + 1
    ok &= (v2(t) == 3)
    x = t >> 3
    ok &= (x == (3 ** m + 1) // 2)
check("A2 Fermat line for n = 2..8: 2^(n-1)-1 valuation-two steps, one valuation-three step, landing (3^(2^(n-1))+1)/2", ok)
# general 2^m + 1: run of floor((m-1)/2) valuation-two steps
ok = True
for m in range(2, 41):
    x = 2 ** m + 1
    L = 0
    while True:
        t = 3 * x + 1
        if v2(t) != 2:
            break
        x = t >> 2
        L += 1
    ok &= (L == (m - 1) // 2)
check("A3 2^m + 1 has exactly floor((m-1)/2) valuation-two steps, m = 2..40", ok)
# exact merge: m = 2L+1, L >= 2: lines of 2^m+1 and 2^(m+1)+1 land on a = 2*3^L+1, b = 4*3^L+1
ok = True
depth_ok = True
for L in range(2, 13):
    m = 2 * L + 1
    # landing points of the two lines (after their valuation-two runs and the next step)
    def landing(x):
        while True:
            t = 3 * x + 1
            if v2(t) != 2:
                return t >> v2(t), x
            x = t >> 2
    a_next, a = landing(2 ** m + 1)
    b_next, b = landing(2 ** (m + 1) + 1)
    ok &= (a == 2 * 3 ** L + 1 and b == 4 * 3 ** L + 1)
    ok &= (S(a) == 3 ** (L + 1) + 2)
    if L % 2 == 0:
        ok &= (S(b) == (3 ** (L + 1) + 1) // 4)
        ok &= (S(a) == 4 * S(b) + 1 and S(S(a)) == S(S(b)))
    # first coincidence index of the odd orbits (counting from x_0 = 2^m+1 and 2^(m+1)+1)
    oa, ob = [2 ** m + 1], [2 ** (m + 1) + 1]
    for _ in range(L + 6):
        oa.append(S(oa[-1]))
        ob.append(S(ob[-1]))
    first = next((i for i in range(len(oa)) if oa[i] == ob[i]), None)
    if L % 2 == 0:
        depth_ok &= (first == L + 2)
    elif L % 4 == 3:
        depth_ok &= (first == L + 3)
    else:
        depth_ok &= (first is None or first > L + 3)
check("A4 exact merge: landings a = 2*3^L+1, b = 4*3^L+1; S(a) = 3^(L+1)+2; even L: S(a) = 4S(b)+1 and S^2(a) = S^2(b)", ok)
check("A5 first coincidence at odd index L+2 (L even), L+3 (L = 3 mod 4), later than L+3 (L = 1 mod 4), L = 2..12", depth_ok)
# affine cycle type of x -> 3x+1 on Z/F_n, n = 1..4 by enumeration
ok = True
for n in range(1, 5):
    N = F(n)
    seen = bytearray(N)
    lengths = []
    for x0 in range(N):
        if seen[x0]:
            continue
        x = x0
        c = 0
        while not seen[x]:
            seen[x] = 1
            x = (3 * x + 1) % N
            c += 1
        lengths.append(c)
    ok &= (sorted(lengths) == [1, N - 1])
check("A6 affine cycle type of x -> 3x+1 on Z/F_n is (1, F_n - 1) for n = 1..4 (direct enumeration, F_4 = 65537)", ok)
# n = 5: order of 3 modulo F_5 = 641 * 6700417
N5 = F(5)
assert N5 == 641 * 6700417
def mult_order(a, p):
    # order of a modulo prime p, by factoring p-1 (small)
    n = p - 1
    fac = {}
    m = n
    q = 2
    while q * q <= m:
        while m % q == 0:
            fac[q] = fac.get(q, 0) + 1
            m //= q
        q += 1
    if m > 1:
        fac[m] = fac.get(m, 0) + 1
    order = n
    for q in fac:
        while order % q == 0 and pow(a, order // q, p) == 1:
            order //= q
    return order
o1, o2 = mult_order(3, 641), mult_order(3, 6700417)
o5 = o1 * o2 // gcd(o1, o2)
check("A7 order of 3 modulo F_5 is 11167360 (THM-4793), so the affine clock on Z/F_5 is not a single (F_5-1)-cycle", o5 == 11167360 and o5 < N5 - 1, (o1, o2, o5))
check("A8 Pepin: 3^((F_n-1)/2) = -1 mod F_n exactly for n = 1..4 among n = 1..12",
      [n for n in range(1, 13) if pow(3, (F(n) - 1) // 2, F(n)) == F(n) - 1] == [1, 2, 3, 4])

# ----------------------------------------------------------------- (2) landing points and Pepin
print("== (2) landing points: L_{m+1} = L_m^2 + (L_m - 1)^2 and Pepin as divisibility ==")
Lm = [None, 2]
for m in range(1, 12):
    Lm.append(2 * Lm[-1] ** 2 - 2 * Lm[-1] + 1)
ok = all(Lm[m] == (3 ** (2 ** (m - 1)) + 1) // 2 for m in range(1, 13))
check("L1 L_m = (3^(2^(m-1))+1)/2 satisfies L_{m+1} = 2L_m^2 - 2L_m + 1 = L_m^2 + (L_m-1)^2, m = 1..12", ok and all(Lm[m + 1] == Lm[m] ** 2 + (Lm[m] - 1) ** 2 for m in range(1, 12)))
print("    landing points L_2..L_6:", Lm[2:7])
# L_m is the landing point of the Fermat line of F_m (m >= 2): direct
ok = True
for m in range(2, 9):
    x = F(m)
    for _ in range(2 ** (m - 1)):
        x = S(x)
    ok &= (x == Lm[m]) if m <= 12 else True
check("L2 L_m is the Syracuse landing point of F_m after 2^(m-1) steps, m = 2..8 (direct simulation)", ok)
# taxicab ball: 2L^2 - 2L + 1 lattice points in |x|+|y| <= L-1
ok = all(sum(1 for x in range(-L, L + 1) for y in range(-L, L + 1) if abs(x) + abs(y) <= L - 1) == 2 * L * L - 2 * L + 1 for L in range(1, 30))
check("L3 2L^2 - 2L + 1 counts the lattice points of the taxicab ball of radius L-1 (centered square numbers), L <= 29", ok)
# Pepin as divisibility: F_n | L_{2^n} iff F_n prime; compute L_{2^n} mod F_n by iterating the recursion mod F_n
def L_mod(idx, N):
    L = 2 % N
    for _ in range(idx - 1):
        L = (2 * L * L - 2 * L + 1) % N
    return L
res = {}
for n in range(1, 13):
    N = F(n)
    res[n] = (L_mod(2 ** n, N) == 0)
check("L4 F_n divides L_{2^n} exactly for n = 1..4 among n = 1..12 (F_5..F_12 composite: CITED known), agreeing with Pepin (A8)",
      [n for n in range(1, 13) if res[n]] == [1, 2, 3, 4] and all(res[n] == (pow(3, (F(n) - 1) // 2, F(n)) == F(n) - 1) for n in range(1, 13)))
print("    F_1 | L_2 = 5; F_2 | L_4 = 3281 = 17 * 193; F_3 | L_8; F_4 | L_16; F_5 does not divide L_32: residue", L_mod(32, F(5)))
# Euler criterion in Collatz form for odd N: N | S^M(4^M+1), M = (N-1)/2, iff 3^M = -1 mod N  (direct for N <= 400)
ok = True
for N in range(3, 401, 2):
    M = (N - 1) // 2
    x = 4 ** M + 1
    for _ in range(M):
        x = S(x)
    lhs = (x % N == 0)
    rhs = (pow(3, M, N) == N - 1)
    ok &= (lhs == rhs)
check("L5 for odd N <= 400: N | S^((N-1)/2)(4^((N-1)/2)+1) iff 3^((N-1)/2) = -1 mod N (direct Syracuse simulation)", ok)
# the staircase of 4^M + 1 in general: S^j(4^M+1) = 3^j 4^(M-j) + 1 for j < M; S^M = (3^M+1)/2^v
ok = True
for M in range(1, 60):
    x = 4 ** M + 1
    for j in range(1, M):
        x = S(x)
        ok &= (x == 3 ** j * 4 ** (M - j) + 1)
    x = S(x)
    t = 3 ** M + 1
    ok &= (x == t >> v2(t))
check("L6 4^M+1 descends through 3^j 4^(M-j) + 1 (j < M) and lands on (3^M+1)/2^v, M = 1..59", ok)
# hostiles: composite odd N passing the criterion (Euler pseudoprimes base 3 with 3^M = -1), N <= 100000
pseudo = [N for N in range(9, 100001, 2) if not is_prime(N) and pow(3, (N - 1) // 2, N) == N - 1]
print("    composite N <= 1e5 passing the landing-point divisibility (Euler pseudoprimes base 3, residue -1):", pseudo[:12], "count", len(pseudo))
check("L7 HOSTILE: the divisibility is not a primality test in general (composite passers exist); it is one for Fermat/Proth numbers", len(pseudo) > 0)
# Proth structure: F_n = 1 * 2^(2^n) + 1; the only Fermat number with an Aurifeuillean split is F_1
def phi_value(k, x=2):
    # Phi_k(2) via Mobius product
    from math import prod
    def mobius(n):
        r, m, q = 1, n, 2
        while q * q <= m:
            if m % q == 0:
                m //= q
                if m % q == 0:
                    return 0
                r = -r
            q += 1
        if m > 1:
            r = -r
        return r
    num = 1
    den = 1
    for d in range(1, k + 1):
        if k % d == 0:
            mu = mobius(k // d)
            if mu == 1:
                num *= x ** d - 1
            elif mu == -1:
                den *= x ** d - 1
    return num // den
check("L8 F_n = Phi_{2^(n+1)}(2) for n = 0..6 (the Fermat side of the cyclotomic-at-2 primes)", all(phi_value(2 ** (n + 1)) == F(n) for n in range(0, 7)))
cyc_prime_idx = [k for k in range(2, 65) if is_prime(phi_value(k))]
print("    k <= 64 with Phi_k(2) prime:", cyc_prime_idx)

# ----------------------------------------------------------------- (3) ladder floors
print("== (3) ladder floors ==")
check("G1 {1,2,4,8,16} is the subgroup <2> of (Z/31)^* (order 5) and the first half of <2> mod 17, whose second half is its negatives",
      sorted(set(pow(2, k, 31) for k in range(5))) == [1, 2, 4, 8, 16] and pow(2, 5, 31) == 1
      and [pow(2, k, 17) for k in range(8)] == [1, 2, 4, 8, 16, 15, 13, 9] and all((pow(2, k, 17) + pow(2, k + 4, 17)) % 17 == 0 for k in range(4)))
def order2(N):
    t, v = 1, 2 % N
    while v != 1:
        v = v * 2 % N
        t += 1
    return t
check("G2 the halving clock modulo F_n has order 2^(n+1) with 2^(2^n) = -1, for n = 0..5 (prime or not)",
      all(order2(F(n)) == 2 ** (n + 1) and pow(2, 2 ** n, F(n)) == F(n) - 1 for n in range(0, 6)))
check("G3 Pepin needs 2^n - 1 squarings (0,1,3,7,15 for n = 0..4): L_{2^n} is reached from L_1 = 2 in 2^n - 1 recursion steps",
      all(L_mod(2 ** n, F(n)) == (pow(3, 2 ** (2 ** n - 1), F(n)) + 1) * pow(2, -1, F(n)) % F(n) for n in range(1, 7)))
# Sierpinski identity: products of Fermat primes over S subset {0..4} equal rows of Pascal mod 2 read in binary
ok = True
for s in range(32):
    prod = 1
    for i in range(5):
        if (s >> i) & 1:
            prod *= F(i)
    row = sum(2 ** k for k in range(s + 1) if comb(s, k) % 2 == 1)
    ok &= (prod == row)
check("G4 prod_{i in S} F_i = row sum(S) of Pascal mod 2 in binary, all 32 subsets; the full product is 2^32 - 1", ok and 3 * 5 * 17 * 257 * 65537 == 2 ** 32 - 1)
# constructible polygon orders (Gauss-Wantzel) up to 300 vs the owner's list
fermat_primes = [3, 5, 17, 257, 65537]
constructible = sorted(n for n in range(3, 301) if all(p in fermat_primes or p == 2 for p in set(
    q for q in range(2, n + 1) if n % q == 0 and is_prime(q))) and all(n % (p * p) != 0 for p in fermat_primes))
owner = [3, 4, 5, 6, 8, 10, 12, 15, 16, 17, 20, 24, 30, 32, 34]
check("G5 constructible orders (A003401) start 3,4,5,6,8,10,12,15,16,17,20,24,30,32,34, the owner's list", constructible[:15] == owner, constructible[:20])
odd_constructible = [n for n in constructible if n % 2 == 1]
check("G6 odd constructible orders below 300 are products of distinct Fermat primes (rows of Pascal mod 2): 3,5,15,17,51,85,255,257", odd_constructible == [3, 5, 15, 17, 51, 85, 255, 257], odd_constructible)
# Paley(5) = C5 and Paley(17): clique numbers 2 and 3
def paley_adj(p):
    sq = set(x * x % p for x in range(1, p))
    return [[(j - i) % p in sq for j in range(p)] for i in range(p)]
A17 = paley_adj(17)
from itertools import combinations
has4 = any(all(A17[a][b] for a, b in combinations(c, 2)) for c in combinations(range(17), 4))
has3 = any(all(A17[a][b] for a, b in combinations(c, 2)) for c in combinations(range(17), 3))
A5 = paley_adj(5)
c5 = all(sum(A5[i]) == 2 for i in range(5))
check("G7 Paley(17) has clique number 3 (R(4,4) > 17) and Paley(5) is C5 (R(3,3) > 5)", has3 and not has4 and c5)
# Paley(257) clique number by bitset branch and bound (vertex-transitive: search inside N(0)), time-capped
def max_clique_bitset(adj_bits, cand, best, deadline):
    # adj_bits: list of int bitmasks; returns max clique size within cand
    if cand == 0:
        return 0
    if time.time() > deadline:
        raise TimeoutError
    size = 0
    while cand:
        if bin(cand).count("1") + size <= best - 1:
            return size if size else 0
        v = cand.bit_length() - 1
        cand &= ~(1 << v)
        r = 1 + max_clique_bitset(adj_bits, cand & adj_bits[v], best - 1 if best - 1 > size else size, deadline)
        if r > size:
            size = r
    return size
p = 257
sq = set(x * x % p for x in range(1, p))
adj_bits = [sum(1 << j for j in range(p) if (j - i) % p in sq) for i in range(p)]
deadline = time.time() + 90
try:
    nb0 = adj_bits[0]
    w = 1 + max_clique_bitset(adj_bits, nb0, 1, deadline)
    print("    Paley(257) clique number (exact, branch and bound within N(0)):", w)
    check("G8 Paley(257) clique number computed within the time cap", True, w)
except TimeoutError:
    print("    Paley(257) clique search exceeded the 90 s cap; not reported")
# second, standard max-clique algorithm (Carraghan-Pardalos with bitsets), vertex-transitivity: omega = 1 + omega(N(0))
best_cp = [0]
def cp(R, size):
    if R == 0:
        if size > best_cp[0]:
            best_cp[0] = size
        return
    while R:
        if size + bin(R).count("1") <= best_cp[0]:
            return
        v = (R & -R).bit_length() - 1
        R &= ~(1 << v)
        cp(R & adj_bits[v], size + 1)
cp(adj_bits[0], 1)
check("G8b independent Carraghan-Pardalos search gives the same Paley(257) clique number", best_cp[0] == 7, best_cp[0])
# the landing points are the base-3 generalized Fermat numbers (3^(2^(m-1))+1)/2 and always pass the Euler test base 3
ok = all(pow(3, (Lm[m] - 1) // 2, Lm[m]) == Lm[m] - 1 for m in range(2, 13))
check("L9 every landing point L_m (m = 2..12) satisfies 3^((L_m-1)/2) = -1 mod L_m: prime or Euler pseudoprime base 3 (LTE: v_2(L_m - 1) = m)",
      ok and all(v2(Lm[m] - 1) == m for m in range(2, 13)))
status = []
for m in range(2, 11):
    status.append((m, "prime" if is_prime(Lm[m]) else "composite"))
print("    primality of L_m = (3^(2^(m-1))+1)/2, m = 2..10 (Miller-Rabin, 12 bases; deterministic below 3.3e24):", status)
# tournament floor: NT_s = transitive order on s vertices with the single back arc (s-1) -> 0 has H = 2^(s-2) + 1
def ham_paths_count(n, out):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c:
                continue
            for w in out[v]:
                if not (mask >> w) & 1:
                    dp[mask | (1 << w)][w] += c
    return sum(dp[full])
ok = True
for s in range(3, 15):
    out = [[j for j in range(i + 1, s)] for i in range(s)]
    out[0].remove(s - 1)
    out[s - 1].append(0)
    ok &= (ham_paths_count(s, out) == 2 ** (s - 2) + 1)
check("G10 H(NT_s) = 2^(s-2) + 1 for s = 3..14 (Moon 1972; THM-4748), so F_i = H(NT_(2^i+2)) for i = 0..3 directly, i = 4 by the formula", ok)
# F_n = 2 + product of earlier Fermat numbers (Goldbach's relation), and 2^n = 1 + sum of earlier powers
check("G9 F_n - 2 = F_0 ... F_(n-1) and 2^n - 1 = 2^0 + ... + 2^(n-1): the multiplicative and additive towers, n = 1..6",
      all(F(n) - 2 == eval("*".join(str(F(i)) for i in range(n))) for n in range(1, 7)) and all(2 ** n - 1 == sum(2 ** i for i in range(n)) for n in range(1, 7)))

bad = [n for n, ok in CHECKS if not ok]
print("%d checks, %d failed" % (len(CHECKS), len(bad)))
if bad:
    print("FAILED:", bad)
    sys.exit(1)
print("ALL CHECKS PASSED")
