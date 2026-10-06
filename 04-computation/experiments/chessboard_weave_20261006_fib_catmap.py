#!/usr/bin/env python3
"""The Fibonacci cat map Q = [[0,1],[1,1]] on the torus chessboard Z_N^2
(chessboard-weave, 2026-10-06).  Q acts on columns: (x,y) -> (y, x+y).

(a) Q^n = F_{n-1} I (mod F_n), Cassini => Q^{2n} = (-1)^n I (mod F_n); Pisano period
    pi(F_n) = order of Q mod F_n, compared with 4n and with the exact rule
    pi(F_n) = 2n (n even >= 4), 4n (n odd >= 5).
(b) N = 8 = F_6: cycle decomposition of Q on the 64 squares, vs colour (i+j mod 2),
    the classes mod 2, and the Chebyshev rings about the board centre.
(c) |Fix(Q^j on Z_N^2)| = prod_i gcd(d_i, N), d_i the invariant factors of
    O_j = Z[phi]/(phi^j - 1) from PingYou Prop. 7.3; N = 2..30, j = 1..36; plus the
    Smith normal form of Q^j - I vs Prop. 7.3 and the group structure of Fix.
(d) Q mod 2 permutes the three nonzero classes cyclically; Fibonacci leapers
    Q^n (0,1)^T = (F_n, F_{n+1}) are colour-preserving iff n = 1 (mod 3); -Q^3 = I mod 2.
Run:  python3 chessboard_weave_20261006_fib_catmap.py
"""
from math import gcd
from collections import Counter, defaultdict


def fib(k):
    if k < 0:
        return (-1) ** (k + 1) * fib(-k)
    a, b = 0, 1
    for _ in range(k):
        a, b = b, a + b
    return a


def lucas(k):
    a, b = 2, 1
    for _ in range(k):
        a, b = b, a + b
    return a


def mmul(A, B, N=None):
    C = ((A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]),
         (A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]))
    if N:
        C = tuple(tuple(x % N for x in row) for row in C)
    return C


def mpow(A, k, N=None):
    R = ((1, 0), (0, 1))
    for _ in range(k):
        R = mmul(R, A, N)
    return R


I2 = ((1, 0), (0, 1))
Q = ((0, 1), (1, 1))


def red(A, N):
    return tuple(tuple(x % N for x in row) for row in A)


def order_mod(A, N):
    if N == 1:
        return 1
    R, k = red(A, N), 1
    while R != red(I2, N):
        R, k = mmul(R, A, N), k + 1
    return k


print("== (a) Q^n mod F_n and the Pisano period of F_n ==")
ok_scalar = ok_cassini = ok_sq = True
for n in range(1, 61):
    Fn = fib(n)
    Qn = mpow(Q, n)
    ok_scalar &= red(Qn, Fn) == red(((fib(n - 1), 0), (0, fib(n - 1))), Fn)
    ok_cassini &= fib(n - 1) * fib(n + 1) - fib(n) ** 2 == (-1) ** n
    ok_sq &= red(mpow(Q, 2 * n), Fn) == red((((-1) ** n, 0), (0, (-1) ** n)), Fn)
print("n = 1..60: Q^n == F_{n-1} I mod F_n: %s; Cassini: %s; Q^{2n} == (-1)^n I mod F_n: %s"
      % (ok_scalar, ok_cassini, ok_sq))
print(" n   F_n  pi(F_n)  4n  pi | 4n  rule(2n even>=4, 4n odd>=5)")
allok = True
for n in range(1, 21):
    p = order_mod(Q, fib(n))
    rule = 2 * n if (n % 2 == 0 and n >= 4) else (4 * n if (n % 2 == 1 and n >= 5) else None)
    allok &= (4 * n) % p == 0 and (rule is None or rule == p)
    print("%2d %5d %6d %5d   %-5s  %s" % (n, fib(n), p, 4 * n, (4 * n) % p == 0, rule if rule else '-'))
print("all n <= 20 consistent:", allok, "(exceptions of the rule: pi(F_1)=pi(F_2)=1, pi(F_3)=pi(2)=3)")

print()
print("== (b) N = 8: cycle decomposition of Q on the 64 squares ==")
N = 8
SQ = [(i, j) for i in range(N) for j in range(N)]


def Qv(v, N=N):
    return (v[1] % N, (v[0] + v[1]) % N)


def ring(v):
    i, j = v
    return int(max(abs(i - 3.5), abs(j - 3.5)) - 0.5)


def colour(v):
    return (v[0] + v[1]) % 2


seen, orbits = set(), []
for v in SQ:
    if v in seen:
        continue
    orb, w = [], v
    while w not in seen:
        seen.add(w)
        orb.append(w)
        w = Qv(w)
    orbits.append(orb)
orbits.sort(key=lambda o: (len(o), min(o)))
print("orbit length multiset:", sorted(Counter(len(o) for o in orbits).items()), " #orbits:", len(orbits),
      " order of Q mod 8:", order_mod(Q, 8))
for o in orbits:
    cc = Counter(colour(v) for v in o)
    m2 = Counter((v[0] % 2, v[1] % 2) for v in o)
    rc = Counter(ring(v) for v in o)
    print("  len %2d start %s: colour(0,1) = (%d,%d); mod-2 classes %s; rings %s"
          % (len(o), o[0], cc[0], cc[1], dict(sorted(m2.items())), dict(sorted(rc.items()))))
    if len(o) <= 12:
        print("      orbit:", o)
# alignment tests
orbit_of = {v: k for k, o in enumerate(orbits) for v in o}
ring_sets = defaultdict(set)
for v in SQ:
    ring_sets[ring(v)].add(v)
maps_rings = {r: sorted(Counter(ring(Qv(v)) for v in ring_sets[r]).items()) for r in range(4)}
print("image of ring r under Q, ring histogram:", maps_rings)
print("ring constant on some orbit of length > 1:", any(len({ring(v) for v in o}) == 1 for o in orbits if len(o) > 1))
# join of the orbit partition and the ring partition (connected components of orbit-ring incidence)
parent = list(range(len(orbits) + 4))
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x
for v in SQ:
    a, b = find(orbit_of[v]), find(len(orbits) + ring(v))
    parent[a] = b
print("blocks of the join (finest partition coarser than both orbits and rings):",
      len({find(k) for k in range(len(orbits))}))
# same for the colour partition
parent = list(range(len(orbits) + 2))
for v in SQ:
    a, b = find(orbit_of[v]), find(len(orbits) + colour(v))
    parent[a] = b
print("blocks of the join of orbits and colours:", len({find(k) for k in range(len(orbits))}))
# affine cat maps v -> Qv + w commuting with the central half-turn h(v) = (7,7) - v
def half_turn(v):
    return ((7 - v[0]) % 8, (7 - v[1]) % 8)
def affine(v, w):
    u = Qv(v)
    return ((u[0] + w[0]) % 8, (u[1] + w[1]) % 8)
comm = [w for w in SQ if all(affine(half_turn(v), w) == half_turn(affine(v, w)) for v in SQ)]
print("translations w with (v -> Qv+w) commuting with the central half-turn v -> (7,7)-v:", comm,
      "(needs 2w = (Q-I)(1,1) = (0,1) mod 8: impossible)")
print("colour of Q(v) equals column parity of v (x mod 2):", all(colour(Qv(v)) == v[0] % 2 for v in SQ))
print("Q^3 and -Q^3 preserve every colour class and every class mod 2 on Z_8^2:",
      all(((v[0] + 2 * v[1]) % 2, (2 * v[0] + 3 * v[1]) % 2) == (v[0] % 2, v[1] % 2) for v in SQ))

print()
print("== (c) |Fix(Q^j on Z_N^2)| = prod gcd(d_i, N) with the invariant factors of O_j (Prop. 7.3) ==")


def paper_factors(j):
    if j % 2 == 1 and j % 3 != 0:
        return (1, lucas(j))
    if j % 2 == 1:
        return (2, lucas(j) // 2)
    if j % 4 == 2:
        return (lucas(j // 2), lucas(j // 2))
    return (fib(j // 2), 5 * fib(j // 2))


def snf2(A):
    d1 = gcd(gcd(A[0][0], A[0][1]), gcd(A[1][0], A[1][1]))
    det = abs(A[0][0] * A[1][1] - A[0][1] * A[1][0])
    return (d1, det // d1)


bad_snf = []
for j in range(1, 61):
    M = mpow(Q, j)
    A = ((M[0][0] - 1, M[0][1]), (M[1][0], M[1][1] - 1))
    det = abs(A[0][0] * A[1][1] - A[0][1] * A[1][0])
    if snf2(A) != paper_factors(j) or det != abs(lucas(j) - 1 - (-1) ** j):
        bad_snf.append(j)
print("SNF(Q^j - I) == Prop. 7.3 factors and |O_j| = |L_j - 1 - (-1)^j| for j = 1..60; failures:", bad_snf)


def exponent_and_order(H, N):
    # the order of (x,y) in Z_N^2 is N / gcd(N, x, y)
    return max(N // gcd(N, gcd(v[0], v[1])) for v in H), len(H)


bad, checked = [], 0
for Nn in range(2, 31):
    for j in range(1, 37):
        M = mpow(Q, j, Nn)
        A = ((M[0][0] - 1) % Nn, M[0][1]), (M[1][0], (M[1][1] - 1) % Nn)
        H = [(x, y) for x in range(Nn) for y in range(Nn)
             if (A[0][0] * x + A[0][1] * y) % Nn == 0 and (A[1][0] * x + A[1][1] * y) % Nn == 0]
        d1, d2 = paper_factors(j)
        g1, g2 = gcd(d1, Nn), gcd(d2, Nn)
        e, o = exponent_and_order(H, Nn)
        checked += 1
        if o != g1 * g2 or e != max(g1, g2) or g2 % g1:
            bad.append((Nn, j, o, g1 * g2, e, (g1, g2)))
print("pairs (N,j) checked:", checked, "; |Fix| and group type Z/gcd(d1,N) + Z/gcd(d2,N) mismatches:", bad)
print("chessboard N = 8: |Fix(Q^j)| for j = 1..24:",
      [gcd(paper_factors(j)[0], 8) * gcd(paper_factors(j)[1], 8) for j in range(1, 25)])

print()
print("== (d) Q mod 2 and the Fibonacci leapers ==")
v = (1, 0)
cyc = [v]
for _ in range(3):
    v = Qv(v, 2)
    cyc.append(v)
print("Q mod 2 on the nonzero classes:", ' -> '.join(map(str, cyc)), "; order of Q mod 2:", order_mod(Q, 2))
names = {(0, 1): 'rook/wazir', (1, 1): 'bishop/ferz', (1, 2): 'knight', (2, 3): 'zebra'}
for n in range(0, 13):
    a, b = fib(n), fib(n + 1)
    assert mpow(Q, n) == ((fib(n - 1), a), (a, b))
    print("  n=%2d Q^n(0,1) = (%3d,%3d) %-12s mod 2 = %s  colour-preserving: %-5s  n = 1 mod 3: %s"
          % (n, a, b, names.get((a, b), ''), (a % 2, b % 2), (a + b) % 2 == 0, n % 3 == 1))
Q3 = mpow(Q, 3)
mQ3 = tuple(tuple(-x for x in row) for row in Q3)
print("Q^3 =", Q3, "(columns: knight (1,2), zebra (2,3));  -Q^3 mod 2 =", red(mQ3, 2), "= I")
