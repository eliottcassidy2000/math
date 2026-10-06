#!/usr/bin/env python3
"""Independent audit (written from scratch) of Proposition 4.3 of
05-knowledge/results/chessboard_weave_20261006.md (cat map Q = [[0,1],[1,1]] on torus boards).
Exact integer arithmetic.
"""
from math import gcd
from collections import Counter

fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)
        print("FAIL:", msg)


def mul(A, B, M=None):
    C = [[A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]],
         [A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]]]
    if M:
        C = [[x % M for x in r] for r in C]
    return C


def mpow(A, n, M=None):
    R = [[1, 0], [0, 1]]
    for _ in range(n):
        R = mul(R, A, M)
    return R


Q = [[0, 1], [1, 1]]
F = [0, 1]
for _ in range(100):
    F.append(F[-1] + F[-2])
I = [[1, 0], [0, 1]]

# Cassini consequences
for n in range(3, 41):
    M = F[n]
    Qn = mpow(Q, n, M)
    check(Qn == [[F[n - 1] % M, 0], [0, F[n - 1] % M]], f"Q^n = F_(n-1) I mod F_n, n={n}")
    Q2n = mpow(Q, 2 * n, M)
    s = (-1) ** n % M
    check(Q2n == [[s, 0], [0, s]], f"Q^2n = (-1)^n I mod F_n, n={n}")
print("Q^n = F_{n-1} I and Q^{2n} = (-1)^n I (mod F_n): checked n=3..40")


def order_mod(M):
    if M == 1:
        return 1
    R = Q
    k = 1
    R = [[x % M for x in r] for r in R]
    while R != [[1 % M, 0], [0, 1 % M]]:
        R = mul(R, Q, M)
        k += 1
    return k


pis = {}
for n in range(3, 31):
    pis[n] = order_mod(F[n])
    if n >= 4 and n % 2 == 0:
        check(pis[n] == 2 * n, f"pi(F_{n}) = 2n")
    if n >= 5 and n % 2 == 1:
        check(pis[n] == 4 * n, f"pi(F_{n}) = 4n")
    check((4 * n) % pis[n] == 0, f"pi(F_{n}) | 4n")
print("Pisano periods pi(F_n), n=3..30:", pis)

# 8x8 board = (Z/8)^2
M = 8
pts = [(i, j) for i in range(M) for j in range(M)]


def Qv(v, M=8):
    return (v[1] % M, (v[0] + v[1]) % M)


check(order_mod(8) == 12, "order of Q mod 8 = 12")
seen = set()
orbits = []
for v in pts:
    if v in seen:
        continue
    orb = [v]
    w = Qv(v)
    while w != v:
        orb.append(w)
        w = Qv(w)
    seen |= set(orb)
    orbits.append(orb)
lens = sorted(len(o) for o in orbits)
print("orbit lengths of Q on (Z/8)^2:", lens)
check(lens == [1, 3, 6, 6, 12, 12, 12, 12], "orbit lengths")
orb01 = [o for o in orbits if (0, 1) in o][0]
k = orb01.index((0, 1))
orb01 = orb01[k:] + orb01[:k]
print("orbit of (0,1):", orb01)
check(mpow(Q, 6, 8) == [[5, 0], [0, 5]], "Q^6 = 5I mod 8")
# x^2 - x - 1 irreducible mod 2
check(all((x * x - x - 1) % 2 != 0 for x in range(2)), "2 inert")


# Fix(Q^j) on Z_N^2 vs invariant factors of Q^j - I
def snf2(A):
    a, b, c, d = A[0][0], A[0][1], A[1][0], A[1][1]
    d1 = gcd(gcd(abs(a), abs(b)), gcd(abs(c), abs(d)))
    det = abs(a * d - b * c)
    return (d1, det // d1)


fix8 = []
for j in range(1, 25):
    Qj = mpow(Q, j)
    cnt = sum(1 for (x, y) in pts if ((Qj[0][0] * x + Qj[0][1] * y - x) % 8 == 0 and (Qj[1][0] * x + Qj[1][1] * y - y) % 8 == 0))
    fix8.append(cnt)
print("|Fix(Q^j) on Z_8^2|, j=1..24:", fix8)
check(fix8[:12] == [1, 1, 4, 1, 1, 16, 1, 1, 4, 1, 1, 64], "Fix sequence j=1..12")
check(fix8[12:] == fix8[:12], "period 12")
Q12mI = [[x - y for x, y in zip(r, s)] for r, s in zip(mpow(Q, 12), I)]
print("Q^12 - I =", Q12mI, " invariant factors:", snf2(Q12mI))
check(snf2(Q12mI) == (8, 40), "O_12 = Z/8 + Z/40")
npairs = 0
for N in range(2, 31):
    P = [(i, j) for i in range(N) for j in range(N)]
    for j in range(1, 37):
        Qj = mpow(Q, j, N)
        cnt = sum(1 for (x, y) in P if ((Qj[0][0] * x + Qj[0][1] * y - x) % N == 0 and (Qj[1][0] * x + Qj[1][1] * y - y) % N == 0))
        QjZ = mpow(Q, j)
        d1, d2 = snf2([[QjZ[0][0] - 1, QjZ[0][1]], [QjZ[1][0], QjZ[1][1] - 1]])
        check(cnt == gcd(d1, N) * gcd(d2, N), f"Fix count N={N} j={j}")
        npairs += 1
print("Fix(Q^j on Z_N^2) = prod gcd(d_i,N): checked", npairs, "pairs (N=2..30, j=1..36)")


# rings on the 8x8 board, squares labelled (i,j) in 0..7
def ring(i, j, N=8):
    return (max(abs(2 * i - (N - 1)), abs(2 * j - (N - 1))) - 1) // 2


one_ring = [o for o in orbits if len(o) > 1 and len(set(ring(*v) for v in o)) == 1]
print("orbits of length >1 inside a single ring:", one_ring)
check(one_ring == [], "no orbit of length >1 stays in one ring")
for kk in range(4):
    imgs = set(ring(*Qv(v)) for v in pts if ring(*v) == kk)
    check(len(imgs) > 1, f"Q(ring {kk}) not in a single ring")
    print(f"Q(ring {kk}) meets rings {sorted(imgs)}")
# join of orbit partition and ring partition
parent = {v: v for v in pts}


def find(v):
    while parent[v] != v:
        parent[v] = parent[parent[v]]
        v = parent[v]
    return v


def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[ra] = rb


for o in orbits:
    for v in o[1:]:
        union(o[0], v)
for kk in range(4):
    rk = [v for v in pts if ring(*v) == kk]
    for v in rk[1:]:
        union(rk[0], v)
blocks = len(set(find(v) for v in pts))
print("blocks of the join (orbits v rings):", blocks)
check(blocks == 1, "join is whole board")
# affine maps v -> Qv + w commuting with the half-turn v -> (7,7) - v
good = []
for w in pts:
    ok = True
    for v in pts:
        a = Qv(((7 - v[0]) % 8, (7 - v[1]) % 8))
        a = ((a[0] + w[0]) % 8, (a[1] + w[1]) % 8)
        b = Qv(v)
        b = ((b[0] + w[0]) % 8, (b[1] + w[1]) % 8)
        b = ((7 - b[0]) % 8, (7 - b[1]) % 8)
        if a != b:
            ok = False
            break
    if ok:
        good.append(w)
print("w with (v -> Qv+w) commuting with the central half-turn:", good)
check(good == [], "no affine Q-map commutes with half-turn")
# colours (scaffolds): does Q preserve colour i+j mod 2? does Q^3?
chg1 = sum(1 for v in pts if (sum(v) - sum(Qv(v))) % 2)
Q3 = mpow(Q, 3, 8)
chg3 = sum(1 for (x, y) in pts if ((x + y) - ((Q3[0][0] * x + Q3[0][1] * y) + (Q3[1][0] * x + Q3[1][1] * y))) % 2)
par3 = sum(1 for (x, y) in pts if ((Q3[0][0] * x + Q3[0][1] * y - x) % 2, (Q3[1][0] * x + Q3[1][1] * y - y) % 2) != (0, 0))
print(f"squares whose colour Q changes: {chg1}/64; colour changed by Q^3: {chg3}/64; parity class moved by Q^3: {par3}/64")
print()
print("TOTAL FAILURES:", len(fails))
