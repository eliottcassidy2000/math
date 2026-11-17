"""Independent subset-DP controls and constructive smaller-block selectors."""
from collections import deque
from functools import cache
from math import gcd, isqrt, lcm
from itertools import product

# (Hamiltonian paths, optional two-cover exponent, outgoing row masks)
TILES = (
    (5, 5, (2, 4, 9, 3)),
    (19, 8, (40, 1, 3, 6, 15, 30)),
    (223, 12, (2, 40, 11, 49, 71, 213, 143, 31)),
    (599, 13, (96, 21, 25, 51, 225, 198, 142, 15)),
)
FIVE_FACTOR = (24, 1, 3, 6, 14)
A = (4, 21, 208, 39, 105, 71, 139, 59)
B = (128, 129, 27, 51, 99, 71, 15, 124)


def need(ok):
    if not ok:
        raise AssertionError("failed exact replay")


def counts(rows, max_q):
    n = len(rows)
    f = (1 << n) - 1
    for i, w in enumerate(rows):
        need(w <= f and not (w >> i & 1))
        for j in range(i):
            need((w >> j & 1) + (rows[j] >> i & 1) == 1)
    dp = [[0] * n for _ in range(f + 1)]
    for i in range(n):
        dp[1 << i][i] = 1
    for m in range(1, f + 1):
        for i, v in enumerate(dp[m]):
            if v == 0:
                continue
            rem = rows[i] & ~m
            while rem:
                bit = rem & -rem
                rem ^= bit
                dp[m ^ bit][bit.bit_length() - 1] += v
    h = [1] + [sum(row) for row in dp[1:]]
    @cache
    def cover(m, q):
        if q == 1:
            return h[m]
        r, out = m, 0
        while True:
            out += h[r] * cover(m ^ r, q - 1)
            if not r:
                return out
            r = (r - 1) & m
    return tuple(cover(f, q) for q in range(1, max_q + 1))


def wrap(rows):
    n = len(rows)
    return tuple(x | (27 << n) for x in rows) + tuple(
        (x << n) | (((1 << n) - 1) if j == 2 else 0)
        for j, x in enumerate((18, 20, 17, 7, 8))
    )


def ordinal(*tiles):
    out, size = [], 0
    for tile in reversed(tiles):
        out = [(r << len(tile)) for r in out]
        out = [r | (((1 << size) - 1) << len(tile)) for r in tile] + out
        size += len(tile)
    return tuple(out)


def exact_dyadic(N):
    if N < 124:
        return None
    e = (N // 62).bit_length() - 1
    u, r = divmod(N - 62 * (1 << e), 15)
    if r or u < 1:
        return None
    residual, cost, ys = u, 0, []
    for a, w, tile in TILES:
        copies = 0
        while residual % a == 0:
            residual //= a
            cost += w
            copies += 1
        ys.append(copies)
    if residual == 1 and cost <= e:
        return e, u, tuple(ys), e - cost
    return None


@cache
def base_wheel():
    M = 7 * 9 * 31
    choices = {1: (0,) * len(TILES)}
    q = deque([1])
    while q:
        x = q.popleft()
        for j, (a, _, _) in enumerate(TILES):
            y = x * a % M
            if y not in choices:
                path = list(choices[x])
                path[j] += 1
                choices[y] = tuple(path)
                q.append(y)
    return choices


def words_for_u(u, a, b, c):
    """Unit u mod 7^a 3^b 31^c, where a>=1,b>=2,c>=1."""
    m = 7**a * 3**b * 31**c
    need(gcd(u, m) == 1 and a >= 1 and b >= 2 and c >= 1)
    y = list(base_wheel()[u % 1953])
    val = 1
    for k, (g, _, _) in zip(y, TILES):
        val = val * pow(g, k, m) % m
    z = u * pow(val, -1, m) % m
    # Separate local discrete logs in the principal kernels.
    stride, log = 1, 0
    for p, depth, height in ((7, 1, a), (3, 2, b), (31, 1, c)):
        step, k = 1, 0
        for j in range(depth, height):
            pj = p**(j+1)
            digit = next(d for d in range(p)
                         if pow(5, 30*(k + d*step), pj) == z % pj)
            k += digit * step
            step *= p
        log += ((k - log)*pow(stride, -1, step) % step)*stride if step > 1 else 0
        stride *= step
    y[0] += 30*log
    v = 1
    for k, (g, _, _) in zip(y, TILES):
        v = v * pow(g, k, m) % m
    need(v == u % m)
    return tuple(y)


def local_for_target(N, a, b, c):
    """Return exact R exponent e and counts of dyadic tiles for a residue of N."""
    m = 7**a * 3**b * 31**c
    need(a >= 1 and b >= 3 and c >= 1)
    need(N % 7 == 0 and gcd(N, 93) == 1)
    d = next(e for e in range(6) if (N-62*2**e) % 9 in (3, 6))
    e = d + 6
    # Raising e by a multiple of lcm(6, ord_m(2)) preserves the demanded unit.
    # Euler's exponent is a safe multiple of ord_m(2).
    ph = (6*7**(a-1))*(2*3**(b-1))*(30*31**(c-1))
    u = ((N - 62*pow(2, e, m)) // 3) * pow(5, -1, m//3) % (m//3)
    y = words_for_u(u, a, b-1, c)
    cost = sum(k*w for k,(_,w,_) in zip(y, TILES))
    step = lcm(6, ph)
    e += max(0, (cost-e+step-1)//step)*step
    v = 1
    for k, (g, _, _) in zip(y, TILES):
        v = v * pow(g, k, m) % m
    need((15*v + 62*pow(2,e,m) - N) % m == 0)
    return e, y, e-cost


def run():
    for h, e, rows in TILES:
        need(counts(rows,2) == (h, 2**e))
    need(counts(FIVE_FACTOR,2) == (11, 100))
    need(counts(A,3) == (345,5640,45225))
    need(counts(B,3) == (345,5640,44505))
    need(counts(wrap(A),2) == (79995,4151100))
    need(counts(wrap(B),2) == (79995,4089180))
    need(len(base_wheel()) == 1080)
    for p,d in ((3,2),(7,1),(31,1)):
        need(pow(5,30,p**d)==1 and pow(5,30,p**(d+1))!=1)
    for a,b,c in ((1,3,1),(3,5,3),(4,7,3)):
        for p in range(1,100,2):
            if gcd(p,93)==1:
                local_for_target(7*p,a,b,c)
    P=2377693121
    need(all(P%d for d in range(2,isqrt(P)+1)))
    need(exact_dyadic(7*P)==(28,56905,(1,1,0,1),2))
    w = [int(x) for x in "100110010"]
    t = [(1 << (len(w)+1)) - (1 << (i+2)) + (not b)
         for i,b in enumerate(w)]
    need(counts((sum(b<<(i+1) for i,b in enumerate(w)),*t),1)==(217,))
    print("literal profiles, exact residue wheel/lifts, feedback fork, and target inverse: PASS")


if __name__ == "__main__":
    run()
