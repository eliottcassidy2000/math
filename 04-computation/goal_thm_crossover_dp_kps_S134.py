# Rigorous UPPER-bound DP for o(m) = max odd-order subgroup of S_m, vs Busch floor f(m).
# Upper bounds are safe: primitive odd => solvable (Feit-Thompson) => affine => degree p^d,
# |G| <= p^d * oddpart(|GL(d,p)|); prime degree = Burnside p*oddpart(p-1).
# Transitive imprimitive: |G| <= O(d)^(n/d) * T(n/d)  (kernel + block-image bounds).
# Intransitive: product over odd orbit sizes.
from functools import lru_cache

def oddpart(x):
    while x % 2 == 0:
        x //= 2
    return x

def isprime(n):
    if n < 2:
        return False
    p = 2
    while p * p <= n:
        if n % p == 0:
            return False
        p += 1
    return True

def primepower(n):
    for p in range(2, n + 1):
        if n % p == 0:
            x, d = n, 0
            while x % p == 0:
                x //= p
                d += 1
            return (p, d) if x == 1 else None
    return None

def glorder(p, d):
    o = 1
    for i in range(d):
        o *= (p ** d - p ** i)
    return o

def primUB(n):
    pp = primepower(n)
    if pp is None:
        return 0
    p, d = pp
    if d == 1:
        return n * oddpart(n - 1)
    return n * oddpart(glorder(p, d))

@lru_cache(None)
def T(n):  # upper bound, odd transitive subgroup of S_n
    if n == 1:
        return 1
    if n % 2 == 0:
        return 0
    best = primUB(n)
    d = 3
    while d * d <= n * n and d < n:
        if n % d == 0 and d > 1:
            q = n // d
            if q > 1:
                cand = O(d) ** q * T(q)
                if cand > best:
                    best = cand
        d += 2
    return best

@lru_cache(None)
def O(m):  # upper bound, any odd subgroup of S_m
    if m == 0:
        return 1
    best = 0
    for j in range(1, m + 1, 2):
        tj = T(j)
        if tj:
            c = tj * O(m - j)
            if c > best:
                best = c
    return best

def f(m):  # Busch floor: min 3^a 5^b, 2a+3b = m-1
    best = None
    for b in range(0, (m - 1) // 3 + 1):
        r = (m - 1) - 3 * b
        if r % 2 == 0:
            v = 3 ** (r // 2) * 5 ** b
            if best is None or v < best:
                best = v
    return best

cross = []
for m in range(3, 62):
    fm, om = f(m), O(m)
    star = " <-- CROSSOVER (o>=f)" if om >= fm else ""
    if om >= fm:
        cross.append(m)
    print(f"m={m:2d}  f={fm:>15d}  oUB={om:>15d}{star}")
print("crossovers 3..61:", cross)

# special checks
# m=27: best multi-orbit (intransitive) product
best_intr27 = max(T(j) * O(27 - j) for j in range(1, 27, 2) if T(j))
print("m=27: max intransitive odd order =", best_intr27, "(3^12 =", 3 ** 12, "), f(27) =", f(27), ", T(27) =", T(27))
# m=54: best partition product excluding 27+27
best54 = 0
def best_excl_2727(m):
    # partitions into odd parts, exclude exactly [27,27]
    best = 0
    def rec(rem, maxp, prod, used27):
        nonlocal best
        if rem == 0:
            if not (used27 == 2 and prod == T(27) ** 2):
                best = max(best, prod)
            return
        j = min(maxp, rem)
        if j % 2 == 0:
            j -= 1
        while j >= 1:
            if T(j):
                rec(rem - j, j, prod * T(j), used27 + (1 if j == 27 else 0))
            j -= 2
    rec(m, m, 1, 0)
    return best
b54 = best_excl_2727(54)
print("m=54: best odd order with orbit shape != (27,27) =", b54, ", f(54) =", f(54), ", 3^26 =", 3 ** 26)
# m=9 window divisors
d9 = [d for d in range(75, 82) if 2835 % d == 0]
print("m=9: odd divisors of oddpart(9!)=2835 in [75,81]:", d9)
# growth check: primitive-affine bound k^(1+log3 k) vs 5^((k-1)/3) for k in 25..80
import math
bad = [k for k in range(25, 200) if (1 + math.log(k, 3)) * math.log(k) >= (k - 1) / 3 * math.log(5)]
print("k in 25..199 where affine bound k^(1+log3 k) NOT < 5^((k-1)/3):", bad)
