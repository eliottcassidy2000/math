from itertools import product
from math import gcd

def mul(z, w):
    a, b = z
    c, d = w
    return (a * c - b * d, a * d + b * c)

def powc(z, n):
    a = (1, 0)
    for _ in range(n):
        a = mul(a, z)
    return a

def zeck(n):
    f = [1, 2]
    while f[-1] <= n:
        f.append(sum(f[-2:]))
    out = []
    for x in reversed(f):
        if x <= n:
            out.append(x)
            n -= x
    return out
for N, factors in [(125, [(2, 1, 3)]), (5525, [(2, 1, 2), (3, 2, 1), (4, 1, 1)])]:
    classes = {}
    r = len(factors)
    for mask in range(1 << r):
        z = (1, 0)
        for j, (a, b, e) in enumerate(factors):
            z = mul(z, powc((a, b if mask & 1 << j else -b), e))
        x, y = map(abs, z)
        m, n = sorted((x, y))
        tri = (n * n - m * m, 2 * m * n, N)
        if not (tri[0] ** 2 + tri[1] ** 2 == N * N and gcd(*tri[:2]) == 1):
            raise AssertionError()
        key = min(mask, mask ^ (1 << r) - 1)
        classes.setdefault(key, tri)
        if not classes[key] == tri:
            raise AssertionError()
    print(N, classes)
    for i, t in classes.items():
        print(f'{i:0{r}b}', t, tuple((zeck(x) for x in t)))
factors = [(2, 1, 2), (3, 2, 1), (8, 7, 1)]
z = (1, 0)
for j, (a, b, e) in enumerate(factors):
    z = mul(z, powc((a, b if j == 2 else -b), e))
x, y = map(abs, z)
m, n = sorted((x, y))
A, B, C = (n * n - m * m, 2 * m * n, 36725)
from math import isqrt
s, t = (isqrt(C + B), isqrt(C - B))
if not (A, B, C, s, t) == (813, 36716, 36725, 271, 3):
    raise AssertionError()
k = 0
num = 3 * s + t
while num % 2 == 0:
    k += 1
    num //= 2
g = gcd(num, t)
sn, tn = (num // g, t // g)
if not (k, g, sn, tn) == (4, 3, 17, 1):
    raise AssertionError()
if not (sn * tn, (sn * sn - tn * tn) // 2, (sn * sn + tn * tn) // 2) == (17, 144, 145):
    raise AssertionError()
print('guarded ternary return verified:', (A, B, C), '->', (17, 144, 145))
