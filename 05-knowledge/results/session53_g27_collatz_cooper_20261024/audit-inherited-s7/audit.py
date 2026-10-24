# Standard-library independent audit. Checks remain active under python -O.
from math import comb
N = 180
one = [1] + [0] * N

def add(a, b):
    return [x + y for x, y in zip(a, b)]

def scale(a, k):
    return [k * x for x in a]

def mul(a, b):
    c = [0] * (N + 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b[:N + 1 - i]):
                c[i + j] += x * y
    return c

def power(a, k):
    c = one
    while k:
        if k & 1:
            c = mul(c, a)
        a = mul(a, a)
        k //= 2
    return c

def inv(a):
    if not a[0] == 1:
        raise AssertionError('audit check failed (original source line 19)')
    c = [1] + [0] * N
    for n in range(1, N + 1):
        c[n] = -sum((a[k] * c[n - k] for k in range(1, n + 1)))
    return c

def dil(a, p):
    return [a[n // p] if n % p == 0 else 0 for n in range(N + 1)]

def chi(n):
    return 0 if n % 7 == 0 else 1 if n % 7 in (1, 2, 4) else -1
th = [1] + [0] * N
for d in range(1, N + 1):
    for n in range(d, N + 1, d):
        th[n] += 2 * chi(d)
euler = one[:]
for d in range(1, N + 1):
    for n in range(N, d - 1, -1):
        euler[n] -= euler[n - d]
h = [0] + power(mul(dil(euler, 7), inv(euler)), 4)[:-1]
den = add(add(one, scale(h, 13)), scale(mul(h, h), 49))
t = mul(h, inv(den))
s = mul(add(one, scale(mul(h, h), -49)), inv(den))
eta = [0] + power(mul(euler, dil(euler, 7)), 3)[:-1]
if not mul(t, power(th, 3)) == eta:
    raise AssertionError('audit check failed (original source line 35)')
if not mul(s, s) == add(add(one, scale(t, -26)), scale(mul(t, t), -27)):
    raise AssertionError('audit check failed (original source line 36)')
if not [n * t[n] for n in range(N + 1)] == mul(mul(t, power(th, 2)), s):
    raise AssertionError('audit check failed (original source line 37)')
print('S4 S5 S6 exact q-series identities checked through q^180')
C7 = mul(power(th, 6), power(add(one, t), 2))
if not all(((v - (i == 0)) % 7 == 0 for i, v in enumerate(C7))):
    raise AssertionError('audit check failed (original source line 40)')
print('C7 mod7 Sturm input through 180, bound4:', C7[:5])
f = mul(power(th, 3), s)
A1 = [0] * (N + 1)
A2 = A1[:]
for d in range(1, N + 1):
    for n in range(d, N + 1, d):
        A1[n] += chi(d) * d * d
        A2[n] += chi(n // d) * d * d
if not scale(f, 8) == add(add(scale(one, 8), scale(A1, -7)), scale(A2, -49)):
    raise AssertionError('audit check failed (original source line 46)')
print('D7 weight3 Eisenstein identity exact through180, bound2')
for p, mod, extra, weight, bound in [(2, 4, [(2, 1), (3, 2)], 12, 24), (7, 49, [(7, 1), (8, 14), (9, 28)], 30, 140)]:
    d = dil(t, p)
    for j, b in extra:
        d = add(d, scale(power(t, j), -b))
    out = mul(mul(power(th, weight - 3), power(dil(th, p), 3)), d)
    if not all((v % mod == 0 for v in out)):
        raise AssertionError('audit check failed (original source line 52)')
    print('D%d modular form check modulo%d through q^180, bound%d PASS' % (p, mod, bound))
c = [1]
for n in range(1, 181):
    F = sum((comb(n, k) ** 2 * comb(n + k, k) * comb(2 * k, n) for k in range((n + 1) // 2, n + 1)))
    num = F - sum((c[j] * c[n - j] for j in range(1, n)))
    if not num % 2 == 0:
        raise AssertionError('audit check failed (original source line 59)')
    c.append(num // 2)
for n in range(180):
    if not (n + 1) ** 2 * c[n + 1] == (26 * n * n + 13 * n + 2) * c[n] + (3 * (3 * n - 1) * (3 * n - 2) * c[n - 1] if n else 0):
        raise AssertionError('audit check failed (original source line 61)')
print('Normalized binomial-sum square root independently matches recurrence through180:', c[:7])
for n in range(180, 20000):
    num = (26 * n * n + 13 * n + 2) * c[n] + 3 * (3 * n - 1) * (3 * n - 2) * c[n - 1]
    if not num % (n + 1) ** 2 == 0:
        raise AssertionError('audit check failed (original source line 66)')
    c.append(num // (n + 1) ** 2)
primes = [p for p in range(2, 501) if all((p % d for d in range(2, int(p ** 0.5) + 1)))]
for p in primes:
    good = p == 7 or chi(p) == 1
    bad = next((n for n in range(1, 20000 // p + 1) if (c[p * n] - c[n]) % (p * p)), None)
    if not (bad is None) == good:
        raise AssertionError('audit check failed (original source line 71)')
    if not good:
        if not (c[p] - 2) % p != 0:
            raise AssertionError('audit check failed (original source line 72)')
print('TheoremD finite range PASS: all primes <=500, all pn<=20000; every inert prime fails n=1 modp')
print('Exact small inert witnesses:', [(p, c[p], (c[p] - 2) % p) for p in (3, 5, 13)])
val = [0] * (N + 1)
for a in reversed(c[:N + 1]):
    val = mul(val, t)
    val[0] += a
if not val == th:
    raise AssertionError('audit check failed (original source line 78)')
print('G(t(q))=theta(q) independently checked through q^180 from binomial F square root')

def digit9(n):
    v = 1
    while n:
        n, d = divmod(n, 9)
        if d > 2:
            return 0
        if d == 1:
            v = -v
    return v % 3
if not all((x % 3 == digit9(n) for n, x in enumerate(c))):
    raise AssertionError('audit check failed (original source line 87)')
if not all(((c[9 * n] - c[n]) % 3 == 0 for n in range(len(c) // 9))):
    raise AssertionError('audit check failed (original source line 88)')
print('New p3 base9 digit product and c(9n)=c(n) independently PASS through n=20000')
